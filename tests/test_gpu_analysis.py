from __future__ import annotations

import subprocess
import sys
import importlib.util
from types import SimpleNamespace

import numpy as np
import pytest

from afdb_integration_kit.gpu import (
    MAX_HEAVY_ATOMS,
    Protein,
    empty_protein,
    result_to_clash_schema,
    result_to_interface_schema,
)
from afdb_integration_kit.gpu._runtime import resolve_device


TORCH_AVAILABLE = importlib.util.find_spec("torch") is not None


def _run_python(script: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )


def _make_test_protein() -> Protein:
    coords = np.zeros((4, MAX_HEAVY_ATOMS, 3), dtype=np.float32)
    mask = np.zeros((4, MAX_HEAVY_ATOMS), dtype=bool)
    atom_names = np.full((4, MAX_HEAVY_ATOMS), "", dtype="<U4")

    res_names = np.array(["ALA", "GLY", "ALA", "GLY"], dtype="<U3")
    chain_ids = np.array(["A", "A", "B", "B"], dtype="<U1")
    res_ids = np.array([10, 11, 20, 21], dtype=np.int32)

    for residue_index in range(4):
        atom_names[residue_index, :4] = ["N", "CA", "C", "O"]
        mask[residue_index, :4] = True

    coords[0, 1] = [0.0, 0.0, 0.0]
    coords[1, 1] = [20.0, 0.0, 0.0]
    coords[2, 1] = [4.0, 0.0, 0.0]
    coords[3, 1] = [0.5, 0.0, 0.0]

    coords[0, 0] = [0.0, -1.0, 0.0]
    coords[0, 2] = [0.0, 1.0, 0.0]
    coords[0, 3] = [0.0, 2.0, 0.0]

    coords[1, 0] = [20.0, -1.0, 0.0]
    coords[1, 2] = [20.0, 1.0, 0.0]
    coords[1, 3] = [20.0, 2.0, 0.0]

    coords[2, 0] = [4.0, -1.0, 0.0]
    coords[2, 2] = [4.0, 1.0, 0.0]
    coords[2, 3] = [4.0, 2.0, 0.0]

    coords[3, 0] = [0.4, 0.0, 0.0]
    coords[3, 2] = [0.5, 1.0, 0.0]
    coords[3, 3] = [0.5, 2.0, 0.0]

    return Protein(
        coords=coords,
        mask=mask,
        atom_names=atom_names,
        res_names=res_names,
        chain_ids=chain_ids,
        res_ids=res_ids,
        path="AF-1234567890123456-test.pdb",
    )


def test_gpu_package_import_is_lightweight() -> None:
    script = """
import sys
import afdb_integration_kit.gpu as gpu
print(gpu.Protein.__name__)
print('torch' in sys.modules)
print('fastpdb' in sys.modules)
print(any(name == 'biotite' or name.startswith('biotite.') for name in sys.modules))
"""
    result = _run_python(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["Protein", "False", "False", "False"]


def test_gpu_parse_protein_reports_missing_optional_dependency() -> None:
    script = """
import importlib
real_import_module = importlib.import_module

def fake_import_module(name, package=None):
    if name == "fastpdb":
        raise ModuleNotFoundError("No module named 'fastpdb'", name="fastpdb")
    return real_import_module(name, package)

importlib.import_module = fake_import_module
import afdb_integration_kit.gpu as gpu
protein = gpu.parse_protein('dummy.pdb')
print(protein.n_residues)
"""
    result = _run_python(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[-1] == "0"
    assert "fastpdb" in result.stdout
    assert "uv pip install '.[production]'" in result.stdout


def test_gpu_analysis_entry_points_report_missing_torch() -> None:
    script = """
import importlib
real_import_module = importlib.import_module

def fake_import_module(name, package=None):
    if name == "torch":
        raise ModuleNotFoundError("No module named 'torch'", name="torch")
    return real_import_module(name, package)

importlib.import_module = fake_import_module
import afdb_integration_kit.gpu as gpu
try:
    gpu.analyze_proteins
except ModuleNotFoundError as exc:
    print(exc)
"""
    result = _run_python(script)
    assert result.returncode == 0, result.stderr
    assert "torch" in result.stdout
    assert "uv pip install '.[production]'" in result.stdout


def test_resolve_device_auto_and_cuda_errors() -> None:
    fake_torch_cuda = SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: True))
    fake_torch_cpu = SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: False))

    assert resolve_device("auto", torch_module=fake_torch_cuda) == "cuda"
    assert resolve_device("auto", torch_module=fake_torch_cpu) == "cpu"
    assert resolve_device("cpu", torch_module=fake_torch_cuda) == "cpu"

    with pytest.raises(RuntimeError, match="CUDA is unavailable"):
        resolve_device("cuda", torch_module=fake_torch_cpu)

    with pytest.raises(ValueError, match="Unsupported device"):
        resolve_device("tpu", torch_module=fake_torch_cuda)


def test_schema_conversion_produces_expected_chain_annotations() -> None:
    interface_contact = SimpleNamespace(
        chain_1="A",
        res_1_id=10,
        aa_1_type="ALA",
        chain_2="B",
        res_2_id=20,
        aa_2_type="GLY",
        distance=5.25,
    )
    clash_contact = SimpleNamespace(
        chain_1="A",
        res_1_id=10,
        aa_1_type="ALA",
        atom_1_name="CA",
        chain_2="B",
        res_2_id=21,
        aa_2_type="GLY",
        atom_2_name="N",
        distance=0.6,
        overlap=1.0,
    )
    result = SimpleNamespace(
        path="AF-1234567890123456-example.pdb",
        interface_contacts=[interface_contact],
        n_backbone_clashes=1,
        n_heavy_clashes=1,
        backbone_clashing_residues=[{"res_id": 10, "chain_id": "A", "aa_type": "ALA"}],
        heavy_clashing_residues=[{"res_id": 21, "chain_id": "B", "aa_type": "GLY"}],
        interface_residues=[],
        backbone_clash_contacts=[clash_contact],
        heavy_clash_contacts=[clash_contact],
    )

    interface_schema = result_to_interface_schema(result, interface_cutoff=8.0)
    clash_schema = result_to_clash_schema(result)

    assert interface_schema is not None
    assert interface_schema["af_id"] == "AF-1234567890123456"
    assert interface_schema["sites"][0]["label"] == "AB_interface"
    assert (
        interface_schema["chains"][0]["residues"][0][
            "additional_residue_annotations"
        ]["partner_chain"]
        == "B"
    )

    assert clash_schema["af_id"] == "AF-1234567890123456"
    assert [site["label"] for site in clash_schema["sites"]] == [
        "backbone_clashes",
        "heavy_atom_clashes",
    ]
    assert (
        clash_schema["chains"][0]["residues"][0]["site_data"][0][
            "confidence_classification"
        ]
        == "high"
    )


@pytest.mark.skipif(
    not TORCH_AVAILABLE,
    reason="PyTorch is not installed",
)
def test_cpu_analysis_and_radius_graph_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from afdb_integration_kit.gpu import (
        analyze_proteins,
        count_clashes,
        get_clash_pairs,
    )
    import afdb_integration_kit.gpu.clashes as clashes_module

    protein = _make_test_protein()

    monkeypatch.setattr(clashes_module, "_RADIUS_GRAPH_USE_FALLBACK", True)
    monkeypatch.setattr(clashes_module, "_radius_graph_ext", None)

    counts = count_clashes([protein], selection="backbone", device="cpu", min_seq_sep=0)
    pairs = get_clash_pairs(protein, selection="backbone", device="cpu", min_seq_sep=0)
    results = analyze_proteins(
        [protein],
        batch_size=1,
        device="cpu",
        min_seq_sep=0,
        progress=False,
    )

    assert counts == [1]
    assert len(pairs) == 1
    assert pairs[0][0] == 10
    assert pairs[0][2] == 21

    assert len(results) == 1
    assert results[0].n_backbone_clashes == 1
    assert results[0].n_heavy_clashes >= 1
    assert {
        res["chain_id"] for res in results[0].interface_residues
    } == {"A", "B"}


def test_empty_protein_is_safe() -> None:
    protein = empty_protein("empty.pdb")
    assert protein.n_residues == 0
    assert protein.n_atoms == 0
