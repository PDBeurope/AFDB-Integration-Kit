"""Donor-independent tests for the synthetic fragment coordinate E2E."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"


def _load_script(name: str) -> ModuleType:
    if str(SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(SCRIPTS_DIR))
    path = SCRIPTS_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def synthesis() -> ModuleType:
    return _load_script("synthesize_fragment_metadata_e2e_assets")


@pytest.fixture(scope="module")
def coordinate_runner() -> ModuleType:
    return _load_script("generate_fragment_coordinate_e2e_example")


def _write_donor_pair(
    directory: Path,
    donor_id: str,
    chain_lengths: list[int],
    *,
    corrupt_dimensions: bool = False,
) -> tuple[Path, Path]:
    lines: list[str] = []
    serial = 1
    scores: list[float] = []
    for chain_index, length in enumerate(chain_lengths):
        chain_id = chr(ord("A") + chain_index)
        for residue_number in range(1, length + 1):
            score = float(10 * (len(scores) + 1))
            scores.append(score)
            for atom_offset, (atom, element) in enumerate(
                (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))
            ):
                x = float(residue_number * 4 + atom_offset)
                y = float(chain_index * 20)
                z = float(atom_offset)
                lines.append(
                    f"ATOM  {serial:5d} {atom:^4s} ALA {chain_id}"
                    f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00{score:6.2f}          {element:>2s}  "
                )
                serial += 1
        lines.append(
            f"TER   {serial:5d}      ALA {chain_id}{length:4d}"
        )
        serial += 1
    lines.append("END")
    pdb_path = directory / f"{donor_id}-model_v1.pdb"
    metadata_path = directory / f"{donor_id}-meta_v1.json"
    pdb_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    size = len(scores) + (1 if corrupt_dimensions else 0)
    pae = [[float(abs(i - j)) for j in range(size)] for i in range(size)]
    payload = {
        "plddt": scores,
        "pae": pae,
        "max_pae": float(size),
        "ptm": 0.5,
    }
    if len(chain_lengths) > 1:
        payload["iptm"] = 0.4
    metadata_path.write_text(json.dumps(payload), encoding="utf-8")
    return pdb_path, metadata_path


def test_crop_tile_alanine_and_score_remapping(
    tmp_path: Path, synthesis: ModuleType
) -> None:
    _write_donor_pair(tmp_path, "DONOR", [3])
    donor = synthesis.load_donor(tmp_path, "DONOR", 1)
    target = synthesis.TargetModel(
        "AF-TEST",
        (synthesis.TargetChain(
            "A", "PTEST", True, 10, 14, 5, "1", "Named fragment"
        ),),
    )

    pdb_lines, scores, provenance = synthesis.synthesize_model(target, donor)
    pdb_path = tmp_path / "AF-TEST-model_v1.pdb"
    metadata_path = tmp_path / "AF-TEST-meta_v1.json"
    pdb_path.write_text("\n".join(pdb_lines) + "\n", encoding="ascii")
    metadata_path.write_text(json.dumps(scores), encoding="utf-8")
    synthesis.validate_asset_pair(pdb_path, metadata_path, target)

    assert scores["plddt"] == [10.0, 20.0, 30.0, 10.0, 20.0]
    assert scores["pae"][3][4] == donor.pae[0][1]
    assert provenance["chains"][0]["tiled"] is True
    assert provenance["chains"][0]["tile_count"] == 2
    atom_lines = [line for line in pdb_lines if line.startswith("ATOM")]
    assert len(atom_lines) == 20
    assert {line[17:20] for line in atom_lines} == {"ALA"}
    first_tile = tuple(float(atom_lines[0][start:end]) for start, end in (
        (30, 38), (38, 46), (46, 54)
    ))
    second_tile = tuple(float(atom_lines[12][start:end]) for start, end in (
        (30, 38), (38, 46), (46, 54)
    ))
    assert first_tile != second_tile


def test_load_donor_rejects_dimension_mismatch(
    tmp_path: Path, synthesis: ModuleType
) -> None:
    _write_donor_pair(tmp_path, "BAD", [2], corrupt_dimensions=True)
    with pytest.raises(ValueError, match="PAE row count|columns"):
        synthesis.load_donor(tmp_path, "BAD", 1)


def test_target_lengths_use_full_sequence_or_inclusive_fragment(
    tmp_path: Path, synthesis: ModuleType
) -> None:
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,"
        "sequence_start,sequence_end,protein_name\n"
        "AF-FULL,1,A,P1,false,,,\n"
        "AF-FRAG,1,A,P1,true,3,7,Fragment\n",
        encoding="utf-8",
    )
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({
        "entries": [{"primary_ac": "P1", "sequence_length": 9}]
    }), encoding="utf-8")

    targets = synthesis.load_targets(manifest, seed)

    assert [target.chains[0].length for target in targets] == [9, 5]
    assert (targets[1].chains[0].sequence_start,
            targets[1].chains[0].sequence_end) == (3, 7)


def test_synthetic_writer_refuses_nonempty_output_before_writing(
    tmp_path: Path, synthesis: ModuleType
) -> None:
    _write_donor_pair(tmp_path, synthesis.MONOMER_DONOR_ID, [3])
    _write_donor_pair(tmp_path, synthesis.DIMER_DONOR_ID, [2, 2])
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,"
        "sequence_start,sequence_end\n"
        "AF-TEST,1,A,P1,false,,\n",
        encoding="utf-8",
    )
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({
        "entries": [{"primary_ac": "P1", "sequence_length": 2}]
    }), encoding="utf-8")
    output = tmp_path / "output"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("keep", encoding="utf-8")

    with pytest.raises(ValueError, match="not empty"):
        synthesis.write_synthetic_assets(
            tmp_path, manifest, seed, output, tmp_path / "provenance.json"
        )

    assert marker.read_text(encoding="utf-8") == "keep"


def test_mapping_validator_distinguishes_same_accession_fragments(
    tmp_path: Path,
    synthesis: ModuleType,
    coordinate_runner: ModuleType,
) -> None:
    mapping = tmp_path / "mapping.csv"
    fieldnames = [
        "model_entity_id", "entity_id", "chain_id", "uniprot_ac",
        "protein_name", "is_fragment", "sequence_start", "sequence_end",
    ]
    rows = []
    expected: dict[str, list[object]] = {}
    for model_index in range(10):
        model_id = f"AF-DUMMY-{model_index}"
        chain = synthesis.TargetChain(
            "A", f"P{model_index}", False, 1, 1, 1, "1", ""
        )
        expected[model_id] = [chain]
        rows.append({
            "model_entity_id": model_id, "entity_id": "1", "chain_id": "A",
            "uniprot_ac": f"P{model_index}", "protein_name": "",
            "is_fragment": "false", "sequence_start": "1",
            "sequence_end": "1",
        })
    model_id = "AF-FRAGMENT-HETERO"
    alpha = synthesis.TargetChain(
        "A", "P27409", True, 1, 46, 46, "1", ""
    )
    omega = synthesis.TargetChain(
        "B", "P27409", True, 961, 1071, 111, "2", ""
    )
    expected[model_id] = [alpha, omega]
    rows.extend([
        {
            "model_entity_id": model_id, "entity_id": "1", "chain_id": "A",
            "uniprot_ac": "P27409", "protein_name": "Development fragment alpha",
            "is_fragment": "true", "sequence_start": "1",
            "sequence_end": "46",
        },
        {
            "model_entity_id": model_id, "entity_id": "2", "chain_id": "B",
            "uniprot_ac": "P27409", "protein_name": "Development fragment omega",
            "is_fragment": "true", "sequence_start": "961",
            "sequence_end": "1071",
        },
    ])
    with mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    coordinate_runner.validate_mapping(mapping, expected)

    rows[-1]["entity_id"] = "1"
    with mapping.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(AssertionError, match="entity_id"):
        coordinate_runner.validate_mapping(mapping, expected)
