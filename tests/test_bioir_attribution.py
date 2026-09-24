"""Truthful, homogeneous BioIR attribution; no inference or network is required."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import duckdb
import gemmi
import pytest

from afdb_integration_kit.modelcif.generate import generate
from afdb_integration_kit.modelcif.provenance import (
    BIOIR_MULTIMER_TOOL,
    BIOIR_PTM_TOOL,
    BIOIR_SOFTWARE_NAME,
    normalize_modelcif_provenance,
)

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "uniprot/templates/bioir_modelcif_metadata.json"
LEGACY_TOOL = "ColabFold v1.6.0 / AlphaFold-Multimer"


def _module(relative: str):
    name = "bioir_test_" + Path(relative).stem
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def pipeline():
    return _module("scripts/production_pipeline.py")


def _config(pipeline, tmp_path, tool=BIOIR_MULTIMER_TOOL, **extra):
    values = dict(
        repo_dir=ROOT,
        input_dir=tmp_path,
        output_dir=tmp_path / "out",
        tool_used=tool,
        homodimer_tool_used=tool,
    )
    values.update(extra)
    return pipeline.Config(**values)


@pytest.mark.parametrize("tool", [BIOIR_PTM_TOOL, BIOIR_MULTIMER_TOOL])
@pytest.mark.parametrize("chain_count", [1, 2])
@pytest.mark.parametrize("version", ["?", "producer-version-fixture"])
def test_normalization_preserves_declared_bioir_method_version_and_groups(
    tool, chain_count, version
):
    payload = json.loads(TEMPLATE.read_text())
    payload["categories"]["_software"]["version"] = [version]
    payload["chains"] = [{"chain_id": chr(65 + i)} for i in range(chain_count)]
    normalize_modelcif_provenance(payload, prediction_tool=tool)
    normalized = copy.deepcopy(payload)
    normalize_modelcif_provenance(payload, allow_default_alphafold_version=True)
    assert payload == normalized
    assert payload["metadata"]["prediction_tool"] == tool
    software = payload["categories"]["_software"]
    assert software["name"][0] == BIOIR_SOFTWARE_NAME
    assert software["version"][0] == version
    assert payload["categories"]["_ma_protocol_step"]["details"][0].endswith(tool)
    assert set(payload["categories"]["_ma_software_group"]["software_id"]) <= set(
        software["pdbx_ordinal"]
    )


@pytest.mark.parametrize("conflict", ["marker", "legacy_row", "missing_marker"])
def test_normalization_rejects_conflicting_or_undeclared_predictor(conflict):
    payload = json.loads(TEMPLATE.read_text())
    tool = BIOIR_MULTIMER_TOOL
    if conflict == "marker":
        payload["metadata"]["prediction_tool"] = BIOIR_PTM_TOOL
    elif conflict == "legacy_row":
        payload["categories"]["_software"]["name"] = ["AlphaFold-Multimer"]
    else:
        tool = None
    with pytest.raises(ValueError):
        normalize_modelcif_provenance(payload, prediction_tool=tool)


@pytest.mark.parametrize("declared", [{}, [], 42, "", "unsupported"])
def test_malformed_prediction_declaration_has_an_actionable_error(declared):
    payload = json.loads(TEMPLATE.read_text())
    payload["metadata"]["prediction_tool"] = declared
    with pytest.raises(ValueError, match="metadata.prediction_tool"):
        normalize_modelcif_provenance(payload)


@pytest.mark.parametrize(
    "tool,other",
    [
        (BIOIR_MULTIMER_TOOL, LEGACY_TOOL),
        (BIOIR_MULTIMER_TOOL, BIOIR_PTM_TOOL),
        (LEGACY_TOOL, BIOIR_PTM_TOOL),
    ],
)
def test_pipeline_requires_matching_explicit_bioir_flags(
    pipeline, tmp_path, tool, other
):
    with pytest.raises(ValueError, match="same explicit BioIR method"):
        _config(pipeline, tmp_path, tool, homodimer_tool_used=other)


@pytest.mark.parametrize("tool", [BIOIR_MULTIMER_TOOL, LEGACY_TOOL])
def test_pipeline_rejects_template_and_dataset_conflicts(pipeline, tmp_path, tool):
    template = json.loads(TEMPLATE.read_text())
    template["metadata"]["prediction_tool"] = BIOIR_PTM_TOOL
    path = tmp_path / "template.json"
    path.write_text(json.dumps(template))
    with pytest.raises(ValueError, match="conflict"):
        _config(pipeline, tmp_path, tool, modelcif_template=path)
    dataset = tmp_path / "dataset.json"
    conflicting = LEGACY_TOOL if tool == BIOIR_MULTIMER_TOOL else BIOIR_PTM_TOOL
    dataset.write_text(
        json.dumps({"toolUsed": conflicting, "homodimerToolUsed": conflicting})
    )
    with pytest.raises(ValueError, match="Dataset|dataset"):
        _config(pipeline, tmp_path, tool, dataset_config=dataset)


@pytest.mark.parametrize(
    "tool,source,valid",
    [
        (BIOIR_PTM_TOOL, "openfold2_ptm_1", True),
        (BIOIR_MULTIMER_TOOL, "alphafold2_multimer_1", True),
        (BIOIR_MULTIMER_TOOL, "alphafold2_multimer_5", True),
        (BIOIR_MULTIMER_TOOL, "openfold2_ptm_1", False),
        (BIOIR_PTM_TOOL, "alphafold2_multimer_1", False),
        (BIOIR_MULTIMER_TOOL, None, False),
        (BIOIR_MULTIMER_TOOL, "unknown", False),
        (BIOIR_MULTIMER_TOOL, [], False),
        (LEGACY_TOOL, "alphafold2_multimer_1", False),
        (LEGACY_TOOL, None, True),
    ],
)
def test_pipeline_checks_actual_score_producer(pipeline, tmp_path, tool, source, valid):
    pdb = tmp_path / "AF-TEST-model_v1.pdb"
    scores = tmp_path / "AF-TEST-meta_v1.json"
    pdb.write_text("END\n")
    payload = {} if source is None else {"bioir_model_source": source}
    scores.write_text(json.dumps(payload))
    original = (pdb.read_bytes(), scores.read_bytes())
    config = _config(pipeline, tmp_path, tool)
    if valid:
        pipeline.validate_prediction_inputs(["AF-TEST"], config)
    else:
        with pytest.raises(ValueError, match="bioir_model_source"):
            pipeline.validate_prediction_inputs(["AF-TEST"], config)
    assert (pdb.read_bytes(), scores.read_bytes()) == original


def test_bioir_cache_binds_template_and_method_but_legacy_hash_is_unchanged(
    pipeline, tmp_path
):
    legacy = _config(pipeline, tmp_path, LEGACY_TOOL)
    legacy_key = (
        f"{legacy.repo_dir}:{legacy.input_dir}:{legacy.mapping_file}:"
        f"{legacy.chain_mapping}:{legacy.workers}:{legacy.model_version}:"
        f"{legacy.batch_size}:{legacy.heterodimers}:{legacy.provider_id}"
    )
    assert legacy.get_hash() == hashlib.sha256(legacy_key.encode()).hexdigest()[:16]
    multimer = _config(pipeline, tmp_path)
    assert multimer.get_hash() != _config(pipeline, tmp_path, BIOIR_PTM_TOOL).get_hash()
    custom = tmp_path / "template.json"
    custom.write_bytes(TEMPLATE.read_bytes())
    versioned = _config(pipeline, tmp_path, modelcif_template=custom)
    before = versioned.get_hash()
    template = json.loads(custom.read_text())
    template["categories"]["_software"]["version"] = ["producer-version-fixture"]
    custom.write_text(json.dumps(template))
    assert versioned.get_hash() != before


def _assets(tmp_path, chain_count):
    model_id = "AF-TEST"
    pdb = tmp_path / f"{model_id}-model_v1.pdb"
    lines = []
    for i in range(chain_count):
        for residue in (1, 2):
            serial = 2 * i + residue
            lines.append(
                f"ATOM  {serial:5d}  CA  ALA {chr(65 + i)}{residue:4d}    "
                f"{float(residue * 4):8.3f}{float(i * 20):8.3f}{0.0:8.3f}  1.00 80.00           C  "
            )
    pdb.write_text("\n".join(lines) + "\nEND\n")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac\n"
        + "".join(f"{model_id},1,{chr(65 + i)},PTEST\n" for i in range(chain_count))
    )
    db = tmp_path / "uniprot.duckdb"
    with duckdb.connect(str(db)) as con:
        con.execute(
            "CREATE TABLE entry(primary_ac VARCHAR, sequence VARCHAR, entry_name VARCHAR, "
            "organism VARCHAR, taxid INTEGER, sequence_version_date VARCHAR)"
        )
        con.execute(
            "INSERT INTO entry VALUES ('PTEST', 'AA', 'TEST', 'Test organism', 1, '2026-01-01')"
        )
    return pdb, manifest, db


@pytest.mark.parametrize("version", ["?", "producer-version-fixture"])
@pytest.mark.parametrize("exporter", ["single", "batch"])
@pytest.mark.parametrize(
    "tool,chains",
    [(BIOIR_PTM_TOOL, 1), (BIOIR_MULTIMER_TOOL, 2), (BIOIR_MULTIMER_TOOL, 1)],
)
def test_real_export_and_cif_conversion_preserve_bioir_attribution(
    tmp_path, exporter, tool, chains, version
):
    pdb, manifest, db = _assets(tmp_path, chains)
    metadata = tmp_path / "metadata" / "AF-TEST.json"
    ids = tmp_path / "ids.txt"
    ids.write_text("AF-TEST\n")
    template_path = tmp_path / "template.json"
    template = json.loads(TEMPLATE.read_text())
    template["categories"]["_software"]["version"] = [version]
    template_path.write_text(json.dumps(template))
    script = (
        "export_modelcif_input.py"
        if exporter == "single"
        else "batch_export_modelcif_input.py"
    )
    command = [
        sys.executable,
        str(ROOT / "uniprot/scripts" / script),
        "--manifest",
        str(manifest),
        "--db",
        str(db),
        "--template",
        str(template_path),
        "--prediction-tool",
        tool,
    ]
    if exporter == "single":
        command += ["--model-id", "AF-TEST", "--out", str(metadata)]
    else:
        command += [
            "--model-ids",
            str(ids),
            "--output-dir",
            str(metadata.parent),
            "--workers",
            "1",
        ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    before = pdb.read_bytes(), metadata.read_bytes()
    output = tmp_path / "model.cif"
    generate(
        str(pdb), str(metadata), str(output), validate_dict_path="", skip_alignment=True
    )
    assert (pdb.read_bytes(), metadata.read_bytes()) == before
    block = gemmi.cif.read_file(str(output)).sole_block()
    software = block.get_mmcif_category("_software.")
    assert software["name"][0] == BIOIR_SOFTWARE_NAME
    if version == "?":
        assert software["version"][0] in (None, "?")
    else:
        assert software["version"][0] == version
    assert all(
        not name.startswith(("AlphaFold", "ColabFold", "OpenFold-TRT"))
        for name in software["name"]
    )
    protocol = block.get_mmcif_category("_ma_protocol_step.")
    assert protocol["details"][0].endswith(tool)
    groups = block.get_mmcif_category("_ma_software_group.")
    assert set(groups["software_id"]) <= set(software["pdbx_ordinal"])
    assert set(protocol["software_group_id"]) <= set(groups["group_id"])
    atoms = block.get_mmcif_category("_atom_site.")
    assert [float(x) for x in atoms["Cartn_x"]] == [4.0, 8.0] * chains
    assert [float(x) for x in atoms["Cartn_y"]] == [
        float(i * 20) for i in range(chains) for _ in range(2)
    ]
    assert [float(x) for x in atoms["B_iso_or_equiv"]] == [80.0] * (2 * chains)


@pytest.mark.parametrize("resume", [False, True])
def test_dry_run_and_resume_do_not_bypass_producer_validation(
    pipeline, tmp_path, resume
):
    pdb, manifest, db = _assets(tmp_path, 1)
    (tmp_path / "AF-TEST-meta_v1.json").write_text(
        json.dumps({"bioir_model_source": "openfold2_ptm_1"})
    )
    ids = tmp_path / "ids.txt"
    ids.write_text("AF-TEST\n")
    provider = tmp_path / "provider.json"
    provider.write_text("{}")
    command = [
        sys.executable,
        str(ROOT / "scripts/production_pipeline.py"),
        "--input-dir",
        str(tmp_path),
        "--output-dir",
        str(tmp_path / "output"),
        "--mapping-file",
        str(ids),
        "--chain-mapping",
        str(manifest),
        "--uniprot-db",
        str(db),
        "--provider-json",
        str(provider),
        "--tool-used",
        BIOIR_MULTIMER_TOOL,
        "--homodimer-tool-used",
        BIOIR_MULTIMER_TOOL,
        "--dry-run",
        "--modelcif-template",
        str(TEMPLATE),
    ]
    if resume:
        command.append("--resume")
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert "bioir_model_source" in result.stderr
    assert not (tmp_path / "output/staging").exists()


@pytest.mark.parametrize("tool", [BIOIR_PTM_TOOL, BIOIR_MULTIMER_TOOL])
def test_valid_cli_dry_run_and_chain_metadata_keep_exact_tool(pipeline, tmp_path, tool):
    pdb, manifest, db = _assets(tmp_path, 2 if tool == BIOIR_MULTIMER_TOOL else 1)
    source = "openfold2_ptm_1" if tool == BIOIR_PTM_TOOL else "alphafold2_multimer_1"
    (tmp_path / "AF-TEST-meta_v1.json").write_text(
        json.dumps({"bioir_model_source": source})
    )
    ids = tmp_path / "ids.txt"
    ids.write_text("AF-TEST\n")
    provider = tmp_path / "provider.json"
    provider.write_text("{}")
    output = tmp_path / "output"
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/production_pipeline.py"),
            "--input-dir",
            str(tmp_path),
            "--output-dir",
            str(output),
            "--mapping-file",
            str(ids),
            "--chain-mapping",
            str(manifest),
            "--uniprot-db",
            str(db),
            "--provider-json",
            str(provider),
            "--tool-used",
            tool,
            "--homodimer-tool-used",
            tool,
            "--dry-run",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    config = json.loads((output / "config/dataset_config.json").read_text())
    assert config["toolUsed"] == config["homodimerToolUsed"] == tool
    exporter = _module("uniprot/scripts/batch_export_metadata.py")
    rows = [
        exporter.ManifestRow(
            model_entity_id="AF-TEST",
            entity_id="1",
            chain_id=chain,
            uniprot_ac="PTEST",
            sequence_start=1,
            sequence_end=2,
            is_fragment=False,
            is_isoform=False,
            entity_type="protein",
            average_plddt=80.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.0,
            fraction_plddt_confident=1.0,
            fraction_plddt_very_high=0.0,
        )
        for chain in ("AB" if tool == BIOIR_MULTIMER_TOOL else "A")
    ]
    records = exporter.build_chain_records(
        "AF-TEST", config, rows, {"PTEST": {"sequence": "AA", "entry_name": "TEST"}}, {}
    )
    assert {record["toolUsed"] for record in records} == {tool}


def test_pipeline_stage09_reaches_actual_bioir_exporter(pipeline, tmp_path):
    pdb, manifest, db = _assets(tmp_path, 1)
    config = _config(
        pipeline,
        tmp_path,
        BIOIR_PTM_TOOL,
        uniprot_db=db,
        workers=1,
        python_cmd=[sys.executable],
    )
    config.output_dir.mkdir()
    (config.output_dir / "model_ids.txt").write_text("AF-TEST\n")
    merged = config.output_dir / "merged_manifests"
    merged.mkdir()
    (merged / "uniprot_afid_mapping.csv").write_bytes(manifest.read_bytes())
    result = pipeline.stage_09_export_modelcif_input(
        ["AF-TEST"],
        config,
        pipeline.PipelineLogger(config.output_dir),
        pipeline.ErrorTracker(),
    )
    assert result["success"], result
    metadata = json.loads(
        (config.output_dir / "modelcif_input/AF-TEST.json").read_text()
    )
    assert metadata["metadata"]["prediction_tool"] == BIOIR_PTM_TOOL
    assert metadata["categories"]["_software"]["name"][0] == BIOIR_SOFTWARE_NAME
