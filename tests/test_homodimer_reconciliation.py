"""Frozen AFCDB homodimer reconciliation contracts.

Fixture provenance: the raw PDB/meta pair, the two selected manifest rows,
one-line model list, example provider/dataset configs, and example UniProt
DuckDB come from PDBe AFDB-Integration-Kit commit 6e0f00f (as inspected from
the public fixture checkout at dcc5e034). Generated CIF/BCIF and output trees
are deliberately excluded.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures" / "homodimer_reconciliation"


def _load_production_pipeline():
    module_path = REPO_ROOT / "scripts" / "production_pipeline.py"
    spec = importlib.util.spec_from_file_location(
        "production_pipeline_homodimer_contract", module_path
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


generate_dataset_config = _load_production_pipeline().generate_dataset_config


def test_generate_dataset_config_preserves_distinct_predictor_names(
    tmp_path: Path,
) -> None:
    path = tmp_path / "dataset_config.json"

    generate_dataset_config(
        path,
        "NVDA",
        tmp_path,
        "OpenFold-TRT / AlphaFold-Multimer",
        homodimer_tool_used="ColabFold v1.6.0 / AlphaFold-Multimer",
    )

    config = json.loads(path.read_text(encoding="utf-8"))
    assert config["toolUsed"] == "OpenFold-TRT / AlphaFold-Multimer"
    assert config["homodimerToolUsed"] == "ColabFold v1.6.0 / AlphaFold-Multimer"


def test_homodimer_chain_manifest_has_two_structure_bound_chains() -> None:
    manifest_path = FIXTURE_ROOT / "config" / "chain_manifest.csv"

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert [row["chain_id"] for row in rows] == ["A", "B"]
    assert {row["uniprot_ac"] for row in rows} == {"Q46806"}
    assert {row["model_entity_id"] for row in rows} == {
        "AF-0000000066074510"
    }
