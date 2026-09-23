"""Focused checks for the fragment-metadata pre-structure E2E fixture."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

from uniprot.scripts.enrich_fragment_manifest import enrich_manifest


REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = REPO_ROOT / "examples" / "fragment_metadata_e2e"


def _runner_module():
    module_path = (
        REPO_ROOT / "scripts" / "generate_fragment_metadata_e2e_example.py"
    )
    spec = importlib.util.spec_from_file_location(
        "fragment_metadata_e2e_runner", module_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_committed_manifest_covers_requested_topologies() -> None:
    rows = _rows(EXAMPLE_DIR / "config" / "canonical_input_manifest.csv")
    by_model: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_model.setdefault(row["model_entity_id"], []).append(row)

    assert {model: len(chains) for model, chains in by_model.items()} == {
        "AF-0000000211971324": 1,
        "AF-0000000212005744": 1,
        "AF-0000000212013504": 2,
        "AF-0000000212013519": 2,
        "AF-0000000212039399": 2,
        "AF-0000000212039401": 2,
        "AF-0000000212039400": 2,
    }
    mixed = by_model["AF-0000000212039400"]
    assert [(row["uniprot_ac"], row["is_fragment"]) for row in mixed] == [
        ("P27409", "true"), ("P28711", "false")
    ]
    source_ids = {
        row["afdb_id"]
        for row in _rows(
            EXAMPLE_DIR / "config" / "source_collaborator_wide.csv"
        )
    }
    assert set(by_model) == source_ids


def test_strict_canonical_enrichment_assigns_fragment_names_and_entities(
    tmp_path: Path,
) -> None:
    module = _runner_module()
    assert module.main(
        [
            "--example-dir", str(EXAMPLE_DIR), "--output-dir",
            str(tmp_path / "generated"), "--python-bin", sys.executable,
        ]
    ) == 0

    # The runner exercises the command-line loader; use its generated database
    # and invoke strict enrichment directly to test the canonical adapter too.
    output = tmp_path / "generated"
    report = enrich_manifest(
        EXAMPLE_DIR / "config" / "canonical_input_manifest.csv",
        tmp_path / "strict.csv", output / "mock_uniprot.duckdb", strict=True,
    )
    assert report["status"] == "success"
    rows = _rows(tmp_path / "strict.csv")
    fragment_rows = [row for row in rows if row["is_fragment"] == "true"]
    assert {row["protein_name"] for row in fragment_rows} == {
        "Development fragment alpha", "Development fragment omega"
    }
    fragment_heterodimer = [
        row
        for row in rows
        if row["model_entity_id"] == "AF-0000000212039401"
    ]
    assert [row["entity_id"] for row in fragment_heterodimer] == ["1", "2"]


def test_prestructure_runner_exports_consumer_inputs(tmp_path: Path) -> None:
    module = _runner_module()
    output = module.run_example(
        argparse.Namespace(
            example_dir=EXAMPLE_DIR,
            output_dir=tmp_path / "generated",
            python_bin=sys.executable,
            template=(
                REPO_ROOT / "uniprot" / "templates"
                / "colabfold_example_modelcif_metadata.json"
            ),
        )
    )

    summary = json.loads(
        (output / "run_summary.json").read_text(encoding="utf-8")
    )
    assert summary["status"] == "success"
    assert summary["model_count"] == 7
    assert summary["chain_row_count"] == 12
    assert (output / "manifests" / "real_cases_from_wide.csv").is_file()
    assert len(list((output / "model_metadata").glob("*.json"))) == 7
    assert len(list((output / "chain_metadata").glob("*.json"))) == 7
    assert len(list((output / "modelcif_input").glob("*.json"))) == 7
    model_batch = json.loads(
        (output / "model_batches" / "models.json").read_text()
    )
    chain_batch = json.loads(
        (output / "chain_batches" / "chains.json").read_text()
    )
    assert len(model_batch) == 7
    assert len(chain_batch) == 12
