from __future__ import annotations

import csv
from pathlib import Path

from scripts.generate_colabfold_e2e_example import (
    DEFAULT_FIXTURE_MANIFEST,
    DEFAULT_FIXTURES_ROOT,
    DEFAULT_EXAMPLE_IDS,
    chain_manifest_rows,
    load_fixture_models,
    merge_csv_files,
)


def test_load_fixture_models_defaults_are_present() -> None:
    models = load_fixture_models(
        DEFAULT_FIXTURE_MANIFEST,
        DEFAULT_FIXTURES_ROOT,
        DEFAULT_EXAMPLE_IDS,
    )

    assert [model.example_id for model in models] == DEFAULT_EXAMPLE_IDS
    assert all(model.category == "monomer" for model in models)
    assert all(model.directory.exists() for model in models)


def test_chain_manifest_rows_use_shared_entity_id_for_same_accession() -> None:
    models = load_fixture_models(
        DEFAULT_FIXTURE_MANIFEST,
        DEFAULT_FIXTURES_ROOT,
        ["AF-0000000065760001"],
    )

    rows = chain_manifest_rows(models)

    assert len(rows) == 2
    assert rows[0]["entity_id"] == "1"
    assert rows[1]["entity_id"] == "1"
    assert rows[0]["uniprot_ac"] == rows[1]["uniprot_ac"] == "Q6GZX4"


def test_merge_csv_files_concatenates_rows_once(tmp_path: Path) -> None:
    input_one = tmp_path / "a.csv"
    input_two = tmp_path / "b.csv"
    output = tmp_path / "merged.csv"
    fieldnames = ["model_entity_id", "chain_id"]

    for path, rows in [
        (input_one, [{"model_entity_id": "AF-1", "chain_id": "A"}]),
        (input_two, [{"model_entity_id": "AF-2", "chain_id": "B"}]),
    ]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    merge_csv_files([input_one, input_two], output)

    with output.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        merged_rows = list(reader)

    assert reader.fieldnames == fieldnames
    assert merged_rows == [
        {"model_entity_id": "AF-1", "chain_id": "A"},
        {"model_entity_id": "AF-2", "chain_id": "B"},
    ]
