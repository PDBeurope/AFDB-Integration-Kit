import json
from pathlib import Path

import pytest

from afdb_integration_kit.manifest.resolver import (
    _deduplicate_accessions,
    build_colabfold_manifest,
    classify_model_ids,
)


def _create_duckdb(path: Path, rows: list[tuple[str, int]]) -> None:
    duckdb = pytest.importorskip("duckdb")
    con = duckdb.connect(str(path))
    try:
        con.execute(
            """
            CREATE TABLE entry (
                primary_ac VARCHAR,
                sequence VARCHAR
            )
            """
        )
        con.executemany(
            "INSERT INTO entry VALUES (?, ?)",
            [(accession, "A" * length) for accession, length in rows],
        )
    finally:
        con.close()


def test_classify_model_ids_accepts_supported_16_digit_formats() -> None:
    model_ids = [
        "AF-0000000000000001",
        "AF_0000000000000002",
        "AF_0000000000000003_AF_0000000000000004",
        "AF-1234",
        "AF_0000000000000005_AF-0000000000000006",
    ]

    classified, unrecognised = classify_model_ids(model_ids)

    assert classified == {
        "AF-0000000000000001": ["AF-0000000000000001"],
        "AF_0000000000000002": ["AF_0000000000000002"],
        "AF_0000000000000003_AF_0000000000000004": [
            "AF_0000000000000003",
            "AF_0000000000000004",
        ],
    }
    assert unrecognised == ["AF-1234", "AF_0000000000000005_AF-0000000000000006"]


def test_build_colabfold_manifest_hyphenates_output_rows() -> None:
    model_ids = [
        "AF_0000000000000001",
        "AF_0000000000000002_AF_0000000000000003",
    ]
    classified, unrecognised = classify_model_ids(model_ids)

    rows, skipped = build_colabfold_manifest(
        model_ids,
        af_to_uniprot={
            "AF_0000000000000001": "P11111",
            "AF_0000000000000002": "Q22222",
            "AF_0000000000000003": "R33333",
        },
        classified=classified,
    )

    assert unrecognised == []
    assert skipped == []
    assert rows == [
        {
            "model_entity_id": "AF-0000000000000001",
            "chain_id": "A",
            "uniprot_ac": "P11111",
        },
        {
            "model_entity_id": "AF-0000000000000001",
            "chain_id": "B",
            "uniprot_ac": "P11111",
        },
        {
            "model_entity_id": "AF-0000000000000002_AF-0000000000000003",
            "chain_id": "A",
            "uniprot_ac": "Q22222",
        },
        {
            "model_entity_id": "AF-0000000000000002_AF-0000000000000003",
            "chain_id": "B",
            "uniprot_ac": "R33333",
        },
    ]


def test_ambiguous_accessions_fail_without_disambiguation_inputs() -> None:
    clean, failed = _deduplicate_accessions(
        {"AF_0000000000000001": {"P11111", "Q22222"}},
    )

    assert clean == {}
    assert failed == ["AF_0000000000000001"]


def test_ambiguous_accessions_can_be_resolved_by_plddt_length(tmp_path: Path) -> None:
    duckdb_path = tmp_path / "entries.duckdb"
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    _create_duckdb(duckdb_path, [("P11111", 3), ("Q22222", 2)])
    (input_dir / "AF-0000000000000001-meta_v1.json").write_text(
        json.dumps({"plddt": [90.0, 80.0, 70.0, 60.0, 50.0, 40.0]}),
        encoding="utf-8",
    )

    clean, failed = _deduplicate_accessions(
        {"AF_0000000000000001": {"P11111", "Q22222"}},
        uniprot_db_path=duckdb_path,
        input_dir=input_dir,
    )

    assert clean == {"AF_0000000000000001": "P11111"}
    assert failed == []
