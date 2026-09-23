from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pytest

from uniprot.scripts import add_fragment_metadata


@pytest.fixture
def db_path(tmp_path: Path) -> Path:
    path = tmp_path / "uniprot.duckdb"
    con = duckdb.connect(str(path))
    try:
        con.execute(
            "CREATE TABLE entry(primary_ac VARCHAR, sequence VARCHAR)"
        )
        con.execute(
            "INSERT INTO entry VALUES ('P11111', 'ACDE'), ('Q22222', 'ABCDE')"
        )
    finally:
        con.close()
    return path


def _fragments(
    tmp_path: Path, value: object, name: str = "fragments.json"
) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _rows(db_path: Path) -> list[tuple[str, int, int, str]]:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return con.execute(
            """
            SELECT uniprot_ac, sequence_start, sequence_end, protein_name
            FROM fragment_metadata
            ORDER BY uniprot_ac, sequence_start, sequence_end
            """
        ).fetchall()
    finally:
        con.close()


def _table_exists(db_path: Path, table: str) -> bool:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return bool(con.execute(
            "SELECT 1 FROM information_schema.tables WHERE table_name = ?",
            [table],
        ).fetchall())
    finally:
        con.close()


def test_insert_update_and_history(db_path: Path, tmp_path: Path) -> None:
    first = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 4,
            "protein_name": "Full example",
        }]},
    )
    applied = add_fragment_metadata.apply_fragments(db_path, first)

    assert applied["counts"]["updated"] == 1
    assert _rows(db_path) == [("P11111", 1, 4, "Full example")]
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        history = con.execute(
            """
            SELECT old_value_json, new_value_json
            FROM fragment_metadata_history WHERE run_id = ?
            """,
            [applied["run_id"]],
        ).fetchone()
    finally:
        con.close()
    assert tuple(json.loads(value) for value in history) == (
        None, "Full example"
    )

    second = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 4,
            "protein_name": "Renamed example",
        }]},
        "second.json",
    )
    report = add_fragment_metadata.apply_fragments(db_path, second)

    assert report["counts"]["updated"] == 1
    assert report["applied_changes"][0]["old_value"] == "Full example"
    assert _rows(db_path) == [("P11111", 1, 4, "Renamed example")]


def test_invalid_records_and_unknown_accessions_do_not_block_valid_ones(
    db_path: Path, tmp_path: Path
) -> None:
    fragments = _fragments(
        tmp_path,
        {
            "P11111": [
                {
                    "sequence_start": 0, "sequence_end": 2,
                    "protein_name": "Bad",
                },
                {
                    "sequence_start": 3, "sequence_end": 2,
                    "protein_name": "Bad",
                },
                {
                    "sequence_start": 2, "sequence_end": 5,
                    "protein_name": "Bad",
                },
                {"sequence_start": 2, "sequence_end": 3, "protein_name": "  "},
                {
                    "sequence_start": 2, "sequence_end": 3,
                    "protein_name": "Valid",
                },
            ],
            "UNKNOWN": [
                {
                    "sequence_start": 1, "sequence_end": 1,
                    "protein_name": "Unknown",
                }
            ],
        },
    )

    report = add_fragment_metadata.apply_fragments(db_path, fragments)

    assert _rows(db_path) == [("P11111", 2, 3, "Valid")]
    assert report["counts"] == {
        "requested": 6, "found": 5, "updated": 1, "unchanged": 0,
        "skipped": 5, "conflicted": 0,
    }
    assert {item["reason_code"] for item in report["skipped_items"]} == {
        "invalid_fragment", "invalid_range", "unknown_accession",
    }


def test_dry_run_creates_no_fragment_or_history_tables(
    db_path: Path, tmp_path: Path
) -> None:
    fragments = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 1, "protein_name": "Dry run",
        }]},
    )

    report = add_fragment_metadata.apply_fragments(
        db_path, fragments, dry_run=True
    )

    assert report["dry_run"] is True
    assert report["run_id"] is None
    assert report["applied_changes"][0]["status"] == "would_update"
    assert not _table_exists(db_path, add_fragment_metadata.FRAGMENT_TABLE)
    assert not _table_exists(db_path, add_fragment_metadata.HISTORY_TABLE)


def test_restore_removes_inserted_fragment_and_restores_update(
    db_path: Path, tmp_path: Path
) -> None:
    initial = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 2, "protein_name": "Original",
        }]},
    )
    add_fragment_metadata.apply_fragments(db_path, initial)
    change = _fragments(
        tmp_path,
        {
            "P11111": [{
                "sequence_start": 1, "sequence_end": 2,
                "protein_name": "Changed",
            }],
            "Q22222": [{
                "sequence_start": 3, "sequence_end": 5,
                "protein_name": "Inserted",
            }],
        },
        "change.json",
    )
    applied = add_fragment_metadata.apply_fragments(db_path, change)

    restored = add_fragment_metadata.restore_run(db_path, applied["run_id"])

    assert restored["counts"]["updated"] == 2
    assert _rows(db_path) == [("P11111", 1, 2, "Original")]


def test_restore_conflict_preserves_later_change(
    db_path: Path, tmp_path: Path
) -> None:
    first = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 2, "protein_name": "First",
        }]},
        "first.json",
    )
    second = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 2, "protein_name": "Second",
        }]},
        "second.json",
    )
    first_report = add_fragment_metadata.apply_fragments(db_path, first)
    add_fragment_metadata.apply_fragments(db_path, second)

    restored = add_fragment_metadata.restore_run(
        db_path, first_report["run_id"]
    )

    assert restored["counts"]["updated"] == 0
    assert restored["counts"]["conflicted"] == 1
    assert restored["conflicts"][0]["reason_code"] == "current_value_conflict"
    assert _rows(db_path) == [("P11111", 1, 2, "Second")]


def test_cli_writes_report_and_prints_summary(
    db_path: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fragments = _fragments(
        tmp_path,
        {"P11111": [{
            "sequence_start": 1, "sequence_end": 1, "protein_name": "CLI",
        }]},
    )
    report_path = tmp_path / "report.json"

    exit_code = add_fragment_metadata.main([
        "--db", str(db_path), "--fragments", str(fragments),
        "--report", str(report_path),
    ])

    assert exit_code == 0
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["input_sha256"]
    assert report["run_id"]
    assert report["counts"]["updated"] == 1
    assert "Run ID:" in capsys.readouterr().out
