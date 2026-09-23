from __future__ import annotations

import hashlib
import json
from pathlib import Path

import duckdb
import pytest

from uniprot.scripts import add_custom_metadata
from uniprot.scripts import export_model_metadata


ENTRY_DDL = """
CREATE TABLE entry(
    primary_ac VARCHAR,
    entry_name VARCHAR,
    reviewed BOOLEAN,
    protein_full_names VARCHAR[],
    protein_short_names VARCHAR[],
    gene_names VARCHAR,
    gene_synonyms VARCHAR[],
    organism VARCHAR,
    organism_common_names VARCHAR[],
    taxid BIGINT,
    length INTEGER,
    sequence_version_date VARCHAR,
    is_uniprot_reference_proteome BOOLEAN,
    md5 VARCHAR,
    sequence VARCHAR,
    release VARCHAR,
    is_isoform BOOLEAN,
    unsupported DOUBLE
)
"""


@pytest.fixture
def db_path(tmp_path: Path) -> Path:
    path = tmp_path / "uniprot.duckdb"
    sequence = "ACDE"
    con = duckdb.connect(str(path))
    try:
        con.execute(ENTRY_DDL)
        con.execute(
            """
            INSERT INTO entry VALUES
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?),
                (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                "P11111",
                "OLD_ENTRY",
                True,
                ["Old protein"],
                ["Old"],
                None,
                ["OLD1"],
                "Example organism",
                [],
                9606,
                len(sequence),
                "2025-01-01",
                True,
                hashlib.md5(sequence.encode("ascii")).hexdigest(),
                sequence,
                "2025_01",
                False,
                1.5,
                "Q22222",
                "SECOND_ENTRY",
                False,
                ["Second protein"],
                [],
                "SECOND",
                [],
                "Other organism",
                ["other"],
                10090,
                len(sequence),
                "2025-01-02",
                False,
                hashlib.md5(sequence.encode("ascii")).hexdigest(),
                sequence,
                "2025_01",
                False,
                2.5,
            ],
        )
    finally:
        con.close()
    return path


def _annotations(
    tmp_path: Path, value: object, name: str = "annotations.json"
) -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _entry_value(db_path: Path, column: str, accession: str = "P11111"):
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return con.execute(
            f'SELECT "{column}" FROM entry WHERE primary_ac = ?',
            [accession],
        ).fetchone()[0]
    finally:
        con.close()


def _table_exists(db_path: Path, table: str) -> bool:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return bool(
            con.execute(
                "SELECT 1 FROM information_schema.tables WHERE table_name = ?",
                [table],
            ).fetchall()
        )
    finally:
        con.close()


def test_scalar_list_null_and_non_blocking_skips(
    db_path: Path, tmp_path: Path
) -> None:
    annotations = _annotations(
        tmp_path,
        {
            "P11111": {
                "protein_full_names": "Custom protein",
                "protein_short_names": ["Custom", "CP"],
                "gene_names": "CUSTOM1",
                "taxid": 123,
                "reviewed": False,
                "gene_synonyms": None,
                "organism_common_names": ["human", 9606],
                "entry_name": ["not", "a", "string"],
                "unsupported": 2.0,
                "uniprotDescription": "export alias",
                "primary_ac": "X99999",
            },
            "Q22222": {"taxid": True, "reviewed": 1},
            "UNKNOWN": {"gene_names": "NOPE"},
        },
    )

    report = add_custom_metadata.apply_annotations(db_path, annotations)

    assert _entry_value(db_path, "protein_full_names") == ["Custom protein"]
    assert _entry_value(db_path, "protein_short_names") == ["Custom", "CP"]
    assert _entry_value(db_path, "gene_names") == "CUSTOM1"
    assert _entry_value(db_path, "taxid") == 123
    assert _entry_value(db_path, "reviewed") is False
    assert _entry_value(db_path, "gene_synonyms") == ["OLD1"]
    assert _entry_value(db_path, "primary_ac") == "P11111"
    assert _entry_value(db_path, "taxid", "Q22222") == 10090
    assert _entry_value(db_path, "reviewed", "Q22222") is False
    assert report["counts"]["updated"] == 5
    reason_codes = {item["reason_code"] for item in report["skipped_items"]}
    assert reason_codes == {
        "invalid_value",
        "unknown_column",
        "immutable_column",
        "unknown_accession",
    }
    unknown_columns = {
        item["column_name"]
        for item in report["skipped_items"]
        if item["reason_code"] == "unknown_column"
    }
    assert "uniprotDescription" in unknown_columns
    assert "unsupported" not in unknown_columns
    unsupported = [
        item
        for item in report["skipped_items"]
        if item.get("duckdb_type") == "DOUBLE"
    ]
    assert len(unsupported) == 1


@pytest.mark.parametrize(
    ("bundle", "reason_fragment"),
    [
        ({"sequence": "ACD"}, "supplied together"),
        (
            {
                "sequence": "ACD",
                "length": 4,
                "md5": hashlib.md5(b"ACD").hexdigest(),
            },
            "does not match sequence length",
        ),
        (
            {"sequence": "ACD", "length": 3, "md5": "0" * 32},
            "md5 does not match",
        ),
        (
            {
                "sequence": "ACd",
                "length": 3,
                "md5": hashlib.md5(b"ACd").hexdigest(),
            },
            "uppercase recognised",
        ),
    ],
)
def test_invalid_sequence_bundle_is_one_skip_and_independent_fields_apply(
    db_path: Path,
    tmp_path: Path,
    bundle: dict[str, object],
    reason_fragment: str,
) -> None:
    bundle["gene_names"] = "INDEPENDENT"
    annotations = _annotations(tmp_path, {"P11111": bundle})

    report = add_custom_metadata.apply_annotations(db_path, annotations)

    assert _entry_value(db_path, "sequence") == "ACDE"
    assert _entry_value(db_path, "gene_names") == "INDEPENDENT"
    sequence_skips = [
        item
        for item in report["skipped_items"]
        if item["reason_code"] == "invalid_sequence_bundle"
    ]
    assert len(sequence_skips) == 1
    assert reason_fragment in sequence_skips[0]["reason"]


def test_valid_sequence_bundle_updates_together(
    db_path: Path, tmp_path: Path
) -> None:
    sequence = "BJOUXZ"
    digest = hashlib.md5(sequence.encode("ascii")).hexdigest()
    annotations = _annotations(
        tmp_path,
        {
            "P11111": {
                "sequence": sequence,
                "length": len(sequence),
                "md5": digest.upper(),
            }
        },
    )

    report = add_custom_metadata.apply_annotations(db_path, annotations)

    assert report["counts"]["updated"] == 3
    assert _entry_value(db_path, "sequence") == sequence
    assert _entry_value(db_path, "length") == len(sequence)
    assert _entry_value(db_path, "md5") == digest


def test_dry_run_changes_neither_entry_nor_history(
    db_path: Path, tmp_path: Path
) -> None:
    annotations = _annotations(
        tmp_path,
        {"P11111": {"gene_names": "DRY", "protein_full_names": "Dry"}},
    )

    report = add_custom_metadata.apply_annotations(
        db_path, annotations, dry_run=True
    )

    assert report["dry_run"] is True
    assert report["run_id"] is None
    assert report["counts"]["updated"] == 2
    assert _entry_value(db_path, "gene_names") is None
    assert _entry_value(db_path, "protein_full_names") == ["Old protein"]
    assert not _table_exists(db_path, add_custom_metadata.HISTORY_TABLE)


def test_history_records_json_values_and_restore_handles_null(
    db_path: Path, tmp_path: Path
) -> None:
    annotations = _annotations(
        tmp_path,
        {
            "P11111": {
                "gene_names": "RESTORE_ME",
                "protein_full_names": "Restored protein",
            }
        },
    )
    apply_report = add_custom_metadata.apply_annotations(db_path, annotations)

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        rows = con.execute(
            """
            SELECT column_name, old_value_json, new_value_json
            FROM custom_annotation_history
            WHERE run_id = ?
            ORDER BY column_name
            """,
            [apply_report["run_id"]],
        ).fetchall()
    finally:
        con.close()

    history = {
        name: (json.loads(old_value), json.loads(new_value))
        for name, old_value, new_value in rows
    }
    assert history["gene_names"] == (None, "RESTORE_ME")
    assert history["protein_full_names"] == (
        ["Old protein"],
        ["Restored protein"],
    )

    restore_report = add_custom_metadata.restore_run(
        db_path, apply_report["run_id"]
    )

    assert restore_report["counts"]["updated"] == 2
    assert restore_report["counts"]["conflicted"] == 0
    assert _entry_value(db_path, "gene_names") is None
    assert _entry_value(db_path, "protein_full_names") == ["Old protein"]


def test_restore_conflict_preserves_later_change(
    db_path: Path, tmp_path: Path
) -> None:
    first = _annotations(
        tmp_path, {"P11111": {"gene_names": "FIRST"}}, "first.json"
    )
    second = _annotations(
        tmp_path, {"P11111": {"gene_names": "SECOND"}}, "second.json"
    )
    first_report = add_custom_metadata.apply_annotations(db_path, first)
    add_custom_metadata.apply_annotations(db_path, second)

    restore_report = add_custom_metadata.restore_run(
        db_path, first_report["run_id"]
    )

    assert restore_report["counts"]["updated"] == 0
    assert restore_report["counts"]["conflicted"] == 1
    assert (
        restore_report["conflicts"][0]["reason_code"]
        == "current_value_conflict"
    )
    assert _entry_value(db_path, "gene_names") == "SECOND"


def test_restore_dry_run_does_not_restore(
    db_path: Path, tmp_path: Path
) -> None:
    annotations = _annotations(
        tmp_path, {"P11111": {"gene_names": "STILL_CUSTOM"}}
    )
    apply_report = add_custom_metadata.apply_annotations(db_path, annotations)

    report = add_custom_metadata.restore_run(
        db_path, apply_report["run_id"], dry_run=True
    )

    assert report["counts"]["updated"] == 1
    assert report["applied_changes"][0]["status"] == "would_restore"
    assert _entry_value(db_path, "gene_names") == "STILL_CUSTOM"


def test_database_write_failure_rolls_back_history_and_entry(
    db_path: Path, tmp_path: Path
) -> None:
    con = duckdb.connect(str(db_path))
    try:
        con.execute(
            """
            CREATE TABLE custom_annotation_history(
                run_id VARCHAR,
                applied_at TIMESTAMP,
                source_sha256 VARCHAR,
                primary_ac VARCHAR,
                column_name VARCHAR,
                old_value_json JSON,
                new_value_json JSON,
                required_extra VARCHAR NOT NULL
            )
            """
        )
    finally:
        con.close()
    annotations = _annotations(
        tmp_path, {"P11111": {"gene_names": "MUST_ROLL_BACK"}}
    )

    with pytest.raises(add_custom_metadata.CustomMetadataError):
        add_custom_metadata.apply_annotations(db_path, annotations)

    assert _entry_value(db_path, "gene_names") is None
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        assert con.execute(
            "SELECT count(*) FROM custom_annotation_history"
        ).fetchone()[0] == 0
    finally:
        con.close()


def test_existing_exporter_observes_custom_protein_name(
    db_path: Path, tmp_path: Path
) -> None:
    annotations = _annotations(
        tmp_path, {"P11111": {"protein_full_names": "Exporter-visible name"}}
    )
    add_custom_metadata.apply_annotations(db_path, annotations)

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        entry = export_model_metadata.fetch_entries(con, ["P11111"])["P11111"]
    finally:
        con.close()

    assert (
        export_model_metadata.derive_description(entry)
        == "Exporter-visible name"
    )


def test_report_file_and_human_cli_output(
    db_path: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    annotations = _annotations(
        tmp_path, {"P11111": {"gene_names": "CLI"}, "UNKNOWN": {"taxid": 1}}
    )
    report_path = tmp_path / "report.json"

    exit_code = add_custom_metadata.main(
        [
            "--db",
            str(db_path),
            "--annotations",
            str(annotations),
            "--report",
            str(report_path),
        ]
    )

    assert exit_code == 0
    machine_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert machine_report["database"] == str(db_path.resolve())
    assert machine_report["input_sha256"]
    assert machine_report["run_id"]
    assert machine_report["counts"]["updated"] == 1
    assert (
        machine_report["skipped_items"][0]["reason_code"]
        == "unknown_accession"
    )
    output = capsys.readouterr().out
    assert "Run ID:" in output
    assert "WARNING skipped UNKNOWN" in output


def test_malformed_json_absent_entry_and_incompatible_schema_fail_cleanly(
    db_path: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    malformed = tmp_path / "bad.json"
    malformed.write_text('{"P11111":', encoding="utf-8")
    assert add_custom_metadata.main(
        ["--db", str(db_path), "--annotations", str(malformed)]
    ) == 1
    assert "not valid JSON" in capsys.readouterr().out

    absent_entry = tmp_path / "absent.duckdb"
    con = duckdb.connect(str(absent_entry))
    con.close()
    valid = _annotations(tmp_path, {"P11111": {"gene_names": "X"}})
    assert add_custom_metadata.main(
        ["--db", str(absent_entry), "--annotations", str(valid)]
    ) == 1
    assert "Required table 'entry'" in capsys.readouterr().out

    incompatible = tmp_path / "incompatible.duckdb"
    con = duckdb.connect(str(incompatible))
    try:
        con.execute(
            "CREATE TABLE entry(primary_ac INTEGER, gene_names VARCHAR)"
        )
    finally:
        con.close()
    assert add_custom_metadata.main(
        ["--db", str(incompatible), "--annotations", str(valid)]
    ) == 1
    assert (
        "'primary_ac' must have DuckDB type VARCHAR"
        in capsys.readouterr().out
    )
