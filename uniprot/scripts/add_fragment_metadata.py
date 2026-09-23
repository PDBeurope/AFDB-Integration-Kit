#!/usr/bin/env python3
"""Load and safely restore named UniProt sequence fragments in DuckDB."""

from __future__ import annotations

import argparse
import hashlib
import json
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import duckdb


FRAGMENT_TABLE = "fragment_metadata"
HISTORY_TABLE = "fragment_metadata_history"
REQUIRED_FRAGMENT_COLUMNS = {
    "uniprot_ac": "VARCHAR",
    "sequence_start": "INTEGER",
    "sequence_end": "INTEGER",
    "protein_name": "VARCHAR",
}
REQUIRED_HISTORY_COLUMNS = {
    "run_id": "VARCHAR",
    "applied_at": "TIMESTAMP",
    "source_sha256": "VARCHAR",
    "uniprot_ac": "VARCHAR",
    "sequence_start": "INTEGER",
    "sequence_end": "INTEGER",
    "old_value_json": "JSON",
    "new_value_json": "JSON",
}


class FragmentMetadataError(RuntimeError):
    """A fatal, user-facing fragment metadata command error."""


class InvalidFragment(ValueError):
    """A non-fatal fragment record validation error."""


@dataclass(frozen=True)
class PlannedChange:
    uniprot_ac: str
    sequence_start: int
    sequence_end: int
    old_value: str | None
    new_value: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load named UniProt sequence fragments into DuckDB, or restore "
            "a previously applied run."
        )
    )
    parser.add_argument(
        "--db", required=True, type=Path,
        help="Existing DuckDB database containing the entry table.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--fragments", type=Path,
        help=(
            "JSON object keyed by UniProt accession containing fragment lists."
        ),
    )
    mode.add_argument(
        "--restore", metavar="RUN_ID",
        help=(
            "Restore all non-conflicting fragment changes recorded for RUN_ID."
        ),
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate and report changes without writing to the database.",
    )
    parser.add_argument(
        "--report", type=Path,
        help="Write the machine-readable report as JSON.",
    )
    return parser.parse_args(argv)


def _resolved(path: Path) -> Path:
    return path.expanduser().resolve()


def _quote_identifier(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def _json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _json_equal(left: Any, right: Any) -> bool:
    return _json_text(left) == _json_text(right)


def load_fragments(path: Path) -> tuple[dict[str, Any], str]:
    """Load the fragment JSON document and return its raw SHA-256."""
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise FragmentMetadataError(
            f"Cannot read fragments file {path}: {exc}"
        ) from exc
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FragmentMetadataError(
            f"Fragments file {path} is not valid JSON: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise FragmentMetadataError(
            "The fragments JSON root must be an object keyed by UniProt "
            "accession."
        )
    return value, hashlib.sha256(raw).hexdigest()


def _describe_table(
    con: duckdb.DuckDBPyConnection, table_name: str
) -> dict[str, str]:
    try:
        rows = con.execute(
            f"DESCRIBE {_quote_identifier(table_name)}"
        ).fetchall()
    except duckdb.Error as exc:
        raise FragmentMetadataError(
            f"Required table {table_name!r} is absent or unreadable: {exc}"
        ) from exc
    return {str(row[0]): str(row[1]).upper() for row in rows}


def inspect_entry_schema(con: duckdb.DuckDBPyConnection) -> None:
    """Validate only the immutable entry fields needed for range checks."""
    columns = _describe_table(con, "entry")
    for name in ("primary_ac", "sequence"):
        if name not in columns:
            raise FragmentMetadataError(
                "Incompatible entry schema: required lookup column "
                f"{name!r} is absent."
            )
    if columns["primary_ac"] != "VARCHAR":
        raise FragmentMetadataError(
            "Incompatible entry schema: 'primary_ac' must have DuckDB type "
            f"VARCHAR, found {columns['primary_ac']}."
        )
    if columns["sequence"] != "VARCHAR":
        raise FragmentMetadataError(
            "Incompatible entry schema: 'sequence' must have DuckDB type "
            f"VARCHAR, found {columns['sequence']}."
        )


def _table_exists(con: duckdb.DuckDBPyConnection, table_name: str) -> bool:
    return bool(con.execute(
        """
        SELECT 1 FROM information_schema.tables
        WHERE table_schema = current_schema() AND table_name = ?
        LIMIT 1
        """,
        [table_name],
    ).fetchall())


def _validate_table_schema(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    required_columns: Mapping[str, str],
) -> None:
    columns = _describe_table(con, table_name)
    for name, expected_type in required_columns.items():
        actual = columns.get(name)
        if actual is None:
            raise FragmentMetadataError(
                f"Incompatible {table_name} schema: required column "
                f"{name!r} is absent."
            )
        if actual != expected_type:
            raise FragmentMetadataError(
                f"Incompatible {table_name} schema: column {name!r} must "
                f"have type {expected_type}, found {actual}."
            )


def _validate_fragment_table(con: duckdb.DuckDBPyConnection) -> None:
    """Require the documented not-null columns and inclusive-range key."""
    _validate_table_schema(con, FRAGMENT_TABLE, REQUIRED_FRAGMENT_COLUMNS)
    rows = con.execute(
        f"DESCRIBE {_quote_identifier(FRAGMENT_TABLE)}"
    ).fetchall()
    nullable = {str(row[0]): str(row[2]).upper() for row in rows}
    for name in REQUIRED_FRAGMENT_COLUMNS:
        if nullable.get(name) != "NO":
            raise FragmentMetadataError(
                f"Incompatible {FRAGMENT_TABLE} schema: column {name!r} "
                "must be NOT NULL."
            )
    primary_keys = con.execute(
        """
        SELECT constraint_column_names
        FROM duckdb_constraints()
        WHERE table_name = ? AND constraint_type = 'PRIMARY KEY'
        """,
        [FRAGMENT_TABLE],
    ).fetchall()
    expected_key = ["uniprot_ac", "sequence_start", "sequence_end"]
    if len(primary_keys) != 1 or list(primary_keys[0][0]) != expected_key:
        raise FragmentMetadataError(
            f"Incompatible {FRAGMENT_TABLE} schema: expected primary key "
            "(uniprot_ac, sequence_start, sequence_end)."
        )


def _create_or_validate_fragment_table(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {_quote_identifier(FRAGMENT_TABLE)} (
            uniprot_ac VARCHAR NOT NULL,
            sequence_start INTEGER NOT NULL,
            sequence_end INTEGER NOT NULL,
            protein_name VARCHAR NOT NULL,
            PRIMARY KEY (uniprot_ac, sequence_start, sequence_end)
        )
        """
    )
    _validate_fragment_table(con)


def _create_or_validate_history(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {_quote_identifier(HISTORY_TABLE)} (
            run_id VARCHAR,
            applied_at TIMESTAMP,
            source_sha256 VARCHAR,
            uniprot_ac VARCHAR,
            sequence_start INTEGER,
            sequence_end INTEGER,
            old_value_json JSON,
            new_value_json JSON
        )
        """
    )
    _validate_table_schema(con, HISTORY_TABLE, REQUIRED_HISTORY_COLUMNS)


def _fetch_entry_lengths(
    con: duckdb.DuckDBPyConnection, accessions: Sequence[str]
) -> dict[str, int]:
    if not accessions:
        return {}
    placeholders = ",".join("?" for _ in accessions)
    rows = con.execute(
        f"""
        SELECT primary_ac, length(sequence)
        FROM entry
        WHERE primary_ac IN ({placeholders})
        """,
        list(accessions),
    ).fetchall()
    result: dict[str, int] = {}
    for accession, length in rows:
        accession = str(accession)
        if accession in result:
            raise FragmentMetadataError(
                "Incompatible entry data: primary_ac must uniquely identify "
                f"a row; found duplicate {accession!r}."
            )
        if length is None:
            raise FragmentMetadataError(
                f"Incompatible entry data: sequence is NULL for {accession!r}."
            )
        result[accession] = int(length)
    return result


def _existing_fragments(
    con: duckdb.DuckDBPyConnection,
) -> dict[tuple[str, int, int], str]:
    if not _table_exists(con, FRAGMENT_TABLE):
        return {}
    _validate_fragment_table(con)
    rows = con.execute(
        f"""
        SELECT uniprot_ac, sequence_start, sequence_end, protein_name
        FROM {_quote_identifier(FRAGMENT_TABLE)}
        """
    ).fetchall()
    return {
        (str(accession), int(start), int(end)): str(name)
        for accession, start, end, name in rows
    }


def _skipped_item(
    uniprot_ac: str,
    reason_code: str,
    reason: str,
    *,
    record_index: int | None = None,
    **details: Any,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "uniprot_ac": uniprot_ac,
        "reason_code": reason_code,
        "reason": reason,
    }
    if record_index is not None:
        item["record_index"] = record_index
    item.update(details)
    return item


def _normalize_record(value: Any) -> tuple[int, int, str]:
    if not isinstance(value, dict):
        raise InvalidFragment("fragment record must be a JSON object")
    required = {"sequence_start", "sequence_end", "protein_name"}
    missing = sorted(required.difference(value))
    extra = sorted(set(value).difference(required))
    if missing:
        raise InvalidFragment(
            "missing required field(s): " + ", ".join(missing)
        )
    if extra:
        raise InvalidFragment("unknown field(s): " + ", ".join(extra))
    start = value["sequence_start"]
    end = value["sequence_end"]
    name = value["protein_name"]
    if type(start) is not int or type(end) is not int:
        raise InvalidFragment(
            "sequence_start and sequence_end must be JSON integers "
            "(booleans are not integers)"
        )
    if start <= 0 or end <= 0:
        raise InvalidFragment("sequence positions must be positive")
    if start > end:
        raise InvalidFragment(
            "sequence_start must be less than or equal to sequence_end"
        )
    if not isinstance(name, str):
        raise InvalidFragment("protein_name must be a JSON string")
    normalized_name = name.strip()
    if not normalized_name:
        raise InvalidFragment("protein_name must not be empty or whitespace")
    return start, end, normalized_name


def plan_fragment_changes(
    con: duckdb.DuckDBPyConnection,
    fragments: Mapping[str, Any],
) -> tuple[list[PlannedChange], list[dict[str, Any]], dict[str, int]]:
    """Build all valid changes before any write transaction is opened."""
    lengths = _fetch_entry_lengths(con, list(fragments))
    existing = _existing_fragments(con)
    changes: list[PlannedChange] = []
    skipped: list[dict[str, Any]] = []
    unchanged = 0
    requested = 0
    found = 0
    seen_keys: set[tuple[str, int, int]] = set()

    for accession, records in fragments.items():
        if not isinstance(records, list):
            requested += 1
            skipped.append(_skipped_item(
                accession,
                "invalid_accession_value",
                "fragment value for an accession must be a JSON list",
            ))
            continue
        requested += len(records)
        sequence_length = lengths.get(accession)
        if sequence_length is None:
            for index, _record in enumerate(records):
                skipped.append(_skipped_item(
                    accession,
                    "unknown_accession",
                    "accession was not found in entry.primary_ac",
                    record_index=index,
                ))
            continue
        found += len(records)
        for index, record in enumerate(records):
            try:
                start, end, name = _normalize_record(record)
            except InvalidFragment as exc:
                skipped.append(_skipped_item(
                    accession, "invalid_fragment", str(exc), record_index=index
                ))
                continue
            if end > sequence_length:
                skipped.append(_skipped_item(
                    accession,
                    "invalid_range",
                    f"sequence_end {end} exceeds entry.sequence length "
                    f"{sequence_length}",
                    record_index=index,
                    sequence_start=start,
                    sequence_end=end,
                ))
                continue
            key = (accession, start, end)
            if key in seen_keys:
                skipped.append(_skipped_item(
                    accession,
                    "duplicate_fragment_key",
                    "the input contains more than one record for this "
                    "accession and inclusive range",
                    record_index=index,
                    sequence_start=start,
                    sequence_end=end,
                ))
                continue
            seen_keys.add(key)
            old_value = existing.get(key)
            if old_value == name:
                unchanged += 1
                continue
            changes.append(
                PlannedChange(accession, start, end, old_value, name)
            )

    return changes, skipped, {
        "requested": requested,
        "found": found,
        "updated": len(changes),
        "unchanged": unchanged,
        "skipped": len(skipped),
        "conflicted": 0,
    }


def _change_report(
    change: PlannedChange, *, action: str, dry_run: bool
) -> dict[str, Any]:
    return {
        "uniprot_ac": change.uniprot_ac,
        "sequence_start": change.sequence_start,
        "sequence_end": change.sequence_end,
        "old_value": change.old_value,
        "new_value": change.new_value,
        "action": action,
        "status": f"would_{action}" if dry_run else action,
    }


def _apply_changes(
    con: duckdb.DuckDBPyConnection,
    changes: Sequence[PlannedChange],
    *,
    run_id: str,
    source_sha256: str,
) -> None:
    applied_at = datetime.now(timezone.utc).replace(tzinfo=None)
    try:
        con.execute("BEGIN TRANSACTION")
        _create_or_validate_fragment_table(con)
        _create_or_validate_history(con)
        for change in changes:
            con.execute(
                f"""
                INSERT INTO {_quote_identifier(HISTORY_TABLE)}
                    (run_id, applied_at, source_sha256, uniprot_ac,
                     sequence_start, sequence_end, old_value_json,
                     new_value_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    run_id, applied_at, source_sha256, change.uniprot_ac,
                    change.sequence_start, change.sequence_end,
                    _json_text(change.old_value), _json_text(change.new_value),
                ],
            )
        for change in changes:
            if change.old_value is None:
                updated = con.execute(
                    f"""
                    INSERT INTO {_quote_identifier(FRAGMENT_TABLE)}
                        (uniprot_ac, sequence_start, sequence_end,
                         protein_name)
                    SELECT ?, ?, ?, ?
                    WHERE NOT EXISTS (
                        SELECT 1 FROM {_quote_identifier(FRAGMENT_TABLE)}
                        WHERE uniprot_ac = ? AND sequence_start = ?
                          AND sequence_end = ?
                    )
                    RETURNING uniprot_ac
                    """,
                    [
                        change.uniprot_ac, change.sequence_start,
                        change.sequence_end, change.new_value,
                        change.uniprot_ac, change.sequence_start,
                        change.sequence_end,
                    ],
                ).fetchall()
            else:
                updated = con.execute(
                    f"""
                    UPDATE {_quote_identifier(FRAGMENT_TABLE)}
                    SET protein_name = ?
                    WHERE uniprot_ac = ? AND sequence_start = ?
                      AND sequence_end = ?
                      AND protein_name IS NOT DISTINCT FROM ?
                    RETURNING uniprot_ac
                    """,
                    [
                        change.new_value, change.uniprot_ac,
                        change.sequence_start, change.sequence_end,
                        change.old_value,
                    ],
                ).fetchall()
            if len(updated) != 1:
                raise FragmentMetadataError(
                    "Fragment metadata changed after validation for "
                    f"{change.uniprot_ac}:{change.sequence_start}-"
                    f"{change.sequence_end}; no changes were applied."
                )
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except duckdb.Error:
            pass
        raise


def apply_fragments(
    db_path: Path, fragments_path: Path, *, dry_run: bool = False
) -> dict[str, Any]:
    """Validate and optionally apply one fragment metadata document."""
    db_path = _resolved(db_path)
    fragments_path = _resolved(fragments_path)
    if not db_path.is_file():
        raise FragmentMetadataError(
            f"DuckDB database does not exist: {db_path}"
        )
    fragments, source_sha256 = load_fragments(fragments_path)
    try:
        con = duckdb.connect(str(db_path), read_only=dry_run)
    except duckdb.Error as exc:
        raise FragmentMetadataError(
            f"Cannot open DuckDB database {db_path}: {exc}"
        ) from exc
    try:
        inspect_entry_schema(con)
        changes, skipped, counts = plan_fragment_changes(con, fragments)
        run_id = None if dry_run else str(uuid.uuid4())
        if not dry_run:
            _apply_changes(
                con, changes, run_id=run_id, source_sha256=source_sha256
            )
    except duckdb.Error as exc:
        raise FragmentMetadataError(f"DuckDB operation failed: {exc}") from exc
    finally:
        con.close()
    return {
        "operation": "apply",
        "database": str(db_path),
        "fragments": str(fragments_path),
        "input_sha256": source_sha256,
        "run_id": run_id,
        "dry_run": dry_run,
        "counts": counts,
        "applied_changes": [
            _change_report(change, action="update", dry_run=dry_run)
            for change in changes
        ],
        "skipped_items": skipped,
        "conflicts": [],
    }


def _read_json_value(value: Any, *, label: str) -> Any:
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError as exc:
            raise FragmentMetadataError(
                f"History contains invalid JSON for {label}: {exc}"
            ) from exc
    return value


def _load_history(
    con: duckdb.DuckDBPyConnection, run_id: str
) -> list[dict[str, Any]]:
    relation = con.execute(
        f"""
        SELECT run_id, applied_at, source_sha256, uniprot_ac, sequence_start,
               sequence_end, old_value_json, new_value_json
        FROM {_quote_identifier(HISTORY_TABLE)}
        WHERE run_id = ?
        ORDER BY applied_at, uniprot_ac, sequence_start, sequence_end
        """,
        [run_id],
    )
    names = [item[0] for item in relation.description]
    return [dict(zip(names, values)) for values in relation.fetchall()]


def _plan_restore(
    con: duckdb.DuckDBPyConnection, history: Sequence[Mapping[str, Any]]
) -> tuple[list[PlannedChange], list[dict[str, Any]], int, int]:
    changes: list[PlannedChange] = []
    conflicts: list[dict[str, Any]] = []
    found = 0
    unchanged = 0
    for row in history:
        accession = str(row["uniprot_ac"])
        start = int(row["sequence_start"])
        end = int(row["sequence_end"])
        expected = _read_json_value(
            row["new_value_json"], label=f"{accession}:{start}-{end} new value"
        )
        old_value = _read_json_value(
            row["old_value_json"], label=f"{accession}:{start}-{end} old value"
        )
        current_rows = con.execute(
            f"""
            SELECT protein_name FROM {_quote_identifier(FRAGMENT_TABLE)}
            WHERE uniprot_ac = ? AND sequence_start = ? AND sequence_end = ?
            """,
            [accession, start, end],
        ).fetchall()
        if not current_rows:
            conflicts.append({
                "uniprot_ac": accession, "sequence_start": start,
                "sequence_end": end, "reason_code": "fragment_unavailable",
                "reason": "recorded fragment is no longer present",
                "expected_current_value": expected,
            })
            continue
        found += 1
        current = current_rows[0][0]
        if not _json_equal(current, expected):
            conflicts.append({
                "uniprot_ac": accession, "sequence_start": start,
                "sequence_end": end, "reason_code": "current_value_conflict",
                "reason": (
                    "current value does not match the value written by this "
                    "run"
                ),
                "expected_current_value": expected,
                "actual_current_value": current,
            })
            continue
        if _json_equal(current, old_value):
            unchanged += 1
            continue
        changes.append(
            PlannedChange(accession, start, end, current, old_value)
        )
    return changes, conflicts, found, unchanged


def _restore_changes(
    con: duckdb.DuckDBPyConnection, changes: Sequence[PlannedChange]
) -> None:
    for change in changes:
        if change.new_value is None:
            updated = con.execute(
                f"""
                DELETE FROM {_quote_identifier(FRAGMENT_TABLE)}
                WHERE uniprot_ac = ? AND sequence_start = ?
                  AND sequence_end = ?
                  AND protein_name IS NOT DISTINCT FROM ?
                RETURNING uniprot_ac
                """,
                [
                    change.uniprot_ac, change.sequence_start,
                    change.sequence_end, change.old_value,
                ],
            ).fetchall()
        else:
            updated = con.execute(
                f"""
                UPDATE {_quote_identifier(FRAGMENT_TABLE)}
                SET protein_name = ?
                WHERE uniprot_ac = ? AND sequence_start = ?
                  AND sequence_end = ?
                  AND protein_name IS NOT DISTINCT FROM ?
                RETURNING uniprot_ac
                """,
                [
                    change.new_value, change.uniprot_ac,
                    change.sequence_start, change.sequence_end,
                    change.old_value,
                ],
            ).fetchall()
        if len(updated) != 1:
            raise FragmentMetadataError(
                "Fragment metadata changed while restoring "
                f"{change.uniprot_ac}:{change.sequence_start}-"
                f"{change.sequence_end}; "
                "the restore was rolled back."
            )


def restore_run(
    db_path: Path, run_id: str, *, dry_run: bool = False
) -> dict[str, Any]:
    """Restore one run only where no later change conflicts with it."""
    db_path = _resolved(db_path)
    if not db_path.is_file():
        raise FragmentMetadataError(
            f"DuckDB database does not exist: {db_path}"
        )
    try:
        con = duckdb.connect(str(db_path), read_only=dry_run)
    except duckdb.Error as exc:
        raise FragmentMetadataError(
            f"Cannot open DuckDB database {db_path}: {exc}"
        ) from exc
    history: list[dict[str, Any]] = []
    changes: list[PlannedChange] = []
    conflicts: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    source_sha256: str | None = None
    found = 0
    unchanged = 0
    transaction_open = False
    try:
        inspect_entry_schema(con)
        if not _table_exists(con, HISTORY_TABLE):
            skipped.append(_skipped_item(
                "", "history_unavailable",
                f"table {HISTORY_TABLE!r} does not exist; "
                "no run can be restored",
                run_id=run_id,
            ))
        elif not _table_exists(con, FRAGMENT_TABLE):
            _validate_table_schema(
                con, HISTORY_TABLE, REQUIRED_HISTORY_COLUMNS
            )
            skipped.append(_skipped_item(
                "", "fragment_table_unavailable",
                f"table {FRAGMENT_TABLE!r} does not exist; "
                "no run can be restored",
                run_id=run_id,
            ))
        else:
            _validate_table_schema(
                con, HISTORY_TABLE, REQUIRED_HISTORY_COLUMNS
            )
            _validate_fragment_table(con)
            if not dry_run:
                con.execute("BEGIN TRANSACTION")
                transaction_open = True
            history = _load_history(con, run_id)
            if not history:
                skipped.append(_skipped_item(
                    "", "unknown_run_id",
                    f"no history rows were found for run ID {run_id}",
                    run_id=run_id,
                ))
            else:
                hashes = {str(row["source_sha256"]) for row in history}
                source_sha256 = (
                    next(iter(hashes)) if len(hashes) == 1 else None
                )
                changes, conflicts, found, unchanged = _plan_restore(
                    con, history
                )
                if not dry_run:
                    _restore_changes(con, changes)
            if transaction_open:
                con.execute("COMMIT")
                transaction_open = False
    except duckdb.Error as exc:
        if transaction_open:
            con.execute("ROLLBACK")
        raise FragmentMetadataError(f"DuckDB operation failed: {exc}") from exc
    except Exception:
        if transaction_open:
            try:
                con.execute("ROLLBACK")
            except duckdb.Error:
                pass
        raise
    finally:
        con.close()
    counts = {
        "requested": len(history), "found": found, "updated": len(changes),
        "unchanged": unchanged, "skipped": len(skipped),
        "conflicted": len(conflicts),
    }
    return {
        "operation": "restore", "database": str(db_path), "fragments": None,
        "input_sha256": source_sha256, "run_id": run_id, "dry_run": dry_run,
        "counts": counts,
        "applied_changes": [
            _change_report(change, action="restore", dry_run=dry_run)
            for change in changes
        ],
        "skipped_items": skipped, "conflicts": conflicts,
    }


def write_report(path: Path, report: Mapping[str, Any]) -> None:
    path = _resolved(path)
    try:
        path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:
        raise FragmentMetadataError(
            f"Cannot write report file {path}: {exc}"
        ) from exc


def print_human_report(report: Mapping[str, Any]) -> None:
    prefix = "DRY RUN " if report.get("dry_run") else ""
    counts = report["counts"]
    print(f"{prefix}{report['operation']} report for {report['database']}")
    if report.get("input_sha256"):
        print(f"Input SHA-256: {report['input_sha256']}")
    if report.get("run_id"):
        label = (
            "Restored run ID" if report["operation"] == "restore" else "Run ID"
        )
        print(f"{label}: {report['run_id']}")
    print(
        "Counts: " + " ".join(
            f"{key}={counts[key]}" for key in (
                "requested", "found", "updated", "unchanged", "skipped",
                "conflicted",
            )
        )
    )
    for item in report.get("skipped_items", []):
        target = item.get("uniprot_ac") or item.get("run_id") or "<run>"
        if "sequence_start" in item:
            target += f":{item['sequence_start']}-{item['sequence_end']}"
        print(f"WARNING skipped {target}: {item['reason']}")
    for item in report.get("conflicts", []):
        print(
            "WARNING conflict "
            f"{item['uniprot_ac']}:{item['sequence_start']}-"
            f"{item['sequence_end']}: "
            f"{item['reason']}"
        )


def _fatal_report(args: argparse.Namespace, message: str) -> dict[str, Any]:
    return {
        "operation": "restore" if args.restore else "apply",
        "database": str(_resolved(args.db)),
        "fragments": (
            str(_resolved(args.fragments)) if args.fragments else None
        ),
        "input_sha256": None,
        "run_id": args.restore,
        "dry_run": bool(args.dry_run),
        "fatal_error": message,
        "counts": {
            "requested": 0, "found": 0, "updated": 0, "unchanged": 0,
            "skipped": 0, "conflicted": 0,
        },
        "applied_changes": [], "skipped_items": [], "conflicts": [],
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.fragments is not None:
            report = apply_fragments(
                args.db, args.fragments, dry_run=args.dry_run
            )
        else:
            report = restore_run(args.db, args.restore, dry_run=args.dry_run)
        if args.report:
            write_report(args.report, report)
        print_human_report(report)
        return 0
    except (FragmentMetadataError, OSError) as exc:
        message = str(exc)
        print(f"ERROR: {message}")
        if args.report:
            try:
                write_report(args.report, _fatal_report(args, message))
            except FragmentMetadataError as report_exc:
                print(f"ERROR: {report_exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
