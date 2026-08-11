#!/usr/bin/env python3
"""Apply and restore custom metadata in a UniProt DuckDB database."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import duckdb


HISTORY_TABLE = "custom_annotation_history"
SEQUENCE_COLUMNS = frozenset({"sequence", "length", "md5"})
AMINO_ACID_SEQUENCE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBJOUXZ]+$")
INTEGER_RANGES = {
    "TINYINT": (-(2**7), 2**7 - 1),
    "SMALLINT": (-(2**15), 2**15 - 1),
    "INTEGER": (-(2**31), 2**31 - 1),
    "BIGINT": (-(2**63), 2**63 - 1),
    "HUGEINT": (-(2**127), 2**127 - 1),
    "UTINYINT": (0, 2**8 - 1),
    "USMALLINT": (0, 2**16 - 1),
    "UINTEGER": (0, 2**32 - 1),
    "UBIGINT": (0, 2**64 - 1),
    "UHUGEINT": (0, 2**128 - 1),
}
REQUIRED_HISTORY_COLUMNS = {
    "run_id": "VARCHAR",
    "applied_at": "TIMESTAMP",
    "source_sha256": "VARCHAR",
    "primary_ac": "VARCHAR",
    "column_name": "VARCHAR",
    "old_value_json": "JSON",
    "new_value_json": "JSON",
}


class CustomMetadataError(RuntimeError):
    """A fatal, user-facing custom metadata command error."""


class InvalidAnnotation(ValueError):
    """A non-fatal value validation error."""


@dataclass(frozen=True)
class EntryColumn:
    name: str
    duckdb_type: str


@dataclass(frozen=True)
class PlannedChange:
    primary_ac: str
    column_name: str
    old_value: Any
    new_value: Any


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply custom annotations to existing UniProt entry columns, "
            "or restore a previously applied run."
        )
    )
    parser.add_argument(
        "--db",
        required=True,
        type=Path,
        help="Existing DuckDB database containing the entry table.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--annotations",
        type=Path,
        help="JSON object keyed by entry.primary_ac.",
    )
    mode.add_argument(
        "--restore",
        metavar="RUN_ID",
        help="Restore all non-conflicting field updates recorded for RUN_ID.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and report changes without writing to the database.",
    )
    parser.add_argument(
        "--report",
        type=Path,
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


def load_annotations(path: Path) -> tuple[dict[str, Any], str]:
    """Load an annotation object and return it with its source SHA-256."""
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise CustomMetadataError(
            f"Cannot read annotations file {path}: {exc}"
        ) from exc
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CustomMetadataError(
            f"Annotations file {path} is not valid JSON: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise CustomMetadataError(
            "The annotations JSON root must be an object keyed by "
            "UniProt accession."
        )
    return value, hashlib.sha256(raw).hexdigest()


def _describe_table(
    con: duckdb.DuckDBPyConnection, table_name: str
) -> dict[str, EntryColumn]:
    try:
        rows = con.execute(
            f"DESCRIBE {_quote_identifier(table_name)}"
        ).fetchall()
    except duckdb.Error as exc:
        raise CustomMetadataError(
            f"Required table {table_name!r} is absent or unreadable: {exc}"
        ) from exc
    return {
        str(row[0]): EntryColumn(
            name=str(row[0]), duckdb_type=str(row[1]).upper()
        )
        for row in rows
    }


def inspect_entry_schema(
    con: duckdb.DuckDBPyConnection,
) -> dict[str, EntryColumn]:
    """Inspect and minimally validate the target entry schema at runtime."""
    columns = _describe_table(con, "entry")
    primary_ac = columns.get("primary_ac")
    if primary_ac is None:
        raise CustomMetadataError(
            "Incompatible entry schema: required lookup column "
            "'primary_ac' is absent."
        )
    if primary_ac.duckdb_type != "VARCHAR":
        raise CustomMetadataError(
            "Incompatible entry schema: 'primary_ac' must have DuckDB type "
            f"VARCHAR, found {primary_ac.duckdb_type}."
        )
    return columns


def _table_exists(con: duckdb.DuckDBPyConnection, table_name: str) -> bool:
    rows = con.execute(
        """
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema = current_schema() AND table_name = ?
        LIMIT 1
        """,
        [table_name],
    ).fetchall()
    return bool(rows)


def _validate_history_schema(con: duckdb.DuckDBPyConnection) -> None:
    columns = _describe_table(con, HISTORY_TABLE)
    for name, expected_type in REQUIRED_HISTORY_COLUMNS.items():
        actual = columns.get(name)
        if actual is None:
            raise CustomMetadataError(
                f"Incompatible {HISTORY_TABLE} schema: required column "
                f"{name!r} is absent."
            )
        if actual.duckdb_type != expected_type:
            raise CustomMetadataError(
                f"Incompatible {HISTORY_TABLE} schema: column {name!r} must "
                f"have type {expected_type}, found {actual.duckdb_type}."
            )


def _create_or_validate_history(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {_quote_identifier(HISTORY_TABLE)}(
            run_id VARCHAR,
            applied_at TIMESTAMP,
            source_sha256 VARCHAR,
            primary_ac VARCHAR,
            column_name VARCHAR,
            old_value_json JSON,
            new_value_json JSON
        )
        """
    )
    _validate_history_schema(con)


def _type_kind(duckdb_type: str) -> str | None:
    if duckdb_type == "VARCHAR":
        return "string"
    if duckdb_type == "BOOLEAN":
        return "boolean"
    if duckdb_type in INTEGER_RANGES:
        return "integer"
    if duckdb_type == "VARCHAR[]":
        return "string_list"
    return None


def normalize_value(value: Any, column: EntryColumn) -> Any:
    """Strictly validate one JSON value against its inspected DuckDB type."""
    if value is None:
        raise InvalidAnnotation("JSON null is not an annotation value")

    kind = _type_kind(column.duckdb_type)
    if kind is None:
        raise InvalidAnnotation(
            f"DuckDB type {column.duckdb_type} is not supported"
        )
    if kind == "string":
        if not isinstance(value, str):
            raise InvalidAnnotation("expected a JSON string")
        return value
    if kind == "boolean":
        if type(value) is not bool:
            raise InvalidAnnotation("expected a JSON boolean")
        return value
    if kind == "integer":
        if type(value) is not int:
            raise InvalidAnnotation(
                "expected a JSON integer (booleans are not integers)"
            )
        lower, upper = INTEGER_RANGES[column.duckdb_type]
        if not lower <= value <= upper:
            raise InvalidAnnotation(
                f"integer is outside the range of {column.duckdb_type}"
            )
        return value
    if isinstance(value, str):
        return [value]
    if not isinstance(value, list) or any(
        not isinstance(item, str) for item in value
    ):
        raise InvalidAnnotation(
            "expected a JSON string or a JSON list containing only strings"
        )
    return value


def _validate_sequence_bundle(
    requested: Mapping[str, Any],
    columns: Mapping[str, EntryColumn],
) -> dict[str, Any]:
    missing_values = sorted(SEQUENCE_COLUMNS.difference(requested))
    if missing_values:
        raise InvalidAnnotation(
            "sequence, length, and md5 must be supplied together; missing "
            + ", ".join(missing_values)
        )
    missing_columns = sorted(SEQUENCE_COLUMNS.difference(columns))
    if missing_columns:
        raise InvalidAnnotation(
            "target entry table has no column(s): "
            + ", ".join(missing_columns)
        )

    try:
        sequence = normalize_value(requested["sequence"], columns["sequence"])
        length = normalize_value(requested["length"], columns["length"])
        md5 = normalize_value(requested["md5"], columns["md5"])
    except InvalidAnnotation as exc:
        raise InvalidAnnotation(
            f"invalid sequence bundle type: {exc}"
        ) from exc

    if not isinstance(sequence, str):
        raise InvalidAnnotation("sequence column must be a VARCHAR")
    if type(length) is not int:
        raise InvalidAnnotation("length column must be an integer type")
    if not isinstance(md5, str):
        raise InvalidAnnotation("md5 column must be a VARCHAR")
    if not sequence:
        raise InvalidAnnotation("sequence must be non-empty")
    if AMINO_ACID_SEQUENCE.fullmatch(sequence) is None:
        raise InvalidAnnotation(
            "sequence must contain only uppercase recognised one-letter "
            "residue codes"
        )
    if length != len(sequence):
        raise InvalidAnnotation(
            f"length {length} does not match sequence length {len(sequence)}"
        )
    expected_md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
    if md5.lower() != expected_md5:
        raise InvalidAnnotation(
            f"md5 does not match sequence (expected {expected_md5})"
        )
    return {"sequence": sequence, "length": length, "md5": expected_md5}


def _fetch_requested_rows(
    con: duckdb.DuckDBPyConnection, accessions: Sequence[str]
) -> dict[str, dict[str, Any]]:
    if not accessions:
        return {}
    placeholders = ",".join("?" for _ in accessions)
    relation = con.execute(
        f"SELECT * FROM entry WHERE primary_ac IN ({placeholders})",
        list(accessions),
    )
    names = [item[0] for item in relation.description]
    fetched: dict[str, dict[str, Any]] = {}
    for values in relation.fetchall():
        row = dict(zip(names, values))
        accession = str(row["primary_ac"])
        if accession in fetched:
            raise CustomMetadataError(
                "Incompatible entry data: primary_ac must uniquely identify "
                "a row; "
                f"found duplicate {accession!r}."
            )
        fetched[accession] = row
    return fetched


def _skipped_item(
    primary_ac: str,
    reason_code: str,
    reason: str,
    column_name: str | None = None,
    **details: Any,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "primary_ac": primary_ac,
        "column_name": column_name,
        "reason_code": reason_code,
        "reason": reason,
    }
    item.update(details)
    return item


def plan_annotation_changes(
    con: duckdb.DuckDBPyConnection,
    annotations: Mapping[str, Any],
    columns: Mapping[str, EntryColumn],
) -> tuple[list[PlannedChange], list[dict[str, Any]], dict[str, int]]:
    """Build a write plan and retain independent validation issues."""
    rows = _fetch_requested_rows(con, list(annotations))
    changes: list[PlannedChange] = []
    skipped: list[dict[str, Any]] = []
    unchanged = 0
    requested_fields = 0
    valid_fields = 0

    for accession, requested in annotations.items():
        if not isinstance(requested, dict):
            skipped.append(
                _skipped_item(
                    accession,
                    "invalid_accession_value",
                    "annotation value for an accession must be a JSON object",
                )
            )
            continue
        requested_fields += len(requested)
        row = rows.get(accession)
        if row is None:
            skipped.append(
                _skipped_item(
                    accession,
                    "unknown_accession",
                    "accession was not found in entry.primary_ac",
                    requested_columns=list(requested),
                )
            )
            continue

        requested_sequence_columns = SEQUENCE_COLUMNS.intersection(requested)
        if requested_sequence_columns:
            try:
                bundle = _validate_sequence_bundle(requested, columns)
            except InvalidAnnotation as exc:
                skipped.append(
                    _skipped_item(
                        accession,
                        "invalid_sequence_bundle",
                        str(exc),
                        column_name="sequence,length,md5",
                        supplied_columns=sorted(requested_sequence_columns),
                    )
                )
            else:
                for column_name in ("sequence", "length", "md5"):
                    new_value = bundle[column_name]
                    valid_fields += 1
                    if _json_equal(row[column_name], new_value):
                        unchanged += 1
                    else:
                        changes.append(
                            PlannedChange(
                                primary_ac=accession,
                                column_name=column_name,
                                old_value=row[column_name],
                                new_value=new_value,
                            )
                        )

        for column_name, raw_value in requested.items():
            if column_name in SEQUENCE_COLUMNS:
                continue
            if column_name == "primary_ac":
                skipped.append(
                    _skipped_item(
                        accession,
                        "immutable_column",
                        "primary_ac is the immutable lookup key",
                        column_name=column_name,
                    )
                )
                continue
            column = columns.get(column_name)
            if column is None:
                skipped.append(
                    _skipped_item(
                        accession,
                        "unknown_column",
                        "column is not present in the target entry table",
                        column_name=column_name,
                    )
                )
                continue
            try:
                new_value = normalize_value(raw_value, column)
            except InvalidAnnotation as exc:
                skipped.append(
                    _skipped_item(
                        accession,
                        "invalid_value",
                        str(exc),
                        column_name=column_name,
                        duckdb_type=column.duckdb_type,
                    )
                )
                continue
            valid_fields += 1
            if _json_equal(row[column_name], new_value):
                unchanged += 1
            else:
                changes.append(
                    PlannedChange(
                        primary_ac=accession,
                        column_name=column_name,
                        old_value=row[column_name],
                        new_value=new_value,
                    )
                )

    counts = {
        "requested": len(annotations),
        "found": sum(accession in rows for accession in annotations),
        "updated": len(changes),
        "unchanged": unchanged,
        "skipped": len(skipped),
        "conflicted": 0,
        "requested_fields": requested_fields,
        "valid_fields": valid_fields,
    }
    return changes, skipped, counts


def _change_report(
    change: PlannedChange, *, action: str, dry_run: bool
) -> dict[str, Any]:
    return {
        "primary_ac": change.primary_ac,
        "column_name": change.column_name,
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
        _create_or_validate_history(con)
        for change in changes:
            con.execute(
                f"""
                INSERT INTO {_quote_identifier(HISTORY_TABLE)}
                    (run_id, applied_at, source_sha256, primary_ac,
                     column_name,
                     old_value_json, new_value_json)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    run_id,
                    applied_at,
                    source_sha256,
                    change.primary_ac,
                    change.column_name,
                    _json_text(change.old_value),
                    _json_text(change.new_value),
                ],
            )
        for change in changes:
            column = _quote_identifier(change.column_name)
            updated = con.execute(
                f"""
                UPDATE entry
                SET {column} = ?
                WHERE primary_ac = ?
                  AND {column} IS NOT DISTINCT FROM ?
                RETURNING primary_ac
                """,
                [change.new_value, change.primary_ac, change.old_value],
            ).fetchall()
            if len(updated) != 1:
                raise CustomMetadataError(
                    "Entry data changed after validation for "
                    f"{change.primary_ac}.{change.column_name}; "
                    "no updates were applied."
                )
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except duckdb.Error:
            pass
        raise


def apply_annotations(
    db_path: Path,
    annotations_path: Path,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Validate and optionally apply one annotation document."""
    db_path = _resolved(db_path)
    annotations_path = _resolved(annotations_path)
    if not db_path.is_file():
        raise CustomMetadataError(f"DuckDB database does not exist: {db_path}")
    annotations, source_sha256 = load_annotations(annotations_path)

    try:
        con = duckdb.connect(str(db_path), read_only=dry_run)
    except duckdb.Error as exc:
        raise CustomMetadataError(
            f"Cannot open DuckDB database {db_path}: {exc}"
        ) from exc
    try:
        columns = inspect_entry_schema(con)
        changes, skipped, counts = plan_annotation_changes(
            con, annotations, columns
        )
        run_id = None if dry_run else str(uuid.uuid4())
        if not dry_run:
            _apply_changes(
                con,
                changes,
                run_id=run_id,
                source_sha256=source_sha256,
            )
    except duckdb.Error as exc:
        raise CustomMetadataError(f"DuckDB operation failed: {exc}") from exc
    finally:
        con.close()

    return {
        "operation": "apply",
        "database": str(db_path),
        "annotations": str(annotations_path),
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
            raise CustomMetadataError(
                f"History contains invalid JSON for {label}: {exc}"
            ) from exc
    return value


def _load_history(
    con: duckdb.DuckDBPyConnection, run_id: str
) -> list[dict[str, Any]]:
    relation = con.execute(
        f"""
        SELECT run_id, applied_at, source_sha256, primary_ac, column_name,
               old_value_json, new_value_json
        FROM {_quote_identifier(HISTORY_TABLE)}
        WHERE run_id = ?
        ORDER BY applied_at, primary_ac, column_name
        """,
        [run_id],
    )
    names = [item[0] for item in relation.description]
    return [dict(zip(names, values)) for values in relation.fetchall()]


def _plan_restore(
    con: duckdb.DuckDBPyConnection,
    history: Sequence[Mapping[str, Any]],
    columns: Mapping[str, EntryColumn],
) -> tuple[list[PlannedChange], list[dict[str, Any]], int, int]:
    changes: list[PlannedChange] = []
    conflicts: list[dict[str, Any]] = []
    found = 0
    unchanged = 0

    for history_row in history:
        accession = str(history_row["primary_ac"])
        column_name = str(history_row["column_name"])
        expected = _read_json_value(
            history_row["new_value_json"],
            label=f"{accession}.{column_name} new value",
        )
        restore_value = _read_json_value(
            history_row["old_value_json"],
            label=f"{accession}.{column_name} old value",
        )
        if column_name == "primary_ac" or column_name not in columns:
            conflicts.append(
                {
                    "primary_ac": accession,
                    "column_name": column_name,
                    "reason_code": "column_unavailable",
                    "reason": (
                        "recorded entry column is no longer available "
                        "or mutable"
                    ),
                    "expected_current_value": expected,
                }
            )
            continue

        column = _quote_identifier(column_name)
        current_rows = con.execute(
            f"SELECT {column} FROM entry WHERE primary_ac = ?",
            [accession],
        ).fetchall()
        if len(current_rows) > 1:
            raise CustomMetadataError(
                "Incompatible entry data: primary_ac must uniquely identify "
                "a row; "
                f"found duplicate {accession!r}."
            )
        if not current_rows:
            conflicts.append(
                {
                    "primary_ac": accession,
                    "column_name": column_name,
                    "reason_code": "accession_unavailable",
                    "reason": (
                        "recorded accession is no longer present in entry"
                    ),
                    "expected_current_value": expected,
                }
            )
            continue

        found += 1
        current = current_rows[0][0]
        if not _json_equal(current, expected):
            conflicts.append(
                {
                    "primary_ac": accession,
                    "column_name": column_name,
                    "reason_code": "current_value_conflict",
                    "reason": (
                        "current value does not match the value written by "
                        "this run"
                    ),
                    "expected_current_value": expected,
                    "actual_current_value": current,
                }
            )
            continue
        if _json_equal(current, restore_value):
            unchanged += 1
            continue
        changes.append(
            PlannedChange(
                primary_ac=accession,
                column_name=column_name,
                old_value=current,
                new_value=restore_value,
            )
        )
    return changes, conflicts, found, unchanged


def _restore_changes(
    con: duckdb.DuckDBPyConnection,
    changes: Sequence[PlannedChange],
) -> None:
    for change in changes:
        column = _quote_identifier(change.column_name)
        updated = con.execute(
            f"""
            UPDATE entry
            SET {column} = ?
            WHERE primary_ac = ?
              AND {column} IS NOT DISTINCT FROM ?
            RETURNING primary_ac
            """,
            [change.new_value, change.primary_ac, change.old_value],
        ).fetchall()
        if len(updated) != 1:
            raise CustomMetadataError(
                "Entry data changed while restoring "
                f"{change.primary_ac}.{change.column_name}; "
                "the restore was rolled back."
            )


def restore_run(
    db_path: Path,
    run_id: str,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Restore fields from a run when their current values remain unchanged."""
    db_path = _resolved(db_path)
    if not db_path.is_file():
        raise CustomMetadataError(f"DuckDB database does not exist: {db_path}")
    try:
        con = duckdb.connect(str(db_path), read_only=dry_run)
    except duckdb.Error as exc:
        raise CustomMetadataError(
            f"Cannot open DuckDB database {db_path}: {exc}"
        ) from exc

    history: list[dict[str, Any]] = []
    changes: list[PlannedChange] = []
    conflicts: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    found = 0
    unchanged = 0
    source_sha256: str | None = None
    transaction_open = False
    try:
        columns = inspect_entry_schema(con)
        if _table_exists(con, HISTORY_TABLE):
            _validate_history_schema(con)
            if not dry_run:
                con.execute("BEGIN TRANSACTION")
                transaction_open = True
            history = _load_history(con, run_id)
            if history:
                hashes = {str(row["source_sha256"]) for row in history}
                source_sha256 = (
                    next(iter(hashes)) if len(hashes) == 1 else None
                )
                changes, conflicts, found, unchanged = _plan_restore(
                    con, history, columns
                )
                if not dry_run:
                    _restore_changes(con, changes)
            else:
                skipped.append(
                    _skipped_item(
                        "",
                        "unknown_run_id",
                        f"no history rows were found for run ID {run_id}",
                        run_id=run_id,
                    )
                )
            if transaction_open:
                con.execute("COMMIT")
                transaction_open = False
        else:
            skipped.append(
                _skipped_item(
                    "",
                    "history_unavailable",
                    f"table {HISTORY_TABLE!r} does not exist; "
                    "no run can be restored",
                    run_id=run_id,
                )
            )
    except duckdb.Error as exc:
        if transaction_open:
            try:
                con.execute("ROLLBACK")
            except duckdb.Error:
                pass
        raise CustomMetadataError(f"DuckDB operation failed: {exc}") from exc
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
        "requested": len(history),
        "found": found,
        "updated": len(changes),
        "unchanged": unchanged,
        "skipped": len(skipped),
        "conflicted": len(conflicts),
        "requested_fields": len(history),
        "valid_fields": found,
    }
    return {
        "operation": "restore",
        "database": str(db_path),
        "annotations": None,
        "input_sha256": source_sha256,
        "run_id": run_id,
        "dry_run": dry_run,
        "counts": counts,
        "applied_changes": [
            _change_report(change, action="restore", dry_run=dry_run)
            for change in changes
        ],
        "skipped_items": skipped,
        "conflicts": conflicts,
    }


def write_report(path: Path, report: Mapping[str, Any]) -> None:
    path = _resolved(path)
    try:
        path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:
        raise CustomMetadataError(
            f"Cannot write report file {path}: {exc}"
        ) from exc


def print_human_report(report: Mapping[str, Any]) -> None:
    operation = str(report["operation"])
    dry_run = bool(report.get("dry_run"))
    prefix = "DRY RUN " if dry_run else ""
    counts = report["counts"]
    print(f"{prefix}{operation} report for {report['database']}")
    if report.get("input_sha256"):
        print(f"Input SHA-256: {report['input_sha256']}")
    if report.get("run_id"):
        label = "Restored run ID" if operation == "restore" else "Run ID"
        print(f"{label}: {report['run_id']}")
    print(
        "Counts: "
        f"requested={counts['requested']} "
        f"found={counts['found']} "
        f"updated={counts['updated']} "
        f"unchanged={counts['unchanged']} "
        f"skipped={counts['skipped']} "
        f"conflicted={counts['conflicted']}"
    )
    for item in report.get("skipped_items", []):
        target = item.get("primary_ac") or item.get("run_id") or "<run>"
        if item.get("column_name"):
            target += f".{item['column_name']}"
        print(f"WARNING skipped {target}: {item['reason']}")
    for item in report.get("conflicts", []):
        target = f"{item['primary_ac']}.{item['column_name']}"
        print(f"WARNING conflict {target}: {item['reason']}")


def _fatal_report(args: argparse.Namespace, message: str) -> dict[str, Any]:
    return {
        "operation": "restore" if args.restore else "apply",
        "database": str(_resolved(args.db)),
        "annotations": (
            str(_resolved(args.annotations)) if args.annotations else None
        ),
        "input_sha256": None,
        "run_id": args.restore,
        "dry_run": bool(args.dry_run),
        "fatal_error": message,
        "counts": {
            "requested": 0,
            "found": 0,
            "updated": 0,
            "unchanged": 0,
            "skipped": 0,
            "conflicted": 0,
            "requested_fields": 0,
            "valid_fields": 0,
        },
        "applied_changes": [],
        "skipped_items": [],
        "conflicts": [],
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.annotations is not None:
            report = apply_annotations(
                args.db,
                args.annotations,
                dry_run=args.dry_run,
            )
        else:
            report = restore_run(
                args.db,
                args.restore,
                dry_run=args.dry_run,
            )
        if args.report:
            write_report(args.report, report)
        print_human_report(report)
        return 0
    except (CustomMetadataError, OSError) as exc:
        message = str(exc)
        print(f"ERROR: {message}")
        if args.report:
            try:
                write_report(args.report, _fatal_report(args, message))
            except CustomMetadataError as report_exc:
                print(f"ERROR: {report_exc}")
        return 1
    except Exception as exc:  # pragma: no cover - final CLI guard
        message = f"Unexpected custom metadata failure: {exc}"
        print(f"ERROR: {message}")
        if args.report:
            try:
                write_report(args.report, _fatal_report(args, message))
            except CustomMetadataError as report_exc:
                print(f"ERROR: {report_exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
