#!/usr/bin/env python3
"""Stream-normalise chain manifests and add fragment names from DuckDB.

The collaborator's wide manifest contains chain A/B information on one row.
This command turns it into the canonical one-row-per-chain form without ever
accumulating the normalised manifest in memory.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import duckdb


CANONICAL_COLUMNS = (
    "model_entity_id",
    "entity_id",
    "chain_id",
    "uniprot_ac",
    "is_fragment",
    "sequence_start",
    "sequence_end",
    "protein_name",
)
WIDE_CHAIN_COLUMNS = frozenset(
    {
        "afdb_id",
        "chain_a",
        "chain_b",
        "chain_a_id",
        "chain_b_id",
        "chain_a_uniprot",
        "chain_b_uniprot",
        "chain_a_is_polyprotein",
        "chain_b_is_polyprotein",
        "chain_a_is_polyprotein?",
        "chain_b_is_polyprotein?",
        "chain_a_start",
        "chain_b_start",
        "chain_a_end",
        "chain_b_end",
        "chain_a_length",
        "chain_b_length",
    }
)
TRUE_VALUES = frozenset({"1", "true", "t", "yes", "y"})
FALSE_VALUES = frozenset({"0", "false", "f", "no", "n", ""})
REPORT_SAMPLE_LIMIT = 1_000


class ManifestEnrichmentError(RuntimeError):
    """A fatal, user-facing manifest enrichment error."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the manifest enricher."""
    parser = argparse.ArgumentParser(
        description=(
            "Stream-normalise a collaborator-wide or canonical chain "
            "manifest and enrich fragment protein names from DuckDB."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Collaborator-wide or canonical CSV/TSV manifest.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Destination canonical CSV/TSV manifest.",
    )
    parser.add_argument(
        "--db",
        required=True,
        type=Path,
        help="DuckDB database containing fragment_metadata.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Return a non-zero status and leave no output file when a "
            "fragment has no unambiguous protein name."
        ),
    )
    parser.add_argument(
        "--report",
        type=Path,
        help="Optional path for a JSON report.",
    )
    return parser.parse_args(argv)


def _resolved(path: Path) -> Path:
    return path.expanduser().resolve()


def _read_dialect(path: Path) -> csv.Dialect:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            sample = handle.read(16_384)
    except OSError as exc:
        raise ManifestEnrichmentError(
            f"Cannot read input manifest {path}: {exc}"
        ) from exc
    if not sample:
        raise ManifestEnrichmentError(f"Input manifest {path} is empty.")
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        # CSV is the public default; a one-column header is still valid CSV.
        return csv.excel


def _headers(path: Path, dialect: csv.Dialect) -> list[str]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, dialect=dialect)
            headers = next(reader, None)
    except OSError as exc:
        raise ManifestEnrichmentError(
            f"Cannot read input manifest {path}: {exc}"
        ) from exc
    if not headers or any(not name for name in headers):
        raise ManifestEnrichmentError(
            "Input manifest must have a non-empty header row."
        )
    if len(headers) != len(set(headers)):
        raise ManifestEnrichmentError("Input manifest has duplicate columns.")
    return headers


def detect_schema(headers: Iterable[str]) -> str:
    """Return ``wide`` or ``canonical`` for supported input headers."""
    names = set(headers)
    if {"model_entity_id", "chain_id", "uniprot_ac"} <= names:
        return "canonical"
    has_a_occupancy = "chain_a_id" in names or "chain_a" in names
    has_a_accession = "chain_a_uniprot" in names
    if "afdb_id" in names and has_a_occupancy and has_a_accession:
        return "wide"
    raise ManifestEnrichmentError(
        "Unsupported manifest schema. Expected canonical columns "
        "model_entity_id, chain_id, uniprot_ac, or collaborator-wide "
        "columns afdb_id, chain_a_id, chain_a_uniprot."
    )


def _normalise_boolean(value: str, *, field: str, row_number: int) -> bool:
    text = value.strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ManifestEnrichmentError(
        f"Row {row_number}: {field} must be a boolean value, found {value!r}."
    )


def _normalise_range_value(
    value: str, *, field: str, row_number: int, required: bool
) -> str:
    text = value.strip()
    if not text:
        if required:
            raise ManifestEnrichmentError(
                f"Row {row_number}: {field} is required for a fragment."
            )
        return ""
    try:
        number = int(text)
    except ValueError as exc:
        raise ManifestEnrichmentError(
            f"Row {row_number}: {field} must be an integer, found {value!r}."
        ) from exc
    if number < 1:
        raise ManifestEnrichmentError(
            f"Row {row_number}: {field} must be positive, found {value!r}."
        )
    return str(number)


def _component_key(
    accession: str, is_fragment: bool, start: str, end: str
) -> tuple[str, ...]:
    if is_fragment:
        if int(start) > int(end):
            raise ManifestEnrichmentError(
                f"Fragment range for {accession} has start {start} after "
                f"end {end}."
            )
        return accession, start, end
    return (accession,)


class _EntityAssigner:
    """Assign sequential entity IDs for one contiguous model at a time."""

    def __init__(self) -> None:
        self._model_id: str | None = None
        self._ids: dict[tuple[str, ...], str] = {}

    def entity_id(self, model_id: str, component: tuple[str, ...]) -> str:
        if model_id != self._model_id:
            self._model_id = model_id
            self._ids = {}
        if component not in self._ids:
            self._ids[component] = str(len(self._ids) + 1)
        return self._ids[component]


def _fragment_lookup(db_path: Path) -> dict[tuple[str, str, str], list[str]]:
    """Load the small fragment-name table, preserving malformed duplicates."""
    if not db_path.is_file():
        raise ManifestEnrichmentError(
            f"DuckDB database does not exist: {db_path}"
        )
    try:
        con = duckdb.connect(str(db_path), read_only=True)
    except duckdb.Error as exc:
        raise ManifestEnrichmentError(
            f"Cannot open DuckDB database {db_path}: {exc}"
        ) from exc
    try:
        exists = con.execute(
            """
            SELECT 1 FROM information_schema.tables
            WHERE table_schema = current_schema()
              AND table_name = 'fragment_metadata'
            LIMIT 1
            """
        ).fetchone()
        if not exists:
            return {}
        columns = {
            row[0]
            for row in con.execute("DESCRIBE fragment_metadata").fetchall()
        }
        expected = {
            "uniprot_ac", "sequence_start", "sequence_end", "protein_name"
        }
        if not expected <= columns:
            raise ManifestEnrichmentError(
                "fragment_metadata is missing required columns: "
                + ", ".join(sorted(expected - columns))
            )
        rows = con.execute(
            """
            SELECT uniprot_ac, sequence_start, sequence_end, protein_name
            FROM fragment_metadata
            """
        ).fetchall()
    except duckdb.Error as exc:
        raise ManifestEnrichmentError(
            f"Cannot read fragment_metadata from {db_path}: {exc}"
        ) from exc
    finally:
        con.close()
    lookup: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for accession, start, end, name in rows:
        key = (str(accession), str(start), str(end))
        lookup[key].append(str(name))
    return dict(lookup)


def _base_row(
    raw: Mapping[str, str], schema: str, row_number: int
) -> tuple[dict[str, str], dict[str, str]]:
    """Normalise one source row, returning base values and preserved extras."""
    if None in raw:
        raise ManifestEnrichmentError(
            f"Row {row_number}: contains more fields than the header."
        )
    if schema == "canonical":
        model_id = (raw.get("model_entity_id") or "").strip()
        chain_id = (raw.get("chain_id") or "").strip()
        accession = (raw.get("uniprot_ac") or "").strip()
        if not model_id or not chain_id or not accession:
            raise ManifestEnrichmentError(
                f"Row {row_number}: model_entity_id, chain_id, and "
                "uniprot_ac must all be populated."
            )
        fragment_value = raw.get("is_fragment") or ""
        is_fragment = _normalise_boolean(
            fragment_value, field="is_fragment", row_number=row_number
        )
        start = _normalise_range_value(
            raw.get("sequence_start") or "", field="sequence_start",
            row_number=row_number, required=is_fragment,
        )
        end = _normalise_range_value(
            raw.get("sequence_end") or "", field="sequence_end",
            row_number=row_number, required=is_fragment,
        )
        if is_fragment:
            _component_key(accession, True, start, end)
        base = {
            "model_entity_id": model_id,
            "chain_id": chain_id,
            "uniprot_ac": accession,
            "is_fragment": str(is_fragment).lower(),
            "sequence_start": start,
            "sequence_end": end,
            "protein_name": (raw.get("protein_name") or "").strip(),
        }
        extras = {
            key: value or "" for key, value in raw.items()
            if key not in CANONICAL_COLUMNS
        }
        return base, extras
    raise AssertionError("_base_row only handles canonical rows")


def _wide_rows(
    raw: Mapping[str, str], row_number: int
) -> Iterable[tuple[dict[str, str], dict[str, str]]]:
    """Turn one collaborator-wide row into zero, one, or two chain rows."""
    if None in raw:
        raise ManifestEnrichmentError(
            f"Row {row_number}: contains more fields than the header."
        )
    model_id = (raw.get("afdb_id") or "").strip()
    if not model_id:
        raise ManifestEnrichmentError(f"Row {row_number}: afdb_id is empty.")
    extras = {
        key: value or "" for key, value in raw.items()
        if not _is_wide_chain_column(key) and key not in CANONICAL_COLUMNS
    }
    for chain_label in ("a", "b"):
        occupancy_column = (
            f"chain_{chain_label}_id"
            if f"chain_{chain_label}_id" in raw
            else f"chain_{chain_label}"
        )
        occupancy = (raw.get(occupancy_column) or "").strip()
        if not occupancy:
            continue
        accession = (raw.get(f"chain_{chain_label}_uniprot") or "").strip()
        if not accession:
            raise ManifestEnrichmentError(
                f"Row {row_number}: chain_{chain_label}_uniprot is empty for "
                f"populated {occupancy_column!r}."
            )
        fragment_column = f"chain_{chain_label}_is_polyprotein?"
        if fragment_column not in raw:
            fragment_column = f"chain_{chain_label}_is_polyprotein"
        is_fragment = _normalise_boolean(
            raw.get(fragment_column) or "", field=fragment_column,
            row_number=row_number,
        )
        start = _normalise_range_value(
            raw.get(f"chain_{chain_label}_start") or "",
            field=f"chain_{chain_label}_start", row_number=row_number,
            required=is_fragment,
        )
        end = _normalise_range_value(
            raw.get(f"chain_{chain_label}_end") or "",
            field=f"chain_{chain_label}_end", row_number=row_number,
            required=is_fragment,
        )
        if is_fragment:
            _component_key(accession, True, start, end)
        yield {
            "model_entity_id": model_id,
            # The collaborator's source ID is an occupancy/source identifier,
            # never an authoritative protein name.  The canonical chain ID is
            # the literal A or B selected by that populated source column.
            "chain_id": chain_label.upper(),
            "uniprot_ac": accession,
            "is_fragment": str(is_fragment).lower(),
            "sequence_start": start,
            "sequence_end": end,
            "protein_name": "",
        }, extras


def _record_lookup(
    report: dict[str, Any], key: tuple[str, str, str], kind: str
) -> None:
    report["counts"][f"{kind}_fragment_lookups"] += 1
    records = report[f"{kind}_fragment_lookups"]
    if len(records) < REPORT_SAMPLE_LIMIT:
        records.append(
            {
                "uniprot_ac": key[0],
                "sequence_start": int(key[1]),
                "sequence_end": int(key[2]),
            }
        )
    else:
        report[f"{kind}_fragment_lookups_truncated"] += 1


def _enrich_name(
    row: dict[str, str], lookup: Mapping[tuple[str, str, str], list[str]],
    report: dict[str, Any],
) -> None:
    if row["is_fragment"] != "true":
        return
    report["counts"]["fragment_rows"] += 1
    if row["protein_name"]:
        report["counts"]["fragment_names_preserved"] += 1
        return
    key = (row["uniprot_ac"], row["sequence_start"], row["sequence_end"])
    names = lookup.get(key, [])
    if len(names) == 1:
        row["protein_name"] = names[0]
        report["counts"]["fragment_names_enriched"] += 1
    elif not names:
        _record_lookup(report, key, "unmatched")
    else:
        _record_lookup(report, key, "ambiguous")


def _output_fieldnames(headers: Sequence[str], schema: str) -> list[str]:
    source_extras = (
        (name for name in headers if name not in CANONICAL_COLUMNS)
        if schema == "canonical"
        else (
            name
            for name in headers
            if not _is_wide_chain_column(name)
            and name not in CANONICAL_COLUMNS
        )
    )
    return list(CANONICAL_COLUMNS) + list(source_extras)


def _is_wide_chain_column(name: str) -> bool:
    """Whether a wide-source field describes an individual A/B chain."""
    return (
        name in WIDE_CHAIN_COLUMNS
        or name.startswith("chain_a_")
        or name.startswith("chain_b_")
    )


def enrich_manifest(
    input_path: Path,
    output_path: Path,
    db_path: Path,
    *,
    strict: bool = False,
) -> dict[str, Any]:
    """Create an enriched canonical manifest and return its JSON report.

    A strict failure is represented in the returned report.  The destination
    is only replaced after a successful strict check, so consumers never see
    a partially valid strict output.
    """
    input_path = _resolved(input_path)
    output_path = _resolved(output_path)
    db_path = _resolved(db_path)
    if input_path == output_path:
        raise ManifestEnrichmentError(
            "Input and output manifest paths must differ."
        )
    dialect = _read_dialect(input_path)
    headers = _headers(input_path, dialect)
    schema = detect_schema(headers)
    lookup = _fragment_lookup(db_path)
    output_fields = _output_fieldnames(headers, schema)
    report: dict[str, Any] = {
        "input": str(input_path),
        "output": str(output_path),
        "database": str(db_path),
        "schema": schema,
        "strict": strict,
        "status": "success",
        "counts": {
            "input_rows": 0,
            "output_rows": 0,
            "fragment_rows": 0,
            "fragment_names_enriched": 0,
            "fragment_names_preserved": 0,
            "unmatched_fragment_lookups": 0,
            "ambiguous_fragment_lookups": 0,
        },
        "unmatched_fragment_lookups": [],
        "ambiguous_fragment_lookups": [],
        "unmatched_fragment_lookups_truncated": 0,
        "ambiguous_fragment_lookups_truncated": 0,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
        )
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(
                destination, fieldnames=output_fields, extrasaction="raise",
                delimiter=dialect.delimiter,
            )
            writer.writeheader()
            assigner = _EntityAssigner()
            with input_path.open(
                "r", encoding="utf-8-sig", newline=""
            ) as source:
                reader = csv.DictReader(source, dialect=dialect)
                for row_number, raw in enumerate(reader, start=2):
                    report["counts"]["input_rows"] += 1
                    row_iter = (
                        _wide_rows(raw, row_number)
                        if schema == "wide"
                        else (_base_row(raw, schema, row_number),)
                    )
                    for base, extras in row_iter:
                        component = _component_key(
                            base["uniprot_ac"],
                            base["is_fragment"] == "true",
                            base["sequence_start"], base["sequence_end"],
                        )
                        base["entity_id"] = assigner.entity_id(
                            base["model_entity_id"], component
                        )
                        _enrich_name(base, lookup, report)
                        writer.writerow({**base, **extras})
                        report["counts"]["output_rows"] += 1
        missing = (
            report["counts"]["unmatched_fragment_lookups"]
            + report["counts"]["ambiguous_fragment_lookups"]
        )
        if strict and missing:
            report["status"] = "strict_failure"
            Path(temporary_name).unlink()
            temporary_name = None
        else:
            Path(temporary_name).replace(output_path)
            temporary_name = None
    except OSError as exc:
        raise ManifestEnrichmentError(
            f"Cannot write output manifest {output_path}: {exc}"
        ) from exc
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return report


def write_report(path: Path, report: Mapping[str, Any]) -> None:
    """Write a structured enrichment report."""
    path = _resolved(path)
    try:
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    except OSError as exc:
        raise ManifestEnrichmentError(
            f"Cannot write report {path}: {exc}"
        ) from exc


def print_human_report(report: Mapping[str, Any]) -> None:
    """Print the concise human report required for every invocation."""
    counts = report["counts"]
    print(
        f"{report['status']} {report['schema']} manifest: "
        f"{counts['input_rows']} input rows -> "
        f"{counts['output_rows']} chain rows"
    )
    print(
        "Fragment names: "
        f"{counts['fragment_names_enriched']} enriched, "
        f"{counts['fragment_names_preserved']} preserved, "
        f"{counts['unmatched_fragment_lookups']} unmatched, "
        f"{counts['ambiguous_fragment_lookups']} ambiguous"
    )
    if report["status"] == "strict_failure":
        print(
            "Strict mode left the destination unchanged due to missing names."
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface."""
    args = parse_args(argv)
    try:
        report = enrich_manifest(
            args.input, args.output, args.db, strict=args.strict
        )
        if args.report:
            write_report(args.report, report)
        print_human_report(report)
        return 1 if report["status"] == "strict_failure" else 0
    except ManifestEnrichmentError as exc:
        print(f"error: {exc}")
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
