#!/usr/bin/env python3
"""
Materialise UniProt subset Parquet outputs into a DuckDB database with helpful indexes.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import duckdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a DuckDB database from UniProt subset Parquet tables."
    )
    parser.add_argument(
        "--parquet-dir",
        required=True,
        type=Path,
        help="Directory containing entry.parquet.",
    )
    parser.add_argument(
        "--db",
        required=True,
        type=Path,
        help="Path to the DuckDB database file to create.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing DuckDB file if present.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print build summaries and sample rows.",
    )
    parser.add_argument(
        "--model-metadata",
        type=Path,
        help="Optional path to model metadata (CSV, Parquet, or DuckDB) to load into the database.",
    )
    parser.add_argument(
        "--model-table",
        help="When --model-metadata is a DuckDB file, the table name to import (defaults to the first table).",
    )
    return parser.parse_args()


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s %(message)s")


def validate_inputs(parquet_dir: Path) -> Path:
    entry_path = parquet_dir / "entry.parquet"
    if not entry_path.exists():
        raise FileNotFoundError(
            f"Required Parquet file missing under {parquet_dir}: {entry_path}"
        )
    return entry_path


def prepare_db_path(db_path: Path, force: bool) -> None:
    if db_path.exists():
        if not force:
            raise FileExistsError(
                f"{db_path} already exists. Use --force to overwrite."
            )
        logging.info("Removing existing database %s", db_path)
        db_path.unlink()
    db_path.parent.mkdir(parents=True, exist_ok=True)


def load_model_metadata(
    con: duckdb.DuckDBPyConnection, model_path: Path, table_name: str | None
) -> None:
    logging.info("Loading model metadata from %s", model_path)
    suffix = model_path.suffix.lower()
    con.execute("DROP TABLE IF EXISTS model_metadata")
    if suffix == ".csv":
        con.execute(
            "CREATE TABLE model_metadata AS SELECT * FROM read_csv_auto(?, header=True)",
            [str(model_path)],
        )
    elif suffix in {".parquet", ".pq"}:
        con.execute(
            "CREATE TABLE model_metadata AS SELECT * FROM read_parquet(?)",
            [str(model_path)],
        )
    elif suffix == ".duckdb":
        logging.info("Attaching DuckDB database %s for metadata import.", model_path)
        con.execute("ATTACH DATABASE ? AS metadata_db", [str(model_path)])
        tables = [row[0] for row in con.execute("PRAGMA metadata_db.show_tables").fetchall()]
        if not tables:
            raise ValueError(f"No tables found in metadata database {model_path}")
        chosen_table = table_name or tables[0]
        if chosen_table not in tables:
            raise ValueError(
                f"Table {chosen_table!r} not found in metadata database {model_path}. "
                f"Available tables: {tables}"
            )
        logging.info("Importing table %s from metadata database.", chosen_table)
        con.execute(
            f"CREATE TABLE model_metadata AS SELECT * FROM metadata_db.\"{chosen_table}\""
        )
        con.execute("DETACH DATABASE metadata_db")
    else:
        raise ValueError(f"Unsupported model metadata format: {model_path}")
    logging.info(
        "Loaded %d model metadata rows.",
        con.execute("SELECT COUNT(*) FROM model_metadata").fetchone()[0],
    )


def create_tables_and_indexes(
    con: duckdb.DuckDBPyConnection,
    entry_path: Path,
) -> None:
    con.execute("DROP TABLE IF EXISTS entry")
    con.execute("DROP TABLE IF EXISTS sequence")
    con.execute("DROP TABLE IF EXISTS accession")

    logging.debug("Creating entry table from %s", entry_path)
    con.execute(
        "CREATE TABLE entry AS SELECT * FROM read_parquet(?)",
        [str(entry_path)],
    )

    logging.debug("Creating index idx_entry_primary on entry(primary_ac)")
    con.execute("CREATE INDEX idx_entry_primary ON entry(primary_ac)")
    logging.debug("Creating index idx_entry_release ON entry(release)")
    con.execute("CREATE INDEX idx_entry_release ON entry(release)")


def _table_exists(con: duckdb.DuckDBPyConnection, name: str) -> bool:
    tables = [row[0] for row in con.execute("PRAGMA show_tables").fetchall()]
    return name in tables


def run_sanity_checks(con: duckdb.DuckDBPyConnection) -> None:
    entry_count = con.execute("SELECT COUNT(*) FROM entry").fetchone()[0]
    logging.info("entry rows: %s", entry_count)

    sample = con.execute(
        "SELECT primary_ac, length(sequence) AS seqlen "
        "FROM entry "
        "LIMIT 3"
    ).fetchall()
    if sample:
        logging.info("Sample entries (primary_ac, seqlen): %s", sample)
    else:
        logging.info("No sample rows available (empty table).")

    if _table_exists(con, "model_metadata"):
        meta_count = con.execute("SELECT COUNT(*) FROM model_metadata").fetchone()[0]
        logging.info("model_metadata rows: %s", meta_count)


def build_database(args: argparse.Namespace) -> int:
    parquet_dir = args.parquet_dir.resolve()
    db_path = args.db.resolve()

    if args.model_table and not args.model_metadata:
        raise ValueError("--model-table requires --model-metadata pointing to a DuckDB file.")

    logging.debug("Using Parquet directory %s", parquet_dir)
    logging.debug("Target DuckDB path %s", db_path)

    entry_path = validate_inputs(parquet_dir)
    prepare_db_path(db_path, args.force)

    con = duckdb.connect(str(db_path))
    try:
        create_tables_and_indexes(
            con, entry_path
        )
        if args.model_metadata:
            load_model_metadata(con, args.model_metadata.resolve(), args.model_table)
        con.commit()

        if args.verbose:
            run_sanity_checks(con)

        logging.info("Build complete.")
    finally:
        con.close()
    return 0


def main() -> int:
    args = parse_args()
    configure_logging(args.verbose)
    try:
        return build_database(args)
    except Exception as exc:  # pragma: no cover - CLI guard
        logging.error("Build failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
