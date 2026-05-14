#!/usr/bin/env python3
"""
Fetch UniProt entries via REST API and build a DuckDB database.
Useful for creating small test databases without needing full UniProt flat files.

This script delegates to :mod:`afdb_integration_kit.uniprot.api` for all
fetching, parsing, and database-building logic.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Ensure the repo root is importable when running as a standalone script
_repo_root = Path(__file__).resolve().parent.parent.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from afdb_integration_kit.uniprot.api import (  # noqa: E402
    ENTRY_SCHEMA,
    build_duckdb,
    entries_to_parquet,
    fetch_entries,
    parse_entry,
)

# Re-export for any callers that still import from this module
get_entry_schema = lambda: ENTRY_SCHEMA  # noqa: E731
fetch_uniprot_entries = fetch_entries
parse_uniprot_entry = parse_entry

__all__ = [
    "get_entry_schema",
    "fetch_uniprot_entries",
    "parse_uniprot_entry",
    "entries_to_parquet",
    "build_duckdb",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch UniProt entries and build DuckDB",
    )
    parser.add_argument(
        "--accessions", required=True, type=Path,
        help="File with UniProt accessions, one per line",
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output DuckDB file path",
    )
    parser.add_argument(
        "--batch-size", type=int, default=25,
        help="Accessions per API request",
    )
    parser.add_argument(
        "--release", default="2025_01",
        help="Release tag to store with entries",
    )
    args = parser.parse_args()

    accessions = [
        line.strip()
        for line in args.accessions.read_text().splitlines()
        if line.strip()
    ]
    logger.info("Loaded %d accessions from %s", len(accessions), args.accessions)

    entries = fetch_entries(accessions, args.batch_size, release=args.release)
    logger.info("Fetched %d entries from UniProt API", len(entries))

    if not entries:
        logger.error("No entries fetched!")
        return 1

    parquet_path = args.output.parent / "entry.parquet"
    entries_to_parquet(entries, parquet_path)
    build_duckdb(parquet_path, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
