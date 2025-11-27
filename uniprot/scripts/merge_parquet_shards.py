#!/usr/bin/env python3
"""
Merge shard-level Parquet outputs into a single file.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge multiple Parquet shard files into one Parquet output."
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Path to the merged Parquet file.",
    )
    parser.add_argument(
        "--compression",
        default="zstd",
        help="Compression codec for the output Parquet file (default: zstd).",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Parquet files or glob patterns to merge (e.g., uniprot/outputs/parquet/sprot-shard-*/entry.parquet).",
    )
    return parser.parse_args()


def merge_parquet(inputs: List[str], output: Path, compression: str) -> None:
    dataset = ds.dataset(inputs, format="parquet")
    writer: pq.ParquetWriter | None = None
    try:
        for batch in dataset.to_batches():
            if writer is None:
                writer = pq.ParquetWriter(output, dataset.schema, compression=compression)
            writer.write_table(pa.Table.from_batches([batch]))
    finally:
        if writer is not None:
            writer.close()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    logging.info("Merging %d parquet inputs into %s", len(args.inputs), args.output)
    merge_parquet(args.inputs, args.output, args.compression)
    logging.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
