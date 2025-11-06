#!/usr/bin/env python3
"""
Combine per-accession AF metadata JSON files into chunked batches.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine single-entry AF metadata JSON files into batched outputs."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing per-accession JSON files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where combined metadata files will be written.",
    )
    parser.add_argument(
        "--output-prefix",
        default="AF-metadata",
        help="Prefix for combined JSON filenames (default: AF-metadata).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Number of records per combined JSON (default: 10,000).",
    )
    parser.add_argument(
        "--pattern",
        default="*.json",
        help="Glob pattern for per-accession files inside input directory.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s %(message)s",
    )


def iter_input_files(directory: Path, pattern: str) -> List[Path]:
    files = sorted(directory.glob(pattern))
    filtered = [path for path in files if path.is_file()]
    if not filtered:
        raise FileNotFoundError(f"No files matching pattern {pattern!r} in {directory}")
    return filtered


def chunked(iterable: Iterable[Path], size: int) -> Iterable[List[Path]]:
    batch: List[Path] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def load_records(paths: List[Path]) -> List[dict]:
    records: List[dict] = []
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        records.append(data)
    return records


def write_batch(
    output_dir: Path,
    prefix: str,
    index: int,
    total: int,
    records: List[dict],
) -> Path:
    filename = f"{prefix}-{index}-of-{total}.json"
    path = output_dir / filename
    with path.open("w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return path


def combine_metadata(args: argparse.Namespace) -> int:
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.debug("Scanning input directory %s", input_dir)
    files = iter_input_files(input_dir, args.pattern)
    logging.info("Found %d per-accession JSON files.", len(files))

    total_batches = (len(files) + args.chunk_size - 1) // args.chunk_size
    logging.info(
        "Combining into %d batch%s (chunk size %d).",
        total_batches,
        "" if total_batches == 1 else "es",
        args.chunk_size,
    )

    for idx, group in enumerate(chunked(files, args.chunk_size), start=1):
        logging.debug("Processing batch %d containing %d files.", idx, len(group))
        records = load_records(group)
        output_path = write_batch(output_dir, args.output_prefix, idx, total_batches, records)
        logging.info(
            "Wrote %d records to %s",
            len(records),
            output_path,
        )
    return 0


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)
    try:
        return combine_metadata(args)
    except Exception as exc:  # pragma: no cover - CLI guard
        logging.error("Combining metadata failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
