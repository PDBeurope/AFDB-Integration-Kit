#!/usr/bin/env python3
"""
Shard UniProt flat-files by accession for faster parallel extraction.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import sys
from pathlib import Path
from typing import List, Sequence

from tqdm import tqdm

# Allow running as a standalone script without installing the package.
try:
    from uniprot.scripts.extract_subset import extract_accessions, iter_records
except ModuleNotFoundError:
    repo_root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(repo_root))
    from uniprot.scripts.extract_subset import extract_accessions, iter_records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split UniProt flat-files into stable shards keyed by accession."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Root directory for shard outputs.",
    )
    parser.add_argument(
        "-r",
        "--release",
        required=True,
        help="Release tag used to nest shard outputs (e.g., 2025_03).",
    )
    parser.add_argument(
        "--shard-count",
        type=int,
        default=8,
        help="Number of shards to produce. Defaults to 8.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=5,
        help="Compression level for shard outputs (1-9). Defaults to 5.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more UniProt *.dat.gz files to shard.",
    )
    return parser.parse_args()


def stable_shard_for_accession(accession: str, shard_count: int) -> int:
    digest = hashlib.md5(accession.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % shard_count


def detect_source_label(path: Path) -> str:
    name = path.name.lower()
    if "sprot" in name:
        return "sprot"
    if "trembl" in name:
        return "trembl"
    return Path(path.stem).stem.lower()


def open_shard_writers(
    shard_dir: Path,
    source_label: str,
    shard_count: int,
    gzip_level: int,
) -> List[gzip.GzipFile]:
    shard_dir.mkdir(parents=True, exist_ok=True)
    writers: List[gzip.GzipFile] = []
    for idx in range(shard_count):
        shard_path = shard_dir / f"{source_label}-shard-{idx:02d}.dat.gz"
        writer = gzip.open(shard_path, mode="wt", encoding="utf-8", compresslevel=gzip_level)
        writers.append(writer)
    return writers


def write_record(handle: gzip.GzipFile, lines: Sequence[str]) -> None:
    for line in lines:
        handle.write(line)
        handle.write("\n")
    handle.write("//\n")


def shard_file(
    input_path: Path,
    outdir: Path,
    release: str,
    shard_count: int,
    gzip_level: int,
) -> None:
    source_label = detect_source_label(input_path)
    shard_dir = outdir / release / source_label
    logging.info(
        "Sharding %s into %d parts under %s", input_path, shard_count, shard_dir
    )
    writers = open_shard_writers(shard_dir, source_label, shard_count, gzip_level)
    total_records = 0
    try:
        for record_lines in tqdm(iter_records(input_path), desc=f"{source_label} records", unit="rec", leave=False):
            total_records += 1
            accessions = extract_accessions(record_lines)
            if not accessions:
                continue
            shard_idx = stable_shard_for_accession(accessions[0], shard_count)
            write_record(writers[shard_idx], record_lines)
    finally:
        for writer in writers:
            writer.close()
    logging.info("Finished %s (%d records).", input_path, total_records)


def main() -> int:
    args = parse_args()
    if not 1 <= args.gzip_level <= 9:
        raise SystemExit("--gzip-level must be between 1 and 9.")
    if args.shard_count < 1:
        raise SystemExit("--shard-count must be at least 1.")

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(message)s",
    )
    for input_path in args.inputs:
        shard_file(
            input_path=input_path,
            outdir=args.outdir,
            release=args.release,
            shard_count=args.shard_count,
            gzip_level=args.gzip_level,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
