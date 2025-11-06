#!/usr/bin/env python3
"""
Stream a UniProt flat-file release, keep a target subset, and materialize Parquet tables.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

AC_PREFIX = "AC"
DE_PREFIX = "DE"
GN_PREFIX = "GN"
ID_PREFIX = "ID"
DT_PREFIX = "DT"
OS_PREFIX = "OS"
OX_PREFIX = "OX"
SQ_PREFIX = "SQ"
KW_PREFIX = "KW"

RE_DE_FIELD = re.compile(r"(Full|Short)=([^;]+)")
RE_TAXID = re.compile(r"NCBI_TaxID=(\d+)")
RE_SQ_LEN = re.compile(r"(\d+)\s+AA;")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter UniProt flat-files by accession list and emit Parquet tables."
    )
    parser.add_argument(
        "-t",
        "--targets",
        type=Path,
        help="Path to targets.txt (one accession per line).",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        type=Path,
        help="CSV/TSV mapping file that contains a 'uniprot_ac' column to select accessions.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for parquet files.",
    )
    parser.add_argument(
        "-r",
        "--release",
        required=True,
        help="Release tag stored with each entry row (e.g., 2025_03).",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more UniProt *.dat.gz files to stream.",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=1000,
        help="Number of rows to buffer before flushing to Parquet.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    return parser.parse_args()


def load_targets(path: Path) -> set[str]:
    targets: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            token = line.strip()
            if token:
                targets.add(token)
    return targets


def load_targets_from_mapping(path: Path) -> set[str]:
    import csv

    delimiter = "," if path.suffix.lower() != ".tsv" else "\t"
    targets: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None or "uniprot_ac" not in [name.strip() for name in reader.fieldnames]:
            raise ValueError(
                f"Mapping file {path} must contain a 'uniprot_ac' column."
            )
        for row in reader:
            accession = (row.get("uniprot_ac") or "").strip()
            if accession:
                targets.add(accession)
    return targets


def iter_records(path: Path) -> Iterable[List[str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        current: List[str] = []
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "//":
                if current:
                    yield current
                current = []
                continue
            current.append(line)
        if current:
            logging.warning("File %s ended without record terminator.", path)


def extract_accessions(lines: Sequence[str]) -> List[str]:
    accessions: List[str] = []
    for line in lines:
        if not line.startswith(AC_PREFIX):
            continue
        content = line[5:].strip()
        if not content:
            continue
        tokens = [token.strip() for token in content.split(";") if token.strip()]
        accessions.extend(tokens)
    return accessions


def strip_annotations(value: str) -> str:
    cleaned = value.strip()
    cleaned = re.split(r"\s*[\{\[]", cleaned, 1)[0].strip()
    return cleaned.rstrip(";").strip()


def parse_de_sections(lines: Sequence[str]) -> tuple[List[str], List[str]]:
    full_names: List[str] = []
    short_names: List[str] = []
    current_section: Optional[str] = None
    recname_captured_full = False
    for line in lines:
        if not line.startswith(DE_PREFIX):
            continue
        content = line[5:].rstrip()
        if not content:
            continue
        content = content.lstrip()
        section_match = re.match(r"(RecName|AltName)\s*:\s*(.*)", content)
        if section_match:
            current_section = section_match.group(1)
            content = section_match.group(2).strip()
        else:
            other_match = re.match(r"[A-Za-z]+Name\s*:\s*(.*)", content)
            if other_match:
                current_section = None
                content = other_match.group(1).strip()
        for match in RE_DE_FIELD.finditer(content):
            kind = match.group(1)
            raw_val = match.group(2)
            value = strip_annotations(raw_val)
            if not value:
                continue
            if current_section == "RecName":
                if kind == "Full":
                    if not recname_captured_full:
                        recname_captured_full = True
                    full_names.append(value)
                elif kind == "Short":
                    short_names.append(value)
            elif current_section == "AltName":
                if kind == "Full":
                    full_names.append(value)
                elif kind == "Short":
                    short_names.append(value)
    return full_names, short_names


def parse_organism(raw: str) -> tuple[Optional[str], List[str]]:
    value = raw.strip()
    if not value:
        return None, []
    value = value.rstrip(".")
    common_names = [match.strip().rstrip(".") for match in re.findall(r"\(([^)]+)\)", value)]
    if "(" in value:
        main = value.split("(", 1)[0].strip()
    else:
        main = value.strip()
    return main or None, [name for name in common_names if name]


class ParquetBufferWriter:
    def __init__(self, path: Path, schema: pa.Schema, batch_size: int) -> None:
        self.path = path
        self.schema = schema
        self.batch_size = batch_size
        self.buffer: List[Dict[str, object]] = []
        self.writer: Optional[pq.ParquetWriter] = None

    def append(self, row: Dict[str, object]) -> None:
        self.buffer.append(row)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.buffer:
            return
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)
        if self.writer is None:
            self.writer = pq.ParquetWriter(
                self.path, self.schema, compression="zstd"
            )
        self.writer.write_table(table)
        self.buffer.clear()

    def close(self) -> None:
        self.flush()
        if self.writer is not None:
            self.writer.close()
        else:
            empty_table = pa.Table.from_pylist([], schema=self.schema)
            pq.write_table(empty_table, self.path, compression="zstd")


def build_entry_payload(
    lines: Sequence[str],
    accessions: Sequence[str],
    release: str,
    reviewed: bool,
) -> Dict[str, object]:
    entry_name: Optional[str] = None
    protein_full_names, protein_short_names = parse_de_sections(lines)
    gene_name: Optional[str] = None
    gene_synonyms: List[str] = []
    organism: Optional[str] = None
    organism_common_names: List[str] = []
    taxid: Optional[int] = None
    length: Optional[int] = None
    sequence_version_date: Optional[str] = None
    is_reference_proteome = False
    seq_chunks: List[str] = []
    in_sequence = False

    for line in lines:
        prefix = line[:2]
        if prefix == ID_PREFIX and entry_name is None:
            remainder = line[5:].strip()
            if remainder:
                entry_name = remainder.split()[0]
        elif prefix == GN_PREFIX:
            content = line[5:].strip()
            if not content:
                continue
            tokens = [token.strip() for token in content.split(";") if token.strip()]
            for token in tokens:
                if token.startswith("Name=") and gene_name is None:
                    gene_name_candidate = strip_annotations(token.split("=", 1)[1])
                    if gene_name_candidate:
                        gene_name = gene_name_candidate
                elif token.startswith("Synonyms="):
                    synonym_blob = token.split("=", 1)[1]
                    for raw_synonym in synonym_blob.split(","):
                        synonym = strip_annotations(raw_synonym)
                        if synonym and synonym not in gene_synonyms:
                            gene_synonyms.append(synonym)
        elif prefix == OS_PREFIX and organism is None:
            org_value, common_names = parse_organism(line[5:])
            organism = org_value
            organism_common_names = common_names
        elif prefix == OX_PREFIX and taxid is None:
            match = RE_TAXID.search(line)
            if match:
                taxid = int(match.group(1))
        elif prefix == DT_PREFIX and sequence_version_date is None and "sequence version" in line:
            date_token = line[5:].split(",", 1)[0].strip()
            try:
                sequence_version_date = datetime.strptime(date_token, "%d-%b-%Y").date().isoformat()
            except ValueError:
                sequence_version_date = date_token
        elif prefix == KW_PREFIX and not is_reference_proteome:
            if "Reference proteome" in line:
                is_reference_proteome = True
        elif prefix == SQ_PREFIX:
            in_sequence = True
            if length is None:
                len_match = RE_SQ_LEN.search(line)
                if len_match:
                    length = int(len_match.group(1))
            continue
        elif in_sequence:
            aa_only = "".join(ch for ch in line if ch.isalpha())
            if aa_only:
                seq_chunks.append(aa_only)

    seq_str = "".join(seq_chunks) if seq_chunks else None
    seq_md5 = hashlib.md5(seq_str.encode("utf-8")).hexdigest() if seq_str else None
    primary_ac = accessions[0]
    entry_row: Dict[str, object] = {
        "primary_ac": primary_ac,
        "entry_name": entry_name,
        "reviewed": reviewed,
        "protein_full_names": protein_full_names or None,
        "protein_short_names": protein_short_names or None,
        "gene_names": gene_name,
        "gene_synonyms": gene_synonyms or None,
        "organism": organism,
        "organisme_common_names": organism_common_names or None,
        "taxid": taxid,
        "length": length,
        "sequence_version_date": sequence_version_date,
        "is_uniprot_reference_proteome": is_reference_proteome,
        "md5": seq_md5,
        "sequence": seq_str,
        "release": release,
    }
    return entry_row


def process_inputs(
    inputs: Sequence[Path],
    targets: set[str],
    release: str,
    entry_writer: ParquetBufferWriter,
) -> None:
    total = kept = 0
    remaining = set(targets)
    total_targets = len(remaining)
    pbar = tqdm(total=total_targets, desc="Matched targets", unit="acc", leave=False) if total_targets else None
    try:
        for input_path in inputs:
            reviewed = "sprot" in input_path.name.lower()
            logging.info("Scanning %s (reviewed=%s)", input_path, reviewed)
            for record_lines in iter_records(input_path):
                total += 1
                accessions = extract_accessions(record_lines)
                if not accessions:
                    continue
                matching = [ac for ac in accessions if ac in remaining]
                if not matching:
                    continue
                kept += 1
                entry_row = build_entry_payload(record_lines, accessions, release, reviewed)
                entry_writer.append(entry_row)
                for ac in matching:
                    remaining.discard(ac)
                if pbar:
                    pbar.update(len(matching))
                if not remaining:
                    logging.info("All target accessions processed; stopping early.")
                    break
            if not remaining:
                break
    finally:
        if pbar:
            pbar.close()
    logging.info("Processed %d records, kept %d matching targets.", total, kept)


def ensure_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def get_entry_schema() -> pa.Schema:
    return pa.schema(
        [
            ("primary_ac", pa.string()),
            ("entry_name", pa.string()),
            ("reviewed", pa.bool_()),
            ("protein_full_names", pa.list_(pa.string())),
            ("protein_short_names", pa.list_(pa.string())),
            ("gene_names", pa.string()),
            ("gene_synonyms", pa.list_(pa.string())),
            ("organism", pa.string()),
            ("organisme_common_names", pa.list_(pa.string())),
            ("taxid", pa.int64()),
            ("length", pa.int32()),
            ("sequence_version_date", pa.string()),
            ("is_uniprot_reference_proteome", pa.bool_()),
            ("md5", pa.string()),
            ("sequence", pa.string()),
            ("release", pa.string()),
        ]
    )


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(message)s",
    )
    targets: set[str] = set()
    if args.targets:
        targets |= load_targets(args.targets)
    if args.mapping:
        targets |= load_targets_from_mapping(args.mapping)

    if not targets:
        logging.error("No accessions provided. Use --targets and/or --mapping.")
        return 1

    logging.info("Loaded %d unique target accessions.", len(targets))
    ensure_output_dir(args.outdir)
    entry_schema = get_entry_schema()

    entry_writer = ParquetBufferWriter(args.outdir / "entry.parquet", entry_schema, args.batch_size)

    try:
        process_inputs(
            args.inputs,
            targets,
            args.release,
            entry_writer,
        )
    finally:
        entry_writer.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
