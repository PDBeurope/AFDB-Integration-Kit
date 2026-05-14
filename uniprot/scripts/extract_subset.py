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
RE_ISOFORM_ID = re.compile(r"IsoId=([^;]+)")
RE_ISOFORM_SEQUENCE = re.compile(r"Sequence=([^;]+)")


def stable_shard_for_accession(accession: str, shard_count: int) -> int:
    digest = hashlib.md5(accession.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % shard_count


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
    parser.add_argument(
        "--shard-count",
        type=int,
        help="Total shard count used to partition accessions (pairs with --shard-index).",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        help="Shard index (0-based) to process when running over sharded inputs.",
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


def parse_alt_products(lines: Sequence[str]) -> Dict[str, List[str]]:
    """
    Extract isoform definitions from the CC ALTERNATIVE PRODUCTS block.

    Returns a mapping of isoform accession -> list of VSP ids (or ["Displayed"]).
    """

    isoforms: Dict[str, List[str]] = {}
    in_block = False
    for raw in lines:
        if not raw.startswith("CC"):
            continue
        content = raw[5:].strip()
        if content.startswith("-!- ALTERNATIVE PRODUCTS:"):
            in_block = True
            continue
        if in_block and content.startswith("-!- "):
            # next CC subsection begins; stop parsing alternative products
            break
        if not in_block:
            continue
        iso_match = RE_ISOFORM_ID.search(content)
        seq_match = RE_ISOFORM_SEQUENCE.search(content)
        if not iso_match or not seq_match:
            continue
        iso_id = iso_match.group(1).strip()
        seq_field = seq_match.group(1).strip()
        if not iso_id or not seq_field:
            continue
        tokens = [token.strip() for token in seq_field.split(",") if token.strip()]
        isoforms[iso_id] = tokens or []
    return isoforms


def parse_var_seq(lines: Sequence[str]) -> Dict[str, tuple[int, int, str]]:
    """
    Parse FT VAR_SEQ features into a mapping of VSP id -> (start, end, replacement).
    """

    varseqs: Dict[str, tuple[int, int, str]] = {}
    total_lines = len(lines)
    idx = 0
    while idx < total_lines:
        line = lines[idx]
        if not line.startswith("FT"):
            idx += 1
            continue
        content = line[5:].strip()
        if not content.startswith("VAR_SEQ"):
            idx += 1
            continue
        parts = content.split()
        if len(parts) < 2:
            idx += 1
            continue
        range_token = parts[1]
        try:
            start_str, end_str = range_token.split("..")
            start, end = int(start_str), int(end_str)
        except ValueError:
            idx += 1
            continue

        replacement = ""
        vsp_id: Optional[str] = None
        idx += 1
        while idx < total_lines:
            cont = lines[idx]
            if not cont.startswith("FT"):
                break
            cont_content = cont[5:].strip()
            # New feature starts when column 6 is not space (continuation lines have a space here)
            if cont_content and not cont[5].isspace():
                break
            if cont_content.startswith("/note="):
                note = cont_content[len("/note="):].strip().strip('"')
                if "->" in note:
                    replacement = re.split(r"\s*\(", note.split("->", 1)[1].strip(), 1)[0].strip()
                elif note.lower().startswith("missing"):
                    replacement = ""
                else:
                    replacement = re.split(r"\s*\(", note, 1)[0].strip()
            elif cont_content.startswith("/id="):
                vsp_id = cont_content[len("/id="):].strip().strip('"')
            idx += 1
        if vsp_id:
            varseqs[vsp_id] = (start, end, replacement)
    return varseqs


def apply_var_seq_edits(sequence: str, edits: List[tuple[int, int, str]]) -> str:
    """
    Apply a list of VAR_SEQ edits (start, end, replacement) to a sequence.

    Edits are applied from highest start coordinate to lowest to avoid offset shifts.
    """

    if not edits:
        return sequence
    patched = sequence
    for start, end, replacement in sorted(edits, key=lambda item: item[0], reverse=True):
        # UniProt coordinates are 1-based inclusive
        left = patched[: start - 1]
        right = patched[end:]
        patched = f"{left}{replacement}{right}"
    return patched


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


def parse_organism(raw: str) -> tuple[Optional[str], List[str], List[str]]:
    """Parse OS lines into scientific name, common names, and synonyms."""

    value = raw.strip()
    if not value:
        return None, [], []
    value = value.rstrip(".")
    paren_values = [match.strip().rstrip(".") for match in re.findall(r"\(([^)]+)\)", value)]
    main = value.split("(", 1)[0].strip() if "(" in value else value.strip()
    if paren_values:
        common_names = [paren_values[0]] if paren_values[0] else []
        synonyms = [name for name in paren_values[1:] if name]
    else:
        common_names = []
        synonyms = []
    return main or None, common_names, synonyms


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
) -> List[Dict[str, object]]:
    entry_name: Optional[str] = None
    protein_full_names, protein_short_names = parse_de_sections(lines)
    gene_name: Optional[str] = None
    gene_synonyms: List[str] = []
    gene_ordered_locus_names: List[str] = []
    gene_orf_names: List[str] = []
    organism: Optional[str] = None
    organism_common_names: List[str] = []
    organism_synonyms: List[str] = []
    taxid: Optional[int] = None
    length: Optional[int] = None
    sequence_version_date: Optional[str] = None
    is_reference_proteome = False
    seq_chunks: List[str] = []
    in_sequence = False
    gn_lines: List[str] = []
    os_lines: List[str] = []

    def flush_gene_block(
        name: Optional[str],
        synonyms: List[str],
        locus_names: List[str],
        orf_names: List[str],
    ) -> tuple[Optional[str], List[str], List[str], List[str]]:
        primary = name or (synonyms[0] if synonyms else None) or (
            locus_names[0] if locus_names else (orf_names[0] if orf_names else None)
        )
        return primary, synonyms, locus_names, orf_names

    for line in lines:
        prefix = line[:2]
        if prefix == ID_PREFIX and entry_name is None:
            remainder = line[5:].strip()
            if remainder:
                entry_name = remainder.split()[0]
        elif prefix == GN_PREFIX:
            gn_lines.append(line[5:].strip())
        elif prefix == OS_PREFIX:
            os_lines.append(line[5:].strip())
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

    if gn_lines:
        current_name: Optional[str] = None
        current_synonyms: List[str] = []
        current_locus: List[str] = []
        current_orf: List[str] = []
        tokens = []
        for raw in gn_lines:
            raw_clean = raw.rstrip(".")
            tokens.extend([token.strip() for token in raw_clean.split(";") if token.strip()])
        for token in tokens:
            if token.lower() == "and":
                gene_name, gene_synonyms, gene_ordered_locus_names, gene_orf_names = flush_gene_block(
                    current_name, current_synonyms, current_locus, current_orf
                )
                current_name, current_synonyms, current_locus, current_orf = None, [], [], []
                continue
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            key = key.strip()
            values = [strip_annotations(val) for val in value.split(",")]
            values = [val for val in values if val]
            if key == "Name":
                if current_name is None and values:
                    current_name = values[0]
            elif key == "Synonyms":
                for synonym in values:
                    if synonym and synonym not in current_synonyms:
                        current_synonyms.append(synonym)
            elif key == "OrderedLocusNames":
                for locus in values:
                    if locus and locus not in current_locus:
                        current_locus.append(locus)
            elif key == "ORFNames":
                for orf in values:
                    if orf and orf not in current_orf:
                        current_orf.append(orf)
        if gene_name is None:
            gene_name, gene_synonyms, gene_ordered_locus_names, gene_orf_names = flush_gene_block(
                current_name, current_synonyms, current_locus, current_orf
            )

    if os_lines:
        organism_value, common_names, synonyms = parse_organism(" ".join(os_lines))
        organism = organism or organism_value
        organism_common_names = common_names
        organism_synonyms = synonyms

    seq_str = "".join(seq_chunks) if seq_chunks else None
    seq_md5 = hashlib.md5(seq_str.encode("utf-8")).hexdigest() if seq_str else None
    primary_ac = accessions[0]
    base_row: Dict[str, object] = {
        "primary_ac": primary_ac,
        "entry_name": entry_name,
        "reviewed": reviewed,
        "protein_full_names": protein_full_names or None,
        "protein_short_names": protein_short_names or None,
        "gene_names": gene_name,
        "gene_synonyms": gene_synonyms or None,
        "gene_ordered_locus_names": gene_ordered_locus_names or None,
        "gene_orf_names": gene_orf_names or None,
        "organism": organism,
        "organism_common_names": organism_common_names or None,
        "organism_synonyms": organism_synonyms or None,
        "taxid": taxid,
        "length": length,
        "sequence_version_date": sequence_version_date,
        "is_uniprot_reference_proteome": is_reference_proteome,
        "md5": seq_md5,
        "sequence": seq_str,
        "release": release,
        "is_isoform": False,
    }

    rows: List[Dict[str, object]] = [base_row]

    if not seq_str:
        return rows

    isoform_map = parse_alt_products(lines)
    varseqs = parse_var_seq(lines)

    for isoform_id, tokens in isoform_map.items():
        if isoform_id == primary_ac:
            continue
        if not tokens:
            continue
        if len(tokens) == 1 and tokens[0].lower() == "displayed":
            iso_seq = seq_str
        else:
            edits: List[tuple[int, int, str]] = []
            missing_vsp: List[str] = []
            for token in tokens:
                vsp = token.strip()
                if vsp not in varseqs:
                    missing_vsp.append(vsp)
                    continue
                edits.append(varseqs[vsp])
            if missing_vsp:
                logging.debug(
                    "Skipping isoform %s due to missing VSP ids: %s",
                    isoform_id,
                    ", ".join(missing_vsp),
                )
                continue
            iso_seq = apply_var_seq_edits(seq_str, edits)
        iso_row = dict(base_row)
        iso_row["primary_ac"] = isoform_id
        iso_row["sequence"] = iso_seq
        iso_row["length"] = len(iso_seq)
        iso_row["md5"] = hashlib.md5(iso_seq.encode("utf-8")).hexdigest()
        iso_row["is_isoform"] = True
        rows.append(iso_row)

    return rows


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
                isoform_ids = list(parse_alt_products(record_lines).keys())
                candidate_ids = accessions + isoform_ids
                matching = [ac for ac in candidate_ids if ac in remaining]
                if not matching:
                    continue
                kept += 1
                entry_rows = build_entry_payload(record_lines, accessions, release, reviewed)
                for row in entry_rows:
                    entry_writer.append(row)
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
            ("gene_ordered_locus_names", pa.list_(pa.string())),
            ("gene_orf_names", pa.list_(pa.string())),
            ("organism", pa.string()),
            ("organism_common_names", pa.list_(pa.string())),
            ("organism_synonyms", pa.list_(pa.string())),
            ("taxid", pa.int64()),
            ("length", pa.int32()),
            ("sequence_version_date", pa.string()),
            ("is_uniprot_reference_proteome", pa.bool_()),
            ("md5", pa.string()),
            ("sequence", pa.string()),
            ("release", pa.string()),
            ("is_isoform", pa.bool_()),
        ]
    )


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(message)s",
    )
    if (args.shard_count is None) ^ (args.shard_index is None):
        logging.error("Provide both --shard-count and --shard-index to enable shard filtering.")
        return 1
    if args.shard_count is not None and args.shard_count < 1:
        logging.error("--shard-count must be at least 1.")
        return 1
    if args.shard_index is not None and args.shard_index < 0:
        logging.error("--shard-index must be non-negative.")
        return 1

    targets: set[str] = set()
    if args.targets:
        targets |= load_targets(args.targets)
    if args.mapping:
        targets |= load_targets_from_mapping(args.mapping)
    initial_target_count = len(targets)
    if initial_target_count == 0:
        logging.error("No accessions provided. Use --targets and/or --mapping.")
        return 1

    if args.shard_count is not None and args.shard_index is not None:
        if args.shard_index >= args.shard_count:
            logging.error("--shard-index must be less than --shard-count.")
            return 1
        targets = {
            ac for ac in targets if stable_shard_for_accession(ac, args.shard_count) == args.shard_index
        }
        logging.info(
            "Filtered targets to shard %d of %d (%d of %d accessions).",
            args.shard_index,
            args.shard_count,
            len(targets),
            initial_target_count,
        )
        if not targets:
            ensure_output_dir(args.outdir)
            entry_writer = ParquetBufferWriter(args.outdir / "entry.parquet", get_entry_schema(), args.batch_size)
            entry_writer.close()
            logging.info(
                "No targets fall into shard %d of %d; wrote empty parquet and exiting.",
                args.shard_index,
                args.shard_count,
            )
            return 0

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
