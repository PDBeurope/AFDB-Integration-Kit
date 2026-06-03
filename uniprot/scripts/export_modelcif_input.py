#!/usr/bin/env python3
"""
Generate ModelCIF generator input JSON by combining UniProt-derived data with a
base metadata template and a manifest describing chain/entity assignments.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os

import orjson
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import duckdb


LOG = logging.getLogger("uniprot.export_modelcif_input")


@dataclass(frozen=True)
class ManifestEntry:
    model_entity_id: str
    entity_id: str
    chain_id: str
    uniprot_ac: str
    sequence_start: int | None = None
    sequence_end: int | None = None


@dataclass
class EntityAssignment:
    entity_id: str
    uniprot_ac: str
    chain_ids: List[str]
    sequence_start: int | None = None
    sequence_end: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate ModelCIF generator metadata JSON using UniProt/DuckDB data "
            "and a manifest describing chain/entity assignments."
        )
    )
    parser.add_argument(
        "--model-id",
        required=True,
        help="Model entity identifier (e.g., AF-0000000000000004).",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="CSV manifest containing model_entity_id, entity_id, chain_id, uniprot_ac columns.",
    )
    parser.add_argument(
        "--db",
        required=True,
        type=Path,
        help="DuckDB database produced by build_duckdb.py.",
    )
    parser.add_argument(
        "--template",
        required=True,
        type=Path,
        help="JSON template containing static metadata/categories blocks.",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Destination JSON file for ModelCIF generator input.",
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


def load_manifest(path: Path, model_id: str) -> List[ManifestEntry]:
    if not path.exists():
        raise FileNotFoundError(f"Manifest file {path} does not exist.")
    entries: List[ManifestEntry] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = {"model_entity_id", "entity_id", "chain_id", "uniprot_ac"}
        if reader.fieldnames is None or expected - set(reader.fieldnames):
            raise ValueError(f"Manifest {path} must contain columns {expected}, found {reader.fieldnames}")
        for row in reader:
            if (row.get("model_entity_id") or "").strip() != model_id:
                continue
            entry = ManifestEntry(
                model_entity_id=row["model_entity_id"].strip(),
                entity_id=row["entity_id"].strip(),
                chain_id=row["chain_id"].strip(),
                uniprot_ac=(row.get("uniprot_ac") or "").strip(),
                sequence_start=int(row["sequence_start"]) if (row.get("sequence_start") or "").strip() else None,
                sequence_end=int(row["sequence_end"]) if (row.get("sequence_end") or "").strip() else None,
            )
            if not entry.entity_id or not entry.chain_id or not entry.uniprot_ac:
                raise ValueError(f"Incomplete manifest row for model {model_id}: {row}")
            entries.append(entry)
    if not entries:
        raise ValueError(f"No manifest entries found for model {model_id} in {path}.")
    return entries


def group_entities(entries: Sequence[ManifestEntry]) -> List[EntityAssignment]:
    grouped: "OrderedDict[str, EntityAssignment]" = OrderedDict()
    for entry in entries:
        key = entry.entity_id
        current = grouped.get(key)
        if current is None:
            grouped[key] = EntityAssignment(
                entity_id=entry.entity_id,
                uniprot_ac=entry.uniprot_ac,
                chain_ids=[entry.chain_id],
                sequence_start=entry.sequence_start,
                sequence_end=entry.sequence_end,
            )
        else:
            if entry.uniprot_ac != current.uniprot_ac:
                raise ValueError(
                    f"Entity {entry.entity_id} has conflicting UniProt accessions: "
                    f"{current.uniprot_ac} vs {entry.uniprot_ac}"
                )
            if entry.chain_id in current.chain_ids:
                raise ValueError(
                    f"Duplicate chain {entry.chain_id!r} listed for entity {entry.entity_id}."
                )
            if (
                entry.sequence_start != current.sequence_start
                or entry.sequence_end != current.sequence_end
            ):
                raise ValueError(
                    f"Entity {entry.entity_id} has inconsistent fragment ranges across chains."
                )
            current.chain_ids.append(entry.chain_id)
    return list(grouped.values())


_CRC64_TABLE: List[int] = []
_CRC64_POLY = 0xC96C5795D7870F42
for byte in range(256):
    crc = byte << 56
    for _ in range(8):
        if crc & (1 << 63):
            crc = (crc << 1) ^ _CRC64_POLY
        else:
            crc <<= 1
        crc &= 0xFFFFFFFFFFFFFFFF
    _CRC64_TABLE.append(crc)


def crc64_ecma(data: str) -> str:
    """Compute CRC64-ECMA hash, matching UniProt's checksum."""
    crc = 0xFFFFFFFFFFFFFFFF
    for ch in data.encode("utf-8"):
        idx = ((crc >> 56) ^ ch) & 0xFF
        crc = _CRC64_TABLE[idx] ^ ((crc << 8) & 0xFFFFFFFFFFFFFFFF)
    crc ^= 0xFFFFFFFFFFFFFFFF
    return f"{crc:016X}"


def fetch_entry(con: duckdb.DuckDBPyConnection, accession: str) -> Dict[str, object]:
    query = "SELECT * FROM entry WHERE primary_ac = ?"
    relation = con.execute(query, [accession])
    table = relation.fetch_arrow_table()
    if table.num_rows == 0:
        raise ValueError(f"Accession {accession} not found in entry table.")
    return table.to_pylist()[0]


def ensure_category(template: Dict[str, object], name: str) -> Dict[str, List[str]]:
    categories = template.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("Template 'categories' must be a dictionary.")
    category = categories.setdefault(name, {})
    if not isinstance(category, dict):
        raise TypeError(f"Category {name} must be a dictionary.")
    return category  # type: ignore[return-value]


def serialise_list(values: Iterable[object]) -> List[object]:
    return list(values)


def normalise_optional_text(value: object, placeholder: str = "?") -> str:
    """Return a string value, falling back to '?' when no data is available."""

    def _flatten(obj: object) -> List[str]:
        if obj is None:
            return []
        if isinstance(obj, str):
            stripped = obj.strip()
            return [stripped] if stripped else []
        if isinstance(obj, (list, tuple, set)):
            flattened: List[str] = []
            for item in obj:
                flattened.extend(_flatten(item))
            return flattened
        return [str(obj)]

    parts = _flatten(value)
    return ", ".join(parts) if parts else placeholder


def populate_categories(
    template: Dict[str, object],
    entities: Sequence[EntityAssignment],
    entries_by_entity: Dict[str, Dict[str, object]],
    model_id: str,
) -> None:
    categories = template.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("Template must contain a 'categories' dictionary.")

    # _ma_target_ref_db_details
    target_ref = ensure_category(template, "_ma_target_ref_db_details")
    ref_fields: Dict[str, List[object]] = defaultdict(list)
    for ordinal, assignment in enumerate(entities, start=1):
        entry = entries_by_entity[assignment.entity_id]
        sequence: str = entry.get("sequence") or ""
        if not sequence:
            raise ValueError(f"No sequence stored for accession {assignment.uniprot_ac}.")
        seq_start = assignment.sequence_start or 1
        seq_end = assignment.sequence_end or len(sequence)
        if seq_start < 1 or seq_end < seq_start or seq_end > len(sequence):
            raise ValueError(
                f"Invalid fragment range {seq_start}-{seq_end} for accession {assignment.uniprot_ac}."
            )
        crc64 = crc64_ecma(sequence)
        ref_fields["target_entity_id"].append(assignment.entity_id)
        ref_fields["db_name"].append("UNP")
        ref_fields["db_accession"].append(assignment.uniprot_ac)
        ref_fields["db_code"].append(normalise_optional_text(entry.get("entry_name")))
        ref_fields["gene_name"].append(normalise_optional_text(entry.get("gene_names")))
        taxid = entry.get("taxid")
        ref_fields["ncbi_taxonomy_id"].append(normalise_optional_text(taxid))
        ref_fields["organism_scientific"].append(normalise_optional_text(entry.get("organism")))
        ref_fields["seq_db_align_begin"].append(seq_start)
        ref_fields["seq_db_align_end"].append(seq_end)
        ref_fields["seq_db_isoform"].append("?")
        ref_fields["seq_db_sequence_checksum"].append(crc64)
        ref_fields["seq_db_sequence_version_date"].append(
            normalise_optional_text(entry.get("sequence_version_date"))
        )
    target_ref.clear()
    for field, values in ref_fields.items():
        target_ref[field] = serialise_list(values)

    # _ma_target_entity
    target_entity = ensure_category(template, "_ma_target_entity")
    target_entity["data_id"] = serialise_list("1" for _ in entities)
    target_entity["entity_id"] = serialise_list(ent.entity_id for ent in entities)
    target_entity["origin"] = serialise_list("reference database" for _ in entities)

    # _ma_target_entity_instance
    target_entity_instance = ensure_category(template, "_ma_target_entity_instance")
    asym_ids: List[str] = []
    details: List[str] = []
    entity_refs: List[str] = []
    for assignment in entities:
        for chain_id in assignment.chain_ids:
            asym_ids.append(chain_id)
            details.append(f"Chain {chain_id} from UniProt {assignment.uniprot_ac}")
            entity_refs.append(assignment.entity_id)
    target_entity_instance["asym_id"] = asym_ids
    target_entity_instance["details"] = details
    target_entity_instance["entity_id"] = entity_refs

    # _entity
    entity_cat = ensure_category(template, "_entity")
    entity_ids: List[str] = []
    entity_types: List[str] = []
    src_method: List[str] = []
    descriptions: List[str] = []
    for assignment in entities:
        entry = entries_by_entity[assignment.entity_id]
        full_names = entry.get("protein_full_names") or []
        if isinstance(full_names, str):
            full_names = [full_names]
        description = full_names[0] if full_names else entry.get("entry_name") or assignment.uniprot_ac
        entity_ids.append(assignment.entity_id)
        entity_types.append("polymer")
        src_method.append("man")
        descriptions.append(description)
    entity_cat["id"] = entity_ids
    entity_cat["type"] = entity_types
    entity_cat["src_method"] = src_method
    entity_cat["pdbx_description"] = descriptions

    # _entity_poly
    entity_poly = ensure_category(template, "_entity_poly")
    entity_poly["entity_id"] = serialise_list(ent.entity_id for ent in entities)
    entity_poly["nstd_linkage"] = serialise_list("no" for _ in entities)
    entity_poly["nstd_monomer"] = serialise_list("no" for _ in entities)
    sequences: List[str] = []
    poly_types: List[str] = []
    for assignment in entities:
        entry = entries_by_entity[assignment.entity_id]
        sequence: str = entry.get("sequence") or ""
        if not sequence:
            raise ValueError(f"No sequence stored for accession {assignment.uniprot_ac}.")
        seq_start = assignment.sequence_start or 1
        seq_end = assignment.sequence_end or len(sequence)
        if seq_start < 1 or seq_end < seq_start or seq_end > len(sequence):
            raise ValueError(
                f"Invalid fragment range {seq_start}-{seq_end} for accession {assignment.uniprot_ac}."
            )
        sequences.append(sequence[seq_start - 1:seq_end])
        poly_types.append("polypeptide(L)")
    entity_poly["pdbx_seq_one_letter_code"] = sequences
    entity_poly["type"] = poly_types

    # _struct_ref
    struct_ref = ensure_category(template, "_struct_ref")
    struct_ref_ids: List[int] = []
    db_codes: List[str] = []
    db_names: List[str] = []
    struct_entity_ids: List[str] = []
    pdbx_db_accession: List[str] = []
    for idx, assignment in enumerate(entities, start=1):
        entry = entries_by_entity[assignment.entity_id]
        struct_ref_ids.append(idx)
        db_codes.append(normalise_optional_text(entry.get("entry_name")))
        db_names.append("UNP")
        struct_entity_ids.append(assignment.entity_id)
        pdbx_db_accession.append(assignment.uniprot_ac)
    struct_ref["id"] = struct_ref_ids
    struct_ref["db_code"] = db_codes
    struct_ref["db_name"] = db_names
    struct_ref["entity_id"] = struct_entity_ids
    struct_ref["pdbx_db_accession"] = pdbx_db_accession

    # _struct_ref_seq
    struct_ref_seq = ensure_category(template, "_struct_ref_seq")
    align_ids: List[str] = []
    ref_ids: List[str] = []
    pdb_codes: List[str] = []
    strand_ids: List[str] = []
    seq_align_begin: List[int] = []
    seq_align_end: List[int] = []
    db_align_begin: List[int] = []
    db_align_end: List[int] = []
    align_counter = 1
    for idx, assignment in enumerate(entities, start=1):
        entry = entries_by_entity[assignment.entity_id]
        sequence: str = entry.get("sequence") or ""
        if not sequence:
            raise ValueError(f"No sequence stored for accession {assignment.uniprot_ac}.")
        seq_start = assignment.sequence_start or 1
        seq_end = assignment.sequence_end or len(sequence)
        fragment_length = seq_end - seq_start + 1
        for chain_id in assignment.chain_ids:
            align_ids.append(align_counter)
            ref_ids.append(str(idx))
            pdb_codes.append(model_id)
            strand_ids.append(chain_id)
            seq_align_begin.append(1)
            seq_align_end.append(fragment_length)
            db_align_begin.append(seq_start)
            db_align_end.append(seq_end)
            align_counter += 1
    struct_ref_seq["align_id"] = align_ids
    struct_ref_seq["ref_id"] = ref_ids
    struct_ref_seq["pdbx_PDB_id_code"] = pdb_codes
    struct_ref_seq["pdbx_strand_id"] = strand_ids
    struct_ref_seq["seq_align_beg"] = seq_align_begin
    struct_ref_seq["seq_align_end"] = seq_align_end
    struct_ref_seq["db_align_beg"] = db_align_begin
    struct_ref_seq["db_align_end"] = db_align_end


def build_chains(entries: Sequence[ManifestEntry]) -> List[Dict[str, str]]:
    chains = []
    for entry in entries:
        chains.append(
            {
                "chain_id": entry.chain_id,
                "entity_id": entry.entity_id,
                "uniprot_accession": entry.uniprot_ac,
            }
        )
    return chains


def update_model_identifiers(template: Dict[str, object], model_id: str) -> None:
    categories = template.get("categories")
    if not isinstance(categories, dict):
        return
    entry = categories.get("_entry")
    if isinstance(entry, dict):
        entry["id"] = [model_id]
    database_2 = categories.get("_database_2")
    if isinstance(database_2, dict):
        database_2["database_code"] = [model_id]
    status = categories.get("_pdbx_database_status")
    if isinstance(status, dict):
        status["entry_id"] = [model_id]


def generate_input(args: argparse.Namespace) -> None:
    manifest_entries = load_manifest(args.manifest, args.model_id)
    entity_assignments = group_entities(manifest_entries)

    template_path = args.template
    if not template_path.exists():
        raise FileNotFoundError(f"Template JSON {template_path} does not exist.")
    template = orjson.loads(template_path.read_bytes())
    if "metadata" not in template or "categories" not in template:
        raise ValueError("Template JSON must contain 'metadata' and 'categories' keys.")

    if not args.db.exists():
        raise FileNotFoundError(f"DuckDB database {args.db} does not exist.")

    entries_by_entity: Dict[str, Dict[str, object]] = {}
    with duckdb.connect(str(args.db), read_only=True) as con:
        duckdb_mem = os.environ.get("DUCKDB_MEMORY_LIMIT", "512MB")
        con.execute(f"SET memory_limit = '{duckdb_mem}'")
        for assignment in entity_assignments:
            if assignment.entity_id in entries_by_entity:
                continue
            entry = fetch_entry(con, assignment.uniprot_ac)
            entries_by_entity[assignment.entity_id] = entry

    populate_categories(template, entity_assignments, entries_by_entity, args.model_id)
    update_model_identifiers(template, args.model_id)
    template["chains"] = build_chains(manifest_entries)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("wb") as handle:
        handle.write(orjson.dumps(template, option=orjson.OPT_INDENT_2))
        handle.write(b"\n")
    LOG.info("Wrote ModelCIF metadata JSON for %s to %s", args.model_id, args.out)


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)
    try:
        generate_input(args)
    except Exception as exc:  # pragma: no cover - CLI guard
        LOG.error("Failed to export ModelCIF metadata: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
