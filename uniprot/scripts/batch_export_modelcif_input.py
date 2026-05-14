#!/usr/bin/env python3
"""
Batch export ModelCIF generator input JSON for multiple models.

This script processes multiple models in a single invocation, avoiding the overhead
of repeatedly loading the manifest CSV, DuckDB connections, and template JSON.
"""

from __future__ import annotations

import argparse
import copy
import csv
import logging
import os
import sys
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Set, Tuple

import duckdb
import orjson


LOG = logging.getLogger("uniprot.batch_export_modelcif_input")


@dataclass(frozen=True)
class ManifestEntry:
    model_entity_id: str
    entity_id: str
    chain_id: str
    uniprot_ac: str


@dataclass
class EntityAssignment:
    entity_id: str
    uniprot_ac: str
    chain_ids: List[str]


@dataclass
class ManifestData:
    by_model: Dict[str, List[ManifestEntry]]
    all_accessions: Set[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch export ModelCIF generator metadata JSON for multiple models."
    )
    parser.add_argument(
        "--model-ids",
        required=True,
        type=Path,
        help="File containing model entity IDs (one per line).",
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
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write output JSON files.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 8,
        help="Number of parallel workers (default: all available CPUs).",
    )
    parser.add_argument(
        "--failed-ids-file",
        help="If set, append failed model IDs and errors to this TSV file.",
    )
    parser.add_argument(
        "--stage-name",
        default="stage_09_export_modelcif_input",
        help="Stage label for the failed-ids file.",
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s %(message)s",
    )


def load_model_ids(path: Path) -> List[str]:
    """Load model IDs from a file (one per line)."""
    model_ids = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            model_id = line.strip()
            if model_id and not model_id.startswith("#"):
                model_ids.append(model_id)
    return model_ids


def load_manifest(path: Path, model_ids: List[str]) -> ManifestData:
    """Load manifest entries for specified models into memory."""
    if not path.exists():
        raise FileNotFoundError(f"Manifest file {path} does not exist.")
    
    model_ids_set = set(model_ids)
    by_model: Dict[str, List[ManifestEntry]] = defaultdict(list)
    all_accessions: Set[str] = set()
    
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = {"model_entity_id", "entity_id", "chain_id", "uniprot_ac"}
        if reader.fieldnames is None or expected - set(reader.fieldnames):
            raise ValueError(f"Manifest {path} must contain columns {expected}, found {reader.fieldnames}")
        
        for row in reader:
            model_id = (row.get("model_entity_id") or "").strip()
            if model_id not in model_ids_set:
                continue
            
            entity_id = (row.get("entity_id") or "").strip()
            chain_id = (row.get("chain_id") or "").strip()
            uniprot_ac = (row.get("uniprot_ac") or "").strip()
            
            if not entity_id or not chain_id or not uniprot_ac:
                raise ValueError(f"Incomplete manifest row for model {model_id}: {row}")
            
            entry = ManifestEntry(
                model_entity_id=model_id,
                entity_id=entity_id,
                chain_id=chain_id,
                uniprot_ac=uniprot_ac,
            )
            by_model[model_id].append(entry)
            all_accessions.add(uniprot_ac)
    
    return ManifestData(by_model=dict(by_model), all_accessions=all_accessions)


def group_entities(entries: Sequence[ManifestEntry]) -> List[EntityAssignment]:
    """Group manifest entries by entity_id."""
    grouped: "OrderedDict[str, EntityAssignment]" = OrderedDict()
    for entry in entries:
        key = entry.entity_id
        current = grouped.get(key)
        if current is None:
            grouped[key] = EntityAssignment(
                entity_id=entry.entity_id,
                uniprot_ac=entry.uniprot_ac,
                chain_ids=[entry.chain_id],
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
            current.chain_ids.append(entry.chain_id)
    return list(grouped.values())


def fetch_entries_batch(con: duckdb.DuckDBPyConnection, accessions: Set[str]) -> Dict[str, Dict[str, Any]]:
    """Fetch all entries in a single query."""
    if not accessions:
        return {}
    accession_list = list(accessions)
    placeholders = ",".join("?" for _ in accession_list)
    query = f"SELECT * FROM entry WHERE primary_ac IN ({placeholders})"
    relation = con.execute(query, accession_list)
    table = relation.fetch_arrow_table()
    rows = table.to_pylist()
    return {row["primary_ac"]: row for row in rows}


# CRC64 table for checksum computation
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


def ensure_category(template: Dict[str, object], name: str) -> Dict[str, List[str]]:
    categories = template.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("Template 'categories' must be a dictionary.")
    category = categories.setdefault(name, {})
    if not isinstance(category, dict):
        raise TypeError(f"Category {name} must be a dictionary.")
    return category  # type: ignore[return-value]


def populate_categories(
    template: Dict[str, object],
    entities: Sequence[EntityAssignment],
    entries_by_accession: Dict[str, Dict[str, object]],
    model_id: str,
) -> None:
    """Populate template categories for a single model."""
    categories = template.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("Template must contain a 'categories' dictionary.")

    # Build entity -> entry mapping
    entries_by_entity: Dict[str, Dict[str, object]] = {}
    for assignment in entities:
        if assignment.entity_id not in entries_by_entity:
            entries_by_entity[assignment.entity_id] = entries_by_accession[assignment.uniprot_ac]

    # _ma_target_ref_db_details
    target_ref = ensure_category(template, "_ma_target_ref_db_details")
    ref_fields: Dict[str, List[object]] = defaultdict(list)
    for ordinal, assignment in enumerate(entities, start=1):
        entry = entries_by_entity[assignment.entity_id]
        sequence: str = entry.get("sequence") or ""
        if not sequence:
            raise ValueError(f"No sequence stored for accession {assignment.uniprot_ac}.")
        crc64 = crc64_ecma(sequence)
        ref_fields["target_entity_id"].append(assignment.entity_id)
        ref_fields["db_name"].append("UNP")
        ref_fields["db_accession"].append(assignment.uniprot_ac)
        ref_fields["db_code"].append(normalise_optional_text(entry.get("entry_name")))
        ref_fields["gene_name"].append(normalise_optional_text(entry.get("gene_names")))
        taxid = entry.get("taxid")
        ref_fields["ncbi_taxonomy_id"].append(normalise_optional_text(taxid))
        ref_fields["organism_scientific"].append(normalise_optional_text(entry.get("organism")))
        ref_fields["seq_db_align_begin"].append(1)
        ref_fields["seq_db_align_end"].append(len(sequence))
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
        sequences.append(sequence)
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
        for chain_id in assignment.chain_ids:
            align_ids.append(align_counter)
            ref_ids.append(str(idx))
            pdb_codes.append(model_id)
            strand_ids.append(chain_id)
            length = len(sequence)
            seq_align_begin.append(1)
            seq_align_end.append(length)
            db_align_begin.append(1)
            db_align_end.append(length)
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
    """Build chains list for template."""
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
    """Update model identifiers in template categories."""
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


def process_model(
    model_id: str,
    manifest_entries: List[ManifestEntry],
    entries_by_accession: Dict[str, Dict[str, Any]],
    base_template: Dict[str, Any],
    output_dir: Path,
) -> bool:
    """Process a single model and write output JSON."""
    try:
        entity_assignments = group_entities(manifest_entries)
        
        # Deep copy template to avoid mutation
        template = copy.deepcopy(base_template)
        
        populate_categories(template, entity_assignments, entries_by_accession, model_id)
        update_model_identifiers(template, model_id)
        template["chains"] = build_chains(manifest_entries)
        
        output_path = output_dir / f"{model_id}.json"
        with output_path.open("wb") as handle:
            handle.write(orjson.dumps(template, option=orjson.OPT_INDENT_2))
            handle.write(b"\n")
        
        return True
    except Exception as exc:
        LOG.error("Error processing model %s: %s", model_id, exc)
        return False


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    if not args.db.exists():
        LOG.error("DuckDB database %s does not exist.", args.db)
        return 1

    if not args.template.exists():
        LOG.error("Template JSON %s does not exist.", args.template)
        return 1

    LOG.info("Loading model IDs from %s", args.model_ids)
    model_ids = load_model_ids(args.model_ids)
    LOG.info("Found %d model IDs to process.", len(model_ids))

    if not model_ids:
        LOG.warning("No model IDs to process.")
        return 0

    LOG.info("Loading template from %s", args.template)
    base_template = orjson.loads(args.template.read_bytes())
    if "metadata" not in base_template or "categories" not in base_template:
        LOG.error("Template JSON must contain 'metadata' and 'categories' keys.")
        return 1

    LOG.info("Loading manifest from %s", args.manifest)
    manifest = load_manifest(args.manifest, model_ids)
    LOG.info("Found manifest entries for %d models.", len(manifest.by_model))
    LOG.info("Found %d unique accessions to fetch.", len(manifest.all_accessions))

    LOG.info("Connecting to DuckDB at %s", args.db)
    con = duckdb.connect(str(args.db), read_only=True)
    duckdb_mem = os.environ.get("DUCKDB_MEMORY_LIMIT", "512MB")
    con.execute(f"SET memory_limit = '{duckdb_mem}'")
    try:
        LOG.info("Fetching all entries from DuckDB...")
        entries_by_accession = fetch_entries_batch(con, manifest.all_accessions)
        LOG.info("Fetched %d entries.", len(entries_by_accession))
    finally:
        con.close()

    # Validate we have all needed accessions
    missing_accessions = manifest.all_accessions - set(entries_by_accession.keys())
    if missing_accessions:
        LOG.error("Missing UniProt entries for accessions: %s", ", ".join(missing_accessions))
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def process_single_model(model_id: str) -> Tuple[bool, str]:
        """Process a single model. Returns (success, model_id)."""
        manifest_entries = manifest.by_model.get(model_id)
        if not manifest_entries:
            LOG.warning("Model %s not found in manifest, skipping.", model_id)
            return False, model_id
        success = process_model(model_id, manifest_entries, entries_by_accession, base_template, args.output_dir)
        return success, model_id

    # Process models in parallel
    success_count = 0
    error_count = 0
    failed_ids: list[str] = []
    LOG.info("Processing %d models with %d workers...", len(model_ids), args.workers)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_single_model, model_id): model_id for model_id in model_ids}
        for future in as_completed(futures):
            success, model_id = future.result()
            if success:
                success_count += 1
            else:
                error_count += 1
                failed_ids.append(model_id)

    LOG.info("Batch export complete: %d succeeded, %d failed.", success_count, error_count)

    if failed_ids and args.failed_ids_file:
        failed_path = Path(args.failed_ids_file)
        failed_path.parent.mkdir(parents=True, exist_ok=True)
        with failed_path.open("a") as fh:
            for mid in sorted(failed_ids):
                fh.write(f"{mid}\t{args.stage_name}\texport failed\n")
        LOG.info("Appended %d failed IDs to %s", len(failed_ids), failed_path)

    if success_count == 0 and error_count > 0:
        return 1
    if error_count > 0:
        LOG.warning("Partial success: %d/%d models skipped (missing manifests or metadata).",
                     error_count, success_count + error_count)
    return 0


if __name__ == "__main__":
    sys.exit(main())
