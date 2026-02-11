#!/usr/bin/env python3
"""
Export a single AlphaFold metadata JSON record for a UniProt accession.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import duckdb


@dataclass
class ManifestRow:
    model_entity_id: str
    entity_id: str
    chain_id: str
    uniprot_ac: Optional[str]
    sequence_start: Optional[int]
    sequence_end: Optional[int]
    is_fragment: Optional[bool]
    is_isoform: Optional[bool]
    entity_type: Optional[str]
    average_plddt: Optional[float]
    fraction_plddt_very_low: Optional[float]
    fraction_plddt_low: Optional[float]
    fraction_plddt_confident: Optional[float]
    fraction_plddt_very_high: Optional[float]


@dataclass
class ModelMetadataRow:
    iptm: Optional[float]
    average_plddt: Optional[float]
    complex_name: Optional[str]
    is_am_data: Optional[bool]


@dataclass
class ManifestData:
    by_model: Dict[str, List[ManifestRow]]
    accession_to_models: Dict[str, List[str]]
    model_metadata: Dict[str, ModelMetadataRow]


@dataclass
class ComponentPayload:
    accession: str
    uniprot_id: Optional[str]
    uniprot_description: Optional[str]
    tax_id: Optional[int]
    organism: Optional[str]
    organism_common_names: List[str]
    organism_synonyms: List[str]
    gene: Optional[str]
    gene_synonyms: List[str]
    sequence: str
    checksum: str
    sequence_version_date: Optional[str]
    sequence_start: int
    sequence_end: int
    entity_type: str
    is_isoform: bool
    is_fragment: bool


def parse_bool_field(value: Optional[str], field_name: str) -> Optional[bool]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"1", "true", "t", "yes", "y"}:
        return True
    if lowered in {"0", "false", "f", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean value {value!r} for column '{field_name}'.")


def parse_int_field(value: Optional[str], field_name: str) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Invalid integer value {value!r} for column '{field_name}'.") from exc


def parse_float_field(value: Optional[str], field_name: str) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f"Invalid float value {value!r} for column '{field_name}'.") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a single AF metadata JSON entry for a UniProt accession."
    )
    parser.add_argument(
        "--model-entity-id",
        required=True,
        help="Model entity identifier to export.",
    )
    parser.add_argument(
        "--db",
        required=True,
        type=Path,
        help="DuckDB database containing entry (and optional model_metadata) tables.",
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to dataset-level JSON config (toolUsed, providerId, etc.).",
    )
    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="CSV/TSV mapping between model_entity_id and UniProt accession.",
    )
    parser.add_argument(
        "--model-manifest",
        type=Path,
        help="Optional CSV/TSV model-level manifest (model_entity_id, ipTM, average_plddt, complexName, isAMdata).",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Destination JSON file for the single record.",
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


def load_config(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def normalise_versions(latest: Any, all_versions: Any) -> tuple[int, List[int]]:
    def parse_int(value: Any, default: int) -> int:
        if value is None:
            return default
        if isinstance(value, int):
            return value
        try:
            return int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid latestVersion value: {value!r}") from exc

    latest_version = parse_int(latest, 1)
    if all_versions is None:
        versions = [latest_version]
    elif isinstance(all_versions, list):
        versions = [parse_int(v, latest_version) for v in all_versions]
    elif isinstance(all_versions, str):
        try:
            parsed = json.loads(all_versions)
        except json.JSONDecodeError:
            parsed = [v.strip() for v in all_versions.strip("[]").split(",") if v.strip()]
        if isinstance(parsed, list):
            versions = [parse_int(v, latest_version) for v in parsed]
        else:
            versions = [latest_version]
    else:
        versions = [parse_int(all_versions, latest_version)]
    return latest_version, sorted(set(versions))


def ensure_iso_date(value: Any) -> Optional[str]:
    if value is None:
        return None
    if hasattr(value, "isoformat"):
        iso = value.isoformat()
        if iso.endswith("+00:00"):
            iso = iso[:-6] + "Z"
        if "T" not in iso:
            iso = f"{iso}T00:00:00Z"
        return iso
    date_str = str(value)
    if not date_str:
        return None
    if "T" in date_str:
        return date_str
    return f"{date_str}T00:00:00Z"


def load_manifest(path: Path) -> ManifestData:
    delimiter = "," if path.suffix.lower() != ".tsv" else "\t"
    by_model: Dict[str, List[ManifestRow]] = defaultdict(list)
    accession_to_models: Dict[str, List[str]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        headers = [name.strip() for name in (reader.fieldnames or [])]
        required = {"model_entity_id", "entity_id", "chain_id"}
        if not required.issubset(headers):
            raise ValueError(
                f"Mapping file {path} must contain columns {required}, found {reader.fieldnames}"
            )
        for row in reader:
            model_id = (row.get("model_entity_id") or "").strip()
            entity_id = (row.get("entity_id") or "").strip()
            chain_id = (row.get("chain_id") or "").strip()
            if not model_id or not entity_id or not chain_id:
                continue
            uniprot_ac = (row.get("uniprot_ac") or "").strip() or None
            manifest_row = ManifestRow(
                model_entity_id=model_id,
                entity_id=entity_id,
                chain_id=chain_id,
                uniprot_ac=uniprot_ac,
                sequence_start=parse_int_field(row.get("sequence_start"), "sequence_start"),
                sequence_end=parse_int_field(row.get("sequence_end"), "sequence_end"),
                is_fragment=parse_bool_field(row.get("is_fragment"), "is_fragment"),
                is_isoform=parse_bool_field(row.get("is_isoform"), "is_isoform"),
                entity_type=(row.get("entity_type") or "").strip() or None,
                average_plddt=parse_float_field(row.get("average_plddt"), "average_plddt"),
                fraction_plddt_very_low=parse_float_field(
                    row.get("fraction_plddt_very_low"), "fraction_plddt_very_low"
                ),
                fraction_plddt_low=parse_float_field(row.get("fraction_plddt_low"), "fraction_plddt_low"),
                fraction_plddt_confident=parse_float_field(
                    row.get("fraction_plddt_confident"), "fraction_plddt_confident"
                ),
                fraction_plddt_very_high=parse_float_field(
                    row.get("fraction_plddt_very_high"), "fraction_plddt_very_high"
                ),
            )
            has_start = manifest_row.sequence_start is not None
            has_end = manifest_row.sequence_end is not None
            if has_start != has_end:
                raise ValueError(
                    f"Manifest row {model_id}:{chain_id} must provide both sequence_start and sequence_end or leave both blank."
                )
            if has_start and has_end:
                if manifest_row.sequence_start < 1 or manifest_row.sequence_end < 1:
                    raise ValueError(
                        f"Manifest row {model_id}:{chain_id} has invalid residue range ({manifest_row.sequence_start}, {manifest_row.sequence_end})."
                    )
                if manifest_row.sequence_start > manifest_row.sequence_end:
                    raise ValueError(
                        f"Manifest row {model_id}:{chain_id} has sequence_start greater than sequence_end."
                    )
            if manifest_row.is_fragment:
                if not (has_start and has_end):
                    raise ValueError(
                        f"Manifest row {model_id}:{chain_id} is marked as fragment but missing sequence_start/sequence_end."
                    )
            by_model[model_id].append(manifest_row)
            if uniprot_ac:
                models = accession_to_models[uniprot_ac]
                if model_id not in models:
                    models.append(model_id)
    if not by_model:
        raise ValueError(f"No mappings found in {path}")
    return ManifestData(by_model=dict(by_model), accession_to_models=dict(accession_to_models), model_metadata={})


def load_model_manifest(path: Optional[Path]) -> Dict[str, ModelMetadataRow]:
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"Model manifest {path} does not exist.")
    delimiter = "," if path.suffix.lower() != ".tsv" else "\t"
    model_rows: Dict[str, ModelMetadataRow] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        headers = [name.strip() for name in (reader.fieldnames or [])]
        if "model_entity_id" not in headers:
            raise ValueError(f"Model manifest {path} must contain 'model_entity_id' column, found {headers}")
        for row in reader:
            model_id = (row.get("model_entity_id") or "").strip()
            if not model_id:
                continue
            iptm = parse_float_field(row.get("ipTM") or row.get("complexPredictionAccuracy_ipTM"), "ipTM")
            avg_plddt = parse_float_field(row.get("average_plddt"), "average_plddt")
            complex_name = (row.get("complexName") or "").strip() or None
            is_am_data = parse_bool_field(row.get("isAMdata"), "isAMdata") if "isAMdata" in headers else None
            model_rows[model_id] = ModelMetadataRow(
                iptm=iptm,
                average_plddt=avg_plddt,
                complex_name=complex_name,
                is_am_data=is_am_data,
            )
    return model_rows


def fetch_entries(con: duckdb.DuckDBPyConnection, accessions: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    if not accessions:
        return {}
    placeholders = ",".join("?" for _ in accessions)
    query = f"SELECT * FROM entry WHERE primary_ac IN ({placeholders})"
    relation = con.execute(query, list(accessions))
    table = relation.fetch_arrow_table()
    rows = table.to_pylist()
    return {row["primary_ac"]: row for row in rows}


def as_string_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        items: List[str] = []
        for item in value:
            if item is None:
                continue
            text = str(item).strip()
            if text:
                items.append(text)
        return items
    if isinstance(value, str):
        cleaned = value.strip()
        return [cleaned] if cleaned else []
    return []


def derive_description(entry: Dict[str, Any]) -> Optional[str]:
    full_names = as_string_list(entry.get("protein_full_names"))
    if full_names:
        return full_names[0]
    short_names = as_string_list(entry.get("protein_short_names"))
    if short_names:
        return short_names[0]
    return None


def derive_gene(entry: Dict[str, Any]) -> tuple[Optional[str], List[str]]:
    gene_names = as_string_list(entry.get("gene_names"))
    gene_synonyms = as_string_list(entry.get("gene_synonyms"))
    locus_names = as_string_list(entry.get("gene_ordered_locus_names"))
    orf_names = as_string_list(entry.get("gene_orf_names"))
    candidates = gene_names + locus_names + orf_names
    gene = candidates[0] if candidates else None
    return gene, gene_synonyms


def flatten_unique(values: Sequence[Any]) -> List[str]:
    """Flatten a nested list structure into a unique, order-preserving list of strings."""

    flat: List[str] = []
    seen: set[str] = set()

    def _iter(value: Any) -> List[str]:
        if value is None:
            return []
        if isinstance(value, (list, tuple)):
            items: List[str] = []
            for item in value:
                items.extend(_iter(item))
            return items
        text = str(value).strip()
        return [text] if text else []

    for item in values:
        for text in _iter(item):
            if text not in seen:
                seen.add(text)
                flat.append(text)
    return flat


def summarise_tri_state(flags: Sequence[bool]) -> str:
    if not flags:
        return "none"
    if all(flags):
        return "all"
    if any(flags):
        return "mixed"
    return "none"


OLIGOMERIC_STATE_MAP = {
    1: "monomer",
    2: "dimer",
    3: "trimer",
    4: "tetramer",
    5: "pentamer",
    6: "hexamer",
    7: "heptamer",
    8: "octamer",
    9: "nonamer",
    10: "decamer",
}


def derive_oligomeric_state(subunit_count: int) -> str:
    if subunit_count < 1:
        raise ValueError(f"Oligomeric state requires at least one subunit, received {subunit_count}.")
    return OLIGOMERIC_STATE_MAP.get(subunit_count, "oligomer")


def aggregate_metric(rows: Sequence[ManifestRow], attr: str) -> Optional[float]:
    values = [getattr(row, attr) for row in rows if getattr(row, attr) is not None]
    if not values:
        return None
    return round(sum(values) / len(values), 2)


def aggregate_model_metrics(rows: Sequence[ManifestRow]) -> Dict[str, Optional[float]]:
    return {
        "globalMetricValue": aggregate_metric(rows, "average_plddt"),
        "fractionPlddtVeryLow": aggregate_metric(rows, "fraction_plddt_very_low"),
        "fractionPlddtLow": aggregate_metric(rows, "fraction_plddt_low"),
        "fractionPlddtConfident": aggregate_metric(rows, "fraction_plddt_confident"),
        "fractionPlddtVeryHigh": aggregate_metric(rows, "fraction_plddt_very_high"),
    }


def resolve_consistent_value(
    rows: Sequence[ManifestRow],
    attr: str,
    label: str,
    accession: str,
) -> Optional[Any]:
    values = [getattr(row, attr) for row in rows if getattr(row, attr) is not None]
    if not values:
        return None
    first = values[0]
    for value in values[1:]:
        if value != first:
            model_id = rows[0].model_entity_id
            raise ValueError(
                f"Inconsistent {label} values for accession {accession} in model {model_id}."
            )
    return first


def resolve_model_level_value(
    rows: Sequence[ManifestRow],
    attr: str,
    label: str,
) -> Optional[Any]:
    values = [getattr(row, attr) for row in rows if getattr(row, attr) is not None]
    if not values:
        return None
    first = values[0]
    for value in values[1:]:
        if value != first:
            model_id = rows[0].model_entity_id
            raise ValueError(
                f"Inconsistent {label} values for model {model_id}."
            )
    return first


def build_component_payload(
    accession: str,
    entry: Dict[str, Any],
    rows: Sequence[ManifestRow],
    default_entity_type: Optional[str],
) -> ComponentPayload:
    sequence = entry.get("sequence") or ""
    if not sequence:
        raise ValueError(f"No sequence stored for accession {accession}.")
    seq_len = len(sequence)
    start_override = resolve_consistent_value(rows, "sequence_start", "sequence_start", accession)
    end_override = resolve_consistent_value(rows, "sequence_end", "sequence_end", accession)
    seq_start = start_override or 1
    seq_end = end_override or seq_len
    if seq_start < 1 or seq_start > seq_len:
        raise ValueError(f"sequence_start {seq_start} is out of bounds for accession {accession}.")
    if seq_end < seq_start or seq_end > seq_len:
        raise ValueError(f"sequence_end {seq_end} is out of bounds for accession {accession}.")
    subsequence = sequence[seq_start - 1 : seq_end]
    checksum = hashlib.md5(subsequence.encode("utf-8")).hexdigest()

    fragment_override = resolve_consistent_value(rows, "is_fragment", "is_fragment", accession)
    is_fragment = bool(fragment_override) if fragment_override is not None else (seq_start != 1 or seq_end != seq_len)

    entity_type_override = resolve_consistent_value(rows, "entity_type", "entity_type", accession)
    entity_type = entity_type_override or default_entity_type or "protein"

    isoform_override = resolve_consistent_value(rows, "is_isoform", "is_isoform", accession)
    is_isoform = bool(isoform_override) if isoform_override is not None else False

    uniprot_id = entry.get("entry_name") or accession
    uniprot_description = derive_description(entry)
    if uniprot_description == uniprot_id:
        # Guardrail: avoid ID-as-description duplication when descriptive DE names are absent.
        uniprot_description = None
    gene, gene_synonyms = derive_gene(entry)

    sequence_version_date = ensure_iso_date(entry.get("sequence_version_date"))
    tax_id = entry.get("taxid")
    if tax_id is not None:
        try:
            tax_id = int(tax_id)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid taxId value {tax_id!r} for accession {accession}.") from None

    organism_common_names = as_string_list(entry.get("organism_common_names"))
    organism_synonyms = as_string_list(entry.get("organism_synonyms"))

    return ComponentPayload(
        accession=accession,
        uniprot_id=uniprot_id,
        uniprot_description=uniprot_description,
        tax_id=tax_id,
        organism=entry.get("organism"),
        organism_common_names=organism_common_names,
        organism_synonyms=organism_synonyms,
        gene=gene,
        gene_synonyms=gene_synonyms,
        sequence=subsequence,
        checksum=checksum,
        sequence_version_date=sequence_version_date,
        sequence_start=seq_start,
        sequence_end=seq_end,
        entity_type=entity_type,
        is_isoform=is_isoform,
        is_fragment=is_fragment,
    )


def build_record(
    model_entity_id: str,
    config: Dict[str, Any],
    manifest_rows: Sequence[ManifestRow],
    entry_map: Dict[str, Dict[str, Any]],
    model_metadata: Dict[str, ModelMetadataRow],
) -> List[Dict[str, Any]]:
    if not manifest_rows:
        raise ValueError(f"No manifest entries available for model {model_entity_id}.")

    default_entity_type = config.get("entityType", "protein")
    component_groups: Dict[tuple, List[ManifestRow]] = defaultdict(list)
    component_order: List[tuple] = []
    for row in manifest_rows:
        accession = row.uniprot_ac
        if not accession:
            raise ValueError(
                f"Chain {row.chain_id} (model {row.model_entity_id}) is missing a UniProt accession."
            )
        normalised_entity_type = row.entity_type or default_entity_type
        isoform_flag = row.is_isoform if row.is_isoform is not None else False
        fragment_flag = row.is_fragment if row.is_fragment is not None else False
        key = (
            accession,
            row.sequence_start,
            row.sequence_end,
            fragment_flag,
            isoform_flag,
            normalised_entity_type,
        )
        component_groups[key].append(row)
        if key not in component_order:
            component_order.append(key)

    missing_entries = [key[0] for key in component_order if key[0] not in entry_map]
    if missing_entries:
        raise ValueError(f"Missing UniProt entries for accessions: {', '.join(missing_entries)}.")

    components: List[ComponentPayload] = []
    component_counts: List[int] = []
    key_to_component: Dict[tuple, ComponentPayload] = {}
    for key in component_order:
        accession = key[0]
        rows_for_key = component_groups[key]
        component_counts.append(len(rows_for_key))
        payload = build_component_payload(accession, entry_map[accession], rows_for_key, default_entity_type)
        components.append(payload)
        key_to_component[key] = payload

    accessions = [component.accession for component in components]
    organism_values: List[str] = []
    for component in components:
        if component.organism:
            organism_values.append(component.organism)
    organism_scientific_names = list(dict.fromkeys(organism_values))

    latest_version, all_versions = normalise_versions(
        config.get("latestVersion"),
        config.get("allVersions"),
    )

    unique_id_template = config.get("uniqueIdTemplate")
    base_unique_id = (
        unique_id_template.format(model_entity_id=model_entity_id)
        if unique_id_template
        else model_entity_id
    )
    version_tag = config.get("versionTag", "v1")

    is_complex = len(manifest_rows) > 1
    complex_composition = (
        [f"{component.accession}_{count}" for component, count in zip(components, component_counts)]
        if is_complex
        else None
    )

    is_uniprot_flag = any(bool(row.uniprot_ac) for row in manifest_rows)

    assembly_type = None
    oligomeric_state = None
    if is_complex:
        unique_accessions = {row.uniprot_ac for row in manifest_rows if row.uniprot_ac}
        assembly_type = "Hetero" if len(unique_accessions) > 1 else "Homo"
        oligomeric_state = derive_oligomeric_state(len(manifest_rows))

    model_meta = model_metadata.get(model_entity_id)
    iptm_value = model_meta.iptm if (is_complex and model_meta) else None

    complex_name = model_meta.complex_name if model_meta else None
    if not complex_name and is_complex and assembly_type == "Homo":
        base_desc = components[0].uniprot_description if components else model_entity_id
        prefix = "Homomer" if not oligomeric_state else oligomeric_state.capitalize()
        if oligomeric_state and oligomeric_state.lower().endswith("mer"):
            prefix = f"Homo{oligomeric_state.lower()}"
        complex_name = f"{prefix.capitalize()} of {base_desc}"

    is_am_data = False
    if model_meta and model_meta.is_am_data is not None:
        is_am_data = bool(model_meta.is_am_data)

    records: List[Dict[str, Any]] = []

    for row in manifest_rows:
        accession = row.uniprot_ac
        if not accession:
            continue
        payload = None
        for comp in components:
            if comp.accession == accession:
                payload = comp
                break
        if payload is None:
            raise ValueError(f"No payload built for accession {accession} in model {model_entity_id}.")

        unique_id = f"{base_unique_id}_{version_tag}_{row.chain_id}"

        isoform_flag = payload.is_isoform or ("-" in accession)

        metrics = {
            "globalMetricValue": row.average_plddt,
            "fractionPlddtVeryLow": row.fraction_plddt_very_low,
            "fractionPlddtLow": row.fraction_plddt_low,
            "fractionPlddtConfident": row.fraction_plddt_confident,
            "fractionPlddtVeryHigh": row.fraction_plddt_very_high,
        }

        def metric_value(field: str) -> Optional[float]:
            value = metrics.get(field)
            if value is not None:
                return round(value, 2)
            return config.get(field)

        entry_row = entry_map[payload.accession]
        is_uniprot_reference_proteome = bool(entry_row.get("is_uniprot_reference_proteome"))
        is_uniprot_reviewed = bool(entry_row.get("reviewed"))

        record = {
            "uniqueId": unique_id,
            "toolUsed": config["toolUsed"],
            "modelCreatedDate": ensure_iso_date(config.get("modelCreatedDate")),
            "modelEntityId": model_entity_id,
            "isComplex": is_complex,
            "complexName": complex_name,
            "uniprotAccession": accession,
            "uniprotId": payload.uniprot_id,
            "uniprotDescription": payload.uniprot_description,
            "gene": payload.gene,
            "geneSynonyms": payload.gene_synonyms,
            "taxId": payload.tax_id,
            "organismScientificName": payload.organism,
            "organismCommonNames": payload.organism_common_names,
            "organismSynonyms": payload.organism_synonyms,
            "assemblyType": assembly_type,
            "oligomericState": oligomeric_state,
            "sequence": payload.sequence,
            "sequenceChecksum": payload.checksum,
            "sequenceVersionDate": payload.sequence_version_date,
            "sequenceStart": payload.sequence_start,
            "sequenceEnd": payload.sequence_end,
            "isIsoform": isoform_flag,
            "isFragment": payload.is_fragment,
            "isUniProt": bool(accession),
            "isUniProtReferenceProteome": is_uniprot_reference_proteome,
            "isUniProtReviewed": is_uniprot_reviewed,
            "globalMetricValue": metric_value("globalMetricValue"),
            "fractionPlddtVeryLow": metric_value("fractionPlddtVeryLow"),
            "fractionPlddtLow": metric_value("fractionPlddtLow"),
            "fractionPlddtConfident": metric_value("fractionPlddtConfident"),
            "fractionPlddtVeryHigh": metric_value("fractionPlddtVeryHigh"),
            "latestVersion": latest_version,
            "allVersions": all_versions,
            "providerId": config["providerId"],
            "entityType": payload.entity_type,
            "isAMdata": is_am_data,
        }

        if complex_composition:
            record["complexComposition"] = ",".join(complex_composition)

        if is_complex and iptm_value is not None:
            record["complexPredictionAccuracy_ipTM"] = round(iptm_value, 2)

        records.append(record)

    if not records:
        raise ValueError(f"No chain records produced for model {model_entity_id}.")

    return records


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    if not args.db.exists():
        logging.error("DuckDB database %s does not exist.", args.db)
        return 1

    config = load_config(args.config)
    manifest = load_manifest(args.mapping)
    model_metadata = load_model_manifest(args.model_manifest)
    manifest.model_metadata = model_metadata

    target_model_id = args.model_entity_id
    manifest_rows = manifest.by_model.get(target_model_id)
    if not manifest_rows:
        logging.error("Model %s not found in manifest %s.", target_model_id, args.mapping)
        return 1

    unique_accessions: List[str] = []
    seen: set[str] = set()
    for row in manifest_rows:
        accession = row.uniprot_ac
        if accession and accession not in seen:
            seen.add(accession)
            unique_accessions.append(accession)

    if not unique_accessions:
        logging.error("Model %s does not reference any UniProt accessions.", target_model_id)
        return 1

    con = duckdb.connect(str(args.db), read_only=True)
    try:
        entry_map = fetch_entries(con, unique_accessions)
    finally:
        con.close()

    missing = [acc for acc in unique_accessions if acc not in entry_map]
    if missing:
        logging.error("Accessions %s not found in entry table.", ", ".join(missing))
        return 1

    try:
        records = build_record(target_model_id, config, manifest_rows, entry_map, model_metadata)
    except ValueError as exc:
        logging.error("%s", exc)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    logging.info("Wrote metadata for model %s to %s", target_model_id, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
