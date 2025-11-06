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
from pathlib import Path
from typing import Any, Dict, List, Optional

import duckdb


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a single AF metadata JSON entry for a UniProt accession."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Primary UniProt accession to export.",
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


def fetch_row(con: duckdb.DuckDBPyConnection, query: str, params: List[Any]) -> Optional[Dict[str, Any]]:
    relation = con.execute(query, params)
    table = relation.fetch_arrow_table()
    if table.num_rows == 0:
        return None
    return table.to_pylist()[0]


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


def load_mapping(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    delimiter = "," if path.suffix.lower() != ".tsv" else "\t"
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        expected = {"model_entity_id", "uniprot_ac"}
        if not expected.issubset({name.strip() for name in reader.fieldnames or []}):
            raise ValueError(
                f"Mapping file {path} must contain columns {expected}, found {reader.fieldnames}"
            )
        for row in reader:
            model_id = row["model_entity_id"].strip()
            accession = row["uniprot_ac"].strip()
            if not model_id or not accession:
                continue
            if accession in mapping and mapping[accession] != model_id:
                raise ValueError(
                    f"Duplicate mapping for accession {accession}: {mapping[accession]} vs {model_id}"
                )
            mapping[accession] = model_id
    if not mapping:
        raise ValueError(f"No mappings found in {path}")
    return mapping


def build_record(
    accession: str,
    config: Dict[str, Any],
    entry: Dict[str, Any],
    model_entity_id: str,
) -> Dict[str, Any]:
    sequence = entry.get("sequence") or ""
    if not sequence:
        raise ValueError(f"No sequence stored for accession {accession}.")
    seq_start = 1
    seq_end = len(sequence)
    subsequence = sequence
    checksum = entry.get("md5")
    if not checksum:
        checksum = hashlib.md5(subsequence.encode("utf-8")).hexdigest()

    protein_full_names_value = entry.get("protein_full_names") or []
    if isinstance(protein_full_names_value, str):
        protein_full_names_value = [protein_full_names_value]
    protein_full_names = list(protein_full_names_value)
    if protein_full_names:
        uniprot_description = protein_full_names[0]
        protein_full_names = protein_full_names[1:]
    else:
        uniprot_description = None
    protein_short_names = entry.get("protein_short_names") or []
    if isinstance(protein_short_names, str):
        protein_short_names = [protein_short_names]

    gene_synonyms_config = config.get("geneSynonyms") or []
    if isinstance(gene_synonyms_config, str) and gene_synonyms_config:
        try:
            parsed_synonyms = json.loads(gene_synonyms_config)
        except json.JSONDecodeError:
            parsed_synonyms = [token.strip() for token in gene_synonyms_config.split(",") if token.strip()]
        gene_synonyms_list = parsed_synonyms if isinstance(parsed_synonyms, list) else [parsed_synonyms]
    elif isinstance(gene_synonyms_config, list):
        gene_synonyms_list = gene_synonyms_config
    else:
        gene_synonyms_list = []

    entry_gene_synonyms = entry.get("gene_synonyms") or []
    if isinstance(entry_gene_synonyms, str):
        entry_gene_synonyms = [entry_gene_synonyms]

    def dedupe_synonyms(values: List[str]) -> List[str]:
        seen: set[str] = set()
        ordered: List[str] = []
        for value in values:
            if not value:
                continue
            cleaned = value.strip()
            if cleaned and cleaned not in seen:
                seen.add(cleaned)
                ordered.append(cleaned)
        return ordered

    gene_synonyms = dedupe_synonyms(gene_synonyms_list + list(entry_gene_synonyms))

    latest_version, all_versions = normalise_versions(
        config.get("latestVersion"),
        config.get("allVersions"),
    )

    version_tag = config.get("versionTag", "v1")
    ordinal = config.get("ordinal", 1)
    unique_id = config.get("uniqueIdTemplate")
    if unique_id:
        unique_id = unique_id.format(model_entity_id=model_entity_id)
    else:
        unique_id = f"{model_entity_id}_{version_tag}_{ordinal}"

    record = {
        "toolUsed": config["toolUsed"],
        "providerId": config["providerId"],
        "entityType": config.get("entityType", "protein"),
        "isUniProt": bool(config.get("isUniProt", True)),
        "uniqueId": unique_id,
        "modelEntityId": model_entity_id,
        "modelCreatedDate": ensure_iso_date(config.get("modelCreatedDate")),
        "sequenceVersionDate": ensure_iso_date(entry.get("sequence_version_date")),
        "uniprotAccession": accession,
        "uniprotId": entry.get("entry_name"),
        "gene": entry.get("gene_names"),
        "organismScientificName": entry.get("organism"),
        "sequenceChecksum": checksum,
        "globalMetricValue": config.get("globalMetricValue"),
        "fractionPlddtVeryLow": config.get("fractionPlddtVeryLow"),
        "fractionPlddtLow": config.get("fractionPlddtLow"),
        "fractionPlddtConfident": config.get("fractionPlddtConfident"),
        "fractionPlddtVeryHigh": config.get("fractionPlddtVeryHigh"),
        "latestVersion": latest_version,
        "allVersions": all_versions,
        "uniprotDescription": uniprot_description,
        "proteinFullNames": protein_full_names,
        "proteinShortNames": protein_short_names,
        "geneSynonyms": gene_synonyms,
        "sequence": subsequence,
        "sequenceStart": seq_start,
        "sequenceEnd": seq_end,
        "isUniProtReferenceProteome": entry.get("is_uniprot_reference_proteome"),
        "isUniProtReviewed": entry.get("reviewed"),
        "taxId": entry.get("taxid"),
        "stoichiometry": config.get("stoichiometry"),
        "organismCommonNames": entry.get("organisme_common_names") or [],
        "organismScientificNameT": entry.get("organism"),
    }

    missing = [field for field in ("uniqueId", "modelEntityId", "modelCreatedDate") if record.get(field) is None]
    if missing:
        raise ValueError(
            f"Missing required model metadata fields {missing} for accession {accession}."
        )

    return record


def main() -> int:
    args = parse_args()
    configure_logging(args.log_level)

    if not args.db.exists():
        logging.error("DuckDB database %s does not exist.", args.db)
        return 1

    config = load_config(args.config)
    mapping = load_mapping(args.mapping)
    accession = args.accession

    # Use a shared read-only connection so multiple tasks can access the same DB concurrently.
    con = duckdb.connect(str(args.db), read_only=True)
    try:
        entry = fetch_row(
            con,
            "SELECT * FROM entry WHERE primary_ac = ?",
            [accession],
        )
        if entry is None:
            logging.error("Accession %s not found in entry table.", accession)
            return 1

        if accession not in mapping:
            logging.error("Accession %s not found in mapping file %s.", accession, args.mapping)
            return 1

        record = build_record(accession, config, entry, mapping[accession])
    finally:
        con.close()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    logging.info("Wrote metadata for %s to %s", accession, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
