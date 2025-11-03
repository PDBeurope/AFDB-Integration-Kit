from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

METADATA_PATTERN = re.compile(r"^AF-metadata-\d+-of-\d+\.json$")
ISO_DATETIME_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
MD5_PATTERN = re.compile(r"^[a-f0-9]{32}$", re.IGNORECASE)
UNIQUE_ID_PATTERN = re.compile(r"^(AF-\d{16})_v(\d+)_([1-9]\d*)$")
COMPLEX_COMPONENT_PATTERN = re.compile(r"^[0-9A-Fa-f]{16}_[1-9]\d*(?:,[0-9A-Fa-f]{16}_[1-9]\d*)*$")

FLOAT_FIELDS_0_1 = [
    "fractionPlddtVeryLow",
    "fractionPlddtLow",
    "fractionPlddtConfident",
    "fractionPlddtVeryHigh",
]

REQUIRED_FIELDS = {
    "uniqueId": str,
    "toolUsed": str,
    "modelCreatedDate": str,
    "modelEntityId": str,
    "providerId": str,
    "entityType": str,
    "sequence": str,
    "sequenceChecksum": str,
    "sequenceStart": int,
    "sequenceEnd": int,
    "isUniProt": bool,
    "globalMetricValue": (int, float),
    "fractionPlddtVeryLow": (int, float),
    "fractionPlddtLow": (int, float),
    "fractionPlddtConfident": (int, float),
    "fractionPlddtVeryHigh": (int, float),
    "latestVersion": int,
    "allVersions": list,
    "stoichiometry": int,
}

OPTIONAL_FIELD_TYPES = {
    "complexName": str,
    "complexComposition": str,
    "uniprotAccession": str,
    "uniprotId": str,
    "uniprotDescription": str,
    "geneSynonyms": list,
    "gene": str,
    "isUniProtReferenceProteome": bool,
    "isUniProtReviewed": bool,
    "taxId": int,
    "organismScientificName": str,
    "otherTaxIds": list,
    "otherOrganismScientificNames": list,
    "ipTM": (int, float),
    "ipSAE": (int, float),
    "sequenceVersionDate": str,
    "organismCommonNames": list,
    "organismSynonyms": list,
    "proteinFullNames": list,
    "proteinShortNames": list,
    "keywords": list,
    "taxonomyLineage": list,
    "functions": list,
    "alternativeNames": list,
    "catalyticActivities": list,
    "organismScientificNameT": str,
}

UNIPROT_REQUIRED_FIELDS = {
    field: OPTIONAL_FIELD_TYPES[field]
    for field in (
        "uniprotId",
        "uniprotDescription",
        "geneSynonyms",
        "gene",
        "isUniProtReviewed",
        "taxId",
        "organismScientificName",
        "sequenceVersionDate",
        "organismCommonNames",
        "proteinFullNames",
        "proteinShortNames",
        "keywords",
        "taxonomyLineage",
        "functions",
        "alternativeNames",
        "catalyticActivities",
    )
}


@register_check("metadata")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    metadata_files = [p for p in files if METADATA_PATTERN.match(p.name)]

    for metadata_path in sorted(metadata_files):
        results.extend(_validate_metadata_file(metadata_path))

    return results


def _validate_metadata_file(path: Path) -> List[ValidationResult]:
    results: List[ValidationResult] = []

    if not path.exists():
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_missing_file",
                message="Metadata batch JSON file is missing.",
                suggested_fix="Provide the AF-metadata-*-of-*.json file generated for this dataset.",
            )
        )
        return results

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_json_parse_error",
                message=f"Failed to parse metadata JSON: {exc}",
                suggested_fix="Ensure the metadata batch file contains valid JSON.",
            )
        )
        return results

    if not isinstance(payload, list) or not payload:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_structure",
                message="Metadata batch must be a non-empty JSON array.",
                suggested_fix="Wrap entity metadata objects in a JSON array.",
            )
        )
        return results

    seen_unique_ids: Set[str] = set()
    entry_count = 0

    for index, entry in enumerate(payload, start=1):
        if not isinstance(entry, dict):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_entry_not_object",
                    message=f"Entry #{index} is not a JSON object.",
                    suggested_fix="Ensure each metadata entry is a JSON object with field/value pairs.",
                )
            )
            continue

        entry_results = _validate_entry(entry, path, index, seen_unique_ids)
        results.extend(entry_results)
        if not any(res.level is Level.ERROR for res in entry_results):
            entry_count += 1

    if not any(res.level is Level.ERROR for res in results):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.INFO,
                code="metadata_summary",
                message=f"Validated metadata batch with {entry_count} entries.",
                metrics={"entry_count": float(entry_count)},
            )
        )

    return results


def _validate_entry(
    entry: dict,
    path: Path,
    index: int,
    seen_unique_ids: Set[str],
) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    location = entry.get("uniqueId") or f"entry_{index}"

    # Required fields presence and type
    for field, expected_type in REQUIRED_FIELDS.items():
        if field not in entry:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_missing_field",
                    message=f"{location}: required field '{field}' is missing.",
                    suggested_fix=f"Populate the '{field}' field for every metadata entry.",
                )
            )
            continue
        if not isinstance(entry[field], expected_type):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_type",
                    message=f"{location}: field '{field}' must be of type {type_name(expected_type)}.",
                    suggested_fix=f"Ensure '{field}' is stored as {type_name(expected_type)}.",
                )
            )

    # Optional types
    for field, expected_type in OPTIONAL_FIELD_TYPES.items():
        if field in entry and not isinstance(entry[field], expected_type):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_type",
                    message=f"{location}: field '{field}' must be of type {type_name(expected_type)}.",
                    suggested_fix=f"Ensure '{field}' is stored as {type_name(expected_type)} when present.",
                )
            )

    if entry.get("isUniProt") is True:
        for field in UNIPROT_REQUIRED_FIELDS:
            if field not in entry:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_missing_uniprot_field",
                        message=f"{location}: field '{field}' is required when isUniProt is true.",
                        suggested_fix=f"Populate the '{field}' field for UniProt-sourced entries.",
                    )
                )

    # Unique ID format
    unique_id = entry.get("uniqueId")
    model_entity_id = entry.get("modelEntityId")
    if isinstance(unique_id, str) and isinstance(model_entity_id, str):
        match = UNIQUE_ID_PATTERN.match(unique_id)
        if not match or match.group(1) != model_entity_id:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_unique_id",
                    message=f"{location}: uniqueId must follow '<modelEntityId>_v<version>_<ordinal>' format.",
                    suggested_fix="Construct uniqueId as modelEntityId + '_v' + version + '_' + ordinal.",
                )
            )
        elif unique_id in seen_unique_ids:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_duplicate_unique_id",
                    message=f"{location}: duplicate uniqueId detected.",
                    suggested_fix="Ensure each metadata entry has a unique uniqueId.",
                )
            )
        else:
            seen_unique_ids.add(unique_id)
    else:
        unique_id = None

    # Date formats
    for field in ("modelCreatedDate", "sequenceVersionDate"):
        if field in entry and entry[field] is not None:
            value = entry[field]
            if not isinstance(value, str) or not ISO_DATETIME_PATTERN.match(value):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_invalid_datetime",
                        message=f"{location}: field '{field}' must be an ISO 8601 date (YYYY-MM-DDT00:00:00Z).",
                        suggested_fix=f"Format '{field}' as YYYY-MM-DDT00:00:00Z.",
                    )
                )

    # MD5 checksum
    checksum = entry.get("sequenceChecksum")
    if isinstance(checksum, str) and not MD5_PATTERN.match(checksum):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_checksum",
                message=f"{location}: sequenceChecksum must be a 32-character hexadecimal MD5.",
                suggested_fix="Compute the MD5 checksum of the sequence and store it as lowercase hexadecimal.",
            )
        )

    # Sequence length and residue bounds
    sequence = entry.get("sequence")
    seq_start = entry.get("sequenceStart")
    seq_end = entry.get("sequenceEnd")
    if isinstance(sequence, str) and isinstance(seq_start, int) and isinstance(seq_end, int):
        if seq_start < 1 or seq_end < seq_start:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_sequence_bounds",
                    message=f"{location}: sequenceStart must be >=1 and <= sequenceEnd.",
                    suggested_fix="Check the reported residue range for this sequence.",
                )
            )
        else:
            expected_length = seq_end - seq_start + 1
            if len(sequence) != expected_length:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_sequence_length_mismatch",
                        message=f"{location}: sequence length ({len(sequence)}) does not match sequenceStart/sequenceEnd ({expected_length}).",
                        suggested_fix="Ensure sequenceStart, sequenceEnd, and sequence length are consistent.",
                    )
                )

    # Fractions between 0 and 1
    for field in FLOAT_FIELDS_0_1:
        value = entry.get(field)
        if value is not None:
            if not isinstance(value, (int, float)) or not (0.0 <= float(value) <= 1.0):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_fraction_out_of_range",
                        message=f"{location}: field '{field}' must be between 0 and 1.",
                        suggested_fix=f"Clamp '{field}' to the [0, 1] interval.",
                    )
                )

    global_metric = entry.get("globalMetricValue")
    if global_metric is not None:
        if not isinstance(global_metric, (int, float)) or not (0.0 <= float(global_metric) <= 100.0):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_metric_out_of_range",
                    message=f"{location}: globalMetricValue must be between 0 and 100.",
                    suggested_fix="Clamp globalMetricValue to the [0, 100] interval.",
                )
            )

    # allVersions consistency
    latest_version = entry.get("latestVersion")
    all_versions = entry.get("allVersions")
    if isinstance(latest_version, int) and latest_version < 1:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_latest_version",
                message=f"{location}: latestVersion must be a positive integer.",
                suggested_fix="Ensure latestVersion is >= 1.",
            )
        )
    if isinstance(all_versions, list):
        if not all(isinstance(val, int) for val in all_versions):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_all_versions",
                    message=f"{location}: allVersions must be a list of integers.",
                    suggested_fix="Populate allVersions with integer version numbers.",
                )
            )
        elif not all_versions:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_all_versions",
                    message=f"{location}: allVersions must not be empty.",
                    suggested_fix="Populate allVersions with at least one version number.",
                )
            )
        elif latest_version not in all_versions:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_latest_not_in_all_versions",
                    message=f"{location}: latestVersion must appear in allVersions.",
                    suggested_fix="Include latestVersion in the allVersions list.",
                )
            )

    # optional lists should contain strings or ints
    _validate_string_list(entry, "geneSynonyms", path, location, results)
    _validate_string_list(entry, "organismCommonNames", path, location, results)
    _validate_string_list(entry, "organismSynonyms", path, location, results)
    _validate_string_list(entry, "proteinFullNames", path, location, results)
    _validate_string_list(entry, "proteinShortNames", path, location, results)
    _validate_string_list(entry, "keywords", path, location, results)
    _validate_string_list(entry, "taxonomyLineage", path, location, results)
    _validate_string_list(entry, "functions", path, location, results)
    _validate_string_list(entry, "alternativeNames", path, location, results)
    _validate_string_list(entry, "catalyticActivities", path, location, results)
    _validate_int_list(entry, "otherTaxIds", path, location, results)
    _validate_string_list(entry, "otherOrganismScientificNames", path, location, results)

    if "otherTaxIds" in entry and "otherOrganismScientificNames" in entry:
        tax_ids = entry.get("otherTaxIds")
        other_names = entry.get("otherOrganismScientificNames")
        if isinstance(tax_ids, list) and isinstance(other_names, list) and len(tax_ids) != len(other_names):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_other_lists_mismatch",
                    message=f"{location}: otherTaxIds and otherOrganismScientificNames must have the same number of entries.",
                    suggested_fix="Align the otherOrganismScientificNames list with otherTaxIds.",
                )
            )

    complex_composition = entry.get("complexComposition")
    if isinstance(complex_composition, str) and complex_composition and not COMPLEX_COMPONENT_PATTERN.match(complex_composition):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_complex_composition",
                message=f"{location}: complexComposition must be comma-separated '<checksum>_<stoichiometry>' pairs.",
                suggested_fix="Format complexComposition as comma-separated hex checksum plus underscore stoichiometry values.",
            )
        )

    stoichiometry = entry.get("stoichiometry")
    if isinstance(stoichiometry, int) and stoichiometry < 1:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_stoichiometry",
                message=f"{location}: stoichiometry must be a positive integer.",
                suggested_fix="Ensure stoichiometry reflects the number of copies (>= 1).",
            )
        )

    return results


def _validate_string_list(entry: dict, field: str, path: Path, location: str, results: List[ValidationResult]) -> None:
    value = entry.get(field)
    if value is None:
        return
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_type",
                message=f"{location}: field '{field}' must be a list of strings.",
                suggested_fix=f"Ensure '{field}' is an array of strings.",
            )
        )


def _validate_int_list(entry: dict, field: str, path: Path, location: str, results: List[ValidationResult]) -> None:
    value = entry.get(field)
    if value is None:
        return
    if not isinstance(value, list) or not all(isinstance(item, int) for item in value):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_type",
                message=f"{location}: field '{field}' must be a list of integers.",
                suggested_fix=f"Ensure '{field}' is an array of integers.",
            )
        )


def type_name(expected: object) -> str:
    if isinstance(expected, tuple):
        return " or ".join(type_name(t) for t in expected)
    if expected is str:
        return "string"
    if expected is int:
        return "integer"
    if expected is bool:
        return "boolean"
    if expected is list:
        return "list"
    if expected is float:
        return "number"
    return getattr(expected, "__name__", str(expected))


__all__ = ["run"]
