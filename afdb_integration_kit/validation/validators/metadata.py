from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Set

import orjson

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

METADATA_PATTERN = re.compile(r"^AF-metadata-\d+-of-\d+\.json$")
ISO_DATETIME_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
MD5_PATTERN = re.compile(r"^[a-f0-9]{32}$", re.IGNORECASE)
UNIQUE_ID_PATTERN = re.compile(r"^(AF-\d{16})_v(\d+)$")
COMPLEX_COMPONENT_PATTERN = re.compile(r"^[A-Z0-9-]+_[1-9]\d*$")

FLOAT_FIELDS_0_1 = [
    "fractionPlddtVeryLow",
    "fractionPlddtLow",
    "fractionPlddtConfident",
    "fractionPlddtVeryHigh",
]

TRI_STATE_VALUES = {"all", "none", "mixed"}
ASSEMBLY_TYPES = {"Homo", "Hetero"}
OLIGOMERIC_STATES = {
    "monomer",
    "dimer",
    "trimer",
    "tetramer",
    "pentamer",
    "hexamer",
    "heptamer",
    "octamer",
    "nonamer",
    "decamer",
    "oligomer",
}

REQUIRED_SCALAR_FIELDS = {
    "uniqueId": str,
    "toolUsed": str,
    "modelCreatedDate": str,
    "modelEntityId": str,
    "providerId": str,
    "isComplex": bool,
    "organismScientificName": str,
    "isFragment": str,
    "isUniProt": str,
    "globalMetricValue": (int, float, type(None)),
    "fractionPlddtVeryLow": (int, float, type(None)),
    "fractionPlddtLow": (int, float, type(None)),
    "fractionPlddtConfident": (int, float, type(None)),
    "fractionPlddtVeryHigh": (int, float, type(None)),
    "latestVersion": int,
    "allVersions": list,
}

LIST_FIELD_TYPES = {
    "accession": (str,),
    "uniprotId": (str,),
    "description": (str,),
    "taxId": (int, type(None)),
    "sequence": (str,),
    "sequenceChecksum": (str,),
    "sequenceVersionDate": (str, type(None)),
    "sequenceStart": (int,),
    "sequenceEnd": (int,),
    "entityType": (str,),
    "isIsoform": (bool,),
}


@register_check("metadata")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    cfg = ctx.config.get("metadata", {})
    allow_single_file = bool(cfg.get("allow_single_file"))

    metadata_files: List[Path] = []
    for path in files:
        if METADATA_PATTERN.match(path.name):
            metadata_files.append(path)
        elif allow_single_file and path.suffix.lower() == ".json":
            metadata_files.append(path)

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
        payload = orjson.loads(path.read_bytes())
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

    if isinstance(payload, dict):
        payload_list: List[dict] = [payload]
    else:
        payload_list = payload

    if not isinstance(payload_list, list) or not payload_list:
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

    for index, entry in enumerate(payload_list, start=1):
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

    list_field_values: Dict[str, List[object]] = {}

    for field, expected_type in REQUIRED_SCALAR_FIELDS.items():
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
        value = entry[field]
        if not isinstance(value, expected_type):
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

    for field, element_types in LIST_FIELD_TYPES.items():
        value = entry.get(field)
        if not isinstance(value, list) or not value:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_type",
                    message=f"{location}: field '{field}' must be a non-empty list.",
                    suggested_fix=f"Populate '{field}' with one entry per complex component.",
                )
            )
            continue
        list_field_values[field] = value
        for item in value:
            if item is None and type(None) in element_types:
                continue
            if not isinstance(item, element_types):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_invalid_type",
                        message=f"{location}: field '{field}' must contain values of type {type_name(element_types)}.",
                        suggested_fix=f"Ensure every '{field}' entry matches {type_name(element_types)}.",
                    )
                )
                break

    is_complex = entry.get("isComplex")
    assembly_value = entry.get("assemblyType")
    oligomer_value = entry.get("oligomericState")
    if is_complex is True:
        if not isinstance(assembly_value, str) or assembly_value not in ASSEMBLY_TYPES:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_assembly_type",
                    message=f"{location}: assemblyType must be one of {sorted(ASSEMBLY_TYPES)} when isComplex is true.",
                    suggested_fix="Set assemblyType to 'Homo' or 'Hetero'.",
                )
            )
        if not isinstance(oligomer_value, str) or oligomer_value not in OLIGOMERIC_STATES:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_oligomeric_state",
                    message=f"{location}: oligomericState must be one of {sorted(OLIGOMERIC_STATES)} when isComplex is true.",
                    suggested_fix="Use the Latin oligomeric names (e.g., dimer, trimer, ...).",
                )
            )
    else:
        if assembly_value is not None and (not isinstance(assembly_value, str) or assembly_value not in ASSEMBLY_TYPES):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_assembly_type",
                    message=f"{location}: assemblyType must be one of {sorted(ASSEMBLY_TYPES)} when provided.",
                    suggested_fix="Set assemblyType to 'Homo' or 'Hetero'.",
                )
            )
        if oligomer_value is not None and (not isinstance(oligomer_value, str) or oligomer_value not in OLIGOMERIC_STATES):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_oligomeric_state",
                    message=f"{location}: oligomericState must be one of {sorted(OLIGOMERIC_STATES)} when provided.",
                    suggested_fix="Use the Latin oligomeric names (e.g., dimer, trimer, ...).",
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
                    message=f"{location}: uniqueId must follow '<modelEntityId>_v<version>' format.",
                    suggested_fix="Construct uniqueId as modelEntityId + '_v' + version.",
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

    # Date formats
    model_created = entry.get("modelCreatedDate")
    if isinstance(model_created, str) and not ISO_DATETIME_PATTERN.match(model_created):
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_datetime",
                message=f"{location}: field 'modelCreatedDate' must be an ISO 8601 date (YYYY-MM-DDT00:00:00Z).",
                suggested_fix="Format 'modelCreatedDate' as YYYY-MM-DDT00:00:00Z.",
            )
        )

    seq_dates = list_field_values.get("sequenceVersionDate")
    if seq_dates:
        for idx, value in enumerate(seq_dates, start=1):
            if value is None:
                continue
            if not isinstance(value, str) or not ISO_DATETIME_PATTERN.match(value):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_invalid_datetime",
                        message=f"{location}: sequenceVersionDate entry #{idx} must be an ISO 8601 date (YYYY-MM-DDT00:00:00Z).",
                        suggested_fix="Format each sequenceVersionDate as YYYY-MM-DDT00:00:00Z or use null when unavailable.",
                    )
                )

    # MD5 checksums
    checksums = list_field_values.get("sequenceChecksum")
    if checksums:
        for idx, checksum in enumerate(checksums, start=1):
            if not isinstance(checksum, str) or not MD5_PATTERN.match(checksum):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_invalid_checksum",
                        message=f"{location}: sequenceChecksum entry #{idx} must be a 32-character hexadecimal MD5.",
                        suggested_fix="Compute the MD5 checksum of each sequence and store it as lowercase hexadecimal.",
                    )
                )

    # Sequence value checks
    sequences = list_field_values.get("sequence")
    if sequences:
        for idx, seq in enumerate(sequences, start=1):
            if not isinstance(seq, str) or not seq:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_invalid_sequence",
                        message=f"{location}: sequence entry #{idx} must be a non-empty string.",
                        suggested_fix="Populate each sequence with the residues for that component.",
                    )
                )

    accession_values = list_field_values.get("accession")
    if accession_values:
        expected_len = len(accession_values)
        for field, values in list_field_values.items():
            if field == "accession":
                continue
            if len(values) != expected_len:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_list_length_mismatch",
                        message=f"{location}: field '{field}' must contain {expected_len} entries (found {len(values)}).",
                        suggested_fix="Ensure every per-component list aligns with the accession list.",
                    )
                )
        sequence_starts = list_field_values.get("sequenceStart")
        sequence_ends = list_field_values.get("sequenceEnd")
        if sequence_starts and sequence_ends and len(sequence_starts) == len(sequence_ends) == expected_len:
            for idx, (start, end) in enumerate(zip(sequence_starts, sequence_ends), start=1):
                if not isinstance(start, int) or not isinstance(end, int):
                    continue
                if start < 1 or end < start:
                    results.append(
                        ValidationResult(
                            check="metadata",
                            file=path,
                            level=Level.ERROR,
                            code="metadata_invalid_sequence_bounds",
                            message=f"{location}: sequenceStart/sequenceEnd for component #{idx} are inconsistent.",
                            suggested_fix="Ensure sequenceStart >= 1 and sequenceEnd >= sequenceStart for every component.",
                        )
                    )
                elif sequences and len(sequences) == expected_len:
                    seq_value = sequences[idx - 1]
                    if isinstance(seq_value, str) and len(seq_value) != end - start + 1:
                        results.append(
                            ValidationResult(
                                check="metadata",
                                file=path,
                                level=Level.ERROR,
                                code="metadata_sequence_length_mismatch",
                                message=f"{location}: sequence #{idx} length does not match the reported residue range.",
                                suggested_fix="Align each sequence with its sequenceStart/sequenceEnd window.",
                            )
                        )

    for field in ("isFragment", "isUniProt"):
        value = entry.get(field)
        if isinstance(value, str) and value not in TRI_STATE_VALUES:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_tri_state",
                    message=f"{location}: field '{field}' must be one of {sorted(TRI_STATE_VALUES)}.",
                    suggested_fix=f"Set '{field}' to 'all', 'none', or 'mixed'.",
                )
            )

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
        if not all_versions:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_all_versions",
                    message=f"{location}: allVersions must not be empty.",
                    suggested_fix="Provide at least one version number.",
                )
            )
        elif not all(isinstance(val, int) and val >= 1 for val in all_versions):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_all_versions",
                    message=f"{location}: allVersions must contain positive integers.",
                    suggested_fix="Populate allVersions with positive integer version numbers.",
                )
            )
        elif isinstance(latest_version, int) and latest_version not in all_versions:
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

    complex_composition = entry.get("complexComposition")
    is_complex = entry.get("isComplex")
    if is_complex is True:
        if not isinstance(complex_composition, list) or not complex_composition:
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_missing_complex_composition",
                    message=f"{location}: complexComposition is required for complex entries.",
                    suggested_fix="Provide a list of '<accession>_<stoichiometry>' entries.",
                )
            )
        else:
            if accession_values and len(complex_composition) != len(accession_values):
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_complex_length_mismatch",
                        message=f"{location}: complexComposition must describe every component.",
                        suggested_fix="Include one complexComposition token per accession.",
                    )
                )
            for token in complex_composition:
                if not isinstance(token, str) or not COMPLEX_COMPONENT_PATTERN.match(token):
                    results.append(
                        ValidationResult(
                            check="metadata",
                            file=path,
                            level=Level.ERROR,
                            code="metadata_invalid_complex_composition",
                            message=f"{location}: complexComposition entries must follow '<accession>_<stoichiometry>'.",
                            suggested_fix="Format each complexComposition entry as UniProtAccession_Stoichiometry.",
                        )
                    )
    elif complex_composition:
        if not isinstance(complex_composition, list):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_complex_composition",
                    message=f"{location}: complexComposition must be a list of strings.",
                    suggested_fix="Store complexComposition as a list if provided.",
                )
            )
        else:
            for token in complex_composition:
                if not isinstance(token, str) or not COMPLEX_COMPONENT_PATTERN.match(token):
                    results.append(
                        ValidationResult(
                            check="metadata",
                            file=path,
                            level=Level.ERROR,
                            code="metadata_invalid_complex_composition",
                            message=f"{location}: complexComposition entries must follow '<accession>_<stoichiometry>'.",
                            suggested_fix="Format each complexComposition entry as UniProtAccession_Stoichiometry.",
                        )
                    )

    return results


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
