from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Tuple

import jsonschema
import orjson

from afdb_integration_kit.metadata.validator import (
    FORMAT_CHECKER,
    SCHEMA_PATHS,
    SchemaType,
    _load_json_file,
)

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

MODEL_METADATA_PATTERN = re.compile(r"^AF-model-metadata-\d+-of-\d+\.json$")
CHAIN_METADATA_PATTERN = re.compile(r"^AF-chain-metadata-\d+-of-\d+\.json$")


@register_check("metadata")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    cfg = ctx.config.get("metadata", {})
    allow_single_file = bool(cfg.get("allow_single_file"))

    metadata_files: List[Tuple[Path, Optional[str]]] = []
    for path in files:
        if MODEL_METADATA_PATTERN.match(path.name):
            metadata_files.append((path, "model"))
        elif CHAIN_METADATA_PATTERN.match(path.name):
            metadata_files.append((path, "chain"))
        elif allow_single_file and path.suffix.lower() == ".json":
            metadata_files.append((path, None))

    for metadata_path, schema_hint in sorted(metadata_files, key=lambda pair: str(pair[0])):
        results.extend(_validate_metadata_file(metadata_path, schema_hint))

    return results


def _validate_metadata_file(path: Path, schema_hint: Optional[str]) -> List[ValidationResult]:
    results: List[ValidationResult] = []

    if not path.exists():
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_missing_file",
                message="Metadata batch JSON file is missing.",
                suggested_fix="Provide AF-model-metadata-*-of-*.json or AF-chain-metadata-*-of-*.json.",
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
        payload_list: List[object] = [payload]
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
                suggested_fix="Wrap metadata entries in a non-empty JSON array.",
            )
        )
        return results

    schema_type = schema_hint or _infer_schema_type(payload_list)
    if schema_type is None:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_unknown_schema",
                message="Could not infer metadata schema type (model vs chain).",
                suggested_fix=(
                    "Use canonical filename AF-model-metadata-*-of-*.json or "
                    "AF-chain-metadata-*-of-*.json."
                ),
            )
        )
        return results

    shared_schema_type = {
        "model": SchemaType.MODEL_SUMMARY,
        "chain": SchemaType.COLLECTION_DOC,
    }[schema_type]
    validator = jsonschema.Draft202012Validator(
        _load_json_file(SCHEMA_PATHS[shared_schema_type]),
        format_checker=FORMAT_CHECKER,
    )
    valid_entries = 0
    has_errors = False
    seen_unique_ids = set()
    seen_model_ids = set()

    for index, entry in enumerate(payload_list, start=1):
        if not isinstance(entry, dict):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_entry_not_object",
                    message=f"Entry #{index} is not a JSON object.",
                    suggested_fix="Ensure each metadata entry is a JSON object.",
                )
            )
            has_errors = True
            continue

        location = str(entry.get("uniqueId") or entry.get("modelEntityId") or f"entry_{index}")
        entry_errors = sorted(validator.iter_errors(entry), key=_error_sort_key)
        for error in entry_errors:
            path_str = _format_error_path(error.path)
            message = error.message if not path_str else f"{path_str}: {error.message}"
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_schema_validation_error",
                    message=f"{location}: {message}",
                    suggested_fix="Update the entry to exactly match the strict metadata schema.",
                )
            )
        if entry_errors:
            has_errors = True
            continue

        sequence_start = entry.get("sequenceStart")
        sequence_end = entry.get("sequenceEnd")
        if (
            isinstance(sequence_start, int)
            and isinstance(sequence_end, int)
            and (sequence_end - sequence_start) < 1
        ):
            results.append(
                ValidationResult(
                    check="metadata",
                    file=path,
                    level=Level.ERROR,
                    code="metadata_invalid_sequence_bounds",
                    message=(
                        f"{location}: sequence range invalid; sequenceEnd ({sequence_end}) must be greater than "
                        f"sequenceStart ({sequence_start})."
                    ),
                    suggested_fix="Ensure sequenceEnd > sequenceStart.",
                )
            )
            has_errors = True
            continue

        # Dataset-level uniqueness checks that the schema cannot enforce.
        if schema_type == "chain":
            unique_id = entry.get("uniqueId")
            if unique_id in seen_unique_ids:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_duplicate_unique_id",
                        message=f"{location}: duplicate uniqueId '{unique_id}'.",
                        suggested_fix="Ensure each chain metadata entry has a unique uniqueId.",
                    )
                )
                has_errors = True
                continue
            seen_unique_ids.add(unique_id)
        else:
            model_id = entry.get("modelEntityId")
            if model_id in seen_model_ids:
                results.append(
                    ValidationResult(
                        check="metadata",
                        file=path,
                        level=Level.ERROR,
                        code="metadata_duplicate_model_entity_id",
                        message=f"{location}: duplicate modelEntityId '{model_id}'.",
                        suggested_fix="Ensure each model metadata entry has a unique modelEntityId.",
                    )
                )
                has_errors = True
                continue
            seen_model_ids.add(model_id)

        valid_entries += 1

    if not has_errors:
        results.append(
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.INFO,
                code="metadata_summary",
                message=f"Validated {schema_type} metadata batch with {valid_entries} entries.",
                metrics={"entry_count": float(valid_entries)},
            )
        )

    return results


def _infer_schema_type(payload_list: List[object]) -> Optional[str]:
    first_dict = next((item for item in payload_list if isinstance(item, dict)), None)
    if first_dict is None:
        return None
    keys = set(first_dict.keys())
    if {"uniqueId", "toolUsed", "modelCreatedDate"}.issubset(keys):
        return "chain"
    if {"modelEntityId", "latestVersion", "uniprotAccession"}.issubset(keys):
        return "model"
    return None


def _format_error_path(path_parts) -> str:
    parts = [str(part) for part in path_parts]
    return ".".join(parts)


def _error_sort_key(error: jsonschema.ValidationError) -> tuple:
    return (list(error.path), error.message)


__all__ = ["run"]
