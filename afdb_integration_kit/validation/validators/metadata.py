from __future__ import annotations

import re
from pathlib import Path
from typing import List

import jsonschema
import orjson

from afdb_integration_kit.metadata.validator import SchemaType, validate_against_schema

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

METADATA_PATTERN = re.compile(r"^AF-metadata-\d+-of-\d+\.json$")
DEFAULT_SCHEMA_TYPE = SchemaType.MODEL.value


@register_check("metadata")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    cfg = ctx.config.get("metadata", {})
    allow_single_file = bool(cfg.get("allow_single_file"))
    schema_type = str(cfg.get("schema_type", cfg.get("type", DEFAULT_SCHEMA_TYPE)))

    metadata_files: List[Path] = []
    for path in files:
        if _is_metadata_batch(path):
            metadata_files.append(path)
        elif allow_single_file and path.suffix.lower() == ".json":
            metadata_files.append(path)

    for metadata_path in sorted(metadata_files):
        results.extend(_validate_metadata_file(metadata_path, schema_type))

    return results


def _is_metadata_batch(path: Path) -> bool:
    return bool(METADATA_PATTERN.match(path.name))


def _validate_metadata_file(path: Path, schema_type: str) -> List[ValidationResult]:
    if not path.exists():
        return [
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_missing_file",
                message="Metadata JSON file is missing.",
                suggested_fix="Provide the metadata JSON file generated for this dataset.",
            )
        ]

    try:
        schema_enum = SchemaType(schema_type.lower())
    except ValueError:
        expected = ", ".join(schema.value for schema in SchemaType)
        return [
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_invalid_schema_type",
                message=(
                    f"Unknown metadata schema type '{schema_type}'. "
                    f"Expected one of: {expected}."
                ),
                suggested_fix="Set metadata.schema_type in validation config to a supported schema type.",
            )
        ]

    try:
        validate_against_schema(path, schema_enum.value)
    except orjson.JSONDecodeError as exc:
        return [
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_json_parse_error",
                message=f"Failed to parse metadata JSON: {exc}",
                suggested_fix="Ensure the metadata file contains valid JSON.",
            )
        ]
    except jsonschema.ValidationError as exc:
        return [
            ValidationResult(
                check="metadata",
                file=path,
                level=Level.ERROR,
                code="metadata_schema_validation_error",
                message=exc.message,
                location=_format_error_path(exc),
                suggested_fix=(
                    f"Update the metadata file to satisfy the '{schema_enum.value}' JSON schema."
                ),
            )
        ]

    return [
        ValidationResult(
            check="metadata",
            file=path,
            level=Level.INFO,
            code="metadata_schema_valid",
            message=f"Validated metadata file against the '{schema_enum.value}' schema.",
        )
    ]


def _format_error_path(exc: jsonschema.ValidationError) -> str | None:
    if not exc.path:
        return None
    return ".".join(str(part) for part in exc.path)


__all__ = ["run"]
