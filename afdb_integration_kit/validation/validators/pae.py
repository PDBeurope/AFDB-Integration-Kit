from __future__ import annotations

import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, List, Sequence

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

PAE_KEYS = {"predicted_aligned_error", "max_predicted_aligned_error"}


@register_check("pae")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    pattern = ctx.config.get("pae_pattern")

    cfg = ctx.config.get("pae", {})
    enforce_decimal_places = _as_bool(cfg.get("enforce_decimal_places", True), True)

    candidates = [p for p in files if p.name.endswith("-predicted_aligned_error_v1.json")]

    for path in sorted(candidates):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            results.append(
                ValidationResult(
                    check="pae",
                    file=path,
                    level=Level.ERROR,
                    code="pae_json_parse_error",
                    message=f"Failed to parse JSON: {exc}",
                    suggested_fix="Ensure the PAE file is valid JSON following the AFDB schema.",
                )
            )
            continue

        if not isinstance(data, list) or not data:
            results.append(
                ValidationResult(
                    check="pae",
                    file=path,
                    level=Level.ERROR,
                    code="pae_invalid_top_level",
                    message="PAE file must be a non-empty list containing a single object.",
                    suggested_fix="Wrap the PAE payload in a single-element JSON array.",
                )
            )
            continue

        record = data[0]
        if not isinstance(record, dict):
            results.append(
                ValidationResult(
                    check="pae",
                    file=path,
                    level=Level.ERROR,
                    code="pae_invalid_record",
                    message="PAE entry must be a JSON object with predicted_aligned_error and max_predicted_aligned_error.",
                    suggested_fix="Ensure the first element of the array is a JSON object.",
                )
            )
            continue

        missing = PAE_KEYS - record.keys()
        if missing:
            results.append(
                ValidationResult(
                    check="pae",
                    file=path,
                    level=Level.ERROR,
                    code="pae_missing_keys",
                    message=f"PAE record missing required key(s): {', '.join(sorted(missing))}.",
                    suggested_fix="Populate predicted_aligned_error matrix and max_predicted_aligned_error value.",
                )
            )
            continue

        matrix = record.get("predicted_aligned_error")
        max_value = record.get("max_predicted_aligned_error")

        matrix_issue = _validate_matrix(matrix, path, enforce_decimal_places)
        if matrix_issue:
            results.append(matrix_issue)
            continue

        max_issue = _validate_value(max_value, path, "max_predicted_aligned_error", enforce_decimal_places)
        if max_issue:
            results.append(max_issue)
            continue

        dimension = len(matrix) if isinstance(matrix, Sequence) else 0
        results.append(
            ValidationResult(
                check="pae",
                file=path,
                level=Level.INFO,
                code="pae_summary",
                message="PAE matrix validated.",
                metrics={"dimension": float(dimension), "max_value": float(max_value)},
            )
        )

    return results


def _validate_matrix(
    matrix: object,
    path: Path,
    enforce_decimal_places: bool,
) -> ValidationResult | None:
    if not isinstance(matrix, list) or not matrix:
        return ValidationResult(
            check="pae",
            file=path,
            level=Level.ERROR,
            code="pae_invalid_matrix",
            message="predicted_aligned_error must be a non-empty 2D array.",
            suggested_fix="Ensure predicted_aligned_error is an NxN array of floats.",
        )

    size = len(matrix)
    for row in matrix:
        if not isinstance(row, list) or len(row) != size:
            return ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_matrix_not_square",
                message="predicted_aligned_error must be square (NxN).",
                suggested_fix="Ensure every row in predicted_aligned_error has the same length as the number of rows.",
            )

        for idx, value in enumerate(row):
            issue = _validate_value(value, path, f"predicted_aligned_error[{idx}]", enforce_decimal_places)
            if issue:
                return issue

    return None


def _validate_value(
    value: object,
    path: Path,
    field: str,
    enforce_decimal_places: bool,
) -> ValidationResult | None:
    try:
        float(value)
    except (TypeError, ValueError):
        return ValidationResult(
            check="pae",
            file=path,
            level=Level.ERROR,
            code="pae_non_numeric_value",
            message=f"{field} must be numeric.",
            suggested_fix=f"Provide a numeric value for {field}.",
        )

    if enforce_decimal_places and not _has_two_decimal_places(value):
        return ValidationResult(
            check="pae",
            file=path,
            level=Level.ERROR,
            code="pae_decimal_precision",
            message=f"{field} must have at most two decimal places.",
            suggested_fix="Round values in the PAE file to at most two decimal places.",
        )

    return None


def _has_two_decimal_places(value: object) -> bool:
    try:
        dec = Decimal(str(value))
    except (InvalidOperation, ValueError):
        return False
    exp = -dec.as_tuple().exponent if dec.as_tuple().exponent < 0 else 0
    return exp <= 2


def _as_bool(value: object, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


__all__ = ["run"]
