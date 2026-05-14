from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
import orjson

from afdb_integration_kit.quality_assessment.naming import PATTERNS

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

PAE_KEYS = {"predicted_aligned_error", "max_predicted_aligned_error"}


def _validate_single_pae(args: Tuple[Path, bool]) -> List[ValidationResult]:
    """Validate a single PAE file (must be top-level for pickling)."""
    path, enforce_decimal_places = args
    results: List[ValidationResult] = []

    try:
        data = orjson.loads(path.read_bytes())
    except Exception as exc:
        return [
            ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_json_parse_error",
                message=f"Failed to parse JSON: {exc}",
                suggested_fix="Ensure the PAE file is valid JSON following the AFDB schema.",
            )
        ]

    if not isinstance(data, list) or not data:
        return [
            ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_invalid_top_level",
                message="PAE file must be a non-empty list containing a single object.",
                suggested_fix="Wrap the PAE payload in a single-element JSON array.",
            )
        ]

    record = data[0]
    if not isinstance(record, dict):
        return [
            ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_invalid_record",
                message="PAE entry must be a JSON object with predicted_aligned_error and max_predicted_aligned_error.",
                suggested_fix="Ensure the first element of the array is a JSON object.",
            )
        ]

    missing = PAE_KEYS - record.keys()
    if missing:
        return [
            ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_missing_keys",
                message=f"PAE record missing required key(s): {', '.join(sorted(missing))}.",
                suggested_fix="Populate predicted_aligned_error matrix and max_predicted_aligned_error value.",
            )
        ]

    matrix = record.get("predicted_aligned_error")
    max_value = record.get("max_predicted_aligned_error")

    matrix_issue = _validate_matrix(matrix, path, enforce_decimal_places)
    if matrix_issue:
        return [matrix_issue]

    max_issue = _validate_value(max_value, path, "max_predicted_aligned_error", enforce_decimal_places)
    if max_issue:
        return [max_issue]

    dimension = len(matrix) if isinstance(matrix, Sequence) else 0
    return [
        ValidationResult(
            check="pae",
            file=path,
            level=Level.INFO,
            code="pae_summary",
            message="PAE matrix validated.",
            metrics={"dimension": float(dimension), "max_value": float(max_value)},
        )
    ]


@register_check("pae")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    cfg = ctx.config.get("pae", {})
    allow_any_name = bool(cfg.get("allow_any_name"))
    enforce_decimal_places = _as_bool(cfg.get("enforce_decimal_places", True), True)

    pattern = PATTERNS.get("pae")
    candidates: List[Path] = []
    for path in files:
        if allow_any_name:
            candidates.append(path)
        elif pattern and pattern.match(path.name):
            candidates.append(path)

    if not candidates:
        return []

    sorted_candidates = sorted(candidates)
    args = [(p, enforce_decimal_places) for p in sorted_candidates]

    num_workers = min(len(candidates), os.cpu_count() or 4)

    # Use parallel processing for many files (threshold of 10)
    if num_workers > 1 and len(candidates) >= 10:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            batched = list(executor.map(_validate_single_pae, args))
        return [r for batch in batched for r in batch]
    else:
        # Sequential for few files (avoids process spawn overhead)
        return [r for arg in args for r in _validate_single_pae(arg)]


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

    # Convert to numpy array - validates numeric values in one operation
    try:
        arr = np.array(matrix, dtype=np.float64)
    except (ValueError, TypeError):
        return ValidationResult(
            check="pae",
            file=path,
            level=Level.ERROR,
            code="pae_non_numeric_value",
            message="predicted_aligned_error contains non-numeric values.",
            suggested_fix="Ensure all values in the matrix are numeric.",
        )

    # Check matrix is 2D and square
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return ValidationResult(
            check="pae",
            file=path,
            level=Level.ERROR,
            code="pae_matrix_not_square",
            message="predicted_aligned_error must be square (NxN).",
            suggested_fix="Ensure every row in predicted_aligned_error has the same length as the number of rows.",
        )

    # Vectorized decimal place validation - checks ALL values at once
    if enforce_decimal_places:
        rounded = np.round(arr, 2)
        if not np.allclose(arr, rounded, rtol=0, atol=1e-9):
            return ValidationResult(
                check="pae",
                file=path,
                level=Level.ERROR,
                code="pae_decimal_precision",
                message="Some values in predicted_aligned_error exceed two decimal places.",
                suggested_fix="Round all values in the PAE matrix to at most two decimal places.",
            )

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
        f = float(value)
        return f == round(f, 2)
    except (TypeError, ValueError):
        return False


def _as_bool(value: object, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


__all__ = ["run"]
