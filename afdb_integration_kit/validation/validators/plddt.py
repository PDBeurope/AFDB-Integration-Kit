from __future__ import annotations

import json
import math
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, List, Tuple

from afdb_integration_kit.quality_assessment.naming import PATTERNS

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult


def _extract_scores(obj: object) -> Tuple[Iterable[float] | None, str]:
    if isinstance(obj, dict):
        scores = obj.get("confidenceScore")
        if isinstance(scores, list):
            return scores, "object"
        return None, "dict_missing_confidenceScore"
    if isinstance(obj, list):
        return obj, "array"
    return None, "unsupported"


def _normalise_scores(
    raw_scores: Iterable[float],
    *,
    min_score: float,
    max_score: float,
    enforce_decimal_places: bool,
) -> Tuple[List[float], List[int], List[int]]:
    scores: List[float] = []
    invalid_indices: List[int] = []
    decimal_issues: List[int] = []
    for idx, value in enumerate(raw_scores):
        try:
            score = float(value)
        except (TypeError, ValueError):
            invalid_indices.append(idx)
            continue
        if math.isnan(score) or math.isinf(score) or not (min_score <= score <= max_score):
            invalid_indices.append(idx)
            continue
        if enforce_decimal_places and not _has_max_two_decimal_places(value):
            decimal_issues.append(idx)
        scores.append(score)
    return scores, invalid_indices, decimal_issues


def _has_max_two_decimal_places(value: object) -> bool:
    try:
        dec = Decimal(str(value))
    except (InvalidOperation, ValueError):
        return False
    exponent = -dec.as_tuple().exponent if dec.as_tuple().exponent < 0 else 0
    return exponent <= 2


def _as_bool(value: object, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _check_length_alignment(
    data: dict,
    enforce: bool,
    path: Path,
) -> List[ValidationResult]:
    if not enforce:
        return []

    residue_numbers = data.get("residueNumber")
    scores = data.get("confidenceScore")
    categories = data.get("confidenceCategory")

    results: List[ValidationResult] = []

    if not isinstance(scores, list):
        return results

    lengths = []
    for label, values in (
        ("residueNumber", residue_numbers),
        ("confidenceScore", scores),
        ("confidenceCategory", categories),
    ):
        if values is None:
            continue
        if not isinstance(values, list):
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code=f"plddt_invalid_{label}",
                    message=f"Field '{label}' must be a list when present.",
                    suggested_fix=f"Ensure '{label}' is emitted as an array matching confidenceScore length.",
                )
            )
            return results
        lengths.append((label, len(values)))

    if lengths:
        base_len = lengths[0][1]
        mismatched = [(label, length) for label, length in lengths if length != base_len]
        if mismatched:
            details = ", ".join(f"{label}={length}" for label, length in lengths)
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_length_mismatch",
                    message=f"confidenceScore and related arrays differ in length ({details}).",
                    suggested_fix="Ensure residueNumber, confidenceScore, and confidenceCategory lists all align.",
                )
            )
            return results

    if residue_numbers:
        sequential = all(
            isinstance(val, int) and val == idx + 1 for idx, val in enumerate(residue_numbers)
        )
        if not sequential:
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_residue_number_sequence",
                    message="residueNumber must start at 1 and increase sequentially by 1.",
                    suggested_fix="Regenerate confidence file ensuring residueNumber is 1..N.",
                )
            )

    return results


def _check_categories(
    raw_scores: Iterable[object],
    categories: List[str],
    *,
    enforce: bool,
    invalid_indices: Iterable[int],
    path: Path,
) -> List[ValidationResult]:
    if not enforce or not categories:
        return []

    invalid_set = set(invalid_indices)
    allowed = {"V", "H", "M", "L", "D"}
    bad_symbols: List[int] = []
    mismatched: List[int] = []

    for idx, (score_value, cat) in enumerate(zip(raw_scores, categories)):
        if idx in invalid_set:
            continue
        if not isinstance(cat, str):
            bad_symbols.append(idx)
            continue

        cat = cat.strip().upper()
        if cat not in allowed:
            bad_symbols.append(idx)
            continue

        try:
            score = float(score_value)
        except (TypeError, ValueError):
            continue  # Already flagged as invalid elsewhere

        expected = _allowed_categories_for_score(score)
        if cat not in expected:
            mismatched.append(idx)

    results: List[ValidationResult] = []
    if bad_symbols:
        index_list = ", ".join(map(str, bad_symbols[:10]))
        extra = "" if len(bad_symbols) <= 10 else f" (and {len(bad_symbols) - 10} more)"
        results.append(
            ValidationResult(
                check="plddt",
                file=path,
                level=Level.ERROR,
                code="plddt_invalid_category_symbol",
                message=f"{len(bad_symbols)} confidenceCategory entries are missing or invalid at indices {index_list}{extra}.",
                suggested_fix="Ensure confidenceCategory entries use one of V, H, M, L, or D.",
            )
        )

    if mismatched:
        index_list = ", ".join(map(str, mismatched[:10]))
        extra = "" if len(mismatched) <= 10 else f" (and {len(mismatched) - 10} more)"
        results.append(
            ValidationResult(
                check="plddt",
                file=path,
                level=Level.ERROR,
                code="plddt_category_mismatch",
                message=f"{len(mismatched)} confidenceCategory values do not match the score ranges at indices {index_list}{extra}.",
                suggested_fix="Update confidenceCategory to reflect the pLDDT score ranges (V>90, H=70-90, M=50-70, L=30-50, D<30).",
            )
        )

    return results


def _allowed_categories_for_score(score: float) -> set[str]:
    if score > 90.0:
        return {"V", "H"}  # V is optional for very high, allow H as historical behaviour
    if 70.0 <= score <= 90.0:
        return {"H"}
    if 50.0 <= score < 70.0:
        return {"M"}
    if 30.0 <= score < 50.0:
        return {"L"}
    return {"D"}


@register_check("plddt")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []

    pattern = PATTERNS.get("plddt")
    candidates = [p for p in files if pattern and pattern.match(p.name)]

    cfg = ctx.config.get("plddt", {})
    min_score = float(cfg.get("min_score", 0.0))
    max_score = float(cfg.get("max_score", 100.0))
    enforce_length_match = _as_bool(cfg.get("enforce_length_match", True), True)
    enforce_categories = _as_bool(cfg.get("enforce_categories", True), True)
    enforce_decimal_places = _as_bool(cfg.get("enforce_decimal_places", True), True)

    for path in sorted(candidates):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_json_parse_error",
                    message=f"Failed to parse JSON: {exc}",
                    suggested_fix="Ensure the pLDDT JSON follows the AFDB confidence schema.",
                )
            )
            continue

        scores_raw, mode = _extract_scores(data)
        if scores_raw is None:
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_unknown_format",
                    message=f"Unrecognised pLDDT JSON structure ({mode}).",
                    suggested_fix="Ensure the JSON contains a confidenceScore array or is a list of scores.",
                )
            )
            continue

        if isinstance(data, dict):
            results.extend(
                _check_length_alignment(
                    data,
                    enforce_length_match,
                    path,
                )
            )

        raw_scores = list(scores_raw)

        scores, invalid_indices, decimal_issues = _normalise_scores(
            raw_scores,
            min_score=min_score,
            max_score=max_score,
            enforce_decimal_places=enforce_decimal_places,
        )

        if invalid_indices:
            index_list = ", ".join(str(i) for i in invalid_indices[:10])
            extra = "" if len(invalid_indices) <= 10 else f" (and {len(invalid_indices) - 10} more)"
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_invalid_values",
                    message=f"{len(invalid_indices)} invalid pLDDT values at indices {index_list}{extra}.",
                    suggested_fix=f"Ensure all pLDDT scores are numeric values between {min_score} and {max_score}.",
                )
            )

        if decimal_issues:
            index_list = ", ".join(str(i) for i in decimal_issues[:10])
            extra = "" if len(decimal_issues) <= 10 else f" (and {len(decimal_issues) - 10} more)"
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_decimal_precision",
                    message=f"{len(decimal_issues)} pLDDT values exceed two decimal places at indices {index_list}{extra}.",
                    suggested_fix="Format scores with at most two decimal places.",
                )
            )

        categories = []
        if isinstance(data, dict):
            cat_values = data.get("confidenceCategory")
            if isinstance(cat_values, list):
                categories = cat_values

        if categories:
            results.extend(
                _check_categories(
                    raw_scores,
                    categories,
                    enforce=enforce_categories,
                    invalid_indices=invalid_indices,
                    path=path,
                )
            )

        if not scores:
            results.append(
                ValidationResult(
                    check="plddt",
                    file=path,
                    level=Level.ERROR,
                    code="plddt_no_scores",
                    message="No valid pLDDT scores found.",
                    suggested_fix="Populate confidenceScore values with numbers in the range [0, 100].",
                )
            )
            continue

        mean_score = sum(scores) / len(scores)
        pct_high = (sum(1 for score in scores if score >= 70.0) / len(scores)) * 100.0

        results.append(
            ValidationResult(
                check="plddt",
                file=path,
                level=Level.INFO,
                code="plddt_summary",
                message="pLDDT metrics calculated.",
                metrics={
                    "mean": mean_score,
                    "length": float(len(scores)),
                    "pct_ge_70": pct_high,
                },
            )
        )

    return results


__all__ = ["run"]
