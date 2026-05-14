from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

import orjson

from afdb_integration_kit.quality_assessment.naming import PATTERNS

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult


@register_check("relationships")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []

    plddt_pattern = PATTERNS.get("plddt", re.compile(r"^(AF-\d{16})-confidence_(v\d+)\.json$"))
    pae_pattern = PATTERNS.get("pae", re.compile(r"^(AF-\d{16})-predicted_aligned_error_(v\d+)\.json$"))

    plddt_map: Dict[Tuple[str, str], Path] = {}
    pae_map: Dict[Tuple[str, str], Path] = {}

    for path in files:
        name = path.name
        m_plddt = plddt_pattern.match(name)
        if m_plddt:
            key = (m_plddt.group(1), m_plddt.group(2))
            plddt_map[key] = path
            continue

        m_pae = pae_pattern.match(name)
        if m_pae:
            key = (m_pae.group(1), m_pae.group(2))
            pae_map[key] = path

    missing_pae = sorted(key for key in plddt_map if key not in pae_map)
    missing_plddt = sorted(key for key in pae_map if key not in plddt_map)

    for afid, version in missing_pae:
        results.append(
            ValidationResult(
                check="relationships",
                file=plddt_map[(afid, version)],
                level=Level.WARN,
                code="relationship_missing_pae",
                message=f"No PAE file found matching {afid} {version}.",
                suggested_fix="Ensure the PAE JSON is generated alongside the pLDDT file.",
            )
        )

    for afid, version in missing_plddt:
        results.append(
            ValidationResult(
                check="relationships",
                file=pae_map[(afid, version)],
                level=Level.WARN,
                code="relationship_missing_plddt",
                message=f"No pLDDT file found matching {afid} {version}.",
                suggested_fix="Ensure the pLDDT JSON is generated alongside the PAE file.",
            )
        )

    shared = sorted(set(plddt_map) & set(pae_map))
    for afid, version in shared:
        p_path = plddt_map[(afid, version)]
        pae_path = pae_map[(afid, version)]

        plddt_length = _extract_plddt_length(p_path, results)
        pae_dimension = _extract_pae_dimension(pae_path, results)

        if plddt_length is None or pae_dimension is None:
            continue

        if plddt_length != pae_dimension:
            results.append(
                ValidationResult(
                    check="relationships",
                    file=p_path,
                    level=Level.ERROR,
                    code="relationship_length_mismatch",
                    message=(
                        f"Mismatch between pLDDT length ({plddt_length}) and PAE dimension ({pae_dimension}) "
                        f"for {afid} {version}."
                    ),
                    suggested_fix="Regenerate pLDDT and PAE outputs ensuring they cover the same residues.",
                )
            )
        else:
            results.append(
                ValidationResult(
                    check="relationships",
                    file=p_path,
                    level=Level.INFO,
                    code="relationship_summary",
                    message=f"pLDDT and PAE lengths match ({plddt_length}).",
                    metrics={"length": float(plddt_length)},
                )
            )

    return results


def _extract_plddt_length(path: Path, results: List[ValidationResult]) -> int | None:
    try:
        data = orjson.loads(path.read_bytes())
    except Exception as exc:
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_plddt_unreadable",
                message=f"Unable to read pLDDT file: {exc}",
                suggested_fix="Ensure the pLDDT file is valid JSON.",
            )
        )
        return None

    if not isinstance(data, dict):
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_plddt_invalid",
                message="pLDDT payload is not a JSON object.",
                suggested_fix="Ensure the pLDDT file follows the AFDB schema.",
            )
        )
        return None

    scores = data.get("confidenceScore")
    if not isinstance(scores, list):
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_plddt_missing_scores",
                message="pLDDT file does not contain a confidenceScore array.",
                suggested_fix="Generate the pLDDT confidence JSON with confidenceScore entries.",
            )
        )
        return None

    return len(scores)


def _extract_pae_dimension(path: Path, results: List[ValidationResult]) -> int | None:
    try:
        data = orjson.loads(path.read_bytes())
    except Exception as exc:
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_pae_unreadable",
                message=f"Unable to read PAE file: {exc}",
                suggested_fix="Ensure the PAE file is valid JSON.",
            )
        )
        return None

    if not isinstance(data, list) or not data or not isinstance(data[0], dict):
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_pae_invalid",
                message="PAE payload does not follow the expected JSON structure.",
                suggested_fix="Ensure the PAE file contains a single-object list with predicted_aligned_error.",
            )
        )
        return None

    matrix = data[0].get("predicted_aligned_error")
    if not isinstance(matrix, list) or not matrix or not all(isinstance(row, list) for row in matrix):
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_pae_missing_matrix",
                message="PAE file is missing the predicted_aligned_error matrix.",
                suggested_fix="Ensure the PAE file contains the predicted_aligned_error field.",
            )
        )
        return None

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        results.append(
            ValidationResult(
                check="relationships",
                file=path,
                level=Level.WARN,
                code="relationship_pae_unsquare",
                message="PAE matrix is not square, skipping relationship check.",
                suggested_fix="Regenerate the PAE matrix ensuring it is square (NxN).",
            )
        )
        return None

    return size


__all__ = ["run"]
