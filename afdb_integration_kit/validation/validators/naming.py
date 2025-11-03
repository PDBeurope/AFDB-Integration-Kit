from __future__ import annotations

from pathlib import Path
from typing import List

from afdb_integration_kit.quality_assessment.naming import (
    REQUIRED_TYPES,
    validate_dataset_naming,
)

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult


@register_check("naming")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    ok, report = validate_dataset_naming(ctx.root)
    results: List[ValidationResult] = []
    config = ctx.config.get("naming", {})
    required_types = set(config.get("required_types", REQUIRED_TYPES))

    if "error" in report:
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.ERROR,
                code="invalid_dataset_root",
                message=report["error"],
                suggested_fix="Provide a directory containing an AFDB dataset.",
            )
        )
        return results

    # Dataset-level signals
    fasta = report.get("sequences_fasta", {})
    status = fasta.get("status")
    if status == "MISSING":
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root / "sequences.fasta",
                level=Level.ERROR,
                code="missing_sequences_fasta",
                message="Required sequences.fasta file is missing from dataset root.",
                suggested_fix="Add sequences.fasta containing the AFDB sequences for this dataset.",
            )
        )
    elif status == "MULTIPLE":
        examples = ", ".join(fasta.get("paths", [])[:3])
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.WARN,
                code="multiple_sequences_fasta",
                message=f"Multiple sequences.fasta files detected: {examples}.",
                suggested_fix="Keep exactly one sequences.fasta at the dataset root.",
            )
        )

    provider = report.get("provider_metadata", {})
    status = provider.get("status")
    if status == "MISSING":
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root / "provider_metadata.json",
                level=Level.ERROR,
                code="missing_provider_metadata",
                message="Provider metadata JSON is missing.",
                suggested_fix="Place the provider metadata JSON in the dataset root following the naming convention.",
            )
        )
    elif status == "MULTIPLE":
        examples = ", ".join(provider.get("candidates", [])[:3])
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.WARN,
                code="multiple_provider_metadata",
                message=f"Multiple provider metadata JSON files detected: {examples}.",
                suggested_fix="Keep a single provider metadata JSON following the naming convention.",
            )
        )

    foldseek_status = report.get("foldseek_index", {}).get("status")
    if foldseek_status == "PARTIAL":
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.WARN,
                code="foldseek_index_partial",
                message="Foldseek index files (.ffindex/.ffdata) are incomplete.",
                suggested_fix="Ensure both .ffindex and .ffdata files are present if foldseek index is provided.",
            )
        )
    elif foldseek_status == "MISSING":
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.INFO,
                code="foldseek_index_missing",
                message="Foldseek index files (.ffindex/.ffdata) not found.",
            )
        )

    entries = report.get("entries", [])
    entry_errors_detected = False
    for entry in entries:
        afid = entry.get("afid")
        version = entry.get("version")
        missing = entry.get("missing", [])
        for missing_type in missing:
            if missing_type not in required_types:
                continue
            results.append(
                ValidationResult(
                    check="naming",
                    file=ctx.root,
                    level=Level.ERROR,
                    code=f"missing_{missing_type}",
                    message=f"{missing_type} file missing for {afid} {version}.",
                    location=f"{afid}-{version}",
                    suggested_fix=f"Add the {missing_type} file for {afid} {version} using the canonical naming pattern.",
                )
            )
            entry_errors_detected = True

    non_compliant = report.get("non_compliant", [])
    for rel_path in non_compliant:
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root / rel_path,
                level=Level.WARN,
                code="non_compliant_filename",
                message=f"Filename '{rel_path}' includes an AF identifier but does not match a canonical pattern.",
                suggested_fix="Rename the file to match the AFDB naming convention.",
            )
        )

    summary = report.get("summary", {})
    if summary.get("entry_count_selected", 0) == 0:
        results.append(
            ValidationResult(
                check="naming",
                file=ctx.root,
                level=Level.ERROR if summary.get("issues") else Level.WARN,
                code="no_entries_selected",
                message="No AFDB entries were discovered under the dataset root.",
                suggested_fix="Ensure the dataset contains files following the AFDB naming patterns.",
            )
        )

    total = float(summary.get("entry_count_total", 0))
    selected = float(summary.get("entry_count_selected", 0))
    metrics = {"entry_count_total": total, "entry_count_selected": selected}
    has_errors = any(res.level is Level.ERROR for res in results if res.check == "naming")
    results.append(
        ValidationResult(
            check="naming",
            file=ctx.root,
            level=Level.INFO if not has_errors else Level.WARN,
            code="summary",
            message="Naming validation summary.",
            metrics=metrics,
        )
    )

    return results


__all__ = ["run"]
