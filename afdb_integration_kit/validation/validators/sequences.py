from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from ..context import ValidationContext
from ..registry import register_check
from ..results import Level, ValidationResult

HEADER_PATTERN = re.compile(r"^>AFDB:(AF-\d{16})$")
ALLOWED_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@register_check("sequences")
def run(files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    sequences_files = [p for p in files if p.name.lower() in {"sequences.fasta", "sequences.fa"}]

    if not sequences_files:
        return results

    for fasta_path in sequences_files:
        results.extend(_validate_fasta(fasta_path))

    return results


def _validate_fasta(path: Path) -> List[ValidationResult]:
    results: List[ValidationResult] = []
    seen_ids: Dict[str, int] = {}
    current_id: str | None = None
    seq_length = 0
    total_sequences = 0
    empty_sequences: List[str] = []
    seen_any_header = False
    header_valid = False

    if not path.is_file():
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.ERROR,
                code="sequences_missing_file",
                message="Sequences FASTA file is missing.",
                suggested_fix="Provide sequences.fasta alongside the dataset.",
            )
        )
        return results

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\n")
                if not line:
                    continue

                if line.startswith(">"):
                    if header_valid and current_id is not None:
                        if seq_length == 0:
                            empty_sequences.append(current_id)
                        else:
                            total_sequences += 1
                    seq_length = 0
                    seen_any_header = True
                    current_id = _validate_header(line, path, line_number, results)
                    header_valid = current_id is not None
                    if header_valid and current_id:
                        seen_ids[current_id] = seen_ids.get(current_id, 0) + 1
                    continue

                if not line.strip():
                    continue

                if not seen_any_header:
                    results.append(
                        ValidationResult(
                            check="sequences",
                            file=path,
                            level=Level.ERROR,
                            code="sequences_missing_header",
                            message=f"Sequence data encountered before header at line {line_number}.",
                            suggested_fix="Ensure each sequence begins with a header line starting with '>'.",
                        )
                    )
                    continue

                if not header_valid or current_id is None:
                    continue

                invalid_chars = [ch for ch in line if ch.upper() not in ALLOWED_AMINO_ACIDS]
                if invalid_chars:
                    chars = "".join(sorted(set(invalid_chars)))
                    results.append(
                        ValidationResult(
                            check="sequences",
                            file=path,
                            level=Level.ERROR,
                            code="sequences_invalid_characters",
                            message=(
                                f"Sequence for {current_id} contains invalid character(s) '{chars}' at line {line_number}."
                            ),
                            suggested_fix="Restrict FASTA sequence characters to standard amino acid codes.",
                        )
                    )
                seq_length += len(line.strip())

            if header_valid and current_id is not None:
                if seq_length == 0:
                    empty_sequences.append(current_id)
                else:
                    total_sequences += 1
    except UnicodeDecodeError as exc:
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.ERROR,
                code="sequences_encoding_error",
                message=f"Unable to decode FASTA file: {exc}",
                suggested_fix="Ensure the FASTA file is UTF-8 encoded.",
            )
        )
        return results

    duplicates = [seq_id for seq_id, count in seen_ids.items() if count > 1]
    if duplicates:
        dup_list = ", ".join(duplicates[:5])
        extra = "" if len(duplicates) <= 5 else f" (and {len(duplicates) - 5} more)"
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.ERROR,
                code="sequences_duplicate_ids",
                message=f"Duplicate sequence identifiers found: {dup_list}{extra}.",
                suggested_fix="Ensure each FASTA header is unique.",
            )
        )

    if empty_sequences:
        ids = ", ".join(empty_sequences[:5])
        extra = "" if len(empty_sequences) <= 5 else f" (and {len(empty_sequences) - 5} more)"
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.ERROR,
                code="sequences_empty_sequence",
                message=f"Header(s) with no sequence data: {ids}{extra}.",
                suggested_fix="Ensure each FASTA entry includes sequence data following the header.",
            )
        )

    if not results:
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.INFO,
                code="sequences_summary",
                message=f"Validated sequences FASTA with {total_sequences} entries.",
                metrics={"sequence_count": float(total_sequences)},
            )
        )

    return results


def _validate_header(
    header: str,
    path: Path,
    line_number: int,
    results: List[ValidationResult],
) -> str | None:
    match = HEADER_PATTERN.match(header)
    if not match:
        results.append(
            ValidationResult(
                check="sequences",
                file=path,
                level=Level.ERROR,
                code="sequences_invalid_header",
                message=f"Header '{header}' at line {line_number} is not in the format '>AFDB:AF-<16 digits>'.",
                suggested_fix="Rename FASTA headers to follow the AFDB convention, e.g. '>AFDB:AF-0000000000000001'.",
            )
        )
        return None
    return match.group(1)


__all__ = ["run"]
