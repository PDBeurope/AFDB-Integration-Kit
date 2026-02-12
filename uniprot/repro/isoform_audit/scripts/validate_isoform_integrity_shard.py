#!/usr/bin/env python3
"""
Validate isoform parsing integrity for one Swiss-Prot shard.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence

from uniprot.scripts.extract_subset import (
    apply_var_seq_edits,
    build_entry_payload,
    extract_accessions,
    iter_records,
    parse_alt_products,
    parse_var_seq,
)


NON_LOCAL_TOKENS = {"external", "not described"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate isoform parsing for one Swiss-Prot shard.")
    parser.add_argument("--input", required=True, type=Path, help="Swiss-Prot shard path (*.dat.gz).")
    parser.add_argument("--output", required=True, type=Path, help="Validation output JSON path.")
    parser.add_argument("--release", default="2025_04", help="Release label passed to parser.")
    parser.add_argument("--max-issues", type=int, default=100, help="Maximum issue samples to keep.")
    return parser.parse_args()


def append_issue(issues: List[Dict[str, object]], max_issues: int, issue: Dict[str, object]) -> None:
    if len(issues) < max_issues:
        issues.append(issue)


def expected_isoform_length(
    base_len: int,
    tokens: Sequence[str],
    varseqs: Dict[str, tuple[int, int, str]],
) -> int:
    delta = 0
    for token in tokens:
        lowered = token.lower()
        if lowered == "displayed":
            continue
        start, end, replacement = varseqs[token]
        delta += len(replacement) - (end - start + 1)
    return base_len + delta


def main() -> int:
    args = parse_args()

    issues: List[Dict[str, object]] = []
    stats: Dict[str, int] = {
        "records_scanned": 0,
        "records_with_alt_products": 0,
        "declared_isoforms": 0,
        "base_rows_missing": 0,
        "displayed_isoforms": 0,
        "displayed_emitted": 0,
        "displayed_seq_match_base": 0,
        "displayed_seq_mismatch": 0,
        "non_local_isoforms": 0,
        "non_local_emitted_unexpectedly": 0,
        "local_candidate_isoforms": 0,
        "local_emitted": 0,
        "local_missing_unexpectedly": 0,
        "local_length_match_expected": 0,
        "local_length_mismatch": 0,
        "local_seq_match_reconstructed": 0,
        "local_seq_mismatch_reconstructed": 0,
        "unresolved_local_tokens": 0,
    }

    for lines in iter_records(args.input):
        stats["records_scanned"] += 1
        accessions = extract_accessions(lines)
        if not accessions:
            continue
        base_accession = accessions[0]
        iso_map = parse_alt_products(lines)
        if not iso_map:
            continue

        stats["records_with_alt_products"] += 1
        stats["declared_isoforms"] += len(iso_map)

        rows = build_entry_payload(lines, accessions, args.release, reviewed=True)
        emitted = {row["primary_ac"]: row for row in rows}
        base_row = emitted.get(base_accession)
        if base_row is None or not base_row.get("sequence"):
            stats["base_rows_missing"] += 1
            append_issue(
                issues,
                args.max_issues,
                {
                    "type": "base_row_missing",
                    "base_accession": base_accession,
                },
            )
            continue

        base_seq = str(base_row["sequence"])
        base_len = len(base_seq)
        varseqs = parse_var_seq(lines)

        for iso_acc, tokens in iso_map.items():
            lowered_tokens = [token.lower() for token in tokens]
            iso_row = emitted.get(iso_acc)

            is_non_local = any(token in NON_LOCAL_TOKENS for token in lowered_tokens)
            if is_non_local:
                stats["non_local_isoforms"] += 1
                if iso_row is not None:
                    stats["non_local_emitted_unexpectedly"] += 1
                    append_issue(
                        issues,
                        args.max_issues,
                        {
                            "type": "non_local_emitted",
                            "base_accession": base_accession,
                            "isoform_accession": iso_acc,
                            "tokens": tokens,
                        },
                    )
                continue

            is_displayed_only = len(tokens) == 1 and lowered_tokens[0] == "displayed"
            if is_displayed_only:
                stats["displayed_isoforms"] += 1
                if iso_row is None:
                    append_issue(
                        issues,
                        args.max_issues,
                        {
                            "type": "displayed_missing",
                            "base_accession": base_accession,
                            "isoform_accession": iso_acc,
                        },
                    )
                    continue
                stats["displayed_emitted"] += 1
                iso_seq = str(iso_row["sequence"])
                if iso_seq == base_seq:
                    stats["displayed_seq_match_base"] += 1
                else:
                    stats["displayed_seq_mismatch"] += 1
                    append_issue(
                        issues,
                        args.max_issues,
                        {
                            "type": "displayed_seq_mismatch",
                            "base_accession": base_accession,
                            "isoform_accession": iso_acc,
                            "base_length": base_len,
                            "iso_length": len(iso_seq),
                        },
                    )
                continue

            vsp_tokens = [token for token in tokens if token.lower() != "displayed"]
            missing = [token for token in vsp_tokens if token not in varseqs]
            if missing:
                stats["unresolved_local_tokens"] += 1
                if iso_row is not None:
                    append_issue(
                        issues,
                        args.max_issues,
                        {
                            "type": "unresolved_tokens_emitted",
                            "base_accession": base_accession,
                            "isoform_accession": iso_acc,
                            "tokens": tokens,
                            "missing_tokens": missing,
                        },
                    )
                continue

            stats["local_candidate_isoforms"] += 1
            if iso_row is None:
                stats["local_missing_unexpectedly"] += 1
                append_issue(
                    issues,
                    args.max_issues,
                    {
                        "type": "local_missing",
                        "base_accession": base_accession,
                        "isoform_accession": iso_acc,
                        "tokens": tokens,
                    },
                )
                continue

            stats["local_emitted"] += 1
            iso_seq = str(iso_row["sequence"])
            exp_len = expected_isoform_length(base_len, tokens, varseqs)
            if len(iso_seq) == exp_len:
                stats["local_length_match_expected"] += 1
            else:
                stats["local_length_mismatch"] += 1
                append_issue(
                    issues,
                    args.max_issues,
                    {
                        "type": "length_mismatch",
                        "base_accession": base_accession,
                        "isoform_accession": iso_acc,
                        "tokens": tokens,
                        "expected_length": exp_len,
                        "observed_length": len(iso_seq),
                    },
                )

            edits = [varseqs[token] for token in vsp_tokens]
            reconstructed = apply_var_seq_edits(base_seq, edits)
            if reconstructed == iso_seq:
                stats["local_seq_match_reconstructed"] += 1
            else:
                stats["local_seq_mismatch_reconstructed"] += 1
                append_issue(
                    issues,
                    args.max_issues,
                    {
                        "type": "sequence_reconstruction_mismatch",
                        "base_accession": base_accession,
                        "isoform_accession": iso_acc,
                        "tokens": tokens,
                    },
                )

    output = {
        "input_file": str(args.input),
        "stats": stats,
        "issues_sample": issues,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
