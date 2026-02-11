#!/usr/bin/env python3
"""
Scan one Swiss-Prot shard for ALTERNATIVE PRODUCTS isoform coverage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence

from uniprot.scripts.extract_subset import (
    build_entry_payload,
    extract_accessions,
    iter_records,
    parse_alt_products,
    parse_var_seq,
)


NON_LOCAL_SEQUENCE_SOURCES = {"external", "not described"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit isoform parsing for one Swiss-Prot shard.")
    parser.add_argument("--input", required=True, type=Path, help="Swiss-Prot shard (*.dat.gz).")
    parser.add_argument("--output", required=True, type=Path, help="Output JSON summary path.")
    parser.add_argument(
        "--focus",
        nargs="*",
        default=[],
        help="Optional list of base accessions to capture detailed examples for.",
    )
    parser.add_argument("--release", default="2025_04", help="Release label for parser call.")
    return parser.parse_args()


def is_emittable_isoform(tokens: Sequence[str], varseqs: Dict[str, tuple[int, int, str]]) -> bool:
    if not tokens:
        return False
    lowered = [token.lower() for token in tokens]
    if len(tokens) == 1 and lowered[0] == "displayed":
        return True
    if any(token in NON_LOCAL_SEQUENCE_SOURCES for token in lowered):
        return False
    vsp_tokens = [token for token in tokens if token.lower() != "displayed"]
    return bool(vsp_tokens) and all(token in varseqs for token in vsp_tokens)


def main() -> int:
    args = parse_args()
    focus_accessions = set(args.focus)

    stats = {
        "input_file": str(args.input),
        "records_scanned": 0,
        "records_with_alt_products": 0,
        "declared_isoforms": 0,
        "emittable_isoforms": 0,
        "emitted_isoforms": 0,
        "missing_emittable_isoforms": 0,
        "non_local_isoforms": 0,
        "displayed_not_dash1_isoforms": 0,
        "focus_examples": [],
    }
    missing_examples: List[Dict[str, object]] = []

    for lines in iter_records(args.input):
        stats["records_scanned"] += 1
        accessions = extract_accessions(lines)
        if not accessions:
            continue
        primary_ac = accessions[0]
        isoform_map = parse_alt_products(lines)
        if not isoform_map:
            continue

        stats["records_with_alt_products"] += 1
        stats["declared_isoforms"] += len(isoform_map)
        varseqs = parse_var_seq(lines)
        payload_rows = build_entry_payload(lines, accessions, args.release, reviewed=True)
        emitted_isoforms = {row["primary_ac"] for row in payload_rows if row.get("is_isoform")}
        stats["emitted_isoforms"] += len(emitted_isoforms)

        local_focus_details = {
            "primary_ac": primary_ac,
            "isoforms": [],
        }
        has_focus = primary_ac in focus_accessions

        for isoform_id, tokens in isoform_map.items():
            lowered = [token.lower() for token in tokens]
            emittable = is_emittable_isoform(tokens, varseqs)
            emitted = isoform_id in emitted_isoforms
            if emittable:
                stats["emittable_isoforms"] += 1
                if not emitted:
                    stats["missing_emittable_isoforms"] += 1
                    if len(missing_examples) < 50:
                        missing_examples.append(
                            {
                                "primary_ac": primary_ac,
                                "isoform_id": isoform_id,
                                "tokens": tokens,
                            }
                        )
            elif any(token in NON_LOCAL_SEQUENCE_SOURCES for token in lowered):
                stats["non_local_isoforms"] += 1

            if "displayed" in lowered and not isoform_id.endswith("-1"):
                stats["displayed_not_dash1_isoforms"] += 1

            if has_focus:
                local_focus_details["isoforms"].append(
                    {
                        "isoform_id": isoform_id,
                        "tokens": tokens,
                        "emittable": emittable,
                        "emitted": emitted,
                    }
                )

        if has_focus:
            stats["focus_examples"].append(local_focus_details)

    stats["missing_examples"] = missing_examples

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
