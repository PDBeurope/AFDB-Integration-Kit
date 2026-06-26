"""Helpers for enriching complex metadata with iPSAE-derived metrics."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


DEFAULT_COMPLEX_ENRICHMENT_METRICS = [
    "pae_cutoff",
    "dist_cutoff",
    "iptm_af",
    "ipsae",
    "ipsae_d0chn",
    "ipsae_d0dom",
    "iptm_d0chn",
    "pDockQ2",
    "LIS",
    "n0res",
    "n0dom",
    "d0res",
    "d0dom",
    "nres1",
    "nres2",
    "dist_nres1",
    "dist_nres2",
    "pDockQ",
    "n0chn",
    "N_clash_backbone",
    "N_clash_heavyAtom",
]


# Regex: strip _{X}{Y} chain-pair suffix (e.g. "_AB", "_BA") from column name.
_CHAIN_PAIR_SUFFIX_RE = re.compile(r"_([A-Z])([A-Z])$")


def model_id_from_ipsae_pdb_path(pdb_path: str) -> str | None:
    """Extract model ID from an iPSAE CSV pdb_path value."""
    if not pdb_path or not pdb_path.strip():
        return None
    name = Path(pdb_path).stem
    for suffix in ("-model_v1", "-model-v1"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def parse_ipsae_csv(csv_path: Path) -> dict[str, dict[str, Any]]:
    """Read iPSAE summary CSV and return ``{model_id: {column: value}}``."""
    result: dict[str, dict[str, Any]] = {}
    if not csv_path.exists():
        return result

    skip_cols = {"pdb_path", "processing_time_ms"}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            model_id = model_id_from_ipsae_pdb_path(row.get("pdb_path", ""))
            if not model_id:
                continue

            metrics: dict[str, Any] = {}
            for col, val in row.items():
                if col in skip_cols:
                    continue
                try:
                    metrics[col] = float(val)
                except (ValueError, TypeError):
                    metrics[col] = val
            result[model_id] = metrics

    return result


def base_metric_name(col: str) -> str:
    """Strip a chain-pair suffix to get the base metric name."""
    match = _CHAIN_PAIR_SUFFIX_RE.search(col)
    return col[: match.start()] if match else col


def ipsae_json_key(col: str) -> str:
    """Map an iPSAE CSV column name to its ``complexPredictionAccuracy_*`` key."""
    non_ipsae_bases = {
        "iptm_af",
        "pDockQ2",
        "LIS",
        "pDockQ",
        "N_clash_backbone",
        "N_clash_heavyAtom",
    }
    if base_metric_name(col) in non_ipsae_bases:
        return f"complexPredictionAccuracy_{col}"
    if col.startswith("ipsae_"):
        return f"complexPredictionAccuracy_{col}"
    return f"complexPredictionAccuracy_ipsae_{col}"


def build_model_enrichment(
    ipsae_row: dict[str, Any],
    clash_row: dict[str, int],
    metrics_filter: list[str],
) -> dict[str, Any]:
    """Build ``complexPredictionAccuracy_*`` values for one model."""
    filter_set = set(metrics_filter)
    out: dict[str, Any] = {}

    for col, val in ipsae_row.items():
        if base_metric_name(col) in filter_set:
            out[ipsae_json_key(col)] = val

    for key, val in clash_row.items():
        if key in filter_set:
            out[f"complexPredictionAccuracy_{key}"] = val

    return out


def build_chain_enrichment(
    ipsae_row: dict[str, Any],
    chain_id: str,
    metrics_filter: list[str],
) -> dict[str, Any]:
    """Build ``complexPredictionAccuracy_*`` values for one chain.

    General metrics go to all chains. Chain-pair metrics go to the chain whose
    ID matches the first letter of the suffix, e.g. ``ipsae_AB`` goes to chain A.
    """
    filter_set = set(metrics_filter)
    out: dict[str, Any] = {}

    for col, val in ipsae_row.items():
        base = base_metric_name(col)
        if base not in filter_set:
            continue

        match = _CHAIN_PAIR_SUFFIX_RE.search(col)
        if match:
            if match.group(1) == chain_id:
                out[ipsae_json_key(col)] = val
        else:
            out[ipsae_json_key(col)] = val

    return out
