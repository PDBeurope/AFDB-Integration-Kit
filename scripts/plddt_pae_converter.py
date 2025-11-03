from __future__ import annotations

import json
import argparse
import os
from decimal import Decimal, ROUND_HALF_UP
from typing import Any, Dict, List, TypedDict


def _categorise_confidence(score: float) -> str:
    """
    Map a single pLDDT score to category:
        "V" : > 90
        "H" : 70–90  (inclusive of 70 and 90)
        "M" : 50–70
        "L" : 30–50
        "D" : < 30
    """
    if score > 90.0:
        return "V"
    if 70.0 <= score <= 90.0:
        return "H"
    if 50.0 <= score < 70.0:
        return "M"
    if 30.0 <= score < 50.0:
        return "L"
    return "D"


def _round_two_dp(value: float) -> float:
    """Round to two decimal places using bankers-safe decimal arithmetic."""
    return float(Decimal(str(value)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def plddt_to_ingest(plddt: List[float]) -> Dict[str, Any]:
    """Build the AFDB pLDDT payload."""
    n = len(plddt)
    residue_numbers = list(range(1, n + 1))  # 1-based indexing
    categories = [_categorise_confidence(x) for x in plddt]
    return {
        "residueNumber": residue_numbers,
        "confidenceScore": [_round_two_dp(x) for x in plddt],
        "confidenceCategory": categories,
    }


class PAEItem(TypedDict):
    predicted_aligned_error: list[list[float]]
    max_predicted_aligned_error: float


def pae_to_ingest(pae: list[list[float]], max_pae: float) -> list[PAEItem]:
    """Build the AFDB PAE payload with light validation."""
    if not pae or any(len(row) != len(pae) for row in pae):
        raise ValueError("PAE must be a non-empty square matrix (NxN).")
    return [{
        "predicted_aligned_error": [[_round_two_dp(v) for v in row] for row in pae],
        "max_predicted_aligned_error": _round_two_dp(max_pae),
    }]


def convert_file(
    input_path: str,
    out_plddt_path: str | None = None,
    out_pae_path: str | None = None,
    outdir: str | None = None,
) -> Dict[str, str]:
    """Convert a single input JSON to AFDB-format JSONs. Returns written paths."""
    with open(input_path, "r") as f:
        data = json.load(f)

    # Extract required fields
    try:
        plddt = data["plddt"]
        pae = data["pae"]
        max_pae = data["max_pae"]
    except KeyError as e:
        raise KeyError(f"Input JSON is missing required key: {e}")

    plddt_payload = plddt_to_ingest(plddt)
    pae_payload = pae_to_ingest(pae, max_pae)

    # Resolve output paths with AFDB naming convention
    base = os.path.splitext(os.path.basename(input_path))[0]
    default_plddt_name = f"{base}-confidence_v1.json"
    default_pae_name = f"{base}-predicted_aligned_error_v1.json"

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out_plddt_path = out_plddt_path or os.path.join(outdir, default_plddt_name)
        out_pae_path   = out_pae_path   or os.path.join(outdir, default_pae_name)
    else:
        out_plddt_path = out_plddt_path or default_plddt_name
        out_pae_path   = out_pae_path   or default_pae_name

    # Write
    with open(out_plddt_path, "w") as f:
        json.dump(plddt_payload, f, separators=(",", ":"), ensure_ascii=False)
    with open(out_pae_path, "w") as f:
        json.dump(pae_payload, f, separators=(",", ":"), ensure_ascii=False)

    return {"plddt": out_plddt_path, "pae": out_pae_path}


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Convert AlphaFold score JSON to AFDB ingestion format "
                    "(outputs: *-confidence_v1.json and *-predicted_aligned_error_v1.json)."
    )
    p.add_argument("input", help="Path to input JSON with keys: plddt, pae, max_pae")
    p.add_argument("--outdir", help="Directory to write outputs (defaults use AFDB names)")
    p.add_argument("--plddt", help="Explicit output path for pLDDT JSON")
    p.add_argument("--pae", help="Explicit output path for PAE JSON")
    return p


def main(argv: List[str] | None = None) -> None:
    args = _build_arg_parser().parse_args(argv)
    paths = convert_file(args.input, out_plddt_path=args.plddt, out_pae_path=args.pae, outdir=args.outdir)
    print(json.dumps(paths, indent=2))


if __name__ == "__main__":
    main()
