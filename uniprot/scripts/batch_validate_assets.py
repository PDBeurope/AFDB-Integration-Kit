#!/usr/bin/env python3
"""
Batch validation of assets - processes multiple models in a single process.
Uses orjson for fast JSON parsing and gemmi for fast PDB parsing.
Uses ProcessPoolExecutor for CPU-bound parallel execution (bypasses GIL).
"""
import argparse
import logging
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Tuple, Optional
import os

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def load_json(path: Path):
    import orjson
    return orjson.loads(path.read_bytes())


def plddt_len(meta: dict) -> int:
    for key in ("plddt", "confidence", "plddt_scores"):
        if key in meta and isinstance(meta[key], list):
            return len(meta[key])
    return 0


def pae_dim(meta: dict) -> str:
    for key in ("predicted_aligned_error", "pae"):
        val = meta.get(key)
        if isinstance(val, list) and val:
            rows = len(val)
            cols = len(val[0]) if isinstance(val[0], list) else 0
            return str(rows) if rows == cols else f"{rows}x{cols}"
    return "0"


def pdb_residue_count_gemmi(path: Path) -> int:
    import gemmi
    structure = gemmi.read_structure(str(path))
    if not structure:
        return 0
    model = structure[0]
    residues = set()
    for chain in model:
        for residue in chain:
            if residue.name == "HOH":
                continue
            for atom in residue:
                if math.isnan(atom.pos.x) or math.isnan(atom.pos.y) or math.isnan(atom.pos.z):
                    raise ValueError(f"NaN coordinates in {path}")
                break
            residues.add((chain.name, residue.seqid.num))
    return len(residues)


def pdb_residue_count_fallback(path: Path) -> int:
    residues = set()
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                res_name = line[17:20].strip()
                if res_name == "HOH":
                    continue
                chain_id = line[21]
                res_seq = line[22:26].strip()
                residues.add((chain_id, res_seq))
    return len(residues)


def validate_single(model_id: str, meta_path: Path, pdb_path: Path) -> tuple:
    """Validate a single model, return (model_id, status, reason, plddt_len, pdb_residues, pae_dim)"""
    try:
        meta = load_json(meta_path)
        plddt = plddt_len(meta)
        pae = pae_dim(meta)
        
        try:
            pdb_res = pdb_residue_count_gemmi(pdb_path)
        except Exception:
            pdb_res = pdb_residue_count_fallback(pdb_path)
        
        # Check for issues
        if isinstance(pae, str) and "x" in pae:
            return (model_id, "mismatch", f"PAE not square ({pae})", str(plddt), str(pdb_res), pae)
        if pae == "0":
            return (model_id, "mismatch", "PAE missing", str(plddt), str(pdb_res), pae)
        if plddt != pdb_res:
            return (model_id, "mismatch", f"pLDDT={plddt}, PDB residues={pdb_res}, PAE={pae}", str(plddt), str(pdb_res), pae)
        
        return (model_id, "ok", "", str(plddt), str(pdb_res), pae)
    except Exception as e:
        return (model_id, "mismatch", f"json_or_parse_error: {e}", "0", "0", "0")


def _validate_single_worker(args: Tuple[str, Optional[str], Optional[str]]) -> Tuple[str, str, str, str, str, str]:
    """
    Module-level worker function for ProcessPoolExecutor (must be picklable).
    
    Args:
        args: Tuple of (model_id, meta_path_str, pdb_path_str) where paths may be None
    
    Returns:
        Tuple of (model_id, status, reason, plddt_len, pdb_residues, pae_dim)
    """
    model_id, meta_path_str, pdb_path_str = args
    
    if not meta_path_str or not pdb_path_str:
        return (model_id, "mismatch", "missing_files", "0", "0", "0")
    
    return validate_single(model_id, Path(meta_path_str), Path(pdb_path_str))


def main():
    parser = argparse.ArgumentParser(description="Batch validate ColabFold assets")
    parser.add_argument("--ids", required=True, help="File containing model IDs, one per line")
    parser.add_argument("--input-dir", required=True, help="Directory containing input files")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 8, help="Number of parallel workers (default: all available CPUs)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    ids_file = Path(args.ids)
    model_ids = [line.strip() for line in ids_file.read_text().splitlines() if line.strip()]

    if not model_ids:
        logger.error("No model IDs found in %s", ids_file)
        sys.exit(1)

    logger.info("Validating %d models with %d workers", len(model_ids), args.workers)

    # Pre-scan directory to build file index (faster than per-model Path.exists() calls)
    logger.info("Scanning input directory for files...")
    file_index = {}  # model_id -> (meta_path, pdb_path)

    # Scan for all JSON and PDB files at once
    all_json = list(input_dir.glob("**/*-meta_v1.json"))
    all_pdb = list(input_dir.glob("**/*-model_v1.pdb"))

    # Build index from JSON files
    json_by_model = {}
    for p in all_json:
        # Extract model_id from filename (e.g., "AF-xxx-meta_v1.json" -> "AF-xxx")
        model_id = p.name.replace("-meta_v1.json", "")
        json_by_model[model_id] = p

    # Build index from PDB files and match with JSON
    for p in all_pdb:
        model_id = p.name.replace("-model_v1.pdb", "")
        if model_id in json_by_model:
            file_index[model_id] = (json_by_model[model_id], p)

    logger.info("Found %d complete model file pairs", len(file_index))

    # Prepare work items as tuples for the module-level worker function
    # Each item is (model_id, meta_path_str, pdb_path_str)
    work_items = []
    for model_id in model_ids:
        meta_path, pdb_path = file_index.get(model_id, (None, None))
        work_items.append((
            model_id,
            str(meta_path) if meta_path else None,
            str(pdb_path) if pdb_path else None
        ))

    logger.info("Starting ProcessPoolExecutor with %d workers", args.workers)
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = list(executor.map(_validate_single_worker, work_items))

    # Write results
    with open(args.output, "w") as f:
        for row in sorted(results, key=lambda x: x[0]):
            f.write("\t".join(row) + "\n")

    ok_count = sum(1 for r in results if r[1] == "ok")
    logger.info("Validation complete: %d ok, %d issues", ok_count, len(results) - ok_count)


if __name__ == "__main__":
    main()
