#!/usr/bin/env python3
"""
Batch ipSAE calculator - computes interface quality metrics for protein complexes.

Uses the optimized C++ ipsae_optimized binary for high-performance calculations
of pDockQ, pDockQ2, LIS, and ipSAE metrics.
"""
import argparse
import logging
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Binary path relative to this script
IPSAE_CPP_DIR = Path(__file__).parent.parent.parent / "afdb_integration_kit" / "ipsae_cpp"
IPSAE_BINARY = IPSAE_CPP_DIR / "ipsae_optimized"


def ensure_ipsae_binary() -> bool:
    """Ensure the ipSAE binary exists, compiling if necessary.
    
    Returns:
        True if binary is available, False if compilation failed.
    """
    if IPSAE_BINARY.exists():
        return True
    
    logger.info("ipSAE binary not found, attempting to compile...")
    makefile = IPSAE_CPP_DIR / "Makefile"
    if not makefile.exists():
        logger.error(f"Cannot compile: Makefile not found at {makefile}")
        return False
    
    result = subprocess.run(
        ["make"],
        cwd=IPSAE_CPP_DIR,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        logger.error(f"Compilation failed: {result.stderr}")
        logger.error("Make sure you have g++ and OpenMP installed")
        return False
    
    logger.info("Successfully compiled ipSAE binary")
    return IPSAE_BINARY.exists()


def find_pae_file(model_id: str, scores_dir: Path) -> Optional[Path]:
    """Find the PAE JSON file for a model."""
    # Try different naming conventions
    patterns = [
        f"{model_id}-predicted_aligned_error_v1.json",
        f"{model_id}_pae.json",
        f"{model_id}_scores.json",
    ]
    
    for pattern in patterns:
        pae_file = scores_dir / pattern
        if pae_file.exists():
            return pae_file
    
    # Try to find meta file that might contain PAE
    meta_file = scores_dir / f"{model_id}-meta_v1.json"
    if meta_file.exists():
        return meta_file
    
    return None


def find_pdb_file(model_id: str, pdb_dir: Path) -> Optional[Path]:
    """Find the PDB file for a model."""
    patterns = [
        f"{model_id}-model_v1.pdb",
        f"{model_id}.pdb",
    ]
    
    for pattern in patterns:
        pdb_file = pdb_dir / pattern
        if pdb_file.exists():
            return pdb_file
    
    return None


def get_model_ids_from_pdb_dir(pdb_dir: Path) -> list[str]:
    """Extract model IDs from PDB files in directory."""
    model_ids = []
    for pdb_file in pdb_dir.glob("*-model_v1.pdb"):
        # Extract model ID from filename like AF-XXXX-model_v1.pdb
        model_id = pdb_file.stem.replace("-model_v1", "")
        model_ids.append(model_id)
    return sorted(model_ids)


def run_ipsae_batch(
    input_list_file: Path,
    output_dir: Path,
    pae_cutoff: float = 10.0,
    dist_cutoff: float = 10.0,
    num_threads: int = 1
) -> tuple[int, int, int]:
    """
    Run the C++ ipsae_optimized binary in batch mode.
    
    Returns:
        Tuple of (processed, skipped, errors)
    """
    if not ensure_ipsae_binary():
        raise FileNotFoundError(
            f"ipSAE binary not found and compilation failed.\n"
            f"  Expected at: {IPSAE_BINARY}\n"
            f"  To compile manually: cd {IPSAE_CPP_DIR} && make"
        )
    
    # Set OpenMP thread count
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(num_threads)
    
    cmd = [
        str(IPSAE_BINARY),
        "--batch",
        str(input_list_file),
        str(output_dir),
        "--pae-cutoff", str(pae_cutoff),
        "--dist-cutoff", str(dist_cutoff)
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env
    )
    
    if result.returncode != 0:
        logger.error(f"ipSAE binary failed: {result.stderr}")
        return (0, 0, -1)
    
    # Parse output for statistics
    processed = 0
    skipped = 0
    errors = 0
    
    for line in result.stdout.split('\n'):
        if "Processed:" in line:
            processed = int(line.split(':')[1].strip())
        elif "Skipped" in line:
            skipped = int(line.split(':')[1].strip())
        elif "Errors:" in line:
            errors = int(line.split(':')[1].strip())
    
    return (processed, skipped, errors)


def main():
    parser = argparse.ArgumentParser(
        description="Batch ipSAE calculator for protein complex interfaces"
    )
    parser.add_argument(
        "--pae-dir",
        type=Path,
        required=True,
        help="Directory containing PAE JSON files (scores directory)"
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        required=True,
        help="Directory containing PDB files (staging directory)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for ipSAE results"
    )
    parser.add_argument(
        "--model-ids",
        type=Path,
        help="Optional file with model IDs to process (one per line)"
    )
    parser.add_argument(
        "--pae-cutoff",
        type=float,
        default=10.0,
        help="PAE cutoff threshold (default: 10.0)"
    )
    parser.add_argument(
        "--dist-cutoff",
        type=float,
        default=10.0,
        help="Distance cutoff for interface detection (default: 10.0)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of OpenMP threads (default: all available CPUs)"
    )
    
    args = parser.parse_args()
    
    # Validate directories
    if not args.pae_dir.exists():
        logger.error(f"PAE directory does not exist: {args.pae_dir}")
        sys.exit(1)
    
    if not args.pdb_dir.exists():
        logger.error(f"PDB directory does not exist: {args.pdb_dir}")
        sys.exit(1)
    
    args.output_dir.mkdir(exist_ok=True, parents=True)
    
    # Get model IDs
    if args.model_ids and args.model_ids.exists():
        model_ids = [line.strip() for line in args.model_ids.read_text().splitlines() if line.strip()]
        logger.info(f"Loaded {len(model_ids)} model IDs from {args.model_ids}")
    else:
        model_ids = get_model_ids_from_pdb_dir(args.pdb_dir)
        logger.info(f"Found {len(model_ids)} models in PDB directory")
    
    if not model_ids:
        logger.warning("No models found to process")
        sys.exit(0)
    
    # Build input list file for batch mode
    # Format: model_id\tpae_file\tpdb_file
    entries = []
    missing_pae = []
    missing_pdb = []
    
    for model_id in model_ids:
        pae_file = find_pae_file(model_id, args.pae_dir)
        pdb_file = find_pdb_file(model_id, args.pdb_dir)
        
        if not pae_file:
            missing_pae.append(model_id)
            continue
        if not pdb_file:
            missing_pdb.append(model_id)
            continue
        
        entries.append(f"{model_id}\t{pae_file}\t{pdb_file}")
    
    if missing_pae:
        logger.warning(f"Missing PAE files for {len(missing_pae)} models")
    if missing_pdb:
        logger.warning(f"Missing PDB files for {len(missing_pdb)} models")
    
    if not entries:
        logger.error("No valid model entries found")
        sys.exit(1)
    
    logger.info(f"Processing {len(entries)} models with {args.workers} threads")
    
    # Write batch input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write('\n'.join(entries))
        input_list_file = Path(f.name)
    
    try:
        processed, skipped, errors = run_ipsae_batch(
            input_list_file,
            args.output_dir,
            args.pae_cutoff,
            args.dist_cutoff,
            args.workers
        )
        
        logger.info("")
        logger.info("=" * 40)
        logger.info("BATCH ipSAE COMPLETE")
        logger.info("=" * 40)
        logger.info(f"  Total models: {len(entries)}")
        logger.info(f"  Processed (multimer): {processed}")
        logger.info(f"  Skipped (monomer): {skipped}")
        logger.info(f"  Errors: {errors}")
        logger.info(f"  Output: {args.output_dir}")
        
        if errors > 0:
            sys.exit(1)
    
    finally:
        # Cleanup temp file
        input_list_file.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
