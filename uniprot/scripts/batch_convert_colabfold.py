#!/usr/bin/env python3
"""
Batch ColabFold converter - processes multiple models across worker processes.
Uses ProcessPoolExecutor with fork so all workers inherit the parent's
DuckDB prefetch cache via copy-on-write (zero-cost, read-only after population).
"""
import argparse
import csv
import logging
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Module-level variables shared across all workers (set before fork)
_WORKER_INPUT_DIR = None
_WORKER_OUTPUT_DIR = None
_WORKER_MANIFEST = None
_WORKER_DUCKDB = None
_WORKER_CHAIN_MANIFEST_DIR = None
_WORKER_MODEL_MANIFEST_DIR = None


def _find_input_files(model_id: str):
    """Find input files - check multiple possible locations."""
    input_dir = _WORKER_INPUT_DIR

    direct_meta = input_dir / f"{model_id}-meta_v1.json"
    direct_pdb = input_dir / f"{model_id}-model_v1.pdb"

    if direct_meta.exists() and direct_pdb.exists():
        return direct_meta, direct_pdb

    digits = model_id.replace("AF-", "")[:16]
    shard_parts = [digits[i:i+4] for i in range(0, 16, 4)]
    shard_dir = input_dir / "/".join(shard_parts)
    shard_meta = shard_dir / f"{model_id}-meta_v1.json"
    shard_pdb = shard_dir / f"{model_id}-model_v1.pdb"

    if shard_meta.exists() and shard_pdb.exists():
        return shard_meta, shard_pdb

    return None, None


def _process_single_model(model_id: str) -> tuple:
    """Process a single model."""
    from afdb_integration_kit.colabfold.converter import convert_file

    meta_json, pdb_file = _find_input_files(model_id)

    if not meta_json or not pdb_file:
        return (model_id, False, "Missing input files")

    try:
        convert_file(
            scores_json_path=str(meta_json),
            pdb_path=str(pdb_file),
            outdir=str(_WORKER_OUTPUT_DIR),
            manifest_path=_WORKER_MANIFEST,
            model_entity_id=model_id,
            duckdb_path=_WORKER_DUCKDB,
            chain_manifest_dir=_WORKER_CHAIN_MANIFEST_DIR,
            model_manifest_dir=_WORKER_MODEL_MANIFEST_DIR,
        )
        return (model_id, True, "")
    except Exception as e:
        return (model_id, False, str(e))


def main():
    global _WORKER_INPUT_DIR, _WORKER_OUTPUT_DIR, _WORKER_MANIFEST, _WORKER_DUCKDB
    global _WORKER_CHAIN_MANIFEST_DIR, _WORKER_MODEL_MANIFEST_DIR

    parser = argparse.ArgumentParser(description="Batch convert ColabFold outputs to AFDB format")
    parser.add_argument("--model-ids-file", required=True, help="File containing model IDs, one per line")
    parser.add_argument("--input-dir", required=True, help="Directory containing input files (meta_v1.json, model_v1.pdb)")
    parser.add_argument("--manifest", required=True, help="ColabFold manifest CSV")
    parser.add_argument("--duckdb", help="Optional DuckDB database path")
    parser.add_argument("--output-dir", required=True, help="Output directory for converted files")
    parser.add_argument("--chain-manifest-dir", required=True, help="Directory for chain manifests")
    parser.add_argument("--model-manifest-dir", required=True, help="Directory for model manifests")
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 8, help="Number of parallel workers (default: all available CPUs)")
    parser.add_argument("--failed-ids-file", help="If set, append failed model IDs and errors to this TSV file")
    parser.add_argument("--stage-name", default="stage_03_convert_colabfold", help="Stage label for the failed-ids file")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    chain_manifest_dir = Path(args.chain_manifest_dir)
    model_manifest_dir = Path(args.model_manifest_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    chain_manifest_dir.mkdir(parents=True, exist_ok=True)
    model_manifest_dir.mkdir(parents=True, exist_ok=True)

    model_ids_file = Path(args.model_ids_file)
    model_ids = [line.strip() for line in model_ids_file.read_text().splitlines() if line.strip()]

    if not model_ids:
        logger.error("No model IDs found in %s", model_ids_file)
        sys.exit(1)

    logger.info("Processing %d models with %d workers (ProcessPool/fork)", len(model_ids), args.workers)

    # Ensure repo root is on sys.path for converter imports
    repo_root = str(Path(__file__).resolve().parent.parent.parent)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    # Set module-level config before forking so children inherit via COW
    _WORKER_INPUT_DIR = input_dir
    _WORKER_OUTPUT_DIR = output_dir
    _WORKER_MANIFEST = args.manifest
    _WORKER_DUCKDB = args.duckdb
    _WORKER_CHAIN_MANIFEST_DIR = str(chain_manifest_dir)
    _WORKER_MODEL_MANIFEST_DIR = str(model_manifest_dir)

    # Prefetch DuckDB metadata ONCE in the parent process; forked children
    # inherit the populated cache via copy-on-write (read-only, zero-copy).
    if args.duckdb:
        all_accessions = []
        manifest_path = Path(args.manifest)
        with manifest_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                acc = row.get("uniprot_ac") or row.get("uniprotAccession")
                if acc:
                    all_accessions.append(acc.strip())

        if all_accessions:
            logger.info("Prefetching DuckDB metadata for %d unique accessions...", len(set(all_accessions)))
            from afdb_integration_kit.colabfold.converter import prefetch_duckdb_metadata
            prefetch_duckdb_metadata(args.duckdb, all_accessions)

    # Pre-import converter in parent so forked children skip the import overhead
    from afdb_integration_kit.colabfold import converter as _  # noqa: F811,F401

    success_count = 0
    error_count = 0
    failed_entries: list[tuple[str, str]] = []

    fork_ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=fork_ctx) as executor:
        futures = {executor.submit(_process_single_model, mid): mid for mid in model_ids}
        for future in as_completed(futures):
            model_id, success, error_msg = future.result()
            if success:
                success_count += 1
                if success_count % 100 == 0:
                    logger.info("Processed %d/%d models", success_count, len(model_ids))
            else:
                error_count += 1
                failed_entries.append((model_id, error_msg))
                logger.error("Failed to convert %s: %s", model_id, error_msg)

    logger.info("Batch conversion complete: %d success, %d errors", success_count, error_count)

    if failed_entries and args.failed_ids_file:
        failed_path = Path(args.failed_ids_file)
        failed_path.parent.mkdir(parents=True, exist_ok=True)
        with failed_path.open("a") as fh:
            for model_id, error_msg in sorted(failed_entries):
                fh.write(f"{model_id}\t{args.stage_name}\t{error_msg}\n")
        logger.info("Appended %d failed IDs to %s", len(failed_entries), failed_path)

    if error_count > 0 and success_count == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
