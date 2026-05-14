#!/usr/bin/env python3
"""Prepare ColabFold outputs for the AFDB production pipeline.

Scans a directory of ColabFold PDB + scores-JSON pairs, builds all the config
files the ``production_pipeline.py`` script expects, and symlinks / copies the
input assets into the canonical layout.

By default the ColabFold scores file is symlinked (or copied) directly as the
``*-meta_v1.json`` the pipeline expects -- no JSON parsing involved.  Pass
``--extract-meta`` to explicitly parse and re-write a leaner meta JSON using
``orjson`` (falls back to the stdlib ``json`` when ``orjson`` is unavailable).

Supports **two operating modes**:

Production (fast, 30M-scale)
    Supply ``--chain-mapping`` and ``--uniprot-db`` to skip manifest
    resolution and UniProt fetching entirely.

Dev / small-scale
    Omit those flags and provide ``--build-from-api`` instead.  The script
    streams the AFCDB manifest once to resolve AF IDs to UniProt accessions,
    fetches full UniProt entries via the REST API, and builds a DuckDB.

Works for homodimers, heterodimers, and homomultimers without any PDB
file parsing -- chain-to-accession mapping is derived purely from model-ID
naming conventions and the AFCDB manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import orjson  # 5-10x faster for numerical JSON
except ImportError:
    orjson = None  # type: ignore[assignment]

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Default ColabFold output suffixes
DEFAULT_SCORE_SUFFIX = (
    ".merged_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
)
DEFAULT_PDB_SUFFIX = (
    ".merged_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
)

MAX_IO_WORKERS = 16


# ---------------------------------------------------------------------------
# 1. Scan for matched PDB + JSON pairs  (directory listing only)
# ---------------------------------------------------------------------------

def find_model_pairs(
    source_dir: Path,
    score_suffix: str,
    pdb_suffix: str,
) -> Dict[str, Tuple[Path, Path]]:
    """Find matched score-JSON / PDB pairs.  Returns ``{model_id: (score, pdb)}``."""
    scores: Dict[str, Path] = {}
    pdbs: Dict[str, Path] = {}

    for entry in os.scandir(source_dir):
        if not entry.is_file():
            continue
        name = entry.name
        if name.endswith(score_suffix):
            model_id = name[: -len(score_suffix)]
            scores[model_id] = Path(entry.path)
        elif name.endswith(pdb_suffix):
            model_id = name[: -len(pdb_suffix)]
            pdbs[model_id] = Path(entry.path)

    matched_ids = scores.keys() & pdbs.keys()
    return {mid: (scores[mid], pdbs[mid]) for mid in matched_ids}


# ---------------------------------------------------------------------------
# 2. Meta-JSON extraction  (only used when --extract-meta is set)
# ---------------------------------------------------------------------------

REQUIRED_META_KEYS = ("plddt", "pae", "max_pae")
OPTIONAL_META_KEYS = ("ptm", "iptm")


def extract_meta_json(score_path: Path, output_path: Path) -> None:
    """Parse scores and write only the keys the pipeline needs.

    Uses ``orjson`` when available (5-10x faster for numerical JSON),
    falling back to the stdlib ``json`` module otherwise.
    """
    if orjson is not None:
        raw = score_path.read_bytes()
        data = orjson.loads(raw)
        meta: Dict[str, Any] = {k: data[k] for k in REQUIRED_META_KEYS}
        for k in OPTIONAL_META_KEYS:
            if k in data:
                meta[k] = data[k]
        output_path.write_bytes(orjson.dumps(meta))
        return

    with score_path.open() as fh:
        data = json.load(fh)
    meta = {k: data[k] for k in REQUIRED_META_KEYS}
    for k in OPTIONAL_META_KEYS:
        if k in data:
            meta[k] = data[k]
    with output_path.open("w") as fh:
        json.dump(meta, fh)


# ---------------------------------------------------------------------------
# 3. File preparation  (symlink / copy PDB + meta JSON)
# ---------------------------------------------------------------------------

def _symlink_or_replace(source: Path, dest: Path) -> None:
    """Create a symlink at *dest* pointing to the resolved *source*."""
    if dest.is_symlink() or dest.exists():
        dest.unlink()
    os.symlink(source.resolve(), dest)


def _process_single_model(
    model_id: str,
    score_path: Path,
    pdb_path: Path,
    output_dir: Path,
    use_symlinks: bool,
    extract_meta: bool,
) -> Optional[str]:
    """Prepare a single model's PDB and meta-JSON files.

    Returns *model_id* on failure, ``None`` on success.

    Meta-JSON strategy (fastest first):
      * symlink mode (default): symlink the scores file as meta JSON -- O(1).
      * copy mode (``--copy``): ``shutil.copy2`` the scores file -- no parsing.
      * extract mode (``--extract-meta``): parse + rewrite via orjson/json.
    """
    try:
        dest_meta = output_dir / f"{model_id}-meta_v1.json"
        dest_pdb = output_dir / f"{model_id}-model_v1.pdb"

        if extract_meta:
            extract_meta_json(score_path, dest_meta)
        elif use_symlinks:
            _symlink_or_replace(score_path, dest_meta)
        else:
            shutil.copy2(score_path, dest_meta)

        if use_symlinks:
            _symlink_or_replace(pdb_path, dest_pdb)
        else:
            shutil.copy2(pdb_path, dest_pdb)

        return None
    except Exception as exc:
        logger.warning("Failed to process %s: %s", model_id, exc)
        return model_id


def prepare_input_files(
    model_ids: List[str],
    pairs: Dict[str, Tuple[Path, Path]],
    output_dir: Path,
    use_symlinks: bool = True,
    extract_meta: bool = False,
) -> List[str]:
    """Prepare PDB + meta-JSON for all models.  Returns list of failed model IDs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    failed: List[str] = []

    with ThreadPoolExecutor(max_workers=MAX_IO_WORKERS) as executor:
        futures = {
            executor.submit(
                _process_single_model,
                mid, pairs[mid][0], pairs[mid][1],
                output_dir, use_symlinks, extract_meta,
            ): mid
            for mid in model_ids
        }
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                failed.append(result)

    return failed


# ---------------------------------------------------------------------------
# 4. Config file writers
# ---------------------------------------------------------------------------

def write_manifest_csv(path: Path, rows: List[Dict[str, str]]) -> None:
    """Write the ColabFold manifest CSV."""
    fieldnames = ["model_entity_id", "chain_id", "uniprot_ac"]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_mapping_tsv(path: Path, model_ids: List[str]) -> None:
    """Write the model-ID mapping file (one ID per line)."""
    with path.open("w") as fh:
        for mid in model_ids:
            fh.write(f"{mid}\n")


def write_dataset_config(path: Path, provider_id: str) -> None:
    """Write the dataset configuration JSON."""
    config = {
        "providerId": provider_id,
        "toolUsed": "AlphaFold",
        "latestVersion": 1,
        "allVersions": [1],
        "entityType": "protein",
        "modelCreatedDate": "2026-01-01T00:00:00Z",
        "uniqueIdTemplate": "{model_entity_id}",
        "versionTag": "v1",
    }
    path.write_text(json.dumps(config, indent=2))


def write_provider_json(path: Path, provider_id: str, provider_name: str) -> None:
    """Write the provider configuration JSON."""
    provider = {
        "providerId": provider_id,
        "providerName": provider_name,
        "providerUrl": "https://alphafold.ebi.ac.uk",
        "copyrights": [
            f"Copyright 2026 {provider_name}. All rights reserved."
        ],
    }
    path.write_text(json.dumps(provider, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare ColabFold outputs for the AFDB production pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Production (30M models): pre-built manifest + DuckDB
  %(prog)s \\
    --input-dir /data/colabfold/gpu0 \\
    --output-dir /data/workdir \\
    --chain-mapping /data/prebuilt_manifest.csv \\
    --uniprot-db /data/uniprot_2025_04.duckdb \\
    --provider-id afcdb-heterodimers \\
    --provider-name "AFCDB Heterodimers"

  # Dev (small-scale): auto-resolve from AFCDB manifest + API
  %(prog)s \\
    --input-dir ./gpu0 \\
    --output-dir ./workdir \\
    --build-from-api /data/afdb_toolkit_manifest_file.csv \\
    --provider-id afcdb-heterodimers \\
    --provider-name "AFCDB Heterodimers"
        """,
    )

    # -- Required --
    parser.add_argument(
        "--input-dir", type=Path, required=True,
        help="Directory containing ColabFold PDB and scores-JSON files.",
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory (will contain inputs/ and config/ sub-dirs).",
    )
    parser.add_argument(
        "--provider-id", required=True,
        help="Dataset provider ID (e.g. 'afcdb-heterodimers').",
    )
    parser.add_argument(
        "--provider-name", required=True,
        help="Dataset provider display name (e.g. 'AFCDB Heterodimers').",
    )

    # -- Production mode: pre-built assets --
    prod = parser.add_argument_group("production mode (pre-built assets)")
    prod.add_argument(
        "--chain-mapping", type=Path, default=None,
        help="Pre-built chain-to-accession mapping CSV (model_entity_id, chain_id, uniprot_ac).  Skips accession resolution.",
    )
    prod.add_argument(
        "--uniprot-db", type=Path, default=None,
        help="Pre-built UniProt DuckDB.  Skips REST API fetching.",
    )

    # -- Dev mode: resolve from AFCDB manifest --
    dev = parser.add_argument_group("dev mode (auto-resolve)")
    dev.add_argument(
        "--build-from-api", type=Path, default=None,
        help="Path to afdb_toolkit_manifest_file.csv.  Streams the manifest to "
             "resolve AF-IDs, then fetches full UniProt entries via REST API.",
    )
    dev.add_argument(
        "--uniprot-release", default="2025_01",
        help="UniProt release tag for API-fetched entries (default: %(default)s).",
    )

    # -- Optional tuning --
    parser.add_argument(
        "--pdb-suffix", default=DEFAULT_PDB_SUFFIX,
        help="PDB file suffix pattern (default: ColabFold multimer v3).",
    )
    parser.add_argument(
        "--score-suffix", default=DEFAULT_SCORE_SUFFIX,
        help="Scores JSON file suffix pattern (default: ColabFold multimer v3).",
    )
    parser.add_argument(
        "--copy", action="store_true",
        help="Copy files instead of symlinking (default: symlink).",
    )
    parser.add_argument(
        "--extract-meta", action="store_true",
        help=(
            "Parse each scores JSON and re-write a leaner meta JSON "
            "(uses orjson when available).  By default the scores file is "
            "symlinked/copied directly -- much faster at scale."
        ),
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validate argument combinations."""
    has_manifest = args.chain_mapping is not None
    has_afcdb = args.build_from_api is not None

    if not has_manifest and not has_afcdb:
        print(
            "ERROR: Provide either --chain-mapping (production) or "
            "--build-from-api (dev) to resolve chain-to-accession mapping.",
            file=sys.stderr,
        )
        sys.exit(1)

    if not has_manifest and has_afcdb and not args.build_from_api.exists():
        print(
            f"ERROR: AFCDB manifest not found: {args.build_from_api}",
            file=sys.stderr,
        )
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_arguments()
    validate_args(args)

    inputs_dir = args.output_dir / "inputs"
    config_dir = args.output_dir / "config"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Scan for matched PDB + JSON pairs ----
    logger.info("Scanning %s for matched pairs...", args.input_dir)
    pairs = find_model_pairs(args.input_dir, args.score_suffix, args.pdb_suffix)
    model_ids = sorted(pairs.keys())
    logger.info("Found %d matched pairs", len(pairs))

    if not model_ids:
        logger.error("No matched PDB + JSON pairs found.  Check --pdb-suffix / --score-suffix.")
        return 1

    # ---- 2. Resolve ColabFold manifest ----
    if args.chain_mapping:
        # Production mode: copy the pre-built manifest
        logger.info("Using pre-built ColabFold manifest: %s", args.chain_mapping)
        manifest_dest = config_dir / "colabfold_manifest.csv"
        shutil.copy2(args.chain_mapping, manifest_dest)

        # Read it to filter model_ids to only those present in the manifest
        manifest_model_ids = set()
        with manifest_dest.open(newline="") as fh:
            for row in csv.DictReader(fh):
                manifest_model_ids.add(row["model_entity_id"])
        valid_ids = [m for m in model_ids if m in manifest_model_ids]
        skipped = len(model_ids) - len(valid_ids)
        if skipped:
            logger.warning(
                "%d models not in manifest (skipped); %d remaining",
                skipped, len(valid_ids),
            )
        model_ids = valid_ids
    else:
        # Dev mode: resolve from AFCDB manifest
        # Import here to keep production mode dependency-free
        repo_root = Path(__file__).resolve().parent.parent
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        from afdb_integration_kit.manifest.resolver import resolve_and_build_manifest

        logger.info("Resolving accessions from AFCDB manifest...")
        manifest_rows, skipped, af_to_uniprot = resolve_and_build_manifest(
            model_ids, args.build_from_api,
        )
        if skipped:
            skipped_path = args.output_dir / "skipped_models.txt"
            skipped_path.write_text("\n".join(sorted(skipped)) + "\n")
            logger.warning(
                "%d models skipped (see %s); %d remaining",
                len(skipped), skipped_path, len(model_ids) - len(skipped),
            )
            skipped_set = set(skipped)
            model_ids = [m for m in model_ids if m not in skipped_set]

        write_manifest_csv(config_dir / "colabfold_manifest.csv", manifest_rows)
        logger.info("Wrote %d manifest rows", len(manifest_rows))

    if not model_ids:
        logger.error("No valid models remaining after manifest resolution.")
        return 1

    # ---- 3. Resolve UniProt DuckDB ----
    duckdb_dest = config_dir / "uniprot.duckdb"
    if args.uniprot_db:
        # Production mode: copy / symlink the pre-built DuckDB
        logger.info("Using pre-built UniProt DuckDB: %s", args.uniprot_db)
        if duckdb_dest.exists() or duckdb_dest.is_symlink():
            duckdb_dest.unlink()
        os.symlink(args.uniprot_db.resolve(), duckdb_dest)
    else:
        # Dev mode: fetch via API and build DuckDB
        repo_root = Path(__file__).resolve().parent.parent
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        from afdb_integration_kit.uniprot.api import (
            build_duckdb,
            entries_to_parquet,
            fetch_entries,
        )

        # Collect unique accessions from the manifest
        manifest_path = config_dir / "colabfold_manifest.csv"
        accessions: set = set()
        with manifest_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                accessions.add(row["uniprot_ac"])
        unique_accessions = sorted(accessions)
        logger.info("Fetching %d unique UniProt accessions via API...", len(unique_accessions))

        entries = fetch_entries(unique_accessions, release=args.uniprot_release)
        fetched = {e["primary_ac"] for e in entries}
        missing = set(unique_accessions) - fetched
        if missing:
            logger.warning(
                "%d accessions not found in UniProt: %s",
                len(missing), sorted(missing)[:10],
            )
        logger.info("Fetched %d / %d entries", len(entries), len(unique_accessions))

        parquet_path = config_dir / "entry.parquet"
        entries_to_parquet(entries, parquet_path)
        build_duckdb(parquet_path, duckdb_dest)

    # ---- 4. Write remaining config files ----
    logger.info("Writing config files...")
    write_mapping_tsv(config_dir / "af_mapping.tsv", model_ids)
    write_dataset_config(config_dir / "dataset_config.json", args.provider_id)
    write_provider_json(config_dir / "provider.json", args.provider_id, args.provider_name)

    # ---- 5. Prepare input files (symlink/copy PDB + meta JSON) ----
    use_symlinks = not args.copy
    meta_strategy = (
        "extract (orjson)" if args.extract_meta and orjson is not None
        else "extract (json stdlib)" if args.extract_meta
        else "symlink" if use_symlinks
        else "copy"
    )
    logger.info(
        "Preparing %d input files (meta strategy: %s)...",
        len(model_ids), meta_strategy,
    )
    failed = prepare_input_files(
        model_ids, pairs, inputs_dir, use_symlinks, args.extract_meta,
    )
    if failed:
        logger.warning("%d models failed during file preparation", len(failed))
        failed_set = set(failed)
        model_ids = [m for m in model_ids if m not in failed_set]

    # ---- 6. Summary ----
    summary = {
        "total_matched_pairs": len(pairs),
        "prepared_models": len(model_ids),
        "failed_models": len(failed) if failed else 0,
        "mode": "production" if args.chain_mapping else "dev",
    }
    summary_path = args.output_dir / "build_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    logger.info("Done. %d models prepared.", len(model_ids))
    logger.info("Summary: %s", json.dumps(summary))
    logger.info("")
    logger.info("Next step:")
    logger.info("  python3 scripts/production_pipeline.py \\")
    logger.info("    --input-dir %s \\", inputs_dir)
    logger.info("    --output-dir %s/output \\", args.output_dir)
    logger.info("    --mapping-file %s/af_mapping.tsv \\", config_dir)
    logger.info("    --chain-mapping %s/colabfold_manifest.csv \\", config_dir)
    logger.info("    --dataset-config %s/dataset_config.json \\", config_dir)
    logger.info("    --provider-json %s/provider.json \\", config_dir)
    logger.info("    --uniprot-db %s/uniprot.duckdb", config_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
