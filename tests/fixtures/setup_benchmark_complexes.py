#!/usr/bin/env python3
"""
Set up benchmark environment using the 100 real complexes from examples/complexes/.

This script:
1. Unzips and decompresses the complexes (files are gzipped inside the zips)
2. Copies the ColabFold manifest
3. Copies the pre-built DuckDB
4. Creates Nextflow config files
5. Creates af_mapping.tsv

Usage:
    python tests/fixtures/setup_benchmark_complexes.py --output-dir /tmp/benchmark_complexes --clean
"""

import argparse
import gzip
import json
import logging
import shutil
import zipfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def find_repo_root() -> Path:
    """Find the AFDB-Integration-Kit repository root."""
    current = Path(__file__).resolve()
    # Look for the directory containing this script's expected path
    for parent in [current] + list(current.parents):
        if parent.name == "AFDB-Integration-Kit" and (parent / "examples" / "complexes").exists():
            return parent
    # Fallback: look for .git with workflow directory
    for parent in [current] + list(current.parents):
        if (parent / ".git").exists() and (parent / "workflow").exists():
            return parent
    raise RuntimeError("Could not find AFDB-Integration-Kit repository root")


def unzip_and_decompress_complexes(repo_root: Path, output_dir: Path) -> int:
    """
    Unzip the complexes from examples/complexes/ and decompress gzipped files.

    Returns the number of model files found.
    """
    complexes_dir = repo_root / "examples" / "complexes"
    output_complexes = output_dir / "input" / "complexes"
    output_complexes.mkdir(parents=True, exist_ok=True)

    zip_files = list(complexes_dir.glob("complexes_*.zip"))
    if not zip_files:
        raise FileNotFoundError(f"No complexes zip files found in {complexes_dir}")

    logger.info(f"Found {len(zip_files)} zip files to extract")

    # Extract all zips
    for zip_path in zip_files:
        logger.info(f"Extracting {zip_path.name}...")
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(output_complexes)

    # Decompress gzipped files
    decompressed = 0
    for path in list(output_complexes.glob("*.json")) + list(output_complexes.glob("*.pdb")):
        try:
            content = path.read_bytes()
            if content[:2] == b'\x1f\x8b':  # gzip magic bytes
                with gzip.open(path, 'rt', errors='ignore') as f:
                    text = f.read()
                path.write_text(text)
                decompressed += 1
        except Exception as e:
            logger.warning(f"Failed to decompress {path.name}: {e}")

    logger.info(f"Decompressed {decompressed} gzipped files")

    # Count model files
    model_count = len(list(output_complexes.glob("*-model_v1.pdb")))
    logger.info(f"Found {model_count} model files")

    return model_count


def copy_manifest(repo_root: Path, output_dir: Path) -> Path:
    """Copy the ColabFold manifest to the config directory."""
    manifest_src = repo_root / "colabfold_manifest_100.csv"
    if not manifest_src.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_src}")

    config_dir = output_dir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    manifest_dst = config_dir / "colabfold_manifest.csv"
    shutil.copy(manifest_src, manifest_dst)
    logger.info(f"Copied manifest to {manifest_dst}")

    return manifest_dst


def copy_duckdb(repo_root: Path, output_dir: Path) -> Path:
    """Copy the pre-built DuckDB to the output directory."""
    duckdb_src = repo_root / "uniprot_complexes_100.duckdb"
    if not duckdb_src.exists():
        # Try to unzip if the zip exists
        duckdb_zip = repo_root / "uniprot_complexes_100.zip"
        if duckdb_zip.exists():
            logger.info(f"Unzipping {duckdb_zip}...")
            with zipfile.ZipFile(duckdb_zip, 'r') as zf:
                zf.extractall(repo_root)
        else:
            raise FileNotFoundError(
                f"DuckDB not found: {duckdb_src}\n"
                "Run: python tests/fixtures/fetch_uniprot_subset.py to create it"
            )

    duckdb_dst = output_dir / "uniprot_complexes.duckdb"
    shutil.copy(duckdb_src, duckdb_dst)
    logger.info(f"Copied DuckDB to {duckdb_dst}")

    return duckdb_dst


def create_af_mapping(manifest_path: Path, output_dir: Path) -> Path:
    """Create af_mapping.tsv with unique model IDs from the manifest."""
    config_dir = output_dir / "config"

    # Read manifest and extract unique model IDs
    model_ids = set()
    with open(manifest_path) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(",")
            if parts:
                model_ids.add(parts[0])

    # Write af_mapping.tsv
    af_mapping = config_dir / "af_mapping.tsv"
    with open(af_mapping, "w") as f:
        for model_id in sorted(model_ids):
            f.write(f"{model_id}\n")

    logger.info(f"Created {af_mapping} with {len(model_ids)} model IDs")
    return af_mapping


def create_config_files(output_dir: Path) -> None:
    """Create Nextflow config files."""
    config_dir = output_dir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    # nextflow_local.config
    nextflow_config = config_dir / "nextflow_local.config"
    nextflow_config.write_text("""\
// Local test configuration for AFDB Integration Kit
process {
    executor = 'local'
    cpus = 2
    memory = '4 GB'
}
docker.enabled = false
singularity.enabled = false
""")
    logger.info(f"Created {nextflow_config}")

    # dataset_config.json (required fields from batch_export_metadata.py)
    dataset_config = config_dir / "dataset_config.json"
    dataset_config.write_text(json.dumps({
        "providerId": "colabfold",
        "toolUsed": "ColabFold",
        "latestVersion": 1,
        "allVersions": [1],
        "entityType": "protein",
        "datasetName": "benchmark_complexes_100",
        "description": "100 protein complexes for benchmarking"
    }, indent=2))
    logger.info(f"Created {dataset_config}")

    # provider.json
    provider_json = config_dir / "provider.json"
    provider_json.write_text(json.dumps({
        "providerId": "colabfold",
        "providerName": "ColabFold",
        "providerUrl": "https://colabfold.com",
        "copyrights": ["ColabFold structure predictions"],
        "description": "ColabFold protein structure predictions"
    }, indent=2))
    logger.info(f"Created {provider_json}")


def print_run_command(output_dir: Path, repo_root: Path) -> None:
    """Print the Nextflow command to run the pipeline."""
    print("\n" + "=" * 70)
    print("BENCHMARK SETUP COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print(f"\nTo run the pipeline:\n")
    print(f"""\
nextflow run "{repo_root}/workflow/end_to_end_with_validation_multibatch.nf" \\
    -c "{output_dir}/config/nextflow_local.config" \\
    --posix_root "{output_dir}/input/complexes" \\
    --mapping_file "{output_dir}/config/af_mapping.tsv" \\
    --colabfold_manifest "{output_dir}/config/colabfold_manifest.csv" \\
    --dataset_dir "{output_dir}" \\
    --dataset_config "{output_dir}/config/dataset_config.json" \\
    --provider_json "{output_dir}/config/provider.json" \\
    --archive_root "{output_dir}/input/complexes" \\
    --converter_duckdb "{output_dir}/uniprot_complexes.duckdb" \\
    --uniprot_db "{output_dir}/uniprot_complexes.duckdb" \\
    --toolkit_cmd "python {repo_root}/main.py" \\
    --sample_size 0 \\
    -work-dir "{output_dir}/work" \\
    -with-report "{output_dir}/report.html" \\
    -with-timeline "{output_dir}/timeline.html"
""")


def main():
    parser = argparse.ArgumentParser(
        description="Set up benchmark environment with 100 real complexes"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/benchmark_complexes"),
        help="Output directory for benchmark data (default: /tmp/benchmark_complexes)"
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing output directory before setup"
    )
    args = parser.parse_args()

    repo_root = find_repo_root()
    output_dir = args.output_dir.resolve()

    logger.info(f"Repository root: {repo_root}")
    logger.info(f"Output directory: {output_dir}")

    if args.clean and output_dir.exists():
        logger.info(f"Removing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Unzip and decompress complexes
    logger.info("\n[Step 1/5] Extracting and decompressing complexes...")
    model_count = unzip_and_decompress_complexes(repo_root, output_dir)

    # Step 2: Copy manifest
    logger.info("\n[Step 2/5] Copying ColabFold manifest...")
    manifest_path = copy_manifest(repo_root, output_dir)

    # Step 3: Copy DuckDB
    logger.info("\n[Step 3/5] Copying UniProt DuckDB...")
    copy_duckdb(repo_root, output_dir)

    # Step 4: Create af_mapping.tsv
    logger.info("\n[Step 4/5] Creating af_mapping.tsv...")
    create_af_mapping(manifest_path, output_dir)

    # Step 5: Create config files
    logger.info("\n[Step 5/5] Creating Nextflow config files...")
    create_config_files(output_dir)

    # Print run command
    print_run_command(output_dir, repo_root)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
