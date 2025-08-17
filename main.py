import logging
from pathlib import Path
import json

import typer

from afdb_integration_kit.cif2bcif.convert import (
    run_batch_cif2bcif,
)
from afdb_integration_kit.cif2bcif.convert import run_cif2bcif as cif2bcif_helper
from afdb_integration_kit.dssp.dssp import run_dssp as dssp_helper
from afdb_integration_kit.metadata.validator import validate_against_schema
from afdb_integration_kit.modelcif.generate import generate
from afdb_integration_kit.modelpdb.generate import generate_pdb_headers
from afdb_integration_kit.modelcif_replace.replace import replace_mmcif_with_json as replace_mmcif_with_json
from afdb_integration_kit.quality_assessment.naming import (
    validate_dataset_naming,
    format_human
)

# Set up logger
logger = logging.getLogger("afdb_integration_kit")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)

app = typer.Typer()


@app.command()
def test():
    """
    Runs a series of checks to verify the environment and toolchain.
    """
    import subprocess

    logger.info("--- Verifying Versions ---")
    subprocess.run(["python", "--version"])
    subprocess.run(["node", "--version"])
    subprocess.run(["npm", "--version"])
    logger.info("--- Testing molstar cif2bcif script ---")
    subprocess.run(["cif2bcif", "-h"])
    logger.info("--- Testing DSSP script ---")
    subprocess.run(["mkdssp", "-h"])
    logger.info("--- Testing Gemmi script ---")
    subprocess.run(["gemmi", "--version"])


@app.command()
def run_schema_validation(
    input_file: Path = typer.Option(
        ...,
        "-i",
        "--input",
        help="Input JSON file path to validate against schema.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    type: str = typer.Option(
        ...,
        "-t",
        "--type",
        help="Type of schema to validate against ('model' or 'provider').",
    ),
):
    """
    Validate a JSON file against the specified schema.
    """

    validate_against_schema(input_file, type)


@app.command()
def run_modelcif_gen(
    pdb: Path = typer.Option(
        ...,
        "-p",
        "--pdb",
        help="Input PDB file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    metadata: Path = typer.Option(
        ...,
        "-m",
        "--metadata",
        help="Input metadata JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output mmCIF file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    validate: str = typer.Option(
        None,
        "--validate",
        help="Optionally validate the output CIF file against a ModelCIF dictionary."
        "If used as a flag without a path, it defaults to 'mmcif_ma.dic'.",
    ),
):
    """
    Enrich a PDB file with metadata to produce a feature-rich mmCIF file.
    """
    # Pre-flight checks
    if not pdb.is_file():
        logger.error(f"Input PDB file not found: {pdb}")
        raise typer.Exit(code=1)
    if not metadata.is_file():
        logger.error(f"Input metadata JSON file not found: {metadata}")
        raise typer.Exit(code=1)

    # If validate is passed as a flag (True but no value), default to 'mmcif_ma.dic'
    validate_path = validate if validate else None
    if validate is not None and validate == "":
        validate_path = "mmcif_ma.dic"

    # Call main logic (assuming main is imported or defined elsewhere)
    generate(str(pdb), str(metadata), str(output), validate_path)

@app.command()
def run_modelpdb_gen(
    cif: Path = typer.Option(
        ...,
        "-c",
        "--cif",
        help="Input mmCIF file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    pdb: Path = typer.Option(
        ...,
        "-p",
        "--pdb",
        help="Input PDB file path (with ATOM coordinates).",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    provider: Path = typer.Option(
        ...,
        "-r",
        "--provider",
        help="Input provider JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output PDB file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
):
    """
    Enrich a PDB file with header information from a mmCIF file.
    """
    # Pre-flight checks
    if not cif.is_file():
        logger.error(f"Input mmCIF file not found: {cif}")
        raise typer.Exit(code=1)
    if not pdb.is_file():
        logger.error(f"Input PDB file not found: {pdb}")
        raise typer.Exit(code=1)
    if not provider.is_file():
        logger.error(f"Input provider JSON file not found: {provider}")
        raise typer.Exit(code=1)

    # Call main logic
    generate_pdb_headers(str(cif), str(pdb), str(output), str(provider))


@app.command()
def run_dssp(
    input_file: Path = typer.Option(
        ...,
        "-i",
        "--input",
        help="Input file in CIF format.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
    ),
    output_file: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output file in CIF format.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        readable=False,
        resolve_path=True,
    ),
):
    """
    Run DSSP on a CIF file to generate secondary structure information.
    """
    logger.info(f"Converting {input_file} to {output_file}")
    dssp_helper(input_file, output_file)
    logger.info("Conversion complete.")


@app.command()
def run_cif2bcif(
    input_file: Path = typer.Option(
        ...,
        "-i",
        "--input",
        help="Input file in CIF format.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
    ),
    output_file: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output file in BCIF or BCIF.GZ format.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        readable=False,
        resolve_path=True,
    ),
):
    """
    Convert CIF to BinaryCIF or BinaryCIF.GZ
    """
    logger.info(f"Converting {input_file} to {output_file}")
    cif2bcif_helper(input_file, output_file)
    logger.info("Conversion complete.")


@app.command()
def batch_cif2bcif(
    input_dir: Path = typer.Option(
        ...,
        "--input-dir",
        "-id",
        help="Input directory containing CIF files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-od",
        help="Output directory for BCIF or BCIF.GZ files.",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    workers: int = typer.Option(
        4, "--workers", "-w", help="Number of parallel workers (default: 4)"
    ),
    gzip: bool = typer.Option(
        False, "--gzip", "-gz", help="Output .bcif.gz files instead of .bcif"
    ),
):
    """
    Batch process all CIF files in a directory to BCIF or BCIF.GZ.
    """
    logger.info(
        f"Batch converting CIF files from {input_dir} to {output_dir} "
        f"using {workers} workers. Gzip: {gzip}"
    )
    run_batch_cif2bcif(input_dir, output_dir, workers=workers, gzip=gzip)
    logger.info("Batch conversion complete.")

@app.command()
def run_naming_check(
    root: Path = typer.Option(
        ...,
        "--root",
        help="Path to dataset directory to scan.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    out: Path = typer.Option(
        None,
        "-o",
        "--out",
        help="Optional path to write a human-readable report.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    ids_file: Path = typer.Option(
        None,
        "--ids-file",
        help="Optional file containing AFIDs to check, one per line. Accepts AF-<16>-vN or AF-<16>.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    errors_only: bool = typer.Option(
        False,
        "--errors-only",
        help="Print only failing entries as one-liners.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Print every entry as one-liners.",
    ),
    limit: int = typer.Option(
        5,
        "--limit",
        min=1,
        help="Limit example AFIDs in summary and cap one-liner output when not in --verbose mode.",
    ),
):
    """
    Check AFDB naming conventions and required-file presence for a dataset directory.
    Default output is a HUMAN-READABLE summary. Exit code 0 on PASS, 1 if any issue is found.
    """
    # Read optional ids-file
    ids_pairs = None
    ids_afids = None
    if ids_file is not None:
        from afdb_integration_kit.quality_assessment.naming import parse_ids_file
        ids_pairs, ids_afids = parse_ids_file(ids_file)

    # Validate and build report
    ok, report = validate_dataset_naming(root, ids_pairs=ids_pairs, ids_afids=ids_afids)

    # Human-readable to stdout
    text = format_human(report, errors_only=errors_only, verbose=verbose, limit=limit)
    print(text)

    # Optional save of human-readable
    if out is not None:
        out.write_text(text + "\n")
        logger.info(f"Wrote report: {out}")

    if not ok:
        raise typer.Exit(code=1)


@app.command()
def run_modelcif_replace(
    pdb: Path = typer.Option(
        ...,
        "-i",
        "--cif",
        help="Input CIF file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    metadata: Path = typer.Option(
        ...,
        "-m",
        "--metadata",
        help="Input metadata JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output mmCIF file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
):
    """
    Replace the metadata in a CIF file with metadata.
    """
    replace_mmcif_with_json(str(pdb), str(metadata), str(output))


if __name__ == "__main__":
    app()
