import logging
from pathlib import Path

import typer

from afdb_integration_kit.cif2bcif.convert import (
    run_batch_cif2bcif,
)
from afdb_integration_kit.cif2bcif.convert import run_cif2bcif as cif2bcif_helper
from afdb_integration_kit.dssp.dssp import run_dssp as dssp_helper

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


if __name__ == "__main__":
    app()
