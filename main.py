from pathlib import Path

import typer

from afdb_integration_kit.cif2bcif.convert import run_batch_cif2bcif, run_cif2bcif

app = typer.Typer()


@app.command()
def test():
    """
    Runs a series of checks to verify the environment and toolchain.
    """
    import subprocess

    print("--- Verifying Versions ---")
    subprocess.run(["python", "--version"])
    subprocess.run(["node", "--version"])
    subprocess.run(["npm", "--version"])
    print("\n--- Testing molstar Preprocess Script ---")
    subprocess.run(["node", "molstar/lib/commonjs/servers/model/preprocess", "-h"])


@app.command()
def cif2bcif(
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
    run_cif2bcif(input_file, output_file)


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
    run_batch_cif2bcif(input_dir, output_dir, workers=workers, gzip=gzip)


if __name__ == "__main__":
    app()
