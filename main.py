import concurrent.futures
import subprocess
from pathlib import Path

import typer

app = typer.Typer()


@app.command()
def test():
    """
    Runs a series of checks to verify the environment and toolchain.
    """
    print("--- Verifying Versions ---")
    subprocess.run(["python", "--version"])
    subprocess.run(["node", "--version"])
    subprocess.run(["npm", "--version"])
    print("\n--- Testing molstar Preprocess Script ---")
    subprocess.run(["node", "molstar/lib/commonjs/servers/model/preprocess", "-h"])


def process_file(input_path, output_path):
    command = [
        "node",
        "molstar/lib/commonjs/cli/cif2bcif",
        str(input_path),
        str(output_path),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, input_path, result.stdout, result.stderr


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
    if not (
        str(output_file).endswith(".bcif") or str(output_file).endswith(".bcif.gz")
    ):
        print(
            "[WARNING] The output file extension '"
            f"{output_file.suffix}"
            "' is not '.bcif' or '.bcif.gz'."
        )
    cmd_str = f"node molstar/lib/commonjs/cli/cif2bcif {input_file} {output_file}"
    print(f"Running command: {cmd_str}")
    code, _, out, err = process_file(input_file, output_file)
    if code == 0:
        print("Command executed successfully:")
        print(out)
    else:
        print("Error executing command:")
        print(err)


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
    input_files = list(input_dir.glob("*.cif"))
    if output_dir.exists():
        print(f"Output directory {output_dir} already exists.")
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Created output directory {output_dir}.")

    def process(input_file):
        ext = ".bcif.gz" if gzip else ".bcif"
        output_file = output_dir / (input_file.stem + ext)
        if not (
            str(output_file).endswith(".bcif") or str(output_file).endswith(".bcif.gz")
        ):
            print(
                "[WARNING] The output file extension '"
                f"{output_file.suffix}"
                "' is not '.bcif' or '.bcif.gz'. File: "
            )
            print(f"{output_file}")
        code, _, out, err = process_file(input_file, output_file)
        if code == 0:
            return (input_file.name, True, out)
        else:
            return (input_file.name, False, err)

    print(f"Processing {len(input_files)} files with {workers} workers...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process, f): f for f in input_files}
        for future in concurrent.futures.as_completed(futures):
            fname, ok, msg = future.result()
            if ok:
                print(f"[OK] {fname}")
            else:
                short_msg = msg[:80]
                suffix = "..." if len(msg) > 80 else ""
                error_line = "[ERROR] " + fname + ": " + short_msg + suffix
                print(error_line)


if __name__ == "__main__":
    app()
