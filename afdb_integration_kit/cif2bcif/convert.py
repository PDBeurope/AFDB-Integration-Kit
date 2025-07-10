import concurrent.futures
import subprocess
from pathlib import Path


def process_file(input_path, output_path):
    command = [
        "cif2bcif",
        str(input_path),
        str(output_path),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, input_path, result.stdout, result.stderr


def run_cif2bcif(input_file: Path, output_file: Path):
    if not (
        str(output_file).endswith(".bcif") or str(output_file).endswith(".bcif.gz")
    ):
        print(
            "[WARNING] The output file extension '"
            f"{output_file.suffix}"
            "' is not '.bcif' or '.bcif.gz'."
        )
    print("Running command: cif2bcif")
    print(str(input_file), str(output_file))
    code, _, out, err = process_file(input_file, output_file)
    if code == 0:
        print("Command executed successfully:")
        print(out)
    else:
        print("Error executing command:")
        print(err)


def run_batch_cif2bcif(
    input_dir: Path, output_dir: Path, workers: int = 4, gzip: bool = False
):
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
