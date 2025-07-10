import concurrent.futures
import logging
import subprocess
from pathlib import Path

# Configure logger
logger = logging.getLogger("cif2bcif")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)s] %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


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
        logger.warning(
            "The output file extension '%s' is not '.bcif' or '.bcif.gz'.",
            output_file.suffix,
        )
    logger.info("Running command: cif2bcif")
    logger.info("%s %s", str(input_file), str(output_file))
    code, _, out, err = process_file(input_file, output_file)
    if code == 0:
        logger.info("Command executed successfully:\n%s", out)
    else:
        logger.error("Error executing command:\n%s", err)


def run_batch_cif2bcif(
    input_dir: Path, output_dir: Path, workers: int = 4, gzip: bool = False
):
    input_files = list(input_dir.glob("*.cif"))
    if output_dir.exists():
        logger.info("Output directory %s already exists.", output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Created output directory %s.", output_dir)

    def process(input_file):
        ext = ".bcif.gz" if gzip else ".bcif"
        output_file = output_dir / (input_file.stem + ext)
        if not (
            str(output_file).endswith(".bcif") or str(output_file).endswith(".bcif.gz")
        ):
            logger.warning(
                "The output file extension '%s' is not '.bcif' or '.bcif.gz'. File: %s",
                output_file.suffix,
                output_file,
            )
        code, _, out, err = process_file(input_file, output_file)
        if code == 0:
            return (input_file.name, True, out)
        else:
            return (input_file.name, False, err)

    logger.info("Processing %d files with %d workers...", len(input_files), workers)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(process, f): f for f in input_files}
        for future in concurrent.futures.as_completed(futures):
            fname, ok, msg = future.result()
            if ok:
                logger.info("[OK] %s", fname)
            else:
                short_msg = msg[:80]
                suffix = "..." if len(msg) > 80 else ""
                error_line = "[ERROR] " + fname + ": " + short_msg + suffix
                logger.error(error_line)
