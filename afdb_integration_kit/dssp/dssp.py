import subprocess
import logging
from pathlib import Path

# Configure logger
logger = logging.getLogger("dssp")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)s] %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


def run_command(input_path, output_path):
    command = [
        "mkdssp",
        str(input_path),
        str(output_path),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, input_path, result.stdout, result.stderr


def run_dssp(input_file: Path, output_file: Path):
    # output should be a .cif file
    if not output_file.suffix == ".cif":
        raise ValueError("Output file must have a .cif extension")

    logger.info("Running command: mkdssp")
    code, _, out, err = run_command(input_file, output_file)
    if code == 0:
        logger.info("Command executed successfully:\n%s", out)
    else:
        logger.error("Error executing command:\n%s", err)
