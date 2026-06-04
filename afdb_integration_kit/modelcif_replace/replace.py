import argparse
import logging
from pathlib import Path
from typing import Any, Dict

import gemmi
import orjson

# Configure logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("afdb_integration_kit")


def replace_mmcif_with_json(mmcif_file: str, json_file: str, output_file: str) -> None:
    """
    Replace categories in an mmCIF file with categories from a JSON file
    and save the result.

    Args:
        mmcif_file (str): Path to the input mmCIF file.
        json_file (str): Path to the input JSON file.
        output_file (str): Path to the output mmCIF file.

    Raises:
        ValueError: If the JSON does not contain the expected 'categories' key.
        FileNotFoundError: If input files are not found.
        Exception: For other gemmi or file I/O errors.
    """
    # Load the mmCIF file
    doc = gemmi.cif.read_file(mmcif_file)
    block = doc.sole_block()

    # Load the JSON file
    data: Dict[str, Any] = orjson.loads(Path(json_file).read_bytes())

    if "categories" not in data:
        raise ValueError("Invalid JSON structure: 'categories' key is missing.")

    # Replace categories from JSON into mmCIF
    for category, values in data["categories"].items():
        try:
            item = block.get_mmcif_category(category)
            item.update(values)
            block.set_mmcif_category(category, item)
        except Exception as e:
            logger.warning(f"Failed to replace category '{category}': {e}")

    # Write the updated mmCIF to the output file
    write_options = gemmi.cif.WriteOptions()
    write_options.align_loops = 50
    write_options.align_pairs = 50
    write_options.prefer_pairs = True

    doc.write_file(output_file, write_options)
    logger.info(f"Replaced mmCIF file saved to {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Replace mmCIF categories with JSON data.")
    parser.add_argument("input_mmcif", help="Path to the input mmCIF file")
    parser.add_argument("input_json", help="Path to the input JSON file")
    parser.add_argument("output_mmcif", help="Path to the output mmCIF file")

    args = parser.parse_args()

    try:
        replace_mmcif_with_json(args.input_mmcif, args.input_json, args.output_mmcif)
    except Exception as e:
        logger.error(f"Error during replacement: {e}")
        raise


if __name__ == "__main__":
    main()
