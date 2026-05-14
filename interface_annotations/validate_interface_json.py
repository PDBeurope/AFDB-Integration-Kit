#!/usr/bin/env python3
"""Validate AF interface annotations against the schema."""
import argparse
import json
from pathlib import Path

import jsonschema


def load_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def main():
    parser = argparse.ArgumentParser(description="Validate an AF interface annotation JSON file")
    parser.add_argument("schema", type=Path, help="Path to interface_schema.json")
    parser.add_argument("json_file", type=Path, help="Path to the JSON instance to validate")
    args = parser.parse_args()

    schema = load_json(args.schema)
    instance = load_json(args.json_file)

    jsonschema.validate(instance=instance, schema=schema)
    print("JSON validated - schema requirements satisfied")


if __name__ == "__main__":
    main()
