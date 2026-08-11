#!/usr/bin/env python3
"""Generate and validate the fragment-metadata pre-structure E2E example.

This deliberately stops before coordinate-dependent stages.  The following
stage, once source assets are available, needs ``<AF-ID>-model_v1.pdb`` and
the raw ``<AF-ID>-meta_v1.json`` containing pLDDT/PAE/max-PAE information;
it does not take input mmCIF files.
"""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

import duckdb


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EXAMPLE_DIR = REPO_ROOT / "examples" / "fragment_metadata_e2e"
DEFAULT_TEMPLATE = (
    REPO_ROOT / "uniprot" / "templates"
    / "colabfold_example_modelcif_metadata.json"
)
REAL_MODEL_IDS = (
    "AF-0000000211971324",
    "AF-0000000212005744",
    "AF-0000000212013504",
    "AF-0000000212013519",
    "AF-0000000212039399",
    "AF-0000000212039401",
    "AF-0000000212039400",
)
FRAGMENT_NAMES = {
    ("P27409", "1", "46"): "Development fragment alpha",
    ("P27409", "961", "1071"): "Development fragment omega",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the disposable fragment metadata pre-structure E2E "
            "example."
        )
    )
    parser.add_argument(
        "--example-dir", type=Path, default=DEFAULT_EXAMPLE_DIR
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Generated directory to recreate; its final path component must "
            "be 'generated'."
        ),
    )
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    return parser.parse_args(argv)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def recreate_generated_directory(path: Path) -> Path:
    """Recreate only a deliberately named generated directory."""
    resolved = path.expanduser().resolve()
    if resolved.name != "generated":
        raise ValueError(
            "--output-dir must end in 'generated' to permit recreation."
        )
    if resolved == resolved.parent or resolved == REPO_ROOT:
        raise ValueError("Refusing to recreate an unsafe output directory.")
    if resolved.exists():
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    return resolved


def create_mock_database(seed_path: Path, db_path: Path) -> None:
    """Build the entry schema needed by the current metadata exporters."""
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    entries = seed.get("entries")
    if not isinstance(entries, list):
        raise ValueError("Mock UniProt seed must contain an 'entries' list.")
    con = duckdb.connect(str(db_path))
    try:
        con.execute(
            """
            CREATE TABLE entry (
                primary_ac VARCHAR NOT NULL PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                protein_full_names VARCHAR[],
                protein_short_names VARCHAR[],
                entry_name VARCHAR,
                gene_names VARCHAR[],
                gene_synonyms VARCHAR[],
                gene_ordered_locus_names VARCHAR[],
                gene_orf_names VARCHAR[],
                sequence_version_date VARCHAR,
                taxid INTEGER,
                organism VARCHAR,
                organism_common_names VARCHAR[],
                organism_synonyms VARCHAR[],
                is_uniprot_reference_proteome BOOLEAN,
                reviewed BOOLEAN
            )
            """
        )
        rows: list[tuple[Any, ...]] = []
        for item in entries:
            accession = str(item["primary_ac"])
            length = int(item["sequence_length"])
            protein_name = str(item["protein_name"])
            if length < 1:
                raise ValueError(
                    f"Mock sequence length for {accession} must be positive."
                )
            rows.append(
                (
                    accession,
                    "A" * length,
                    [protein_name],
                    [],
                    f"MOCK_{accession}",
                    [f"mock_{accession.lower()}"],
                    [], [], [], "2026-01-01", 999999,
                    "Mock virus", [], [], True, True,
                )
            )
        con.executemany(
            (
                "INSERT INTO entry VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, "
                "?, ?, ?, ?, ?)"
            ),
            rows,
        )
    finally:
        con.close()


def run_stage(command: list[str], log_path: Path) -> None:
    """Run one real toolkit command and append its complete output to a log."""
    result = subprocess.run(
        command, cwd=REPO_ROOT, text=True, capture_output=True, check=False
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("$ " + shlex.join(command) + "\n")
        handle.write(result.stdout)
        handle.write(result.stderr)
        handle.write("\n")
    if result.returncode:
        raise RuntimeError(
            f"Stage failed ({result.returncode}): {shlex.join(command)}; "
            f"see {log_path}"
        )


def model_ids(rows: Iterable[dict[str, str]]) -> list[str]:
    return list(dict.fromkeys(row["model_entity_id"] for row in rows))


def core_rows(rows: Iterable[dict[str, str]]) -> list[tuple[str, ...]]:
    fields = (
        "model_entity_id", "entity_id", "chain_id", "uniprot_ac",
        "is_fragment", "sequence_start", "sequence_end", "protein_name",
    )
    return sorted(
        tuple(row.get(field, "") for field in fields) for row in rows
    )


def expected_name(row: dict[str, str], names: dict[str, str]) -> str:
    key = (row["uniprot_ac"], row["sequence_start"], row["sequence_end"])
    return FRAGMENT_NAMES.get(key, names[row["uniprot_ac"]])


def validate_generated_outputs(
    output_dir: Path, canonical_source: Path, names: dict[str, str]
) -> None:
    """Assert topology, lookup, entity, and consumer outcomes for all cases."""
    final_rows = read_csv(output_dir / "manifests" / "canonical_enriched.csv")
    source_rows = read_csv(canonical_source)
    if len(final_rows) != 12:
        raise AssertionError(
            f"Expected 12 chain rows, found {len(final_rows)}."
        )
    expected_counts = {
        REAL_MODEL_IDS[0]: 1,
        REAL_MODEL_IDS[1]: 1,
        REAL_MODEL_IDS[2]: 2,
        REAL_MODEL_IDS[3]: 2,
        REAL_MODEL_IDS[4]: 2,
        REAL_MODEL_IDS[5]: 2,
        REAL_MODEL_IDS[6]: 2,
    }
    actual_counts: dict[str, int] = {}
    for row in final_rows:
        model_id = row["model_entity_id"]
        actual_counts[model_id] = actual_counts.get(model_id, 0) + 1
        if row["is_fragment"] == "true":
            expected = FRAGMENT_NAMES.get(
                (row["uniprot_ac"], row["sequence_start"], row["sequence_end"])
            )
            if row["protein_name"] != expected:
                raise AssertionError(f"Fragment name mismatch in {row!r}.")
        elif row["protein_name"]:
            raise AssertionError(
                f"Full-length row unexpectedly has a custom name: {row!r}"
            )
    if actual_counts != expected_counts:
        raise AssertionError(f"Unexpected topology counts: {actual_counts!r}")

    by_model: dict[str, list[dict[str, str]]] = {}
    for row in final_rows:
        by_model.setdefault(row["model_entity_id"], []).append(row)
    for model_id, rows in by_model.items():
        component_to_entity: dict[tuple[str, str, str], str] = {}
        for row in rows:
            component = (
                row["uniprot_ac"], row["sequence_start"], row["sequence_end"]
            ) if row["is_fragment"] == "true" else (
                row["uniprot_ac"], "", ""
            )
            old = component_to_entity.setdefault(component, row["entity_id"])
            if old != row["entity_id"]:
                raise AssertionError(
                    "Repeated component has different entities in "
                    f"{model_id}."
                )
        if len(component_to_entity) != len(set(component_to_entity.values())):
            raise AssertionError(
                f"Distinct components share an entity in {model_id}."
            )

    real_final = [
        row for row in final_rows if row["model_entity_id"] in REAL_MODEL_IDS
    ]
    wide_rows = read_csv(output_dir / "manifests" / "real_cases_from_wide.csv")
    if core_rows(real_final) != core_rows(wide_rows):
        raise AssertionError(
            "Wide-derived rows differ from canonical enrichment for real "
            "cases."
        )

    # Source rows must retain the intended topology and canonical enrichment
    # not silently alter which models/chains are being tested.
    source_pairs = sorted(
        (row["model_entity_id"], row["chain_id"]) for row in source_rows
    )
    final_pairs = sorted(
        (row["model_entity_id"], row["chain_id"]) for row in final_rows
    )
    if source_pairs != final_pairs:
        raise AssertionError(
            "Canonical enrichment changed model or chain coverage."
        )

    for model_id, rows in by_model.items():
        descriptions = {expected_name(row, names) for row in rows}
        model_payload = json.loads(
            (output_dir / "model_metadata" / f"{model_id}.json").read_text()
        )
        if not descriptions <= set(model_payload["uniprotDescription"]):
            raise AssertionError(
                f"Model metadata lacks expected names for {model_id}."
            )
        chain_payload = json.loads(
            (output_dir / "chain_metadata" / f"{model_id}.json").read_text()
        )
        if {item["uniprotDescription"] for item in chain_payload} != {
            expected_name(row, names) for row in rows
        }:
            raise AssertionError(
                f"Chain metadata names mismatch for {model_id}."
            )
        cif_payload = json.loads(
            (output_dir / "modelcif_input" / f"{model_id}.json").read_text()
        )
        entity_descriptions = set(
            cif_payload["categories"]["_entity"]["pdbx_description"]
        )
        if entity_descriptions != descriptions:
            raise AssertionError(
                f"ModelCIF input names mismatch for {model_id}."
            )

    model_batch = json.loads(
        (output_dir / "model_batches" / "models.json").read_text()
    )
    chain_batch = json.loads(
        (output_dir / "chain_batches" / "chains.json").read_text()
    )
    if len(model_batch) != 7 or len(chain_batch) != 12:
        raise AssertionError(
            "Combined metadata batches have unexpected record counts."
        )


def run_example(args: argparse.Namespace) -> Path:
    example_dir = args.example_dir.expanduser().resolve()
    config_dir = example_dir / "config"
    output_dir = recreate_generated_directory(
        args.output_dir or example_dir / "generated"
    )
    for required in (
        config_dir / "source_collaborator_wide.csv",
        config_dir / "canonical_input_manifest.csv",
        config_dir / "fragment_metadata.json",
        config_dir / "mock_uniprot_seed.json",
        args.template,
    ):
        if not required.is_file():
            raise FileNotFoundError(
                f"Required E2E fixture is missing: {required}"
            )

    logs = output_dir / "logs" / "stages.log"
    logs.parent.mkdir(parents=True)
    db_path = output_dir / "mock_uniprot.duckdb"
    create_mock_database(config_dir / "mock_uniprot_seed.json", db_path)
    source_rows = read_csv(config_dir / "canonical_input_manifest.csv")
    ids = model_ids(source_rows)
    ids_path = output_dir / "config" / "model_ids.txt"
    ids_path.parent.mkdir(parents=True)
    ids_path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    dataset_config = {
        "providerId": "FRAGMENT-E2E",
        "toolUsed": "Mock pre-structure E2E",
        "latestVersion": 1,
        "allVersions": [1],
        "modelCreatedDate": "2026-08-10",
        "entityType": "protein",
    }
    write_json(output_dir / "config" / "dataset_config.json", dataset_config)
    write_json(
        output_dir / "config" / "provider.json",
        {
            "providerId": "FRAGMENT-E2E",
            "providerName": "Fragment metadata E2E fixture",
            "providerUrl": (
                "https://github.com/PDBeurope/AFDB-Integration-Kit"
            ),
        },
    )
    python = str(args.python_bin)
    scripts = REPO_ROOT / "uniprot" / "scripts"
    reports = output_dir / "reports"
    manifests = output_dir / "manifests"
    reports.mkdir()
    manifests.mkdir()
    run_stage([
        python, str(scripts / "add_fragment_metadata.py"), "--db",
        str(db_path),
        "--fragments", str(config_dir / "fragment_metadata.json"), "--report",
        str(reports / "add_fragment_metadata.json"),
    ], logs)
    run_stage([
        python, str(scripts / "enrich_fragment_manifest.py"), "--input",
        str(config_dir / "source_collaborator_wide.csv"), "--output",
        str(manifests / "real_cases_from_wide.csv"), "--db", str(db_path),
        "--strict", "--report", str(reports / "enrich_wide.json"),
    ], logs)
    final_manifest = manifests / "canonical_enriched.csv"
    run_stage([
        python, str(scripts / "enrich_fragment_manifest.py"), "--input",
        str(config_dir / "canonical_input_manifest.csv"), "--output",
        str(final_manifest), "--db", str(db_path), "--strict", "--report",
        str(reports / "enrich_canonical.json"),
    ], logs)
    common_export = [
        "--model-ids", str(ids_path), "--db", str(db_path), "--config",
        str(output_dir / "config" / "dataset_config.json"), "--mapping",
        str(final_manifest),
        "--workers", "1",
    ]
    run_stage([
        python, str(scripts / "batch_export_metadata.py"), *common_export,
        "--output-dir", str(output_dir / "model_metadata"), "--export-type",
        "model",
    ], logs)
    run_stage([
        python, str(scripts / "batch_export_metadata.py"), *common_export,
        "--output-dir", str(output_dir / "chain_metadata"), "--export-type",
        "chain",
    ], logs)
    run_stage([
        python, str(scripts / "batch_export_modelcif_input.py"), "--model-ids",
        str(ids_path), "--manifest", str(final_manifest), "--db", str(db_path),
        "--template", str(args.template),
        "--output-dir", str(output_dir / "modelcif_input"), "--workers", "1",
    ], logs)
    for directory, batch_dir, filename in (
        (
            output_dir / "model_metadata",
            output_dir / "model_batches",
            "models.json",
        ),
        (
            output_dir / "chain_metadata",
            output_dir / "chain_batches",
            "chains.json",
        ),
    ):
        run_stage([
            python, str(scripts / "combine_metadata.py"), "--input-dir",
            str(directory), "--output-dir", str(batch_dir),
            "--output-filename",
            filename,
        ], logs)
    seed_entries = json.loads(
        (config_dir / "mock_uniprot_seed.json").read_text()
    )["entries"]
    names = {
        item["primary_ac"]: item["protein_name"]
        for item in seed_entries
    }
    validate_generated_outputs(
        output_dir, config_dir / "canonical_input_manifest.csv", names
    )
    write_json(
        output_dir / "run_summary.json",
        {
            "status": "success",
            "model_count": len(ids),
            "chain_row_count": len(source_rows),
            "stops_before": (
                "coordinate-dependent conversion and ModelCIF generation"
            ),
            "next_inputs": [
                "<AF-ID>-model_v1.pdb",
                "<AF-ID>-meta_v1.json (pLDDT/PAE/max_pAE)",
            ],
        },
    )
    return output_dir


def main(argv: list[str] | None = None) -> int:
    try:
        output_dir = run_example(parse_args(argv))
    except (
        OSError, ValueError, RuntimeError, AssertionError, duckdb.Error
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(
        "Generated and validated fragment metadata E2E example: "
        f"{output_dir}"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
