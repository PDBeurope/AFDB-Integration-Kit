#!/usr/bin/env python3
"""Generate a small ColabFold-like end-to-end AFDB example dataset."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from afdb_integration_kit.complex_metrics import (
    DEFAULT_COMPLEX_ENRICHMENT_METRICS,
    build_chain_enrichment,
    build_model_enrichment,
    parse_ipsae_csv,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FIXTURE_MANIFEST = REPO_ROOT / "tests/fixtures/colabfold_real_examples/manifest.json"
DEFAULT_FIXTURES_ROOT = REPO_ROOT / "tests/fixtures/colabfold_real_examples"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "examples/colabfold_monomer_e2e"
DEFAULT_TEMPLATE = REPO_ROOT / "uniprot/templates/colabfold_example_modelcif_metadata.json"
DEFAULT_EXAMPLE_IDS = [
    "AF-0000000300000001",
    "AF-0000000300000002",
    "AF-0000000300000003",
]


@dataclass(frozen=True)
class ExampleModel:
    category: str
    example_id: str
    chain_spans: list[dict[str, Any]]
    directory: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a small runnable ColabFold-like example under examples/ "
            "using curated local fixtures and the existing AFDB toolkit scripts."
        )
    )
    parser.add_argument(
        "--duckdb",
        required=True,
        type=Path,
        help="DuckDB database with the UniProt entry table.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for the generated example (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--fixture-manifest",
        type=Path,
        default=DEFAULT_FIXTURE_MANIFEST,
        help="Fixture manifest JSON describing the curated ColabFold examples.",
    )
    parser.add_argument(
        "--fixtures-root",
        type=Path,
        default=DEFAULT_FIXTURES_ROOT,
        help="Root directory containing the curated ColabFold fixture files.",
    )
    parser.add_argument(
        "--example-id",
        action="append",
        dest="example_ids",
        default=[],
        help=(
            "Specific fixture example_id to include. Repeat for multiple models. "
            f"Defaults to: {', '.join(DEFAULT_EXAMPLE_IDS)}."
        ),
    )
    parser.add_argument(
        "--python-bin",
        default=sys.executable,
        help="Python interpreter to use for subprocess commands (default: current interpreter).",
    )
    parser.add_argument(
        "--main-script",
        type=Path,
        default=REPO_ROOT / "main.py",
        help="Path to main.py used for run-modelcif-gen/run-dssp/run-modelpdb-gen/run-cif2bcif.",
    )
    parser.add_argument(
        "--modelcif-template",
        type=Path,
        default=DEFAULT_TEMPLATE,
        help="ModelCIF template JSON to use.",
    )
    parser.add_argument(
        "--dssp-algorithm",
        default="mkdssp",
        choices=["mkdssp", "psea", "pydssp", "tmalign"],
        help="DSSP algorithm to use for the example (default: mkdssp).",
    )
    parser.add_argument(
        "--bcif-backend",
        default="biotite",
        choices=["molstar", "biotite", "auto"],
        help="BCIF backend to use for the example (default: biotite).",
    )
    parser.add_argument(
        "--provider-id",
        default="EXAMPLE",
        help="Provider ID recorded in the generated example config.",
    )
    parser.add_argument(
        "--provider-name",
        default="AFDB Integration Kit Example",
        help="Provider name recorded in the generated example config.",
    )
    parser.add_argument(
        "--tool-used",
        default="ColabFold v1.6.0 / AlphaFold-Multimer",
        help="toolUsed value recorded in dataset_config.json.",
    )
    parser.add_argument(
        "--model-created-date",
        default="2026-05-22T00:00:00Z",
        help="modelCreatedDate value recorded in dataset_config.json.",
    )
    parser.add_argument(
        "--pae-cutoff",
        type=float,
        default=10.0,
        help="PAE cutoff passed to iPSAE for complex examples (default: 10.0).",
    )
    parser.add_argument(
        "--dist-cutoff",
        type=float,
        default=15.0,
        help="Distance cutoff passed to iPSAE for complex examples (default: 15.0).",
    )
    return parser.parse_args()


def load_fixture_models(
    manifest_path: Path,
    fixtures_root: Path,
    selected_ids: Iterable[str],
) -> list[ExampleModel]:
    payload = json.loads(manifest_path.read_text())
    selected = list(selected_ids) or list(DEFAULT_EXAMPLE_IDS)
    requested = set(selected)

    by_id: dict[str, ExampleModel] = {}
    for item in payload.get("examples", []):
        example_id = item["example_id"]
        if example_id not in requested:
            continue
        category = item["category"]
        directory = fixtures_root / f"{category}s" / example_id
        by_id[example_id] = ExampleModel(
            category=category,
            example_id=example_id,
            chain_spans=item["chain_spans"],
            directory=directory,
        )

    missing = [example_id for example_id in selected if example_id not in by_id]
    if missing:
        raise ValueError(f"Fixture example_id(s) not found in manifest: {', '.join(missing)}")

    return [by_id[example_id] for example_id in selected]


def chain_manifest_rows(models: Iterable[ExampleModel]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for model in models:
        accession_to_entity: dict[str, str] = {}
        next_entity_id = 1
        for chain in model.chain_spans:
            accession = chain["uniprot_ac"]
            entity_id = accession_to_entity.setdefault(accession, str(next_entity_id))
            if entity_id == str(next_entity_id):
                next_entity_id += 1
            residue_count = chain.get("residue_count")
            if residue_count is None:
                residue_count = chain["sequence_end"] - chain["sequence_start"] + 1
            rows.append(
                {
                    "model_entity_id": model.example_id,
                    "entity_id": entity_id,
                    "chain_id": chain["chain_id"],
                    "uniprot_ac": accession,
                    "sequence_start": 1,
                    "sequence_end": residue_count,
                }
            )
    return rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def merge_csv_files(inputs: list[Path], output: Path) -> None:
    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    for path in sorted(inputs):
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file {path} has no header row.")
            if fieldnames is None:
                fieldnames = list(reader.fieldnames)
            elif list(reader.fieldnames) != fieldnames:
                raise ValueError(
                    f"CSV header mismatch while merging {path}: "
                    f"{reader.fieldnames!r} != {fieldnames!r}"
                )
            rows.extend(reader)

    if fieldnames is None:
        raise ValueError("No CSV files were provided for merging.")

    write_csv(output, fieldnames, rows)


def json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def stage_inputs(models: list[ExampleModel], input_dir: Path) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    for model in models:
        pdb_src = model.directory / f"{model.example_id}-model_v1.pdb"
        scores_src = model.directory / f"{model.example_id}-scores_v1.json"
        if not scores_src.exists():
            scores_src = model.directory / f"{model.example_id}-meta_v1.json"
        if not pdb_src.exists() or not scores_src.exists():
            raise FileNotFoundError(f"Missing fixture files for {model.example_id} in {model.directory}")
        shutil.copy2(pdb_src, input_dir / pdb_src.name)
        shutil.copy2(scores_src, input_dir / f"{model.example_id}-meta_v1.json")


def provider_payload(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "providerId": args.provider_id,
        "providerName": args.provider_name,
        "providerUrl": "https://github.com/PDBeurope/AFDB-Integration-Kit",
        "license": "CC-BY-4.0",
        "copyrights": [
            "Example fixture data for local AFDB toolkit verification only.",
        ],
    }


def dataset_config_payload(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "providerId": args.provider_id,
        "toolUsed": args.tool_used,
        "latestVersion": 1,
        "allVersions": [1],
        "entityType": "protein",
        "modelCreatedDate": args.model_created_date,
        "uniqueIdTemplate": "{model_entity_id}",
        "versionTag": "v1",
    }


def has_complex_models(models: Iterable[ExampleModel]) -> bool:
    return any(len(model.chain_spans) > 1 for model in models)


def enrich_model_jsons(model_json_dir: Path, ipsae_csv: Path) -> int:
    ipsae_data = parse_ipsae_csv(ipsae_csv)
    enriched = 0

    for json_file in sorted(model_json_dir.glob("*.json")):
        enrichment = build_model_enrichment(
            ipsae_data.get(json_file.stem, {}),
            {},
            DEFAULT_COMPLEX_ENRICHMENT_METRICS,
        )
        if not enrichment:
            continue

        data = json.loads(json_file.read_text(encoding="utf-8"))
        data.update(enrichment)
        json_dump(json_file, data)
        enriched += 1

    return enriched


def enrich_chain_jsons(chain_json_dir: Path, ipsae_csv: Path) -> int:
    ipsae_data = parse_ipsae_csv(ipsae_csv)
    enriched = 0

    for json_file in sorted(chain_json_dir.glob("*.json")):
        ipsae_row = ipsae_data.get(json_file.stem, {})
        if not ipsae_row:
            continue

        records = json.loads(json_file.read_text(encoding="utf-8"))
        modified = False
        for record in records:
            chain_id = record.get("uniqueId", "").rsplit("_", 1)[-1]
            enrichment = build_chain_enrichment(
                ipsae_row,
                chain_id,
                DEFAULT_COMPLEX_ENRICHMENT_METRICS,
            )
            if enrichment:
                record.update(enrichment)
                modified = True

        if modified:
            json_dump(json_file, records)
            enriched += 1

    return enriched


def run_command(
    command: list[str],
    cwd: Path,
    commands_log: list[str],
) -> None:
    commands_log.append(shlex.join(command))
    subprocess.run(command, cwd=cwd, check=True)


def try_command(
    command: list[str],
    cwd: Path,
    commands_log: list[str],
) -> bool:
    commands_log.append(shlex.join(command))
    result = subprocess.run(command, cwd=cwd, check=False)
    return result.returncode == 0


def tool_status() -> dict[str, Any]:
    return {
        "mkdssp_on_path": shutil.which("mkdssp") is not None,
        "cif2bcif_on_path": shutil.which("cif2bcif") is not None,
        "pydssp_installed": importlib.util.find_spec("pydssp") is not None,
        "biotite_installed": importlib.util.find_spec("biotite") is not None,
    }


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    readme_backup: str | None = None
    readme_path = output_dir / "README.md"
    if readme_path.exists():
        readme_backup = readme_path.read_text()
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if readme_backup is not None:
        readme_path.write_text(readme_backup)

    models = load_fixture_models(
        args.fixture_manifest.resolve(),
        args.fixtures_root.resolve(),
        args.example_ids,
    )
    model_ids = [model.example_id for model in models]
    commands_log: list[str] = []
    bcif_results: dict[str, dict[str, Any]] = {}

    config_dir = output_dir / "config"
    input_dir = output_dir / "input"
    scores_dir = output_dir / "scores"
    chain_manifest_dir = output_dir / "chain_manifests"
    model_manifest_dir = output_dir / "model_manifests"
    merged_manifest_dir = output_dir / "merged_manifests"
    model_json_dir = output_dir / "model_jsons"
    chain_json_dir = output_dir / "chain_jsons"
    model_batch_dir = output_dir / "model_batches"
    chain_batch_dir = output_dir / "chain_batches"
    ipsae_dir = output_dir / "ipsae"
    modelcif_input_dir = output_dir / "modelcif_input"
    modelcif_dir = output_dir / "modelcif"
    dssp_dir = output_dir / "dssp"
    modelpdb_dir = output_dir / "modelpdb"
    bcif_dir = output_dir / "bcif"

    for directory in [
        config_dir,
        input_dir,
        scores_dir,
        chain_manifest_dir,
        model_manifest_dir,
        merged_manifest_dir,
        model_json_dir,
        chain_json_dir,
        model_batch_dir,
        chain_batch_dir,
        ipsae_dir,
        modelcif_input_dir,
        modelcif_dir,
        dssp_dir,
        modelpdb_dir,
        bcif_dir,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    stage_inputs(models, input_dir)

    ids_file = config_dir / "model_ids.txt"
    ids_file.write_text("\n".join(model_ids) + "\n")
    write_csv(
        config_dir / "chain_manifest.csv",
        ["model_entity_id", "entity_id", "chain_id", "uniprot_ac", "sequence_start", "sequence_end"],
        chain_manifest_rows(models),
    )
    json_dump(config_dir / "selected_examples.json", {
        "example_ids": model_ids,
        "categories": {model.example_id: model.category for model in models},
        "source_directories": {model.example_id: str(model.directory) for model in models},
    })
    json_dump(config_dir / "dataset_config.json", dataset_config_payload(args))
    json_dump(config_dir / "provider.json", provider_payload(args))

    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/batch_validate_assets.py"),
            "--ids",
            str(ids_file),
            "--input-dir",
            str(input_dir),
            "--output",
            str(output_dir / "validation_results.tsv"),
        ],
        REPO_ROOT,
        commands_log,
    )
    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/batch_convert_colabfold.py"),
            "--model-ids-file",
            str(ids_file),
            "--input-dir",
            str(input_dir),
            "--manifest",
            str(config_dir / "chain_manifest.csv"),
            "--duckdb",
            str(args.duckdb.resolve()),
            "--output-dir",
            str(scores_dir),
            "--chain-manifest-dir",
            str(chain_manifest_dir),
            "--model-manifest-dir",
            str(model_manifest_dir),
            "--workers",
            "1",
        ],
        REPO_ROOT,
        commands_log,
    )

    complex_metrics_enabled = has_complex_models(models)
    if complex_metrics_enabled:
        run_command(
            [
                args.python_bin,
                str(REPO_ROOT / "uniprot/scripts/batch_ipsae.py"),
                "--pae-dir",
                str(input_dir),
                "--pdb-dir",
                str(input_dir),
                "--output-dir",
                str(ipsae_dir),
                "--model-ids",
                str(ids_file),
                "--pae-cutoff",
                str(args.pae_cutoff),
                "--dist-cutoff",
                str(args.dist_cutoff),
                "--workers",
                "1",
            ],
            REPO_ROOT,
            commands_log,
        )
        shutil.rmtree(ipsae_dir / "input", ignore_errors=True)

    merge_csv_files(
        list(chain_manifest_dir.glob("*_afid_mapping.csv")),
        merged_manifest_dir / "uniprot_afid_mapping.csv",
    )
    merge_csv_files(
        list(model_manifest_dir.glob("*_model_metadata.csv")),
        merged_manifest_dir / "uniprot_model_metadata.csv",
    )

    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/batch_export_metadata.py"),
            "--model-ids",
            str(ids_file),
            "--db",
            str(args.duckdb.resolve()),
            "--config",
            str(config_dir / "dataset_config.json"),
            "--mapping",
            str(merged_manifest_dir / "uniprot_afid_mapping.csv"),
            "--model-manifest",
            str(merged_manifest_dir / "uniprot_model_metadata.csv"),
            "--output-dir",
            str(model_json_dir),
            "--export-type",
            "model",
            "--workers",
            "1",
        ],
        REPO_ROOT,
        commands_log,
    )
    if complex_metrics_enabled:
        enrich_model_jsons(model_json_dir, ipsae_dir / "ipsae_summary.csv")

    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/batch_export_metadata.py"),
            "--model-ids",
            str(ids_file),
            "--db",
            str(args.duckdb.resolve()),
            "--config",
            str(config_dir / "dataset_config.json"),
            "--mapping",
            str(merged_manifest_dir / "uniprot_afid_mapping.csv"),
            "--model-manifest",
            str(merged_manifest_dir / "uniprot_model_metadata.csv"),
            "--output-dir",
            str(chain_json_dir),
            "--export-type",
            "chain",
            "--workers",
            "1",
        ],
        REPO_ROOT,
        commands_log,
    )
    if complex_metrics_enabled:
        enrich_chain_jsons(chain_json_dir, ipsae_dir / "ipsae_summary.csv")

    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/combine_metadata.py"),
            "--input-dir",
            str(model_json_dir),
            "--output-dir",
            str(model_batch_dir),
            "--output-prefix",
            "AF-metadata",
            "--chunk-size",
            "1000",
        ],
        REPO_ROOT,
        commands_log,
    )
    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/combine_metadata.py"),
            "--input-dir",
            str(chain_json_dir),
            "--output-dir",
            str(chain_batch_dir),
            "--output-prefix",
            "AF-chain-metadata",
            "--chunk-size",
            "1000",
        ],
        REPO_ROOT,
        commands_log,
    )
    run_command(
        [
            args.python_bin,
            str(REPO_ROOT / "uniprot/scripts/batch_export_modelcif_input.py"),
            "--model-ids",
            str(ids_file),
            "--manifest",
            str(merged_manifest_dir / "uniprot_afid_mapping.csv"),
            "--db",
            str(args.duckdb.resolve()),
            "--template",
            str(args.modelcif_template.resolve()),
            "--output-dir",
            str(modelcif_input_dir),
            "--workers",
            "1",
        ],
        REPO_ROOT,
        commands_log,
    )

    batch_modelcif_cmd = [
        args.python_bin,
        str(args.main_script.resolve()),
        "run-batch-modelcif-gen",
        "--input-dir",
        str(input_dir),
        "--metadata-dir",
        str(modelcif_input_dir),
        "--output-dir",
        str(modelcif_dir),
        "--model-version",
        "v1",
        "--skip-validation",
        "--skip-alignment",
        "--workers",
        "1",
    ]
    if complex_metrics_enabled:
        batch_modelcif_cmd.extend(
            [
                "--model-json-dir",
                str(model_json_dir),
                "--cif-qa-metrics",
                "auto",
            ]
        )
    run_command(batch_modelcif_cmd, REPO_ROOT, commands_log)

    for model_id in model_ids:
        run_command(
            [
                args.python_bin,
                str(args.main_script.resolve()),
                "run-dssp",
                "-i",
                str(modelcif_dir / f"{model_id}-model_v1.cif"),
                "-o",
                str(dssp_dir / f"{model_id}-model_v1.cif"),
                "-a",
                args.dssp_algorithm,
            ],
            REPO_ROOT,
            commands_log,
        )
        run_command(
            [
                args.python_bin,
                str(args.main_script.resolve()),
                "run-modelpdb-gen",
                "-c",
                str(dssp_dir / f"{model_id}-model_v1.cif"),
                "-p",
                str(input_dir / f"{model_id}-model_v1.pdb"),
                "-r",
                str(config_dir / "provider.json"),
                "-o",
                str(modelpdb_dir / f"{model_id}-model_v1.pdb"),
            ],
            REPO_ROOT,
            commands_log,
        )

        bcif_from_dssp = [
            args.python_bin,
            str(args.main_script.resolve()),
            "run-cif2bcif",
            "-i",
            str(dssp_dir / f"{model_id}-model_v1.cif"),
            "-o",
            str(bcif_dir / f"{model_id}-model_v1.bcif"),
            "-b",
            args.bcif_backend,
        ]
        if try_command(bcif_from_dssp, REPO_ROOT, commands_log):
            bcif_results[model_id] = {
                "status": "generated",
                "source_cif": f"dssp/{model_id}-model_v1.cif",
                "output_bcif": f"bcif/{model_id}-model_v1.bcif",
            }
            continue

        bcif_from_modelcif = [
            args.python_bin,
            str(args.main_script.resolve()),
            "run-cif2bcif",
            "-i",
            str(modelcif_dir / f"{model_id}-model_v1.cif"),
            "-o",
            str(bcif_dir / f"{model_id}-model_v1.bcif"),
            "-b",
            args.bcif_backend,
        ]
        if try_command(bcif_from_modelcif, REPO_ROOT, commands_log):
            bcif_results[model_id] = {
                "status": "generated_with_fallback",
                "source_cif": f"modelcif/{model_id}-model_v1.cif",
                "output_bcif": f"bcif/{model_id}-model_v1.bcif",
                "reason": (
                    "BCIF generation from the DSSP-enriched CIF failed with the local backend; "
                    "the example BCIF was generated from the pre-DSSP ModelCIF instead."
                ),
            }
            continue

        bcif_results[model_id] = {
            "status": "skipped",
            "reason": (
                "BCIF generation failed for both the DSSP-enriched CIF and the pre-DSSP "
                "ModelCIF with the configured local backend."
            ),
        }

    (config_dir / "commands.txt").write_text("\n".join(commands_log) + "\n")
    json_dump(
        output_dir / "run_summary.json",
        {
            "example_ids": model_ids,
            "duckdb": str(args.duckdb.resolve()),
            "dssp_algorithm": args.dssp_algorithm,
            "bcif_backend": args.bcif_backend,
            "bcif_results": bcif_results,
            "tool_status": tool_status(),
            "notes": {
                "dssp": (
                    "This example uses a local fallback algorithm rather than mkdssp."
                    if args.dssp_algorithm != "mkdssp"
                    else "This example uses mkdssp."
                ),
                "bcif": (
                    "This example uses an explicit BCIF backend instead of relying on Mol* auto-detection."
                    if args.bcif_backend != "molstar"
                    else "This example uses the Mol* backend."
                ),
            },
        },
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
