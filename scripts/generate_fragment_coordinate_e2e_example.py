#!/usr/bin/env python3
"""Run the seven-model synthetic fragment-metadata E2E through ModelCIF."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import duckdb
import gemmi

import generate_fragment_metadata_e2e_example as prestructure
import synthesize_fragment_metadata_e2e_assets as synthesis


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EXAMPLE_DIR = REPO_ROOT / "examples/fragment_metadata_e2e"
DEFAULT_WORKSPACE = Path(
    "/mnt/disks/toolkit-data/viruses/fragment_metadata_synthetic_e2e"
)
DEFAULT_OUTPUT_DIR = DEFAULT_WORKSPACE / "generated"
DEFAULT_DONOR_DIR = Path("/mnt/disks/toolkit-data/viruses/sample_data")
DEFAULT_TEMPLATE = (
    REPO_ROOT / "uniprot/templates/colabfold_example_modelcif_metadata.json"
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--example-dir", type=Path, default=DEFAULT_EXAMPLE_DIR)
    parser.add_argument("--donor-dir", type=Path, default=DEFAULT_DONOR_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    return parser.parse_args(argv)


def run_command(command: list[str], log_path: Path) -> None:
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


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def merge_csv_files(inputs: Iterable[Path], output: Path) -> None:
    paths = sorted(inputs)
    if not paths:
        raise ValueError(f"No CSV files found for {output}.")
    fieldnames: list[str] | None = None
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"CSV has no header: {path}")
            if fieldnames is None:
                fieldnames = list(reader.fieldnames)
            elif fieldnames != list(reader.fieldnames):
                raise ValueError(f"CSV header mismatch while merging {path}.")
            rows.extend(reader)
    assert fieldnames is not None
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _expected_by_model(
    manifest_path: Path, seed_path: Path
) -> dict[str, list[synthesis.TargetChain]]:
    targets = synthesis.load_targets(manifest_path, seed_path)
    return {target.model_id: list(target.chains) for target in targets}


def _expected_name(
    chain: synthesis.TargetChain, seed_names: dict[str, str]
) -> str:
    if chain.is_fragment:
        key = (
            chain.uniprot_ac, str(chain.sequence_start),
            str(chain.sequence_end),
        )
        return prestructure.FRAGMENT_NAMES[key]
    return seed_names[chain.uniprot_ac]


def validate_mapping(
    mapping_path: Path,
    expected: dict[str, list[synthesis.TargetChain]],
) -> None:
    rows = prestructure.read_csv(mapping_path)
    if len(rows) != 12:
        raise AssertionError(f"Expected 12 merged chain rows, found {len(rows)}.")
    observed = {(row["model_entity_id"], row["chain_id"]): row for row in rows}
    for model_id, chains in expected.items():
        for chain in chains:
            row = observed.get((model_id, chain.chain_id))
            if row is None:
                raise AssertionError(f"Missing mapping row {model_id}:{chain.chain_id}.")
            comparisons = {
                "entity_id": chain.entity_id,
                "uniprot_ac": chain.uniprot_ac,
                "is_fragment": str(chain.is_fragment).lower(),
                "sequence_start": str(chain.sequence_start),
                "sequence_end": str(chain.sequence_end),
            }
            for key, value in comparisons.items():
                if row.get(key, "") != value:
                    raise AssertionError(
                        f"Mapping {model_id}:{chain.chain_id} {key} is "
                        f"{row.get(key)!r}, expected {value!r}."
                    )
            expected_fragment_name = ""
            if chain.is_fragment:
                expected_fragment_name = prestructure.FRAGMENT_NAMES[
                    (chain.uniprot_ac, str(chain.sequence_start),
                     str(chain.sequence_end))
                ]
            if row.get("protein_name", "") != expected_fragment_name:
                raise AssertionError(
                    f"Mapping name mismatch for {model_id}:{chain.chain_id}."
                )


def _load_json_records(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        return payload
    if isinstance(payload, dict):
        return [payload]
    raise AssertionError(f"Metadata JSON is not an object/list: {path}")


def validate_metadata(
    output_dir: Path,
    expected: dict[str, list[synthesis.TargetChain]],
    seed_names: dict[str, str],
) -> None:
    model_dir = output_dir / "metadata/model"
    chain_dir = output_dir / "metadata/chain"
    if len(list(model_dir.glob("*.json"))) != 7:
        raise AssertionError("Expected seven model metadata JSON files.")
    if len(list(chain_dir.glob("*.json"))) != 7:
        raise AssertionError("Expected seven chain metadata JSON files.")
    for model_id, chains in expected.items():
        model = _load_json_records(model_dir / f"{model_id}.json")[0]
        names = [_expected_name(chain, seed_names) for chain in chains]
        seen_entities: set[str] = set()
        component_names: list[str] = []
        for chain in chains:
            if chain.entity_id in seen_entities:
                continue
            seen_entities.add(chain.entity_id)
            component_names.append(_expected_name(chain, seed_names))
        if model["uniprotDescription"] != component_names:
            raise AssertionError(f"Model names mismatch for {model_id}.")
        expected_ranges = [
            (chain.sequence_start, chain.sequence_end) for chain in chains
        ]
        # Model-level complex metadata intentionally omits per-chain ranges;
        # those are authoritative in chain records.  Monomers retain a scalar
        # start/end and are checked here as well.
        if len(chains) == 1 and (
            model.get("sequenceStart"), model.get("sequenceEnd")
        ) != expected_ranges[0]:
            raise AssertionError(f"Model range mismatch for {model_id}.")
        chain_records = _load_json_records(chain_dir / f"{model_id}.json")
        if len(chain_records) != len(chains):
            raise AssertionError(f"Chain metadata count mismatch for {model_id}.")
        for record, chain, name in zip(chain_records, chains, names):
            if record["uniprotDescription"] != name:
                raise AssertionError(f"Chain name mismatch for {model_id}.")
            if (
                record["sequenceStart"], record["sequenceEnd"],
                record["isFragment"], len(record["sequence"]),
            ) != (
                chain.sequence_start, chain.sequence_end,
                chain.is_fragment, chain.length,
            ):
                raise AssertionError(f"Chain range mismatch for {model_id}.")
    model_batch_files = list((output_dir / "batches/model").glob("*.json"))
    chain_batch_files = list((output_dir / "batches/chain").glob("*.json"))
    if len(model_batch_files) != 1 or len(chain_batch_files) != 1:
        raise AssertionError("Expected one production-like batch per level.")
    model_batch = _load_json_records(model_batch_files[0])
    chain_batch = _load_json_records(chain_batch_files[0])
    if len(model_batch) != 7 or len(chain_batch) != 12:
        raise AssertionError("Metadata batch record counts are incorrect.")
    by_unique_id = {record["uniqueId"]: record for record in chain_batch}
    for model_id, chains in expected.items():
        for chain in chains:
            unique_id = f"{model_id}_v1_{chain.chain_id}"
            record = by_unique_id.get(unique_id)
            if record is None:
                raise AssertionError(f"Batched chain missing {unique_id}.")
            if record["uniprotDescription"] != _expected_name(
                chain, seed_names
            ):
                raise AssertionError(f"Batched chain name mismatch: {unique_id}.")


def _cif_rows(
    block: gemmi.cif.Block, category: str, columns: list[str]
) -> list[dict[str, str]]:
    table = block.find(f"{category}.", columns)
    if not table:
        raise AssertionError(f"Final CIF lacks category {category}.")
    return [
        {
            column: gemmi.cif.as_string(str(row[index]))
            for index, column in enumerate(columns)
        }
        for row in table
    ]


def _normalize_sequence(value: str) -> str:
    return "".join(value.replace(";", "").split())


def validate_final_cif(
    cif_path: Path,
    chains: list[synthesis.TargetChain],
    seed_names: dict[str, str],
) -> None:
    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    entity_rows = _cif_rows(block, "_entity", ["id", "pdbx_description"])
    descriptions = {row["id"]: row["pdbx_description"] for row in entity_rows}
    expected_entities: dict[str, synthesis.TargetChain] = {}
    for chain in chains:
        expected_entities.setdefault(chain.entity_id, chain)
    if set(descriptions) != set(expected_entities):
        raise AssertionError(f"Entity IDs mismatch in {cif_path.name}.")
    for entity_id, chain in expected_entities.items():
        if descriptions[entity_id] != _expected_name(chain, seed_names):
            raise AssertionError(f"Entity description mismatch in {cif_path.name}.")

    ref_rows = _cif_rows(
        block, "_struct_ref",
        ["id", "entity_id", "pdbx_db_accession", "db_name"],
    )
    ref_entity = {row["id"]: row["entity_id"] for row in ref_rows}
    for row in ref_rows:
        chain = expected_entities[row["entity_id"]]
        if row["pdbx_db_accession"] != chain.uniprot_ac:
            raise AssertionError(f"UniProt reference mismatch in {cif_path.name}.")
    seq_rows = _cif_rows(
        block, "_struct_ref_seq",
        [
            "ref_id", "pdbx_strand_id", "seq_align_beg",
            "seq_align_end", "db_align_beg", "db_align_end",
        ],
    )
    chain_by_id = {chain.chain_id: chain for chain in chains}
    for row in seq_rows:
        chain = chain_by_id[row["pdbx_strand_id"]]
        if ref_entity[row["ref_id"]] != chain.entity_id:
            raise AssertionError(f"Reference entity mismatch in {cif_path.name}.")
        actual = tuple(int(row[key]) for key in (
            "seq_align_beg", "seq_align_end", "db_align_beg", "db_align_end",
        ))
        expected = (1, chain.length, chain.sequence_start, chain.sequence_end)
        if actual != expected:
            raise AssertionError(
                f"Local/UniProt range mismatch in {cif_path.name}: "
                f"{actual} != {expected}."
            )

    target_rows = _cif_rows(
        block, "_ma_target_ref_db_details",
        [
            "target_entity_id", "db_accession", "seq_db_align_begin",
            "seq_db_align_end",
        ],
    )
    for row in target_rows:
        chain = expected_entities[row["target_entity_id"]]
        if (
            row["db_accession"] != chain.uniprot_ac
            or int(row["seq_db_align_begin"]) != chain.sequence_start
            or int(row["seq_db_align_end"]) != chain.sequence_end
        ):
            raise AssertionError(f"Target DB range mismatch in {cif_path.name}.")

    asym_rows = _cif_rows(block, "_struct_asym", ["id", "entity_id"])
    asym_entities = {row["id"]: row["entity_id"] for row in asym_rows}
    if asym_entities != {
        chain.chain_id: chain.entity_id for chain in chains
    }:
        raise AssertionError(f"struct_asym mapping mismatch in {cif_path.name}.")

    poly_rows = _cif_rows(
        block, "_entity_poly", ["entity_id", "pdbx_seq_one_letter_code"]
    )
    poly_sequences = {
        row["entity_id"]: _normalize_sequence(row["pdbx_seq_one_letter_code"])
        for row in poly_rows
    }
    for entity_id, chain in expected_entities.items():
        if poly_sequences[entity_id] != "A" * chain.length:
            raise AssertionError(f"entity_poly sequence mismatch in {cif_path.name}.")

    poly_seq_rows = _cif_rows(
        block, "_entity_poly_seq", ["entity_id", "num", "mon_id"]
    )
    counts: dict[str, set[int]] = defaultdict(set)
    for row in poly_seq_rows:
        if row["mon_id"] != "ALA":
            raise AssertionError(f"Non-ALA entity_poly_seq in {cif_path.name}.")
        counts[row["entity_id"]].add(int(row["num"]))
    for entity_id, chain in expected_entities.items():
        if counts[entity_id] != set(range(1, chain.length + 1)):
            raise AssertionError(f"entity_poly_seq mismatch in {cif_path.name}.")

    atom_rows = _cif_rows(
        block, "_atom_site",
        ["label_asym_id", "label_entity_id", "label_seq_id", "label_comp_id"],
    )
    observed_residues: dict[str, set[int]] = defaultdict(set)
    for row in atom_rows:
        chain = chain_by_id[row["label_asym_id"]]
        if row["label_entity_id"] != chain.entity_id:
            raise AssertionError(f"Atom entity mismatch in {cif_path.name}.")
        if row["label_comp_id"] != "ALA":
            raise AssertionError(f"Non-ALA atom_site in {cif_path.name}.")
        observed_residues[chain.chain_id].add(int(row["label_seq_id"]))
    for chain in chains:
        if observed_residues[chain.chain_id] != set(
            range(1, chain.length + 1)
        ):
            raise AssertionError(f"atom_site residue mismatch in {cif_path.name}.")


def validate_all_outputs(
    output_dir: Path,
    expected: dict[str, list[synthesis.TargetChain]],
    seed_names: dict[str, str],
) -> None:
    source_inputs = output_dir / "source_inputs"
    expected_source_inputs = {
        "canonical_input_manifest.csv",
        "source_collaborator_wide.csv",
        "fragment_metadata.json",
        "mock_uniprot_seed.json",
    }
    if {path.name for path in source_inputs.iterdir() if path.is_file()} != (
        expected_source_inputs
    ):
        raise AssertionError("External source-input snapshot is incomplete.")
    input_dir = output_dir / "input"
    if len(list(input_dir.glob("*-model_v1.pdb"))) != 7:
        raise AssertionError("Expected seven synthesized PDB files.")
    if len(list(input_dir.glob("*-meta_v1.json"))) != 7:
        raise AssertionError("Expected seven synthesized score JSON files.")
    validation_rows = [
        line.split("\t")
        for line in (output_dir / "validation_results.tsv")
        .read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(validation_rows) != 7 or any(row[1] != "ok" for row in validation_rows):
        raise AssertionError("Batch asset validation did not report seven OK rows.")
    scores_dir = output_dir / "converted/scores"
    if len(list(scores_dir.glob("*-confidence_v1.json"))) != 7:
        raise AssertionError("Expected seven converted confidence JSON files.")
    if len(list(scores_dir.glob("*-predicted_aligned_error_v1.json"))) != 7:
        raise AssertionError("Expected seven converted PAE JSON files.")
    mapping = output_dir / "manifests/merged/uniprot_afid_mapping.csv"
    validate_mapping(mapping, expected)
    validate_metadata(output_dir, expected, seed_names)
    cif_files = sorted((output_dir / "modelcif").glob("*-model_v1.cif"))
    if len(cif_files) != 7:
        raise AssertionError(f"Expected seven final ModelCIFs, found {len(cif_files)}.")
    for cif_path in cif_files:
        model_id = cif_path.name.removesuffix("-model_v1.cif")
        validate_final_cif(cif_path, expected[model_id], seed_names)

    # Explicit critical scenarios, kept separate for readable failure output.
    homodimer = expected["AF-0000000212039399"]
    if homodimer[0].entity_id != homodimer[1].entity_id:
        raise AssertionError("Identical fragment homodimer must share an entity.")
    heterodimer = expected["AF-0000000212039401"]
    if (
        heterodimer[0].uniprot_ac != heterodimer[1].uniprot_ac
        or heterodimer[0].entity_id == heterodimer[1].entity_id
        or (heterodimer[1].length, heterodimer[1].sequence_start,
            heterodimer[1].sequence_end) != (111, 961, 1071)
    ):
        raise AssertionError("Fragment heterodimer critical semantics failed.")
    mixed = expected["AF-0000000212039400"]
    if not mixed[0].is_fragment or mixed[1].is_fragment:
        raise AssertionError("Mixed fragment/full topology was not preserved.")


def run_example(args: argparse.Namespace) -> Path:
    example_dir = args.example_dir.expanduser().resolve()
    config_source = example_dir / "config"
    output_dir = prestructure.recreate_generated_directory(args.output_dir)
    donor_dir = args.donor_dir.expanduser().resolve()
    template = args.template.expanduser().resolve()
    required = (
        config_source / "canonical_input_manifest.csv",
        config_source / "source_collaborator_wide.csv",
        config_source / "fragment_metadata.json",
        config_source / "mock_uniprot_seed.json",
        donor_dir / f"{synthesis.MONOMER_DONOR_ID}-model_v1.pdb",
        donor_dir / f"{synthesis.MONOMER_DONOR_ID}-meta_v1.json",
        donor_dir / f"{synthesis.DIMER_DONOR_ID}-model_v1.pdb",
        donor_dir / f"{synthesis.DIMER_DONOR_ID}-meta_v1.json",
        template,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required input(s) missing: " + ", ".join(missing))

    directories = {
        "source_inputs": output_dir / "source_inputs",
        "config": output_dir / "config",
        "input": output_dir / "input",
        "reports": output_dir / "reports",
        "logs": output_dir / "logs",
        "manifests": output_dir / "manifests",
        "scores": output_dir / "converted/scores",
        "chain_manifests": output_dir / "manifests/per_model_chain",
        "model_manifests": output_dir / "manifests/per_model_model",
        "merged": output_dir / "manifests/merged",
        "model_metadata": output_dir / "metadata/model",
        "chain_metadata": output_dir / "metadata/chain",
        "model_batches": output_dir / "batches/model",
        "chain_batches": output_dir / "batches/chain",
        "modelcif_input": output_dir / "modelcif_input",
        "modelcif": output_dir / "modelcif",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    for source_name in (
        "canonical_input_manifest.csv",
        "source_collaborator_wide.csv",
        "fragment_metadata.json",
        "mock_uniprot_seed.json",
    ):
        shutil.copy2(
            config_source / source_name,
            directories["source_inputs"] / source_name,
        )
    log_path = directories["logs"] / "commands.log"
    manifest_source = config_source / "canonical_input_manifest.csv"
    seed_path = config_source / "mock_uniprot_seed.json"
    expected = _expected_by_model(manifest_source, seed_path)
    model_ids = list(expected)
    ids_path = directories["config"] / "model_ids.txt"
    ids_path.write_text("\n".join(model_ids) + "\n", encoding="utf-8")
    dataset_config = {
        "providerId": "FRAGMENT-SYNTHETIC-E2E",
        "toolUsed": "Synthetic software-validation fixture",
        "latestVersion": 1,
        "allVersions": [1],
        "modelCreatedDate": "2026-08-11",
        "entityType": "protein",
    }
    write_json(directories["config"] / "dataset_config.json", dataset_config)
    write_json(directories["config"] / "provider.json", {
        "providerId": "FRAGMENT-SYNTHETIC-E2E",
        "providerName": "Fragment metadata synthetic E2E",
        "providerUrl": "https://github.com/PDBeurope/AFDB-Integration-Kit",
        "warning": synthesis.WARNING,
    })
    db_path = output_dir / "mock_uniprot.duckdb"
    prestructure.create_mock_database(seed_path, db_path)
    python = str(args.python_bin)
    scripts = REPO_ROOT / "uniprot/scripts"
    run_command([
        python, str(scripts / "add_fragment_metadata.py"), "--db", str(db_path),
        "--fragments", str(config_source / "fragment_metadata.json"),
        "--report", str(directories["reports"] / "add_fragment_metadata.json"),
    ], log_path)
    wide_enriched = directories["manifests"] / "real_cases_from_wide.csv"
    run_command([
        python, str(scripts / "enrich_fragment_manifest.py"), "--input",
        str(config_source / "source_collaborator_wide.csv"), "--output",
        str(wide_enriched), "--db", str(db_path), "--strict", "--report",
        str(directories["reports"] / "enrich_wide.json"),
    ], log_path)
    enriched = directories["manifests"] / "canonical_enriched.csv"
    run_command([
        python, str(scripts / "enrich_fragment_manifest.py"), "--input",
        str(manifest_source), "--output", str(enriched), "--db", str(db_path),
        "--strict", "--report",
        str(directories["reports"] / "enrich_canonical.json"),
    ], log_path)
    if prestructure.core_rows(prestructure.read_csv(wide_enriched)) != (
        prestructure.core_rows(prestructure.read_csv(enriched))
    ):
        raise AssertionError("Wide and canonical enrichment core rows differ.")

    run_command([
        python,
        str(REPO_ROOT / "scripts/synthesize_fragment_metadata_e2e_assets.py"),
        "--donor-dir", str(donor_dir), "--manifest", str(manifest_source),
        "--seed", str(seed_path), "--output-dir", str(directories["input"]),
        "--provenance",
        str(directories["reports"] / "synthetic_provenance.json"),
    ], log_path)
    run_command([
        python, str(scripts / "batch_validate_assets.py"), "--ids", str(ids_path),
        "--input-dir", str(directories["input"]), "--output",
        str(output_dir / "validation_results.tsv"), "--workers", "1",
    ], log_path)
    run_command([
        python, str(scripts / "batch_convert_colabfold.py"),
        "--model-ids-file", str(ids_path), "--input-dir",
        str(directories["input"]), "--manifest", str(enriched), "--duckdb",
        str(db_path), "--output-dir", str(directories["scores"]),
        "--chain-manifest-dir", str(directories["chain_manifests"]),
        "--model-manifest-dir", str(directories["model_manifests"]),
        "--workers", "1",
    ], log_path)
    merged_chain = directories["merged"] / "uniprot_afid_mapping.csv"
    merged_model = directories["merged"] / "uniprot_model_metadata.csv"
    merge_csv_files(
        directories["chain_manifests"].glob("*_afid_mapping.csv"),
        merged_chain,
    )
    merge_csv_files(
        directories["model_manifests"].glob("*_model_metadata.csv"),
        merged_model,
    )
    common_export = [
        "--model-ids", str(ids_path), "--db", str(db_path), "--config",
        str(directories["config"] / "dataset_config.json"), "--mapping",
        str(merged_chain), "--model-manifest", str(merged_model),
        "--workers", "1",
    ]
    run_command([
        python, str(scripts / "batch_export_metadata.py"), *common_export,
        "--output-dir", str(directories["model_metadata"]),
        "--export-type", "model",
    ], log_path)
    run_command([
        python, str(scripts / "batch_export_metadata.py"), *common_export,
        "--output-dir", str(directories["chain_metadata"]),
        "--export-type", "chain",
    ], log_path)
    run_command([
        python, str(scripts / "combine_metadata.py"), "--input-dir",
        str(directories["model_metadata"]), "--output-dir",
        str(directories["model_batches"]), "--output-prefix", "AF-metadata",
        "--chunk-size", "1000",
    ], log_path)
    run_command([
        python, str(scripts / "combine_metadata.py"), "--input-dir",
        str(directories["chain_metadata"]), "--output-dir",
        str(directories["chain_batches"]), "--output-prefix",
        "AF-chain-metadata", "--chunk-size", "1000",
    ], log_path)
    run_command([
        python, str(scripts / "batch_export_modelcif_input.py"),
        "--model-ids", str(ids_path), "--manifest", str(merged_chain),
        "--db", str(db_path), "--template", str(template), "--output-dir",
        str(directories["modelcif_input"]), "--workers", "1",
    ], log_path)
    run_command([
        python, str(REPO_ROOT / "main.py"), "run-batch-modelcif-gen",
        "--input-dir", str(directories["input"]), "--metadata-dir",
        str(directories["modelcif_input"]), "--output-dir",
        str(directories["modelcif"]), "--model-version", "v1",
        "--dssp-algorithm", "mkdssp", "--skip-validation",
        "--skip-alignment", "--workers", "1",
    ], log_path)
    seed_entries = json.loads(seed_path.read_text(encoding="utf-8"))["entries"]
    seed_names = {
        str(item["primary_ac"]): str(item["protein_name"])
        for item in seed_entries
    }
    validate_all_outputs(output_dir, expected, seed_names)
    summary = {
        "status": "success",
        "warning": synthesis.WARNING,
        "model_count": 7,
        "chain_count": 12,
        "final_stage": "ModelCIF generation",
        "excluded_stages": ["iPSAE", "DSSP", "modelPDB", "BCIF"],
        "schema_readiness": {
            "metadata_batches": "pre_iPSAE_integration_fixtures",
            "final_production_schema_ready": False,
            "reason": (
                "Complex model and chain schemas require the iPSAE metric "
                "bundle, but iPSAE is intentionally excluded from this run."
            ),
            "required_production_step": (
                "Calculate and enrich iPSAE metrics before final schema "
                "validation and release; do not add placeholder metrics."
            ),
        },
        "commands_log": str(log_path),
        "synthetic_provenance": str(
            directories["reports"] / "synthetic_provenance.json"
        ),
        "validated_outputs": {
            "source_inputs": 4,
            "raw_pdb": 7,
            "raw_score_json": 7,
            "confidence_json": 7,
            "pae_json": 7,
            "model_metadata_json": 7,
            "chain_metadata_json": 7,
            "model_batch_records": 7,
            "chain_batch_records": 12,
            "modelcif": 7,
        },
    }
    write_json(output_dir / "run_summary.json", summary)
    return output_dir


def main(argv: list[str] | None = None) -> int:
    try:
        output = run_example(parse_args(argv))
    except (
        OSError, ValueError, RuntimeError, AssertionError, duckdb.Error,
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Generated and validated fragment coordinate E2E: {output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
