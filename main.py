import json
import logging
from collections import defaultdict
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

import typer
from rich.console import Console

console = Console()

from afdb_integration_kit.cif2bcif.convert import (
    run_batch_cif2bcif,
)
from afdb_integration_kit.cif2bcif.convert import run_cif2bcif as cif2bcif_helper
from afdb_integration_kit.dssp.dssp import run_dssp as dssp_helper
from afdb_integration_kit.dssp.dssp import run_batch_dssp
from afdb_integration_kit.metadata.validator import validate_against_schema
from afdb_integration_kit.modelcif.generate import generate
from afdb_integration_kit.modelcif_replace.replace import replace_mmcif_with_json as replace_mmcif_with_json
from afdb_integration_kit.modelpdb.generate import generate_pdb_headers
from afdb_integration_kit.quality_assessment.naming import (
    parse_ids_file
)
import afdb_integration_kit.validation.validators  # noqa: F401
from afdb_integration_kit.validation import (
    REGISTERED_CHECKS,
    Level,
    ValidationResult,
    run_validation_check,
    run_validations,
    summarise,
    write_results,
)
from afdb_integration_kit.utils.file_guard import require_non_empty_file
from concurrent.futures import ProcessPoolExecutor
from typing import Tuple

# --- Module-level worker functions for ProcessPoolExecutor (must be picklable) ---

def _process_modelcif_single(args: Tuple[str, str, str, str, bool, bool, str, str]) -> Tuple[str, bool, Optional[str]]:
    """
    Process a single model for ModelCIF generation.
    Module-level function required for ProcessPoolExecutor pickling.
    """
    json_file, input_dir, output_dir, model_version, skip_validation, skip_alignment, model_json_dir, cif_qa_metrics = args
    json_path = Path(json_file)
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    model_id = json_path.stem
    pdb_file = input_path / f"{model_id}-model_v1.pdb"

    if not pdb_file.exists():
        pdb_file = input_path / f"{model_id}-model_{model_version}.pdb"

    if not pdb_file.exists():
        return (model_id, False, "PDB file not found")

    output_file = output_path / f"{model_id}-model_{model_version}.cif"

    # Resolve model JSON path for QA metric enrichment
    model_json_path = None
    if model_json_dir:
        candidate = Path(model_json_dir) / f"{model_id}.json"
        if candidate.exists():
            model_json_path = str(candidate)

    try:
        # Import here to ensure module is available in subprocess
        from afdb_integration_kit.modelcif.generate import generate
        generate(
            str(pdb_file), str(json_file), str(output_file), None, False,
            skip_validation, skip_alignment,
            model_json_path, cif_qa_metrics or None,
        )
        return (model_id, True, None)
    except Exception as e:
        return (model_id, False, str(e))


def _process_modelpdb_single(args: Tuple[str, str, str, str, str]) -> Tuple[str, bool, Optional[str]]:
    """
    Process a single model for ModelPDB generation.
    Module-level function required for ProcessPoolExecutor pickling.
    """
    cif_file, pdb_dir, output_dir, provider_str, model_version = args
    cif_path = Path(cif_file)
    pdb_path = Path(pdb_dir)
    output_path = Path(output_dir)
    
    stem = cif_path.stem
    model_id = stem.replace(f"-model_{model_version}", "")
    
    pdb_file = pdb_path / f"{model_id}-model_v1.pdb"
    if not pdb_file.exists():
        pdb_file = pdb_path / f"{model_id}-model_{model_version}.pdb"
    
    if not pdb_file.exists():
        return (model_id, False, "PDB file not found")
    
    output_file = output_path / f"{model_id}-model_{model_version}.pdb"
    
    try:
        # Import here to ensure module is available in subprocess
        from afdb_integration_kit.modelpdb.generate import generate_pdb_headers
        generate_pdb_headers(str(cif_path), str(pdb_file), str(output_file), provider_str)
        return (model_id, True, None)
    except Exception as e:
        return (model_id, False, str(e))


# Set up logger
logger = logging.getLogger("afdb_integration_kit")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)

app = typer.Typer()
def _relative_path(path: Path, root: Optional[Path]) -> str:
    if root is None:
        return str(path)
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _format_validation_results(
    results: Sequence[ValidationResult],
    root: Path,
    *,
    include_info: bool = False,
) -> str:
    grouped: Dict[str, List[ValidationResult]] = defaultdict(list)
    for res in results:
        grouped[res.check].append(res)

    lines: List[str] = []
    for check in sorted(grouped):
        lines.append(f"[{check}]")
        check_results = grouped[check]
        for res in check_results:
            if res.level is Level.INFO and not include_info:
                continue
            location_bits = []
            location_bits.append(_relative_path(res.file, root))
            if res.location:
                location_bits.append(res.location)
            location_text = " | ".join(filter(None, location_bits))
            lines.append(f"- {res.level.value}: {location_text} → {res.message}")
            if res.suggested_fix and res.level is not Level.INFO:
                lines.append(f"  fix: {res.suggested_fix}")
            if res.metrics:
                metrics_text = ", ".join(f"{k}={v:.2f}" for k, v in res.metrics.items())
                lines.append(f"  metrics: {metrics_text}")
        lines.append("")

    return "\n".join(line for line in lines if line).strip()


def _emit_single_validation_results(
    results: Sequence[ValidationResult],
    root: Path,
) -> None:
    text = _format_validation_results(results, root, include_info=True)
    if text:
        typer.echo(text)
    else:
        typer.echo("No validation findings.")
    has_error = any(res.level is Level.ERROR for res in results)
    raise typer.Exit(code=1 if has_error else 0)


def _format_validation_summary(
    results: Sequence[ValidationResult],
    root: Path,
    *,
    include_info: bool = False,
    top_codes: int = 3,
) -> str:
    if not results:
        return ""

    filtered: List[ValidationResult] = [
        res for res in results if include_info or res.level is not Level.INFO
    ]
    if not filtered:
        filtered = list(results)

    by_file: Dict[str, List[ValidationResult]] = defaultdict(list)
    for res in filtered:
        key = _relative_path(res.file, root)
        by_file[key].append(res)

    if not by_file:
        return ""

    level_order = {Level.ERROR.value: 0, Level.WARN.value: 1, Level.INFO.value: 2}
    lines: List[str] = []
    lines.append(f"Files with findings: {len(by_file)}")

    for file_name in sorted(by_file):
        file_results = by_file[file_name]
        counts: Dict[str, int] = defaultdict(int)
        code_counts: Dict[str, int] = defaultdict(int)

        for res in file_results:
            counts[res.level.value] += 1
            if res.level is not Level.INFO:
                code_counts[res.code] += 1

        counts_text = ", ".join(
            f"{level}={counts[level]}"
            for level in sorted(counts, key=lambda lvl: level_order.get(lvl, 99))
        )
        lines.append(f"{file_name} → {counts_text}")

        if code_counts:
            sorted_codes = sorted(code_counts.items(), key=lambda item: (-item[1], item[0]))
            if top_codes > 0:
                sorted_codes = sorted_codes[:top_codes]
            codes_text = ", ".join(f"{code}x{count}" for code, count in sorted_codes)
            lines.append(f"  top codes: {codes_text}")

    return "\n".join(lines)


def _parse_scalar_token(value: str) -> Any:
    lowered = value.lower()
    if lowered in {"", "null", "none", "~"}:
        return None
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        if (value.startswith("\"") and value.endswith("\"")) or (
            value.startswith("'") and value.endswith("'")
        ):
            return value[1:-1]
        return value


def _parse_simple_yaml(text: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    current_key: Optional[str] = None
    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue
        if line[0] not in " \t":
            key, _, value = line.partition(":")
            key = key.strip()
            value = value.strip()
            if value:
                result[key] = _parse_scalar_token(value)
                current_key = None
            else:
                current_key = key
                result[key] = []
            continue

        if current_key is None:
            continue

        stripped = line.lstrip()
        if stripped.startswith("- "):
            container = result.get(current_key)
            if not isinstance(container, list):
                container = []
                result[current_key] = container
            container.append(_parse_scalar_token(stripped[2:].strip()))
        else:
            container = result.get(current_key)
            if not isinstance(container, dict):
                container = {}
                result[current_key] = container
            sub_key, _, sub_value = stripped.partition(":")
            container[sub_key.strip()] = _parse_scalar_token(sub_value.strip())

    return result


def _parse_yaml_text(text: str) -> Dict[str, Any]:
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(text) or {}
        if not isinstance(data, dict):
            raise ValueError("Configuration must be a mapping")
        return data
    except ModuleNotFoundError:
        return _parse_simple_yaml(text)


def _load_default_config() -> Dict[str, Any]:
    defaults_path = resources.files("afdb_integration_kit.validation").joinpath("defaults.yaml")
    text = defaults_path.read_text(encoding="utf-8")
    return _parse_yaml_text(text)


def _load_config_file(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
        if not isinstance(data, dict):
            raise ValueError("JSON config must be an object")
        return data
    return _parse_yaml_text(text)


def _merge_dicts(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged: Dict[str, Any] = dict(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, dict)
        ):
            merged[key] = _merge_dicts(merged[key], value)
        else:
            merged[key] = value
    return merged


@app.command()
def test():
    """
    Runs a series of checks to verify the environment and toolchain.
    """
    import subprocess

    logger.info("--- Verifying Versions ---")
    subprocess.run(["python", "--version"])
    subprocess.run(["node", "--version"])
    subprocess.run(["npm", "--version"])
    logger.info("--- Testing molstar cif2bcif script ---")
    subprocess.run(["cif2bcif", "-h"])
    logger.info("--- Testing DSSP script ---")
    subprocess.run(["mkdssp", "-h"])
    logger.info("--- Testing Gemmi script ---")
    subprocess.run(["gemmi", "--version"])


@app.command(name="list-validations")
def list_validations_cmd() -> None:
    """List the validation checks that are available to run."""
    if not REGISTERED_CHECKS:
        typer.echo("No validations are registered.")
        return

    for name in sorted(REGISTERED_CHECKS):
        typer.echo(name)


@app.command()
def run_schema_validation(
    input_file: Path = typer.Option(
        ...,
        "-i",
        "--input",
        help="Input JSON file path to validate against schema.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    type: str = typer.Option(
        ...,
        "-t",
        "--type",
        help="Type of schema to validate against ('model' or 'provider').",
    ),
):
    """
    Validate a JSON file against the specified schema.
    """

    validate_against_schema(input_file, type)


@app.command("validate-metadata-file")
def validate_metadata_file(
    metadata_file: Path = typer.Option(
        ...,
        "--file",
        "-f",
        help="Metadata JSON file to validate.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
):
    """
    Validate a single metadata JSON file (batch or per-accession).
    """
    require_non_empty_file(metadata_file, description="Metadata JSON file")
    results = run_validation_check(
        "metadata",
        [metadata_file],
        config={"metadata": {"allow_single_file": True}},
    )
    _emit_single_validation_results(results, metadata_file.parent)


@app.command("validate-plddt-file")
def validate_plddt_file(
    plddt_file: Path = typer.Option(
        ...,
        "--file",
        "-f",
        help="pLDDT confidence JSON file to validate.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
):
    """
    Validate a single pLDDT confidence JSON file.
    """
    require_non_empty_file(plddt_file, description="pLDDT confidence JSON file")
    results = run_validation_check(
        "plddt",
        [plddt_file],
        config={"plddt": {"allow_any_name": True}},
    )
    _emit_single_validation_results(results, plddt_file.parent)


@app.command("validate-pae-file")
def validate_pae_file(
    pae_file: Path = typer.Option(
        ...,
        "--file",
        "-f",
        help="PAE JSON file to validate.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
):
    """
    Validate a single PAE JSON file.
    """
    require_non_empty_file(pae_file, description="PAE JSON file")
    results = run_validation_check(
        "pae",
        [pae_file],
        config={"pae": {"allow_any_name": True}},
    )
    _emit_single_validation_results(results, pae_file.parent)


@app.command("validate-relationships-pair")
def validate_relationships_pair(
    plddt_file: Path = typer.Option(
        ...,
        "--plddt-file",
        help="pLDDT confidence JSON file.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    pae_file: Path = typer.Option(
        ...,
        "--pae-file",
        help="PAE JSON file.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
):
    """
    Validate a matching pLDDT/PAE pair for consistency.
    """
    require_non_empty_file(plddt_file, description="pLDDT confidence JSON file")
    require_non_empty_file(pae_file, description="PAE JSON file")
    results = run_validation_check(
        "relationships",
        [plddt_file, pae_file],
    )
    _emit_single_validation_results(results, plddt_file.parent)


@app.command("validate-sequences-file")
def validate_sequences_file(
    fasta_file: Path = typer.Option(
        ...,
        "--file",
        "-f",
        help="Sequences FASTA file to validate.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
):
    """
    Validate a single sequences.fasta/.fa file.
    """
    require_non_empty_file(fasta_file, description="Sequences FASTA file")
    results = run_validation_check(
        "sequences",
        [fasta_file],
        config={"sequences": {"allow_any_name": True}},
    )
    _emit_single_validation_results(results, fasta_file.parent)


@app.command()
def run_modelcif_gen(
    pdb: Path = typer.Option(
        ...,
        "-p",
        "--pdb",
        help="Input PDB file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    metadata: Path = typer.Option(
        ...,
        "-m",
        "--metadata",
        help="Input metadata JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output mmCIF file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    validate: str = typer.Option(
        None,
        "--validate",
        help="Optionally validate the output CIF file against a ModelCIF dictionary."
        "If used as a flag without a path, it defaults to 'mmcif_ma.dic'.",
    ),
    fetch_uniprot: bool = typer.Option(
        False,
        "--fetch-uniprot",
        help="Optionally fetch UniProt data for the chains in the PDB file.",
    ),
):
    """
    Enrich a PDB file with metadata to produce a feature-rich mmCIF file.
    """
    # Pre-flight checks
    require_non_empty_file(pdb, description="Input PDB file")
    require_non_empty_file(metadata, description="Input metadata JSON file")

    # If validate is passed as a flag (True but no value), default to 'mmcif_ma.dic'
    validate_path = validate if validate else None
    if validate is not None and validate == "":
        validate_path = "mmcif_ma.dic"
    if validate_path:
        require_non_empty_file(
            Path(validate_path),
            description="ModelCIF dictionary (--validate)",
        )

    # Call main logic (assuming main is imported or defined elsewhere)
    generate(str(pdb), str(metadata), str(output), validate_path, fetch_uniprot)


@app.command()
def run_batch_modelcif_gen(
    input_dir: Path = typer.Option(
        ...,
        "--input-dir",
        help="Directory containing PDB files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    metadata_dir: Path = typer.Option(
        ...,
        "--metadata-dir",
        help="Directory containing metadata JSON files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        help="Output directory for mmCIF files.",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    model_version: str = typer.Option(
        "v1",
        "--model-version",
        help="Model version string for output filenames.",
    ),
    workers: int = typer.Option(
        8,
        "--workers",
        help="Number of parallel workers.",
    ),
    skip_validation: bool = typer.Option(
        False,
        "--skip-validation",
        help="Skip JSON schema validation (use when input is pre-validated).",
    ),
    skip_alignment: bool = typer.Option(
        False,
        "--skip-alignment",
        help="Skip CIF column alignment for faster writes (intermediate files).",
    ),
    model_json_dir: Optional[Path] = typer.Option(
        None,
        "--model-json-dir",
        help="Directory containing model metadata JSONs for CIF QA metric enrichment.",
    ),
    cif_qa_metrics: Optional[str] = typer.Option(
        None,
        "--cif-qa-metrics",
        help="Comma-separated metric short names to add to CIF (e.g. 'ipsae_AB,iptm_af'). Use 'auto' for all.",
    ),
):
    """
    Batch process multiple PDB files to produce mmCIF files.
    Expects PDB files named {model_id}-model_v1.pdb and JSON files named {model_id}.json.
    Uses ProcessPoolExecutor for CPU-bound parallel execution (bypasses GIL).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all JSON metadata files
    json_files = list(metadata_dir.glob("*.json"))
    if not json_files:
        console.print("[red]No JSON files found in metadata directory[/red]")
        raise typer.Exit(1)

    console.print(f"[blue]Processing {len(json_files)} models with {workers} workers (ProcessPool)...[/blue]")

    # Prepare work items as tuples for the module-level worker function
    work_items = [
        (str(jf), str(input_dir), str(output_dir), model_version, skip_validation, skip_alignment,
         str(model_json_dir) if model_json_dir else "", cif_qa_metrics or "")
        for jf in json_files
    ]
    
    success_count = 0
    error_count = 0
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = executor.map(_process_modelcif_single, work_items)
        for model_id, success, error in results:
            if success:
                success_count += 1
            else:
                error_count += 1
                if error:
                    console.print(f"[yellow]Error {model_id}: {error}[/yellow]")
    
    console.print(f"[green]Batch complete: {success_count} success, {error_count} errors[/green]")


@app.command()
def run_modelpdb_gen(
    cif: Path = typer.Option(
        ...,
        "-c",
        "--cif",
        help="Input mmCIF file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    pdb: Path = typer.Option(
        ...,
        "-p",
        "--pdb",
        help="Input PDB file path (with ATOM coordinates).",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    provider: Path = typer.Option(
        ...,
        "-r",
        "--provider",
        help="Input provider JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output PDB file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
):
    """
    Enrich a PDB file with header information from a mmCIF file.
    """
    # Pre-flight checks
    require_non_empty_file(cif, description="Input mmCIF file")
    require_non_empty_file(pdb, description="Input PDB file")
    require_non_empty_file(provider, description="Input provider JSON file")

    # Call main logic
    generate_pdb_headers(str(cif), str(pdb), str(output), str(provider))


@app.command()
def run_batch_modelpdb_gen(
    cif_dir: Path = typer.Option(
        ...,
        "--cif-dir",
        help="Directory containing mmCIF files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    pdb_dir: Path = typer.Option(
        ...,
        "--pdb-dir",
        help="Directory containing original PDB files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    provider: Path = typer.Option(
        ...,
        "-r",
        "--provider",
        help="Input provider JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        help="Output directory for enriched PDB files.",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    model_version: str = typer.Option(
        "v1",
        "--model-version",
        help="Model version string for output filenames.",
    ),
    workers: int = typer.Option(
        8,
        "--workers",
        help="Number of parallel workers.",
    ),
):
    """
    Batch enrich PDB files with header information from mmCIF files.
    Expects CIF files named {model_id}-model_{version}.cif and PDB files named {model_id}-model_v1.pdb.
    Uses ProcessPoolExecutor for CPU-bound parallel execution (bypasses GIL).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all CIF files
    cif_files = list(cif_dir.glob("*-model_*.cif"))
    if not cif_files:
        console.print("[red]No CIF files found in CIF directory[/red]")
        raise typer.Exit(1)
    
    console.print(f"[blue]Processing {len(cif_files)} models with {workers} workers (ProcessPool)...[/blue]")
    
    provider_str = str(provider)
    
    # Prepare work items as tuples for the module-level worker function
    work_items = [
        (str(cf), str(pdb_dir), str(output_dir), provider_str, model_version)
        for cf in cif_files
    ]
    
    success_count = 0
    error_count = 0
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = executor.map(_process_modelpdb_single, work_items)
        for model_id, success, error in results:
            if success:
                success_count += 1
            else:
                error_count += 1
                if error:
                    console.print(f"[yellow]Error {model_id}: {error}[/yellow]")
    
    console.print(f"[green]Batch complete: {success_count} success, {error_count} errors[/green]")


@app.command()
def run_dssp(
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
        help="Output file in CIF format.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        readable=False,
        resolve_path=True,
    ),
):
    """
    Run DSSP on a CIF file to generate secondary structure information.
    """
    require_non_empty_file(input_file, description="DSSP input CIF file")
    logger.info(f"Converting {input_file} to {output_file}")
    dssp_helper(input_file, output_file)
    logger.info("Conversion complete.")


@app.command(name="batch-dssp")
def batch_dssp(
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
        help="Output directory for DSSP-annotated CIF files.",
        file_okay=False,
        dir_okay=True,
        writable=True,
        resolve_path=True,
    ),
    workers: int = typer.Option(
        8, "--workers", "-w", help="Number of parallel workers (default: 8)"
    ),
    algorithm: str = typer.Option(
        "psea",
        "--algorithm",
        "-a",
        help="Algorithm for secondary structure: 'psea' (geometry), 'pydssp' (H-bond), or 'tmalign' (CA-CA distance)"
    ),
    device: str = typer.Option(
        "cpu",
        "--device",
        "-d",
        help="Device for PyDSSP: 'cpu' or 'cuda'. GPU accelerates H-bond matrix computation."
    ),
):
    """
    Batch process all CIF files in a directory to add secondary structure.
    
    Supports three algorithms:
    - psea: Biotite's P-SEA algorithm (geometry-based, ~95% agreement with DSSP)
    - pydssp: PyDSSP (simplified H-bond DSSP, ~97% agreement with DSSP)
    - tmalign: TM-align make_sec algorithm (CA-CA distance patterns, very fast)
    
    When --device cuda is used with pydssp, the H-bond and SSE computation
    runs on GPU via PyTorch for faster processing.
    """
    if algorithm not in ("psea", "pydssp", "tmalign"):
        console.print(f"[red]Invalid algorithm '{algorithm}'. Use 'psea', 'pydssp', or 'tmalign'.[/red]")
        raise typer.Exit(1)
    
    algo_names = {"psea": "P-SEA (Biotite)", "pydssp": "PyDSSP", "tmalign": "TM-align"}
    algo_name = algo_names[algorithm]
    device_label = f"GPU ({device})" if device != "cpu" else "CPU"
    logger.info(
        f"Batch secondary structure processing from {input_dir} to {output_dir} "
        f"using {workers} workers with {algo_name} on {device_label}"
    )
    success, errors = run_batch_dssp(input_dir, output_dir, workers=workers, algorithm=algorithm, device=device)
    console.print(f"[green]Batch complete: {success} success, {errors} errors[/green]")


@app.command()
def run_cif2bcif(
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
    Convert CIF to BinaryCIF or BinaryCIF.GZ using gemmi.
    """
    require_non_empty_file(input_file, description="cif2bcif input CIF file")
    logger.info(f"Converting {input_file} to {output_file}")
    cif2bcif_helper(input_file, output_file)
    logger.info("Conversion complete.")


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
    logger.info(
        f"Batch converting CIF files from {input_dir} to {output_dir} "
        f"using {workers} workers. Gzip: {gzip}"
    )
    run_batch_cif2bcif(input_dir, output_dir, workers=workers, gzip=gzip)
    logger.info("Batch conversion complete.")


@app.command(name="run-validations")
def run_validation_suite(
    root: Path = typer.Option(
        ...,
        "--root",
        "-r",
        help="Dataset directory to validate.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    checks: List[str] = typer.Option(
        None,
        "--checks",
        "--check",
        "-c",
        help="Validation names to run. Defaults to config enabled list or all registered checks.",
    ),
    config_path: Optional[Path] = typer.Option(
        None,
        "--config",
        help="Optional YAML or JSON configuration file overriding defaults.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    out: Optional[Path] = typer.Option(
        None,
        "--out",
        "--output",
        "-o",
        help="Write machine-readable validation results to this path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    output_format: Optional[str] = typer.Option(
        None,
        "--format",
        help="Output format for --out (json, jsonl, or txt). Defaults to json.",
        case_sensitive=False,
    ),
    fail_on: str = typer.Option(
        "error",
        "--fail-on",
        help="Exit with code 1 on 'error' (default), 'warn', or 'none'.",
        case_sensitive=False,
    ),
    json_output: Optional[Path] = typer.Option(
        None,
        "--json-output",
        help="Compatibility alias for --out when emitting JSON.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    errors_only: bool = typer.Option(
        False,
        "--errors-only",
        help="Show only ERROR findings in the human-readable summary.",
    ),
    details: bool = typer.Option(
        False,
        "--details",
        help="Force detailed findings in the human-readable output.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Include INFO-level findings in the human-readable summary.",
    ),
    summary: bool = typer.Option(
        False,
        "--summary",
        help="Show a condensed per-file summary instead of detailed findings.",
    ),
    summary_limit: int = typer.Option(
        3,
        "--summary-limit",
        min=0,
        help="Maximum number of top error/warn codes to display per file when using --summary (0 for all).",
    ),
):
    """Run one or more validation checks via the shared validation runner."""

    defaults = _load_default_config()
    user_config: Dict[str, Any] = {}
    if config_path is not None:
        try:
            user_config = _load_config_file(config_path)
        except Exception as exc:  # pragma: no cover - defensive
            logger.error("Failed to load config %s: %s", config_path, exc)
            raise typer.Exit(code=1)

    config = _merge_dicts(defaults, user_config)

    requested = [c for c in checks if c] if checks else []
    if not requested:
        default_checks = config.get("enabled_checks")
        if isinstance(default_checks, Sequence):
            requested = [str(c) for c in default_checks]

    checks_arg: Optional[Sequence[str]]
    if not requested or requested == ["all"]:
        checks_arg = None
    else:
        unknown = [name for name in requested if name not in REGISTERED_CHECKS]
        if unknown:
            available = ", ".join(sorted(REGISTERED_CHECKS)) or "<none>"
            raise typer.BadParameter(
                f"Unknown validation(s): {', '.join(unknown)}. Available: {available}",
                param_hint="--checks/-c",
            )
        checks_arg = requested

    try:
        results = run_validations(root, checks=checks_arg, config=config)
    except FileNotFoundError:
        logger.error("Dataset root does not exist: %s", root)
        raise typer.Exit(code=1)
    except ValueError as exc:
        logger.error(str(exc))
        raise typer.Exit(code=1)

    if summary and details:
        raise typer.BadParameter("Cannot combine --summary with --details.", param_hint="--summary/--details")

    display_results = results
    if errors_only:
        display_results = [res for res in results if res.level is Level.ERROR]

    display_block = display_results if display_results else results

    human_text = _format_validation_results(
        display_block,
        root,
        include_info=verbose,
    )
    if not human_text:
        human_text = "No validation findings."

    rendered_text = human_text
    if summary:
        summary_text = _format_validation_summary(
            display_block,
            root,
            include_info=verbose,
            top_codes=summary_limit,
        )
        if not summary_text:
            summary_text = "No validation findings."
        rendered_text = summary_text

    typer.echo(rendered_text)

    summary = summarise(results)
    level_counts = summary.get("levels", {})
    if level_counts:
        counts_line = ", ".join(f"{lvl}={cnt}" for lvl, cnt in level_counts.items())
        typer.echo(f"Level counts: {counts_line}")

    if json_output is not None:
        if out is not None and out != json_output:
            raise typer.BadParameter("--json-output and --out point to different files.")
        out = json_output
        output_format = output_format or "json"

    if out is not None:
        fmt = (output_format or "json").lower()
        if fmt not in {"json", "jsonl", "txt"}:
            raise typer.BadParameter("--format must be 'json', 'jsonl', or 'txt'.")
        if fmt == "txt":
            out.write_text(rendered_text + "\n", encoding="utf-8")
            logger.info("Wrote TXT report to %s", out)
        else:
            write_results(results, out, fmt)  # type: ignore[arg-type]
            logger.info("Wrote %s report to %s", fmt.upper(), out)

    threshold = fail_on.lower()
    if threshold not in {"error", "warn", "none"}:
        raise typer.BadParameter("--fail-on must be one of: error, warn, none")

    exit_code = 0
    if threshold == "error" and any(res.level is Level.ERROR for res in results):
        exit_code = 1
    elif threshold == "warn" and any(res.level in {Level.ERROR, Level.WARN} for res in results):
        exit_code = 1

    raise typer.Exit(code=exit_code)


@app.command()
def run_naming_check(
    root: Path = typer.Option(
        ...,
        "--root",
        help="Path to dataset directory to scan.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    out: Path = typer.Option(
        None,
        "-o",
        "--out",
        help="Optional path to write a human-readable report.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    ids_file: Path = typer.Option(
        None,
        "--ids-file",
        help="Optional file containing AFIDs to check, one per line. Accepts AF-<16>-vN or AF-<16>.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    errors_only: bool = typer.Option(
        False,
        "--errors-only",
        help="Print only failing entries as one-liners.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Print every entry as one-liners.",
    ),
    limit: int = typer.Option(
        5,
        "--limit",
        min=1,
        help="Limit example AFIDs in summary and cap one-liner output when not in --verbose mode.",
    ),
):
    """
    Check AFDB naming conventions and required-file presence for a dataset directory.
    Default output is a HUMAN-READABLE summary. Exit code 0 on PASS, 1 if any issue is found.
    """
    if ids_file is not None:
        ids_pairs, ids_afids = parse_ids_file(ids_file)
        if ids_pairs or ids_afids:
            logger.warning("Filtering by IDs is not currently supported; running all entries.")

    try:
        results = run_validations(root, checks=["naming"], config={})
    except Exception as exc:  # pragma: no cover - defensive
        logger.error("Naming validator failed: %s", exc)
        raise typer.Exit(code=1)

    display = results
    if errors_only:
        display = [res for res in results if res.level is Level.ERROR]

    summary_text = _format_validation_results(
        display if display else results,
        root,
        include_info=verbose,
    )
    if not summary_text:
        summary_text = "No validation findings."

    typer.echo(summary_text)

    if out is not None:
        out.write_text(summary_text + "\n", encoding="utf-8")
        logger.info(f"Wrote report: {out}")

    has_error = any(res.level is Level.ERROR for res in results)
    raise typer.Exit(code=1 if has_error else 0)


@app.command("modelcif-replace")
def run_modelcif_replace(
    pdb: Path = typer.Option(
        ...,
        "-i",
        "--cif",
        help="Input CIF file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    metadata: Path = typer.Option(
        ...,
        "-m",
        "--metadata",
        help="Input metadata JSON file path.",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output mmCIF file path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
):
    """
    Replace the metadata in a CIF file with metadata from a JSON file.
    """
    replace_mmcif_with_json(str(pdb), str(metadata), str(output))
    logger.info("Replaced metadata in %s -> %s", pdb, output)


@app.command("plddt-check")
def run_plddt_check(
    root: Path = typer.Option(
        ...,
        "--root",
        "-r",
        help="Dataset root to scan for AF-*-confidence-*.json files.",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    output: Path = typer.Option(
        None,
        "-o",
        "--output",
        help="Optional path to save the human-readable report.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    errors_only: bool = typer.Option(
        False, "--errors-only", help="Print only failing entries (compact)."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", help="Print all entries (one line each)."
    ),
    limit: int = typer.Option(
        5, "--limit", "-l", help="Max examples to show in summaries."
    ),
    ids: List[str] = typer.Option(
        None,
        "--ids",
        help="Filter by specific AFID/version tokens like 'AF-...-v2'. "
             "You can pass multiple --ids flags.",
    ),
    ids_file: Path = typer.Option(
        None,
        "--ids-file",
        help="Text file of AFID/version, one per line (e.g. AF-...-v2).",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    skip_pae: bool = typer.Option(
        False, "--skip-pae", help="Skip PAE dimension cross-check to save memory/time."
    ),
    bfactor_tolerance: float = typer.Option(
        1.0,
        "--bfactor-tolerance",
        help="Tolerance when comparing structure B-factor min/max to JSON min/max.",
    ),
    with_structure: bool = typer.Option(
        False,
        "--with-structure/--no-with-structure",
        help="Cross-check vs local PDB/mmCIF (residue count and optional B-factor).",
    ),
):
    """
    Validate pLDDT JSON files and summarise results.
    """
    if ids or ids_file:
        logger.warning("Filtering pLDDT checks by IDs is not currently supported; inspecting all files.")
    if skip_pae or with_structure:
        logger.warning("Structure/PAE cross-check flags are not yet implemented in the new validator.")

    try:
        results = run_validations(root, checks=["plddt"], config={})
    except Exception as exc:  # pragma: no cover - defensive
        logger.error("pLDDT validator failed: %s", exc)
        raise typer.Exit(code=1)

    display = results
    if errors_only:
        display = [res for res in results if res.level is Level.ERROR]

    summary_text = _format_validation_results(
        display if display else results,
        root,
        include_info=verbose,
    )
    if not summary_text:
        summary_text = "No validation findings."

    typer.echo(summary_text)
    if output:
        output.write_text(summary_text + "\n", encoding="utf-8")
        logger.info("Wrote report: %s", output)

    has_error = any(res.level is Level.ERROR for res in results)
    raise typer.Exit(code=1 if has_error else 0)


if __name__ == "__main__":
    app()
