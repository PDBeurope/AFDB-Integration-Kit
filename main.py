import json
import logging
from collections import defaultdict
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

import typer

from afdb_integration_kit.cif2bcif.convert import (
    run_batch_cif2bcif,
)
from afdb_integration_kit.cif2bcif.convert import run_cif2bcif as cif2bcif_helper
from afdb_integration_kit.dssp.dssp import run_dssp as dssp_helper
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
    run_validations,
    summarise,
    write_results,
)
from afdb_integration_kit.utils.file_guard import require_non_empty_file

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
    Convert CIF to BinaryCIF or BinaryCIF.GZ
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
