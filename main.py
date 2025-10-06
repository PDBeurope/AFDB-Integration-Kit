import logging
from pathlib import Path
from typing import Dict, List, Optional, Set
import json

import typer

from afdb_integration_kit.cif2bcif.convert import (
    run_batch_cif2bcif,
)
from afdb_integration_kit.cif2bcif.convert import run_cif2bcif as cif2bcif_helper
from afdb_integration_kit.dssp.dssp import run_dssp as dssp_helper
from afdb_integration_kit.metadata.validator import validate_against_schema
from afdb_integration_kit.modelcif.generate import generate
from afdb_integration_kit.modelpdb.generate import generate_pdb_headers
from afdb_integration_kit.modelcif_replace.replace import replace_mmcif_with_json as replace_mmcif_with_json
from afdb_integration_kit.quality_assessment.naming import (
    parse_ids_file
)
from afdb_integration_kit.validation import (
    ensure_default_validators,
    list_validators,
    run_validations,
)
from afdb_integration_kit.validation.runner import ValidationResult

# Set up logger
logger = logging.getLogger("afdb_integration_kit")
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)

app = typer.Typer()


def _format_validation_result(
    result: ValidationResult,
    *,
    errors_only: bool = False,
    verbose: bool = False,
    limit: int = 5,
) -> str:
    if result.formatter is None:
        return json.dumps(result.report, indent=2, sort_keys=True)
    return result.formatter(
        result.report,
        errors_only=errors_only,
        verbose=verbose,
        limit=limit,
    )


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
    ensure_default_validators()
    hooks = list_validators()
    if not hooks:
        typer.echo("No validations registered.")
        return

    for hook in hooks:
        suffix = f" - {hook.description}" if hook.description else ""
        typer.echo(f"{hook.name}{suffix}")


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
):
    """
    Enrich a PDB file with metadata to produce a feature-rich mmCIF file.
    """
    # Pre-flight checks
    if not pdb.is_file():
        logger.error(f"Input PDB file not found: {pdb}")
        raise typer.Exit(code=1)
    if not metadata.is_file():
        logger.error(f"Input metadata JSON file not found: {metadata}")
        raise typer.Exit(code=1)

    # If validate is passed as a flag (True but no value), default to 'mmcif_ma.dic'
    validate_path = validate if validate else None
    if validate is not None and validate == "":
        validate_path = "mmcif_ma.dic"

    # Call main logic (assuming main is imported or defined elsewhere)
    generate(str(pdb), str(metadata), str(output), validate_path)

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
    if not cif.is_file():
        logger.error(f"Input mmCIF file not found: {cif}")
        raise typer.Exit(code=1)
    if not pdb.is_file():
        logger.error(f"Input PDB file not found: {pdb}")
        raise typer.Exit(code=1)
    if not provider.is_file():
        logger.error(f"Input provider JSON file not found: {provider}")
        raise typer.Exit(code=1)

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
    check: List[str] = typer.Option(
        None,
        "--check",
        "-c",
        help="Validation name to run. Repeat for multiple; defaults to all registered.",
    ),
    output: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Optional path to write the combined human-readable report.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    json_output: Path = typer.Option(
        None,
        "--json-output",
        help="Optional path to write an aggregate JSON report.",
        file_okay=True,
        dir_okay=False,
        writable=True,
        resolve_path=True,
    ),
    errors_only: bool = typer.Option(
        False,
        "--errors-only",
        help="Print only failing entries for each validation report.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Print all entries for each validation report.",
    ),
    limit: int = typer.Option(
        5,
        "--limit",
        min=1,
        help="Limit example entries in summaries when not in verbose mode.",
    ),
):
    """Run one or more validation checks via the shared validation runner."""
    ensure_default_validators()
    hooks = {hook.name: hook for hook in list_validators()}
    if not hooks:
        logger.error("No validations are registered.")
        raise typer.Exit(code=1)

    selected = check or list(hooks.keys())
    unknown = [name for name in selected if name not in hooks]
    if unknown:
        available = ", ".join(sorted(hooks))
        raise typer.BadParameter(
            f"Unknown validation(s): {', '.join(unknown)}. Available: {available}",
            param_hint="--check/-c",
        )

    ok, results = run_validations(root, checks=selected)

    human_output = ""
    if results:
        chunks = [
            _format_validation_result(
                result,
                errors_only=errors_only,
                verbose=verbose,
                limit=limit,
            )
            for result in results
        ]
        human_output = "\n\n".join(chunks)
        typer.echo(human_output)
    else:
        typer.echo("No validation results produced.")

    aggregate = {
        "dataset_root": str(Path(root).resolve()),
        "overall_ok": ok,
        "results": {
            result.name: {
                "ok": result.ok,
                "description": result.description,
                "options": result.options,
                "report": result.report,
            }
            for result in results
        },
    }

    if output is not None and human_output:
        output.write_text(human_output + "\n", encoding="utf-8")
        logger.info(f"Wrote report: {output}")

    if json_output is not None:
        json_output.write_text(json.dumps(aggregate, indent=2) + "\n", encoding="utf-8")
        logger.info(f"Wrote JSON report: {json_output}")

    raise typer.Exit(code=0 if ok else 1)


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
    ids_pairs = None
    ids_afids = None
    if ids_file is not None:
        ids_pairs, ids_afids = parse_ids_file(ids_file)

    overrides: Dict[str, Dict[str, object]] = {}
    selection_args: Dict[str, object] = {}
    if ids_pairs:
        selection_args["ids_pairs"] = ids_pairs
    if ids_afids:
        selection_args["ids_afids"] = ids_afids
    if selection_args:
        overrides["naming"] = selection_args

    ok, results = run_validations(root, checks=["naming"], overrides=overrides)
    if not results:
        logger.error("Validator 'naming' is not registered.")
        raise typer.Exit(code=1)

    result = results[0]
    text = _format_validation_result(
        result,
        errors_only=errors_only,
        verbose=verbose,
        limit=limit,
    )
    print(text)

    if out is not None:
        out.write_text(text + "\n", encoding="utf-8")
        logger.info(f"Wrote report: {out}")

    if not ok:
        raise typer.Exit(code=1)


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
    # parse --ids / --ids-file into filters
    ids_pairs: Set[tuple[str, str]] = set()
    ids_afids: Set[str] = set()

    def parse_token(tok: str) -> Optional[tuple[str, str]]:
        # Accept tokens like AF-<16>-vN
        import re
        m = re.match(r"^(AF-\d{16})-(v\d+)$", tok.strip())
        return (m.group(1), m.group(2)) if m else None

    if ids:
        for t in ids:
            pv = parse_token(t)
            if pv:
                ids_pairs.add(pv)
    if ids_file:
        for line in ids_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pv = parse_token(line)
            if pv:
                ids_pairs.add(pv)

    overrides: Dict[str, Dict[str, object]] = {
        "plddt": {
            "skip_pae": skip_pae,
            "bfactor_tolerance": bfactor_tolerance,
            "with_structure": with_structure,
        }
    }
    if ids_pairs:
        overrides["plddt"]["ids_pairs"] = ids_pairs
    if ids_afids:
        overrides["plddt"]["ids_afids"] = ids_afids

    ok, results = run_validations(root, checks=["plddt"], overrides=overrides)
    if not results:
        logger.error("Validator 'plddt' is not registered.")
        raise typer.Exit(code=1)

    result = results[0]
    text = _format_validation_result(
        result,
        errors_only=errors_only,
        verbose=verbose,
        limit=limit,
    )
    print(text)
    if output:
        output.write_text(text + "\n", encoding="utf-8")
        logger.info("Wrote report: %s", output)

    raise typer.Exit(code=0 if ok else 1)


if __name__ == "__main__":
    app()
