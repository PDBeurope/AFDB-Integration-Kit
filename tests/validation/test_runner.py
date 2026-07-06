from __future__ import annotations

import json
from pathlib import Path

import afdb_integration_kit.validation.validators  # noqa: F401
from typer.testing import CliRunner

from afdb_integration_kit.validation import Level, REGISTERED_CHECKS, run_validations
from main import app


FIXTURES_DIR = Path("tests/fixtures/validation")


def test_registry_contains_default_validators() -> None:
    assert "naming" in REGISTERED_CHECKS
    assert "plddt" in REGISTERED_CHECKS
    assert "relationships" in REGISTERED_CHECKS
    assert "sequences" in REGISTERED_CHECKS
    assert "metadata" in REGISTERED_CHECKS


def test_run_validations_respects_checks() -> None:
    dataset = FIXTURES_DIR / "good_dataset"
    results = run_validations(dataset, checks=["naming"], config={})
    assert results, "Expected naming validator to produce results"
    assert all(result.check == "naming" for result in results)


def test_cli_exit_codes_for_plddt() -> None:
    runner = CliRunner()
    good_dataset = FIXTURES_DIR / "good_dataset"
    bad_dataset = FIXTURES_DIR / "bad_dataset"

    ok_result = runner.invoke(
        app,
        [
            "run-validations",
            "--root",
            str(good_dataset),
            "--checks",
            "plddt",
            "--fail-on",
            "error",
        ],
    )
    assert ok_result.exit_code == 0, ok_result.stdout

    bad_result = runner.invoke(
        app,
        [
            "run-validations",
            "--root",
            str(bad_dataset),
            "--checks",
            "plddt",
            "--fail-on",
            "error",
        ],
    )
    assert bad_result.exit_code == 1, bad_result.stdout
    assert "Level counts" in bad_result.stdout


def test_cli_summary_and_details_output(tmp_path) -> None:
    runner = CliRunner()
    bad_dataset = tmp_path / "bad_dataset"
    bad_dataset.mkdir()
    (bad_dataset / "AF-metadata-1-of-1.json").write_text(
        json.dumps([{"uniqueId": 123}]), encoding="utf-8"
    )

    summary_result = runner.invoke(
        app,
        [
            "run-validations",
            "--root",
            str(bad_dataset),
            "--checks",
            "metadata",
            "--fail-on",
            "error",
            "--summary",
            "--summary-limit",
            "2",
        ],
    )
    assert summary_result.exit_code == 1, summary_result.stdout
    assert "Files with findings:" in summary_result.stdout
    assert "AF-metadata-1-of-1.json" in summary_result.stdout
    assert "metadata_schema_validation_error" in summary_result.stdout

    txt_path = tmp_path / "details.txt"
    details_result = runner.invoke(
        app,
        [
            "run-validations",
            "--root",
            str(bad_dataset),
            "--checks",
            "metadata",
            "--fail-on",
            "error",
            "--details",
            "--out",
            str(txt_path),
            "--format",
            "txt",
        ],
    )
    assert details_result.exit_code == 1, details_result.stdout
    assert "ERROR" in details_result.stdout
    assert txt_path.exists()


def test_validate_metadata_file_cli_uses_schema_validator() -> None:
    runner = CliRunner()
    metadata_file = FIXTURES_DIR / "good_dataset" / "AF-metadata-1-of-1.json"

    result = runner.invoke(
        app,
        [
            "validate-metadata-file",
            "--file",
            str(metadata_file),
            "--type",
            "model",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "Validated metadata file against the 'model' schema" in result.stdout


def test_validate_metadata_file_cli_requires_type() -> None:
    runner = CliRunner()
    metadata_file = FIXTURES_DIR / "good_dataset" / "AF-metadata-1-of-1.json"

    result = runner.invoke(
        app,
        [
            "validate-metadata-file",
            "--file",
            str(metadata_file),
        ],
    )

    assert result.exit_code != 0
    assert "--type" in result.output


def test_validate_metadata_file_cli_rejects_unknown_type() -> None:
    runner = CliRunner()
    metadata_file = FIXTURES_DIR / "good_dataset" / "AF-metadata-1-of-1.json"

    result = runner.invoke(
        app,
        [
            "validate-metadata-file",
            "--file",
            str(metadata_file),
            "--type",
            "unknown",
        ],
    )

    assert result.exit_code == 1, result.stdout
    assert "Unknown metadata schema type 'unknown'" in result.stdout
    assert "Expected one of:" in result.stdout


def test_validate_metadata_file_cli_rejects_wrong_schema_type() -> None:
    runner = CliRunner()
    metadata_file = FIXTURES_DIR / "good_dataset" / "AF-metadata-1-of-1.json"

    result = runner.invoke(
        app,
        [
            "validate-metadata-file",
            "--file",
            str(metadata_file),
            "--type",
            "provider",
        ],
    )

    assert result.exit_code == 1, result.stdout
    assert "provider" in result.stdout
    assert "Update the metadata file to satisfy" in result.stdout


def test_run_validations_metadata_uses_schema_validator_for_good_fixture() -> None:
    runner = CliRunner()
    good_dataset = FIXTURES_DIR / "good_dataset"

    result = runner.invoke(
        app,
        [
            "run-validations",
            "--root",
            str(good_dataset),
            "--checks",
            "metadata",
            "--fail-on",
            "error",
            "--verbose",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "Validated metadata file against the 'model' schema" in result.stdout
    assert "metadata_invalid_type" not in result.stdout
    assert "metadata_missing_field" not in result.stdout


def test_plddt_additional_checks(tmp_path) -> None:
    dataset = tmp_path / "dataset"
    dataset.mkdir()
    payload = {
        "residueNumber": [1, 2, 3],
        "confidenceScore": [95.0, 65.123, 25.0],
        "confidenceCategory": ["M", "L", "D"],
    }
    (dataset / "AF-0000000000000001-confidence_v1.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )

    results = run_validations(dataset, checks=["plddt"])
    codes = {res.code for res in results if res.level is Level.ERROR}
    assert "plddt_decimal_precision" in codes
    assert "plddt_category_mismatch" in codes
def test_pae_validator(tmp_path) -> None:
    dataset = tmp_path / "dataset"
    dataset.mkdir()
    payload = [
        {
            "predicted_aligned_error": [
                [0.25, 1.26],
                [1.45, 0.25],
            ],
            "max_predicted_aligned_error": 1.45,
        }
    ]
    pae_path = dataset / "AF-0000000000000001-predicted_aligned_error_v1.json"
    pae_path.write_text(json.dumps(payload), encoding="utf-8")

    results = run_validations(dataset, checks=["pae"])
    assert any(res.code == "pae_summary" for res in results)

    # Test with explicitly non-square matrix (3x2 instead of ragged array)
    broken = [
        {
            "predicted_aligned_error": [
                [0.25, 1.26],
                [1.45, 0.25],
                [0.50, 0.75],
            ],
            "max_predicted_aligned_error": 1.45,
        }
    ]
    pae_path.write_text(json.dumps(broken), encoding="utf-8")
    results = run_validations(dataset, checks=["pae"])
    # Should get either matrix_not_square or non_numeric_value (ragged arrays may fail numpy conversion)
    error_codes = {res.code for res in results if res.level is Level.ERROR}
    assert "pae_matrix_not_square" in error_codes or "pae_non_numeric_value" in error_codes


def test_relationship_validator(tmp_path) -> None:
    dataset = tmp_path / "dataset"
    dataset.mkdir()
    # Write valid pLDDT and PAE
    plddt_payload = {
        "residueNumber": [1, 2, 3],
        "confidenceScore": [90.0, 80.0, 50.0],
        "confidenceCategory": ["H", "H", "M"],
    }
    pae_payload = [
        {
            "predicted_aligned_error": [
                [0.1, 0.2, 0.3],
                [0.2, 0.1, 0.4],
                [0.3, 0.4, 0.1],
            ],
            "max_predicted_aligned_error": 0.4,
        }
    ]
    (dataset / "AF-0000000000000001-confidence_v1.json").write_text(json.dumps(plddt_payload), encoding="utf-8")
    (dataset / "AF-0000000000000001-predicted_aligned_error_v1.json").write_text(json.dumps(pae_payload), encoding="utf-8")

    results = run_validations(dataset, checks=["relationships"])
    assert any(res.code == "relationship_summary" for res in results)

    # Break the PAE dimension
    plddt_payload["confidenceScore"] = [90.0, 80.0]
    plddt_payload["residueNumber"] = [1, 2]
    plddt_payload["confidenceCategory"] = ["H", "H"]
    (dataset / "AF-0000000000000001-confidence_v1.json").write_text(json.dumps(plddt_payload), encoding="utf-8")
    results = run_validations(dataset, checks=["relationships"])
    assert any(res.code == "relationship_length_mismatch" for res in results if res.level is Level.ERROR)


def test_sequences_validator(tmp_path) -> None:
    dataset = tmp_path / "dataset"
    dataset.mkdir()

    fasta_path = dataset / "sequences.fasta"
    fasta_path.write_text(">AFDB:AF-0000000000000001\nACDEFGHIKL\n", encoding="utf-8")

    results = run_validations(dataset, checks=["sequences"])
    assert any(res.code == "sequences_summary" for res in results)

    # Introduce invalid header and character
    fasta_path.write_text(
        ">AFDB:AF-0000000000000001\nACDEFGHIKL\n>BADHEADER\nACDEFGHIKL\n>AFDB:AF-0000000000000002\nACDEFGHIK1\n",
        encoding="utf-8",
    )
    results = run_validations(dataset, checks=["sequences"])
    error_codes = {res.code for res in results if res.level is Level.ERROR}
    assert "sequences_invalid_header" in error_codes
    assert "sequences_invalid_characters" in error_codes


def test_metadata_validator(tmp_path) -> None:
    dataset = tmp_path / "dataset"
    dataset.mkdir()

    fixture_path = FIXTURES_DIR / "good_dataset" / "AF-metadata-1-of-1.json"
    entry = json.loads(fixture_path.read_text(encoding="utf-8"))[0]
    metadata_path = dataset / "AF-metadata-1-of-1.json"
    metadata_path.write_text(json.dumps([entry]), encoding="utf-8")

    results = run_validations(dataset, checks=["metadata"])
    assert any(res.code == "metadata_schema_valid" for res in results)

    # Test missing required field
    bad_entry = dict(entry)
    bad_entry.pop("latestVersion")
    metadata_path.write_text(json.dumps([bad_entry]), encoding="utf-8")
    results = run_validations(dataset, checks=["metadata"])
    assert any(
        res.code == "metadata_schema_validation_error"
        for res in results
        if res.level is Level.ERROR
    )

    # Test invalid type against the shared model schema
    bad_type_entry = dict(entry)
    bad_type_entry["isUniProt"] = "all"
    metadata_path.write_text(json.dumps([bad_type_entry]), encoding="utf-8")
    results = run_validations(dataset, checks=["metadata"])
    assert any(
        res.code == "metadata_schema_validation_error"
        for res in results
        if res.level is Level.ERROR
    )
