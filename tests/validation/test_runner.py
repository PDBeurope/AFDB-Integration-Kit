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
    bad_dataset = FIXTURES_DIR / "bad_dataset"

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
    # The bad_dataset has type mismatches (bool instead of string for isUniProt/isFragment)
    assert "metadata_invalid_type" in summary_result.stdout or "metadata_missing_field" in summary_result.stdout

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

    # Entry with correct types per the validator schema:
    # - isUniProt and isFragment must be strings ("all", "none", "mixed")
    # - isComplex must be boolean
    # - List fields (accession, uniprotId, etc.) must be non-empty lists
    entry = {
        "uniqueId": "AF-0000000000000001_v1",
        "toolUsed": "AlphaFold",
        "modelCreatedDate": "2024-01-01T00:00:00Z",
        "modelEntityId": "AF-0000000000000001",
        "providerId": "DB1",
        "isComplex": False,
        "organismScientificName": "Homo sapiens",
        "isFragment": "none",
        "isUniProt": "all",
        "globalMetricValue": 75.5,
        "fractionPlddtVeryLow": 0.1,
        "fractionPlddtLow": 0.2,
        "fractionPlddtConfident": 0.3,
        "fractionPlddtVeryHigh": 0.4,
        "latestVersion": 1,
        "allVersions": [1],
        # List fields (per-component)
        "accession": ["P01234"],
        "uniprotId": ["P01234"],
        "description": ["Hypothetical protein description."],
        "taxId": [9606],
        "sequence": ["ACDEFGHIKL"],
        "sequenceChecksum": ["0123456789abcdef0123456789abcdef"],
        "sequenceVersionDate": ["2023-01-01T00:00:00Z"],
        "sequenceStart": [1],
        "sequenceEnd": [10],
        "entityType": ["protein"],
        "isIsoform": [False],
    }

    metadata_path = dataset / "AF-metadata-1-of-1.json"
    metadata_path.write_text(json.dumps([entry]), encoding="utf-8")

    results = run_validations(dataset, checks=["metadata"])
    assert any(res.code == "metadata_summary" for res in results)

    # Test missing required field
    bad_entry = dict(entry)
    bad_entry.pop("latestVersion")
    metadata_path.write_text(json.dumps([bad_entry]), encoding="utf-8")
    results = run_validations(dataset, checks=["metadata"])
    assert any(res.code == "metadata_missing_field" for res in results if res.level is Level.ERROR)

    # Test invalid type for list field
    bad_type_entry = dict(entry)
    bad_type_entry["uniprotId"] = "P01234"  # Should be a list, not a string
    metadata_path.write_text(json.dumps([bad_type_entry]), encoding="utf-8")
    results = run_validations(dataset, checks=["metadata"])
    assert any(res.code == "metadata_invalid_type" for res in results if res.level is Level.ERROR)
