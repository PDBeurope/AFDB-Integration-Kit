import json
import tempfile
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest
from jsonschema import ValidationError

from afdb_integration_kit.metadata import validator as schema_validator


# Fixture to temporarily write JSON to a file
@pytest.fixture
def temp_json_file():
    def _write_temp_json(data: dict) -> Path:
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            return Path(f.name)

    return _write_temp_json


@pytest.fixture
def fake_model_schema(tmp_path):
    # Write a minimal model schema for testing
    schema = {
        "type": "object",
        "required": ["name"],
        "properties": {"name": {"type": "string"}},
    }
    schema_path = tmp_path / "test_model_schema.json"
    schema_path.write_text(json.dumps(schema), encoding="utf-8")
    return schema_path


def test_validate_with_overridden_schema(tmp_path, fake_model_schema):
    # Write valid input that matches the test schema
    input_data = {"name": "Test Model"}
    input_file = tmp_path / "valid_input.json"
    input_file.write_text(json.dumps(input_data), encoding="utf-8")

    with patch.dict(
        schema_validator.SCHEMA_PATHS,
        {schema_validator.SchemaType.MODEL: fake_model_schema},
    ):
        schema_validator.validate_against_schema(input_file, "model")


def test_validate_model_batch_with_overridden_schema(tmp_path, fake_model_schema):
    input_data = [{"name": "Test Model"}, {"name": "Second Model"}]
    input_file = tmp_path / "valid_model_batch.json"
    input_file.write_text(json.dumps(input_data), encoding="utf-8")

    with patch.dict(
        schema_validator.SCHEMA_PATHS,
        {schema_validator.SchemaType.MODEL: fake_model_schema},
    ):
        schema_validator.validate_against_schema(input_file, "model")


def test_validate_response_docs_with_overridden_summary_schema(
    tmp_path, fake_model_schema
):
    input_data = {"response": {"docs": [{"name": "Test Model"}]}}
    input_file = tmp_path / "valid_summary_response.json"
    input_file.write_text(json.dumps(input_data), encoding="utf-8")

    with patch.dict(
        schema_validator.SCHEMA_PATHS,
        {schema_validator.SchemaType.MODEL_SUMMARY: fake_model_schema},
    ):
        schema_validator.validate_against_schema(input_file, "model-summary")


def test_invalid_schema_type(temp_json_file):
    input_file = temp_json_file({"foo": "bar"})
    with pytest.raises(ValueError, match="Expected one of"):
        schema_validator.validate_against_schema(input_file, "invalid_type")


def test_invalid_json_file(tmp_path):
    bad_json_path = tmp_path / "bad.json"
    bad_json_path.write_text("{ not: valid json }")

    with pytest.raises(json.JSONDecodeError):
        schema_validator.validate_against_schema(bad_json_path, "model")


def test_schema_validation_error(temp_json_file):
    # Assuming the model schema requires a "name" field
    invalid_data = {"invalid_field": "missing name"}
    input_file = temp_json_file(invalid_data)

    with pytest.raises(ValidationError):
        schema_validator.validate_against_schema(input_file, "model")


def test_load_json_file_success(tmp_path):
    # Create a valid JSON file
    valid_data = {"key": "value"}
    file_path = tmp_path / "valid.json"
    file_path.write_text(json.dumps(valid_data), encoding="utf-8")

    result = schema_validator._load_json_file(file_path)
    assert result == valid_data


def test_load_json_file_not_found():
    non_existent_path = Path("/nonexistent/path/to/file.json")
    with pytest.raises(FileNotFoundError):
        schema_validator._load_json_file(non_existent_path)


def test_load_json_file_invalid_json(tmp_path):
    # Write invalid JSON to a file
    file_path = tmp_path / "invalid.json"
    file_path.write_text("{ invalid json }", encoding="utf-8")

    with pytest.raises(json.JSONDecodeError):
        schema_validator._load_json_file(file_path)


def test_monomer_model_summary_allows_absent_complex_metrics() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    input_file = (
        repo_root
        / "examples/colabfold_monomer_e2e/model_jsons/AF-0000000300000001.json"
    )

    schema_validator.validate_against_schema(input_file, "model-summary")


def test_complex_model_summary_requires_ipsae_metric_block(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    source_file = (
        repo_root
        / "examples/colabfold_complex_e2e/model_jsons/AF-0000000300000101.json"
    )
    payload = json.loads(source_file.read_text(encoding="utf-8"))
    payload.pop("complexPredictionAccuracy_ipsae_BA")

    input_file = tmp_path / "broken_complex_model_summary.json"
    input_file.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValidationError, match="complexPredictionAccuracy_ipsae_BA"):
        schema_validator.validate_against_schema(input_file, "model-summary")


def test_complex_model_summary_batch_validates() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    input_file = (
        repo_root
        / "examples/colabfold_complex_e2e/model_batches/AF-metadata-1-of-1.json"
    )

    schema_validator.validate_against_schema(input_file, "model-summary")


def test_complex_collection_doc_requires_directional_metrics(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    source_file = (
        repo_root
        / "examples/colabfold_complex_e2e/chain_jsons/AF-0000000300000101.json"
    )
    payload = json.loads(source_file.read_text(encoding="utf-8"))
    broken_payload = deepcopy(payload)
    broken_payload[0].pop("complexPredictionAccuracy_ipsae_AB")

    input_file = tmp_path / "broken_complex_chain_doc.json"
    input_file.write_text(json.dumps(broken_payload), encoding="utf-8")

    with pytest.raises(ValidationError, match="complexPredictionAccuracy_ipsae_AB"):
        schema_validator.validate_against_schema(input_file, "collection-doc")


def test_complex_collection_doc_and_batch_validate() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    single_file = (
        repo_root
        / "examples/colabfold_complex_e2e/chain_jsons/AF-0000000300000101.json"
    )
    batch_file = (
        repo_root
        / "examples/colabfold_complex_e2e/chain_batches/AF-chain-metadata-1-of-1.json"
    )

    schema_validator.validate_against_schema(single_file, "collection-doc")
    schema_validator.validate_against_schema(batch_file, "collection-doc")
