import json
import tempfile
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


def test_invalid_schema_type(temp_json_file):
    input_file = temp_json_file({"foo": "bar"})
    with pytest.raises(ValueError, match="Unknown schema type"):
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
