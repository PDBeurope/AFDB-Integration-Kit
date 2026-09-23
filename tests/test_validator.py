import json
import importlib.util
import sys
import tempfile
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest
import jsonschema
from jsonschema import ValidationError
from jsonschema.exceptions import FormatError

from afdb_integration_kit.metadata import validator as schema_validator
from afdb_integration_kit.complex_metrics import (
    DEFAULT_COMPLEX_ENRICHMENT_METRICS,
    build_chain_enrichment,
    build_model_enrichment,
    parse_ipsae_csv,
)


COMPLEX_MODEL_SUMMARY_FIELDS = (
    "complexName", "assemblyType", "oligomericState",
    "oligomericStateDescription", "complexPredictionAccuracy_ipTM",
    "complexPredictionAccuracy_ipsae_pae_cutoff",
    "complexPredictionAccuracy_ipsae_dist_cutoff",
    "complexPredictionAccuracy_iptm_af", "complexPredictionAccuracy_ipsae_AB",
    "complexPredictionAccuracy_ipsae_BA", "complexPredictionAccuracy_ipsae_d0chn_AB",
    "complexPredictionAccuracy_ipsae_d0chn_BA", "complexPredictionAccuracy_ipsae_d0dom_AB",
    "complexPredictionAccuracy_ipsae_d0dom_BA",
    "complexPredictionAccuracy_ipsae_iptm_d0chn_AB",
    "complexPredictionAccuracy_ipsae_iptm_d0chn_BA",
    "complexPredictionAccuracy_pDockQ2_AB", "complexPredictionAccuracy_pDockQ2_BA",
    "complexPredictionAccuracy_LIS_AB", "complexPredictionAccuracy_LIS_BA",
    "complexPredictionAccuracy_ipsae_n0res_AB", "complexPredictionAccuracy_ipsae_n0res_BA",
    "complexPredictionAccuracy_ipsae_n0dom_AB", "complexPredictionAccuracy_ipsae_n0dom_BA",
    "complexPredictionAccuracy_ipsae_d0res_AB", "complexPredictionAccuracy_ipsae_d0res_BA",
    "complexPredictionAccuracy_ipsae_nres1_AB", "complexPredictionAccuracy_ipsae_nres1_BA",
    "complexPredictionAccuracy_ipsae_nres2_AB", "complexPredictionAccuracy_ipsae_nres2_BA",
    "complexPredictionAccuracy_ipsae_dist_nres1_AB",
    "complexPredictionAccuracy_ipsae_dist_nres1_BA",
    "complexPredictionAccuracy_ipsae_dist_nres2_AB",
    "complexPredictionAccuracy_ipsae_dist_nres2_BA",
    "complexPredictionAccuracy_pDockQ", "complexPredictionAccuracy_ipsae_n0chn",
)
COMPLEX_COLLECTION_FIELDS = (
    "complexComposition", "complexPredictionAccuracy_ipTM",
    "complexPredictionAccuracy_ipsae_pae_cutoff",
    "complexPredictionAccuracy_ipsae_dist_cutoff",
    "complexPredictionAccuracy_iptm_af", "complexPredictionAccuracy_pDockQ",
    "complexPredictionAccuracy_ipsae_n0chn",
)
NUMERIC_COMPLEX_FIELDS = tuple(
    field for field in COMPLEX_MODEL_SUMMARY_FIELDS if field.startswith("complexPredictionAccuracy_")
)


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


def test_metadata_schemas_allow_absent_gene_and_reject_placeholders(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    model_payload = json.loads(
        (repo_root / "tests/fixtures/validation/good_dataset/AF-metadata-1-of-1.json")
        .read_text(encoding="utf-8")
    )[0]
    model_payload.setdefault("isComplex", False)
    validation_cases = [
        ("model", model_payload, [None, "", "   "]),
        ("model-summary", None, [[None], [], [""], ["   "]]),
        ("collection-doc", None, [None, "", "   "]),
    ]

    for schema_type, source_entry, invalid_values in validation_cases:
        if source_entry is None:
            schema_path = Path(
                schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]
            ).with_name(f"{schema_type.replace('-', '_')}_schema.json")
            entry = _complex_instance(json.loads(schema_path.read_text(encoding="utf-8")))
        else:
            entry = deepcopy(source_entry)
        entry.pop("gene", None)

        absent_gene_file = tmp_path / f"{schema_type}-absent-gene.json"
        absent_gene_file.write_text(json.dumps(entry), encoding="utf-8")
        schema_validator.validate_against_schema(absent_gene_file, schema_type)

        for index, invalid_value in enumerate(invalid_values):
            invalid_entry = deepcopy(entry)
            invalid_entry["gene"] = invalid_value
            invalid_gene_file = tmp_path / f"{schema_type}-invalid-gene-{index}.json"
            invalid_gene_file.write_text(json.dumps(invalid_entry), encoding="utf-8")

            with pytest.raises(ValidationError):
                schema_validator.validate_against_schema(invalid_gene_file, schema_type)


@pytest.mark.parametrize("schema_type", ["model", "model-summary", "collection-doc", "provider", "modelcif-metadata"])
def test_schema_types_dispatch_to_distinct_schemas(schema_type):
    assert schema_validator.SCHEMA_PATHS[schema_validator.SchemaType(schema_type)].name == {
        "model": "model_schema.json",
        "model-summary": "model_summary_schema.json",
        "collection-doc": "collection_doc_schema.json",
        "provider": "provider_schema.json",
        "modelcif-metadata": "schema.json",
    }[schema_type]


@pytest.mark.parametrize(
    ("schema_name", "required_fields"),
    [
        ("model_summary_schema.json", COMPLEX_MODEL_SUMMARY_FIELDS),
        ("collection_doc_schema.json", COMPLEX_COLLECTION_FIELDS),
    ],
)
def test_complex_schema_contract_requires_exact_fields(schema_name, required_fields):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    conditional = next(
        clause for clause in schema["allOf"]
        if clause["if"]["properties"].get("isComplex") == {"const": True}
    )
    assert tuple(conditional["then"]["required"]) == required_fields


@pytest.mark.parametrize("schema_name", ["model_summary_schema.json", "collection_doc_schema.json"])
@pytest.mark.parametrize("field", NUMERIC_COMPLEX_FIELDS)
def test_complex_metrics_are_numeric_when_present(schema_name, field):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    field_schema = schema["properties"].get(
        field, schema["patternProperties"]["^complexPredictionAccuracy_[A-Za-z0-9_]+$"]
    )
    assert field_schema["type"] in ("number", ["number", "null"])


@pytest.mark.parametrize("schema_name", ["model_summary_schema.json", "collection_doc_schema.json"])
@pytest.mark.parametrize(
    "field",
    [
        "complexPredictionAccuracy_N_clash_backbone",
        "complexPredictionAccuracy_N_clash_heavyAtom",
    ],
)
def test_clash_metrics_remain_optional_numeric(schema_name, field):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    assert schema["patternProperties"]["^complexPredictionAccuracy_[A-Za-z0-9_]+$"]["type"] in (
        "number", ["number", "null"]
    )
    assert field not in schema.get("required", [])
    assert field not in schema["allOf"][0]["then"]["required"]


def _complex_instance(schema):
    conditional_required = schema["allOf"][0]["then"]["required"]
    instance = {}
    for field in [*schema["required"], *conditional_required]:
        field_schema = schema["properties"].get(field)
        if field_schema is None:
            instance[field] = 1.0
            continue
        field_type = field_schema["type"]
        field_type = next(item for item in field_type if item != "null") if isinstance(field_type, list) else field_type
        instance[field] = {
            "string": "value", "integer": 1, "number": 1.0,
            "boolean": False, "array": ["value"],
        }[field_type]
    valid_overrides = {
        "modelEntityId": "AF-0000000000000001",
        "uniqueId": "AF-0000000000000001_v1_A",
        "modelCreatedDate": "2026-01-01T00:00:00Z",
        "sequenceVersionDate": "2026-01-01T00:00:00Z",
        "sequenceChecksum": "0123456789abcdef0123456789abcdef",
        "sequence": "ACDEFG",
        "complexComposition": (
            ["Q46806_2"] if isinstance(instance.get("complexComposition"), list)
            else "Q46806_2"
        ),
        "assemblyType": "Homo",
        "oligomericState": "dimer",
    }
    instance.update(
        {field: value for field, value in valid_overrides.items() if field in instance}
    )
    for field in ("taxId", "taxId_display"):
        if field in instance and schema["properties"][field]["type"] == "array":
            instance[field] = [1]
    if "allVersions" in instance:
        instance["allVersions"] = [1]
    instance["isComplex"] = True
    # The directional iPSAE metrics are a code-level contract enforced by
    # schema_validator._validate_complex_metric_contract, not by the schema, so
    # a schema-derived instance would otherwise omit them.
    unique_id = instance.get("uniqueId")
    chain_id = unique_id.rsplit("_", 1)[-1] if isinstance(unique_id, str) and "_" in unique_id else None
    for field in schema_validator.COMPLEX_COLLECTION_DIRECTIONAL_REQUIRED_FIELDS.get(chain_id, ()):
        instance.setdefault(field, 1.0)
    return instance


@pytest.mark.parametrize(
    ("schema_name", "required_fields"),
    [
        ("model_summary_schema.json", COMPLEX_MODEL_SUMMARY_FIELDS),
        ("collection_doc_schema.json", COMPLEX_COLLECTION_FIELDS),
    ],
)
def test_each_required_complex_field_is_enforced(schema_name, required_fields):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    valid = _complex_instance(schema)
    jsonschema.validate(valid, schema)
    for field in required_fields:
        broken = deepcopy(valid)
        del broken[field]
        with pytest.raises(ValidationError, match=field):
            jsonschema.validate(broken, schema)


@pytest.mark.parametrize("schema_name", ["model_summary_schema.json", "collection_doc_schema.json"])
@pytest.mark.parametrize("field", NUMERIC_COMPLEX_FIELDS)
def test_each_complex_numeric_metric_rejects_string(schema_name, field):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    valid = _complex_instance(schema)
    if field not in valid:
        return
    valid[field] = "not numeric"
    with pytest.raises(ValidationError):
        jsonschema.validate(valid, schema)


@pytest.mark.parametrize(
    ("schema_name", "field", "bad_value"),
    [
        ("model_summary_schema.json", "modelEntityId", "bad-id"),
        ("model_summary_schema.json", "latestVersion", 0),
        ("model_summary_schema.json", "uniprotAccession", []),
        ("model_summary_schema.json", "taxId", ["9606"]),
        ("model_summary_schema.json", "assemblyType", "Mixed"),
        ("model_summary_schema.json", "providerId", ""),
        ("model_summary_schema.json", "globalMetricValue", 101),
        ("model_summary_schema.json", "oligomericState", "undecimer"),
        ("model_summary_schema.json", "taxId", [0]),
        ("collection_doc_schema.json", "uniqueId", "bad-id"),
        ("collection_doc_schema.json", "uniqueId", "AF-0000000000000001_v0_A"),
        ("collection_doc_schema.json", "modelEntityId", "bad-id"),
        ("collection_doc_schema.json", "modelCreatedDate", "not-a-date"),
        ("collection_doc_schema.json", "sequenceVersionDate", "not-a-date"),
        ("collection_doc_schema.json", "sequenceChecksum", "not-md5"),
        ("collection_doc_schema.json", "sequenceStart", 0),
        ("collection_doc_schema.json", "latestVersion", 0),
        ("collection_doc_schema.json", "allVersions", []),
        ("collection_doc_schema.json", "complexComposition", "Q46806_0"),
        ("collection_doc_schema.json", "assemblyType", "Mixed"),
        ("collection_doc_schema.json", "toolUsed", ""),
        ("collection_doc_schema.json", "providerId", ""),
        ("collection_doc_schema.json", "entityType", ""),
        ("collection_doc_schema.json", "sequence", "ACD*"),
        ("collection_doc_schema.json", "globalMetricValue", -1),
        ("collection_doc_schema.json", "fractionPlddtVeryHigh", 1.1),
        ("collection_doc_schema.json", "allVersions", [0]),
        ("collection_doc_schema.json", "oligomericState", "undecimer"),
        ("collection_doc_schema.json", "taxId", 0),
    ],
)
def test_public_schema_constraints_reject_invalid_values(schema_name, field, bad_value):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    instance = _complex_instance(schema)
    instance[field] = bad_value
    with pytest.raises(ValidationError):
        jsonschema.validate(instance, schema, format_checker=schema_validator.FORMAT_CHECKER)


@pytest.mark.parametrize("schema_name", ["model_summary_schema.json", "collection_doc_schema.json"])
def test_monomer_rejects_complex_only_fields(schema_name):
    schema_path = Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL]).with_name(schema_name)
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    instance = _complex_instance(schema)
    instance["isComplex"] = False
    instance["assemblyType"] = None
    instance["oligomericState"] = None
    with pytest.raises(ValidationError):
        jsonschema.validate(instance, schema)


def test_model_summary_requires_sequence_bounds_as_a_pair():
    schema = json.loads(
        Path(schema_validator.SCHEMA_PATHS[schema_validator.SchemaType.MODEL_SUMMARY]).read_text()
    )
    instance = _complex_instance(schema)
    instance["sequenceStart"] = 1
    instance.pop("sequenceEnd", None)
    with pytest.raises(ValidationError):
        jsonschema.validate(instance, schema)


@pytest.mark.parametrize(
    ("value", "valid"),
    [
        ("2026-01-01T00:00:00", False),
        ("2026-01-01T00:00:00Z", True),
        ("2026-01-01T00:00:00+01:00", True),
        ("2026-01-01T00:00:00.123456-05:30", True),
        ("2026-01-01t00:00:00z", True),
        ("1990-12-31T23:59:60Z", True),
        ("1990-12-31T23:59:60+01:00", True),
        ("1990-12-31T23:59:60+25:00", False),
        ("1990-12-31T23:59:61Z", False),
        ("2026-01-01T24:00:00Z", False),
        ("2026-01-01T00:60:00Z", False),
        ("2026-01-01T00:00:00+24:00", False),
        ("2026-01-01T00:00:00+01:60", False),
        ("2026-02-30T00:00:00Z", False),
        ("2026-01-01 00:00:00Z", False),
        ("not-a-timestamp", False),
        ("２０２６-01-01T00:00:00Z", False),
        ("2026-01-01T０0:00:00Z", False),
        ("2026-01-01T00:00:00.١Z", False),
        ("2026-01-01T00:00:00+０1:00", False),
        (123, False),
    ],
)
def test_rfc3339_datetime_requires_timezone(value, valid):
    if valid:
        schema_validator.FORMAT_CHECKER.check(value, "date-time")
    else:
        with pytest.raises(FormatError):
            schema_validator.FORMAT_CHECKER.check(value, "date-time")


def test_metadata_datetime_checker_does_not_mutate_jsonschema_globals():
    fresh = jsonschema.FormatChecker()
    draft = jsonschema.Draft202012Validator.FORMAT_CHECKER
    assert fresh.checkers["date-time"] == draft.checkers["date-time"]
    assert fresh.checkers["date-time"] != schema_validator.FORMAT_CHECKER.checkers["date-time"]


def _load_exporter(name: str):
    repo_root = Path(__file__).resolve().parent.parent
    path = repo_root / "uniprot" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"{name}_task5", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_homodimer_exporters_emit_schema_valid_enriched_artifacts(tmp_path):
    duckdb = pytest.importorskip("duckdb")
    repo_root = Path(__file__).resolve().parent.parent
    fixture = repo_root / "tests/fixtures/homodimer_reconciliation"
    model_id = "AF-0000000066074510"
    config = json.loads((fixture / "config/dataset_config.json").read_text())
    config.update({
        "globalMetricValue": 90.0, "fractionPlddtVeryLow": 0.0,
        "fractionPlddtLow": 0.0, "fractionPlddtConfident": 0.1,
        "fractionPlddtVeryHigh": 0.9,
    })
    ipsae = parse_ipsae_csv(fixture / "input/ipsae_summary.csv")[model_id]

    artifacts = {}
    for exporter_name, schema_type in [
        ("export_model_metadata", "model-summary"),
        ("export_chain_metadata", "collection-doc"),
    ]:
        exporter = _load_exporter(exporter_name)
        manifest = exporter.load_manifest(fixture / "config/chain_manifest.csv")
        rows = manifest.by_model[model_id]
        con = duckdb.connect(str(fixture / "config/uniprot_example_subset.duckdb"), read_only=True)
        try:
            entries = exporter.fetch_entries(con, ["Q46806"])
        finally:
            con.close()
        model_meta = {
            model_id: exporter.ModelMetadataRow(
                iptm=0.83, average_plddt=90.0,
                complex_name="Homodimer of Q46806", is_am_data=False,
            )
        }
        records = exporter.build_record(model_id, config, rows, entries, model_meta)
        records = records if isinstance(records, list) else [records]
        for record in records:
            if schema_type == "model-summary":
                record.update(build_model_enrichment(ipsae, {}, DEFAULT_COMPLEX_ENRICHMENT_METRICS))
            else:
                chain_id = record["uniqueId"].rsplit("_", 1)[-1]
                record.update(build_chain_enrichment(ipsae, chain_id, DEFAULT_COMPLEX_ENRICHMENT_METRICS))
        artifact = tmp_path / f"{schema_type}.json"
        artifact.write_text(json.dumps(records), encoding="utf-8")
        schema_validator.validate_against_schema(artifact, schema_type)
        artifacts[schema_type] = records

    assert all(field in artifacts["model-summary"][0] for field in COMPLEX_MODEL_SUMMARY_FIELDS)
    assert all(
        field in record
        for record in artifacts["collection-doc"]
        for field in COMPLEX_COLLECTION_FIELDS
    )


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
