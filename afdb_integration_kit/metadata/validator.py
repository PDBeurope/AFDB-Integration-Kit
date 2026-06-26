import logging
from enum import Enum
from importlib.resources import files
from pathlib import Path

import jsonschema
import orjson

# --- Logger setup ---
logger = logging.getLogger("schema_validator")
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
logger.setLevel(logging.INFO)


COMPLEX_MODEL_SUMMARY_REQUIRED_FIELDS = (
    "complexPredictionAccuracy_ipTM",
    "complexPredictionAccuracy_ipsae_pae_cutoff",
    "complexPredictionAccuracy_ipsae_dist_cutoff",
    "complexPredictionAccuracy_iptm_af",
    "complexPredictionAccuracy_ipsae_AB",
    "complexPredictionAccuracy_ipsae_BA",
    "complexPredictionAccuracy_ipsae_d0chn_AB",
    "complexPredictionAccuracy_ipsae_d0chn_BA",
    "complexPredictionAccuracy_ipsae_d0dom_AB",
    "complexPredictionAccuracy_ipsae_d0dom_BA",
    "complexPredictionAccuracy_ipsae_iptm_d0chn_AB",
    "complexPredictionAccuracy_ipsae_iptm_d0chn_BA",
    "complexPredictionAccuracy_pDockQ2_AB",
    "complexPredictionAccuracy_pDockQ2_BA",
    "complexPredictionAccuracy_LIS_AB",
    "complexPredictionAccuracy_LIS_BA",
    "complexPredictionAccuracy_ipsae_n0res_AB",
    "complexPredictionAccuracy_ipsae_n0res_BA",
    "complexPredictionAccuracy_ipsae_n0dom_AB",
    "complexPredictionAccuracy_ipsae_n0dom_BA",
    "complexPredictionAccuracy_ipsae_d0res_AB",
    "complexPredictionAccuracy_ipsae_d0res_BA",
    "complexPredictionAccuracy_ipsae_nres1_AB",
    "complexPredictionAccuracy_ipsae_nres1_BA",
    "complexPredictionAccuracy_ipsae_nres2_AB",
    "complexPredictionAccuracy_ipsae_nres2_BA",
    "complexPredictionAccuracy_ipsae_dist_nres1_AB",
    "complexPredictionAccuracy_ipsae_dist_nres1_BA",
    "complexPredictionAccuracy_ipsae_dist_nres2_AB",
    "complexPredictionAccuracy_ipsae_dist_nres2_BA",
    "complexPredictionAccuracy_pDockQ",
    "complexPredictionAccuracy_ipsae_n0chn",
)

COMPLEX_COLLECTION_COMMON_REQUIRED_FIELDS = (
    "complexComposition",
    "complexPredictionAccuracy_ipTM",
    "complexPredictionAccuracy_ipsae_pae_cutoff",
    "complexPredictionAccuracy_ipsae_dist_cutoff",
    "complexPredictionAccuracy_iptm_af",
    "complexPredictionAccuracy_pDockQ",
    "complexPredictionAccuracy_ipsae_n0chn",
)

COMPLEX_COLLECTION_DIRECTIONAL_REQUIRED_FIELDS = {
    "A": (
        "complexPredictionAccuracy_ipsae_AB",
        "complexPredictionAccuracy_ipsae_d0chn_AB",
        "complexPredictionAccuracy_ipsae_d0dom_AB",
        "complexPredictionAccuracy_ipsae_iptm_d0chn_AB",
        "complexPredictionAccuracy_pDockQ2_AB",
        "complexPredictionAccuracy_LIS_AB",
        "complexPredictionAccuracy_ipsae_n0res_AB",
        "complexPredictionAccuracy_ipsae_n0dom_AB",
        "complexPredictionAccuracy_ipsae_d0res_AB",
        "complexPredictionAccuracy_ipsae_nres1_AB",
        "complexPredictionAccuracy_ipsae_nres2_AB",
        "complexPredictionAccuracy_ipsae_dist_nres1_AB",
        "complexPredictionAccuracy_ipsae_dist_nres2_AB",
    ),
    "B": (
        "complexPredictionAccuracy_ipsae_BA",
        "complexPredictionAccuracy_ipsae_d0chn_BA",
        "complexPredictionAccuracy_ipsae_d0dom_BA",
        "complexPredictionAccuracy_ipsae_iptm_d0chn_BA",
        "complexPredictionAccuracy_pDockQ2_BA",
        "complexPredictionAccuracy_LIS_BA",
        "complexPredictionAccuracy_ipsae_n0res_BA",
        "complexPredictionAccuracy_ipsae_n0dom_BA",
        "complexPredictionAccuracy_ipsae_d0res_BA",
        "complexPredictionAccuracy_ipsae_nres1_BA",
        "complexPredictionAccuracy_ipsae_nres2_BA",
        "complexPredictionAccuracy_ipsae_dist_nres1_BA",
        "complexPredictionAccuracy_ipsae_dist_nres2_BA",
    ),
}


# --- Constants ---
class SchemaType(Enum):
    MODEL = "model"
    MODEL_SUMMARY = "model-summary"
    COLLECTION_DOC = "collection-doc"
    PROVIDER = "provider"


SCHEMA_PATHS = {
    SchemaType.MODEL: files("afdb_integration_kit.metadata.resources").joinpath(
        "model_schema.json"
    ),
    SchemaType.MODEL_SUMMARY: files("afdb_integration_kit.metadata.resources").joinpath(
        "model_summary_schema.json"
    ),
    SchemaType.COLLECTION_DOC: files("afdb_integration_kit.metadata.resources").joinpath(
        "collection_doc_schema.json"
    ),
    SchemaType.PROVIDER: files("afdb_integration_kit.metadata.resources").joinpath(
        "provider_schema.json"
    ),
}


# --- Utility Functions ---
def _load_json_file(file_path: Path) -> dict:
    """Load a JSON file and return its content as a dict."""
    try:
        return orjson.loads(file_path.read_bytes())
    except FileNotFoundError:
        logger.error("File not found: %s", file_path)
        raise
    except orjson.JSONDecodeError as e:
        logger.error("Invalid JSON in file %s: %s", file_path, e)
        raise


def _instances_to_validate(data, schema_enum: SchemaType):
    if schema_enum is SchemaType.PROVIDER:
        return [data]
    if (
        schema_enum in {SchemaType.MODEL_SUMMARY, SchemaType.COLLECTION_DOC}
        and isinstance(data, dict)
        and isinstance(data.get("response"), dict)
        and isinstance(data["response"].get("docs"), list)
    ):
        return data["response"]["docs"]
    if isinstance(data, list):
        return data
    return [data]


def _validate_complex_metric_contract(
    entry: dict,
    schema_enum: SchemaType,
    index: int,
) -> None:
    if entry.get("isComplex") is not True:
        return

    if schema_enum is SchemaType.MODEL_SUMMARY:
        missing_fields = [
            field
            for field in COMPLEX_MODEL_SUMMARY_REQUIRED_FIELDS
            if field not in entry
        ]
        if missing_fields:
            raise jsonschema.ValidationError(
                f"entry #{index}: complex model summary is missing required iPSAE metrics: "
                + ", ".join(missing_fields)
            )
        return

    if schema_enum is not SchemaType.COLLECTION_DOC:
        return

    common_missing = [
        field
        for field in COMPLEX_COLLECTION_COMMON_REQUIRED_FIELDS
        if field not in entry
    ]
    if common_missing:
        raise jsonschema.ValidationError(
            f"{entry.get('uniqueId') or f'entry #{index}'}: complex collection doc is missing required "
            f"common iPSAE metrics: {', '.join(common_missing)}"
        )

    unique_id = entry.get("uniqueId")
    chain_id = unique_id.rsplit("_", 1)[-1] if isinstance(unique_id, str) and "_" in unique_id else None
    directional_fields = COMPLEX_COLLECTION_DIRECTIONAL_REQUIRED_FIELDS.get(chain_id)
    if directional_fields is None:
        return

    directional_missing = [field for field in directional_fields if field not in entry]
    if directional_missing:
        raise jsonschema.ValidationError(
            f"{unique_id}: complex collection doc for chain {chain_id} is missing required directional "
            f"iPSAE metrics: {', '.join(directional_missing)}"
        )


def validate_against_schema(input_file: Path, schema_type: str):
    """
    Validate the input JSON file against the specified schema.

    Args:
        input_file (Path): Path to the JSON file to validate.
        schema_type (str): Either 'model' or 'provider'.

    Raises:
        ValueError: If the schema_type is not valid.
        FileNotFoundError, JSONDecodeError,
        jsonschema.ValidationError: For specific errors.
    """
    try:
        schema_enum = SchemaType(schema_type.lower())
    except ValueError:
        expected = ", ".join(schema.value for schema in SchemaType)
        logger.error(
            "Unknown schema type '%s'. Expected one of: %s.", schema_type, expected
        )
        raise ValueError(
            f"Unknown schema type '{schema_type}'. Expected one of: {expected}."
        )

    schema_path = SCHEMA_PATHS[schema_enum]

    logger.debug("Using schema: %s", schema_path)
    logger.debug("Validating file: %s", input_file)

    schema = _load_json_file(schema_path)
    data = _load_json_file(input_file)

    try:
        instances = _instances_to_validate(data, schema_enum)
        if not instances:
            raise jsonschema.ValidationError("No documents found to validate")
        for index, entry in enumerate(instances, start=1):
            jsonschema.validate(instance=entry, schema=schema)
            _validate_complex_metric_contract(entry, schema_enum, index)
        logger.info(
            "Validation successful for '%s' against schema '%s'",
            input_file.name,
            schema_enum.value,
        )
    except jsonschema.ValidationError as e:
        logger.error("Validation error: %s", e.message)
        raise
    except Exception as e:
        logger.exception("Unexpected error during validation: %s", e)
        raise
