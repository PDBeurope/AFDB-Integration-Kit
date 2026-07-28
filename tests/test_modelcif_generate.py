from __future__ import annotations

from pathlib import Path
from typing import Any

import orjson
import pytest

import afdb_integration_kit.modelcif.generate as modelcif_generate
from afdb_integration_kit.modelcif.provenance import normalize_modelcif_provenance


def test_generate_rejects_missing_alphafold_version_without_top_level_fallback(
    monkeypatch, tmp_path: Path
) -> None:
    input_metadata = {
        "metadata": {"version": "3.1"},
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1"],
                "name": ["AlphaFold"],
                "type": ["package"],
                "description": ["Structure prediction"],
                "classification": ["model building"],
            }
        },
        "chains": [],
    }

    monkeypatch.setattr(modelcif_generate, "load_json_file", lambda path: input_metadata)

    with pytest.raises(ValueError, match="_software.version"):
        modelcif_generate.generate(
            "input.pdb",
            "metadata.json",
            str(tmp_path / "output.cif"),
            validate_dict_path="",
            skip_validation=True,
        )


def test_generate_normalizes_monomer_provenance_without_ipsae(
    monkeypatch, tmp_path: Path
) -> None:
    captured_stores: list[Any] = []

    class FakeCifDataStorage:
        def __init__(self) -> None:
            self.data = {
                "_atom_site": {
                    "B_iso_or_equiv": ["10.0", "20.0"],
                    "label_asym_id": ["A", "A"],
                    "label_comp_id": ["ALA", "ALA"],
                    "label_seq_id": ["1", "1"],
                }
            }
            captured_stores.append(self)

        def populate_from_cif_block(self, cif_block: object) -> None:
            return None

        def set_items(self, category_name: str, items_dict: dict[str, list[Any]]) -> None:
            self.data.setdefault(category_name, {}).update(items_dict)

        def set_item(self, category_name: str, item_name: str, item_value: Any) -> None:
            self.data.setdefault(category_name, {})[item_name] = item_value

        def get_data(self) -> dict[str, dict[str, list[Any]]]:
            return self.data

        def write_to_cif(
            self, output_file: str, block_name: str = "model", skip_alignment: bool = False
        ) -> None:
            Path(output_file).write_text("data_test\n", encoding="utf-8")

    input_metadata = {
        "metadata": {"version": "3.1"},
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1"],
                "name": ["AlphaFold"],
                "version": ["2.3.2"],
                "type": ["package"],
                "description": ["Structure prediction"],
                "classification": ["model building"],
            }
        },
        "chains": [],
    }

    monkeypatch.setattr(modelcif_generate, "load_json_file", lambda path: input_metadata)
    monkeypatch.setattr(modelcif_generate, "pdb_to_cif_block", lambda path: object())
    monkeypatch.setattr(modelcif_generate, "CifDataStorage", FakeCifDataStorage)
    monkeypatch.setattr(modelcif_generate, "map_entities_and_chains", lambda *args: None)
    monkeypatch.setattr(modelcif_generate, "add_standard_chem_comp_data", lambda *args: None)
    monkeypatch.setattr(modelcif_generate, "compute_global_plddt", lambda values: -1.0)
    monkeypatch.setattr(modelcif_generate, "compute_local_plddt_metrics", lambda *args: {})
    monkeypatch.setattr(modelcif_generate, "create_polymer_sequence_categories", lambda *args: {})
    monkeypatch.setattr(modelcif_generate, "_clamp_struct_ref_seq_to_entity_poly_seq", lambda *args: None)

    modelcif_generate.generate(
        "input.pdb",
        "metadata.json",
        str(tmp_path / "output.cif"),
        validate_dict_path="",
        skip_validation=True,
    )

    software = captured_stores[0].data["_software"]
    assert software["name"] == ["AlphaFold", "DSSP"]
    assert software["version"] == ["2.3.2", "?"]
    assert "_ma_software_parameter" not in captured_stores[0].data


def test_validate_json_requires_software_version(monkeypatch) -> None:
    input_path = (
        Path(__file__).resolve().parent.parent
        / "examples/colabfold_monomer_e2e/modelcif_input/AF-0000000300000001.json"
    )
    input_metadata = orjson.loads(input_path.read_bytes())
    input_metadata["categories"]["_software"].pop("version")

    monkeypatch.setattr(modelcif_generate, "_SCHEMA_CACHE", None)

    with pytest.raises(SystemExit):
        modelcif_generate.validate_json_with_schema(
            input_metadata, modelcif_generate.JSON_SCHEMA_PATH
        )


def test_normalize_modelcif_provenance_strips_ipsae_for_monomer() -> None:
    payload = {
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1", "2", "3", "4", "5"],
                "name": ["ColabFold", "MMseqs2-GPU", "AlphaFold-Multimer", "Interface QA", "PyDSSP"],
                "version": ["example-1.0", "example-2.0", "2.3.2", "example-1.0", "example-1.0"],
                "type": ["package", "package", "package", "package", "library"],
                "description": ["a", "b", "c", "d", "e"],
                "classification": ["model building", "data collection", "model building", "data processing", "data extraction"],
            }
        },
        "chains": [{"chain_id": "A", "entity_id": "1", "uniprot_accession": "P1"}],
    }

    normalize_modelcif_provenance(payload, dssp_algorithm="mkdssp")

    assert payload["categories"]["_software"]["name"] == ["AlphaFold", "DSSP"]
    assert payload["categories"]["_software"]["version"] == ["2.3.2", "?"]
    assert "_ma_software_parameter" not in payload["categories"]
    assert payload["categories"]["_ma_protocol_step"]["step_name"] == [
        "model inference",
        "secondary structure assignment",
    ]


def test_normalize_modelcif_provenance_missing_alphafold_version_raises() -> None:
    payload = {
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1"],
                "name": ["AlphaFold-Multimer"],
                "type": ["package"],
                "description": ["legacy"],
                "classification": ["model building"],
            }
        },
        "chains": [
            {"chain_id": "A", "entity_id": "1", "uniprot_accession": "P1"},
            {"chain_id": "B", "entity_id": "1", "uniprot_accession": "P1"},
        ],
    }

    with pytest.raises(ValueError, match="_software.version"):
        normalize_modelcif_provenance(payload, dssp_algorithm="mkdssp")


def test_normalize_modelcif_provenance_can_default_version_for_export() -> None:
    payload = {
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1"],
                "name": ["AlphaFold-Multimer"],
                "type": ["package"],
                "description": ["legacy"],
                "classification": ["model building"],
            }
        },
        "chains": [
            {"chain_id": "A", "entity_id": "1", "uniprot_accession": "P1"},
            {"chain_id": "B", "entity_id": "1", "uniprot_accession": "P1"},
        ],
    }

    normalize_modelcif_provenance(
        payload,
        dssp_algorithm="mkdssp",
        allow_default_alphafold_version=True,
    )

    assert payload["categories"]["_software"]["version"] == ["2.3.2", "?", "?"]


def test_normalize_modelcif_provenance_adds_ipsae_for_complex() -> None:
    payload = {
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1"],
                "name": ["AlphaFold-Multimer"],
                "version": ["2.3.2"],
                "type": ["package"],
                "description": ["legacy"],
                "classification": ["model building"],
            },
            "_ma_software_parameter": {
                "parameter_id": ["1", "2"],
                "group_id": ["2", "2"],
                "data_type": ["float", "float"],
                "name": ["pae_cutoff", "dist_cutoff"],
                "value": ["12.5", "18.0"],
                "description": ["old", "old"],
            },
        },
        "chains": [
            {"chain_id": "A", "entity_id": "1", "uniprot_accession": "P1"},
            {"chain_id": "B", "entity_id": "1", "uniprot_accession": "P1"},
        ],
    }

    normalize_modelcif_provenance(payload, dssp_algorithm="mkdssp")

    assert payload["categories"]["_software"]["name"] == [
        "AlphaFold-Multimer",
        "ipSAE",
        "DSSP",
    ]
    assert payload["categories"]["_software"]["version"] == ["2.3.2", "?", "?"]
    assert payload["categories"]["_ma_software_parameter"]["name"] == [
        "pae_cutoff",
        "dist_cutoff",
    ]
    assert payload["categories"]["_ma_software_parameter"]["value"] == ["12.5", "18.0"]
    assert payload["categories"]["_ma_protocol_step"]["step_name"] == [
        "model inference",
        "interface scoring",
        "secondary structure assignment",
    ]
    assert (
        payload["categories"]["_ma_protocol_step"]["details"][1]
        == "Post-processing interface QA metrics computed with ipSAE"
    )


def test_generate_normalizes_complex_provenance_and_injects_metrics(
    monkeypatch, tmp_path: Path
) -> None:
    captured_stores: list[Any] = []

    class FakeCifDataStorage:
        def __init__(self) -> None:
            self.data = {
                "_atom_site": {
                    "B_iso_or_equiv": ["10.0", "20.0"],
                    "label_asym_id": ["A", "A"],
                    "label_comp_id": ["ALA", "ALA"],
                    "label_seq_id": ["1", "1"],
                }
            }
            captured_stores.append(self)

        def populate_from_cif_block(self, cif_block: object) -> None:
            return None

        def set_items(self, category_name: str, items_dict: dict[str, list[Any]]) -> None:
            self.data.setdefault(category_name, {}).update(items_dict)

        def set_item(self, category_name: str, item_name: str, item_value: Any) -> None:
            self.data.setdefault(category_name, {})[item_name] = item_value

        def get_data(self) -> dict[str, dict[str, list[Any]]]:
            return self.data

        def write_to_cif(
            self, output_file: str, block_name: str = "model", skip_alignment: bool = False
        ) -> None:
            Path(output_file).write_text("data_test\n", encoding="utf-8")

    input_metadata = {
        "metadata": {"version": "3.1"},
        "categories": {
            "_software": {
                "pdbx_ordinal": ["1", "2", "3", "4", "5"],
                "name": ["ColabFold", "MMseqs2-GPU", "AlphaFold-Multimer", "Interface QA", "PyDSSP"],
                "version": ["example-1.0", "example-2.0", "2.3.2", "example-1.0", "example-1.0"],
                "type": ["package", "package", "package", "package", "library"],
                "description": ["a", "b", "c", "d", "e"],
                "classification": ["model building", "data collection", "model building", "data processing", "data extraction"],
            }
        },
        "chains": [
            {"chain_id": "A", "entity_id": "1", "uniprot_accession": "P1"},
            {"chain_id": "B", "entity_id": "1", "uniprot_accession": "P1"},
        ],
    }
    model_json = tmp_path / "AF-TEST.json"
    model_json.write_bytes(
        orjson.dumps(
            {
                "complexPredictionAccuracy_ipTM": 0.95,
                "complexPredictionAccuracy_ipsae_AB": 0.91,
                "complexPredictionAccuracy_pDockQ": 0.62,
                "complexPredictionAccuracy_ipsae_pae_cutoff": 10.0,
                "complexPredictionAccuracy_ipsae_dist_cutoff": 15.0,
            }
        )
    )

    monkeypatch.setattr(modelcif_generate, "load_json_file", lambda path: input_metadata)
    monkeypatch.setattr(modelcif_generate, "pdb_to_cif_block", lambda path: object())
    monkeypatch.setattr(modelcif_generate, "CifDataStorage", FakeCifDataStorage)
    monkeypatch.setattr(modelcif_generate, "map_entities_and_chains", lambda *args: None)
    monkeypatch.setattr(modelcif_generate, "add_standard_chem_comp_data", lambda *args: None)
    monkeypatch.setattr(modelcif_generate, "compute_global_plddt", lambda values: -1.0)
    monkeypatch.setattr(modelcif_generate, "compute_local_plddt_metrics", lambda *args: {})
    monkeypatch.setattr(modelcif_generate, "create_polymer_sequence_categories", lambda *args: {})
    monkeypatch.setattr(modelcif_generate, "_clamp_struct_ref_seq_to_entity_poly_seq", lambda *args: None)

    modelcif_generate.generate(
        "input.pdb",
        "metadata.json",
        str(tmp_path / "output.cif"),
        validate_dict_path="",
        skip_validation=True,
        model_json_path=str(model_json),
        cif_qa_metrics="auto",
        dssp_algorithm="mkdssp",
    )

    store = captured_stores[0].data
    assert store["_software"]["name"] == ["AlphaFold-Multimer", "ipSAE", "DSSP"]
    assert "Interface QA" not in store["_software"]["name"]
    assert store["_ma_protocol_step"]["step_name"] == [
        "model inference",
        "interface scoring",
        "secondary structure assignment",
    ]
    assert store["_ma_software_parameter"]["name"] == ["pae_cutoff", "dist_cutoff"]
    assert "ipTM" in store["_ma_qa_metric"]["name"]
    assert "ipsae_AB" in store["_ma_qa_metric"]["name"]
