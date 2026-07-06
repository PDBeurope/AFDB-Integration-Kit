from __future__ import annotations

from pathlib import Path
from typing import Any

import orjson
import pytest

import afdb_integration_kit.modelcif.generate as modelcif_generate


def test_generate_does_not_synthesize_alphafold_software_version(
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
    assert software["name"] == ["AlphaFold"]
    assert "version" not in software


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
