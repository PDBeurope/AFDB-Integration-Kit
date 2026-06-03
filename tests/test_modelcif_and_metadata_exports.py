from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_module(relative_path: str, module_name: str):
    repo_root = Path(__file__).resolve().parent.parent
    module_path = repo_root / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_export_modelcif_input_uses_fragment_ranges() -> None:
    module = _load_module("uniprot/scripts/export_modelcif_input.py", "export_modelcif_input")

    template = {"categories": {}}
    entities = [
        module.EntityAssignment(
            entity_id="1",
            uniprot_ac="P11111",
            chain_ids=["A"],
            sequence_start=5,
            sequence_end=7,
        )
    ]
    entries_by_entity = {
        "1": {
            "entry_name": "TEST_ENTRY",
            "gene_names": "GENE1",
            "taxid": 1234,
            "organism": "Example organism",
            "sequence": "ABCDEFGH",
            "sequence_version_date": "2024-01-01",
            "protein_full_names": ["Example protein"],
        }
    }

    module.populate_categories(template, entities, entries_by_entity, "AF-TEST")
    categories = template["categories"]

    assert categories["_entity_poly"]["pdbx_seq_one_letter_code"] == ["EFG"]
    assert categories["_ma_target_ref_db_details"]["seq_db_align_begin"] == [5]
    assert categories["_ma_target_ref_db_details"]["seq_db_align_end"] == [7]
    assert categories["_struct_ref_seq"]["seq_align_beg"] == [1]
    assert categories["_struct_ref_seq"]["seq_align_end"] == [3]
    assert categories["_struct_ref_seq"]["db_align_beg"] == [5]
    assert categories["_struct_ref_seq"]["db_align_end"] == [7]


def test_batch_export_metadata_generates_hetero_complex_name() -> None:
    module = _load_module("uniprot/scripts/batch_export_metadata.py", "batch_export_metadata")

    config = {
        "providerId": "EXAMPLE",
        "latestVersion": 1,
        "allVersions": [1],
        "toolUsed": "ColabFold",
        "modelCreatedDate": "2026-01-01T00:00:00Z",
        "versionTag": "v1",
        "entityType": "protein",
    }
    manifest_rows = [
        module.ManifestRow(
            model_entity_id="AF-TEST",
            entity_id="1",
            chain_id="A",
            uniprot_ac="P11111",
            sequence_start=1,
            sequence_end=3,
            is_fragment=True,
            is_isoform=False,
            entity_type="protein",
            average_plddt=70.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.0,
            fraction_plddt_confident=1.0,
            fraction_plddt_very_high=0.0,
        ),
        module.ManifestRow(
            model_entity_id="AF-TEST",
            entity_id="2",
            chain_id="B",
            uniprot_ac="Q22222",
            sequence_start=1,
            sequence_end=3,
            is_fragment=True,
            is_isoform=False,
            entity_type="protein",
            average_plddt=80.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.0,
            fraction_plddt_confident=1.0,
            fraction_plddt_very_high=0.0,
        ),
    ]
    entry_map = {
        "P11111": {
            "protein_full_names": ["Protein one"],
            "is_uniprot_reference_proteome": True,
            "reviewed": True,
            "organism_scientific_name": "Example organism",
            "organism_common_names": [],
            "organism_synonyms": [],
            "gene_primary": "GENE1",
            "gene_synonyms": [],
            "sequence": "ABC",
            "sequence_checksum": "checksum1",
            "sequence_version_date": "2024-01-01",
            "entry_name": "PROT1_EXAMPLE",
            "organism_id": 1234,
        },
        "Q22222": {
            "protein_full_names": ["Protein two"],
            "is_uniprot_reference_proteome": True,
            "reviewed": False,
            "organism_scientific_name": "Example organism",
            "organism_common_names": [],
            "organism_synonyms": [],
            "gene_primary": "GENE2",
            "gene_synonyms": [],
            "sequence": "DEF",
            "sequence_checksum": "checksum2",
            "sequence_version_date": "2024-01-01",
            "entry_name": "PROT2_EXAMPLE",
            "organism_id": 1234,
        },
    }
    model_metadata = {
        "AF-TEST": module.ModelMetadataRow(
            iptm=0.42,
            average_plddt=75.0,
            complex_name=None,
            is_am_data=False,
        )
    }

    model_record = module.build_model_record("AF-TEST", config, manifest_rows, entry_map, model_metadata)
    chain_records = module.build_chain_records("AF-TEST", config, manifest_rows, entry_map, model_metadata)

    assert model_record["complexName"] == "Complex of Protein one/Protein two"
    assert chain_records[0]["complexName"] == "Complex of Protein one/Protein two"
    assert chain_records[1]["complexName"] == "Complex of Protein one/Protein two"
