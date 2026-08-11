from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import orjson
import pytest

from afdb_integration_kit.uniprot.naming import protein_description


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


def test_modelcif_template_and_examples_keep_explicit_alphafold_versions() -> None:
    repo_root = Path(__file__).resolve().parent.parent

    template = orjson.loads(
        (
            repo_root / "uniprot/templates/colabfold_example_modelcif_metadata.json"
        ).read_bytes()
    )
    monomer_input = orjson.loads(
        (
            repo_root
            / "examples/colabfold_monomer_e2e/modelcif_input/AF-0000000300000001.json"
        ).read_bytes()
    )
    complex_input = orjson.loads(
        (
            repo_root
            / "examples/colabfold_complex_e2e/modelcif_input/AF-0000000066074510.json"
        ).read_bytes()
    )

    assert template["categories"]["_software"]["version"][0] == "2.3.2"
    assert monomer_input["categories"]["_software"]["name"] == ["AlphaFold", "DSSP"]
    assert monomer_input["categories"]["_software"]["version"] == ["2.3.2", "?"]
    assert complex_input["categories"]["_software"]["name"] == [
        "AlphaFold-Multimer",
        "ipSAE",
        "DSSP",
    ]
    assert complex_input["categories"]["_software"]["version"] == ["2.3.2", "?", "?"]


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


@pytest.mark.parametrize(
    ("relative_path", "module_name"),
    [
        ("uniprot/scripts/export_modelcif_input.py", "single_modelcif_name"),
        (
            "uniprot/scripts/batch_export_modelcif_input.py",
            "batch_modelcif_name",
        ),
    ],
)
def test_modelcif_exporters_prefer_manifest_protein_name(
    relative_path: str,
    module_name: str,
) -> None:
    module = _load_module(relative_path, module_name)
    template = {"categories": {}}
    entities = [
        module.EntityAssignment(
            entity_id="1",
            uniprot_ac="P11111",
            chain_ids=["A"],
            protein_name="Named fragment",
        )
    ]
    entry = {
        "entry_name": "TEST_ENTRY",
        "protein_full_names": ["Whole protein"],
        "sequence": "ABC",
    }
    entries = {"1": entry, "P11111": entry}

    module.populate_categories(template, entities, entries, "AF-TEST")

    assert template["categories"]["_entity"]["pdbx_description"] == [
        "Named fragment"
    ]


def test_protein_description_fallback_precedence() -> None:
    entry = {
        "protein_full_names": ["Full name"],
        "protein_short_names": ["Short name"],
        "entry_name": "ENTRY_NAME",
    }
    assert (
        protein_description("Manifest name", entry, "P11111")
        == "Manifest name"
    )
    assert protein_description(None, entry, "P11111") == "Full name"
    entry["protein_full_names"] = []
    assert protein_description(None, entry, "P11111") == "Short name"
    entry["protein_short_names"] = []
    assert protein_description(None, entry, "P11111") == "ENTRY_NAME"
    entry["entry_name"] = ""
    assert protein_description(None, entry, "P11111") == "P11111"


@pytest.mark.parametrize(
    ("relative_path", "module_name", "loader_kind"),
    [
        (
            "uniprot/scripts/export_model_metadata.py",
            "parse_single_model",
            "metadata",
        ),
        (
            "uniprot/scripts/export_chain_metadata.py",
            "parse_single_chain",
            "metadata",
        ),
        (
            "uniprot/scripts/batch_export_metadata.py",
            "parse_batch_metadata",
            "metadata",
        ),
        (
            "uniprot/scripts/export_modelcif_input.py",
            "parse_single_modelcif",
            "modelcif",
        ),
        (
            "uniprot/scripts/batch_export_modelcif_input.py",
            "parse_batch_modelcif",
            "batch_modelcif",
        ),
    ],
)
def test_manifest_consumers_parse_optional_protein_name(
    tmp_path: Path,
    relative_path: str,
    module_name: str,
    loader_kind: str,
) -> None:
    module = _load_module(relative_path, module_name)
    manifest = tmp_path / f"{module_name}.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,protein_name\n"
        "AF-TEST,1,A,P11111,Named fragment\n",
        encoding="utf-8",
    )

    if loader_kind == "metadata":
        row = module.load_manifest(manifest).by_model["AF-TEST"][0]
    elif loader_kind == "batch_modelcif":
        loaded = module.load_manifest(manifest, ["AF-TEST"])
        row = loaded.by_model["AF-TEST"][0]
    else:
        row = module.load_manifest(manifest, "AF-TEST")[0]

    assert row.protein_name == "Named fragment"


@pytest.mark.parametrize(
    ("relative_path", "module_name", "builder", "expected"),
    [
        (
            "uniprot/scripts/export_model_metadata.py",
            "single_model_metadata_name",
            "build_record",
            "model",
        ),
        (
            "uniprot/scripts/export_chain_metadata.py",
            "single_chain_metadata_name",
            "build_record",
            "chain",
        ),
        (
            "uniprot/scripts/batch_export_metadata.py",
            "batch_model_metadata_name",
            "build_model_record",
            "model",
        ),
        (
            "uniprot/scripts/batch_export_metadata.py",
            "batch_chain_metadata_name",
            "build_chain_records",
            "chain",
        ),
    ],
)
def test_metadata_exporters_prefer_manifest_protein_name(
    relative_path: str,
    module_name: str,
    builder: str,
    expected: str,
) -> None:
    module = _load_module(relative_path, module_name)
    row = module.ManifestRow(
        model_entity_id="AF-TEST",
        entity_id="1",
        chain_id="A",
        uniprot_ac="P11111",
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
        protein_name="Named fragment",
    )
    config = {
        "providerId": "EXAMPLE",
        "latestVersion": 1,
        "allVersions": [1],
        "toolUsed": "ColabFold",
        "modelCreatedDate": "2026-01-01T00:00:00Z",
        "versionTag": "v1",
        "entityType": "protein",
    }
    entry_map = {
        "P11111": {
            "protein_full_names": ["Whole protein"],
            "entry_name": "TEST_ENTRY",
            "sequence": "ABC",
        }
    }

    result = getattr(module, builder)("AF-TEST", config, [row], entry_map, {})

    if expected == "model":
        assert result["uniprotDescription"] == ["Named fragment"]
    else:
        assert result[0]["uniprotDescription"] == "Named fragment"


@pytest.mark.parametrize(
    (
        "relative_path",
        "module_name",
        "builder",
        "record_kind",
        "expected_hetero_name",
    ),
    [
        (
            "uniprot/scripts/export_model_metadata.py",
            "fragment_identity_model_metadata",
            "build_record",
            "model",
            None,
        ),
        (
            "uniprot/scripts/export_chain_metadata.py",
            "fragment_identity_chain_metadata",
            "build_record",
            "chain",
            None,
        ),
        (
            "uniprot/scripts/batch_export_metadata.py",
            "fragment_identity_batch_model_metadata",
            "build_model_record",
            "model",
            "Complex of Example protein/Example protein",
        ),
        (
            "uniprot/scripts/batch_export_metadata.py",
            "fragment_identity_batch_chain_metadata",
            "build_chain_records",
            "chain",
            "Complex of Example protein/Example protein",
        ),
    ],
)
@pytest.mark.parametrize(
    ("fragment_ranges", "expected_assembly", "expected_description", "expected_name"),
    [
        (
            ((1, 46), (961, 1071)),
            "Hetero",
            "Heterodimer",
            None,
        ),
        (
            ((1, 46), (1, 46)),
            "Homo",
            "Homodimer",
            "Homodimer of Example protein",
        ),
    ],
    ids=("distinct-fragments", "repeated-fragment"),
)
def test_metadata_exporters_classify_assembly_by_fragment_identity(
    relative_path: str,
    module_name: str,
    builder: str,
    record_kind: str,
    expected_hetero_name: str | None,
    fragment_ranges: tuple[tuple[int, int], tuple[int, int]],
    expected_assembly: str,
    expected_description: str,
    expected_name: str | None,
) -> None:
    module = _load_module(relative_path, module_name)
    manifest_rows = [
        module.ManifestRow(
            model_entity_id="AF-FRAGMENT-IDENTITY",
            entity_id=str(index),
            chain_id=chain_id,
            uniprot_ac="P27409",
            sequence_start=sequence_start,
            sequence_end=sequence_end,
            is_fragment=True,
            is_isoform=False,
            entity_type="protein",
            average_plddt=80.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.0,
            fraction_plddt_confident=1.0,
            fraction_plddt_very_high=0.0,
        )
        for index, (chain_id, (sequence_start, sequence_end)) in enumerate(
            zip(("A", "B"), fragment_ranges),
            start=1,
        )
    ]
    config = {
        "providerId": "EXAMPLE",
        "latestVersion": 1,
        "allVersions": [1],
        "toolUsed": "ColabFold",
        "modelCreatedDate": "2026-01-01T00:00:00Z",
        "versionTag": "v1",
        "entityType": "protein",
    }
    entry_map = {
        "P27409": {
            "protein_full_names": ["Example protein"],
            "entry_name": "EXAMPLE_ENTRY",
            "sequence": "A" * 1071,
        }
    }

    result = getattr(module, builder)(
        "AF-FRAGMENT-IDENTITY",
        config,
        manifest_rows,
        entry_map,
        {},
    )
    records = [result] if record_kind == "model" else result

    assert {record["assemblyType"] for record in records} == {expected_assembly}
    assert {record["oligomericStateDescription"] for record in records} == {
        expected_description
    }
    if expected_assembly == "Hetero":
        expected_name = expected_hetero_name
    if expected_name is None:
        assert all("complexName" not in record for record in records)
    else:
        assert {record["complexName"] for record in records} == {expected_name}


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

    assert model_record["isUniProt"] is True
    assert model_record["complexName"] == "Complex of Protein one/Protein two"
    assert chain_records[0]["complexName"] == "Complex of Protein one/Protein two"
    assert chain_records[1]["complexName"] == "Complex of Protein one/Protein two"


def test_batch_export_metadata_does_not_generate_hetero_complex_name_for_three_chains() -> None:
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
            model_entity_id="AF-TEST-3",
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
            model_entity_id="AF-TEST-3",
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
        module.ManifestRow(
            model_entity_id="AF-TEST-3",
            entity_id="3",
            chain_id="C",
            uniprot_ac="R33333",
            sequence_start=1,
            sequence_end=3,
            is_fragment=True,
            is_isoform=False,
            entity_type="protein",
            average_plddt=75.0,
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
        "R33333": {
            "protein_full_names": ["Protein three"],
            "is_uniprot_reference_proteome": True,
            "reviewed": False,
            "organism_scientific_name": "Example organism",
            "organism_common_names": [],
            "organism_synonyms": [],
            "gene_primary": "GENE3",
            "gene_synonyms": [],
            "sequence": "GHI",
            "sequence_checksum": "checksum3",
            "sequence_version_date": "2024-01-01",
            "entry_name": "PROT3_EXAMPLE",
            "organism_id": 1234,
        },
    }
    model_metadata = {
        "AF-TEST-3": module.ModelMetadataRow(
            iptm=0.42,
            average_plddt=75.0,
            complex_name=None,
            is_am_data=False,
        )
    }

    model_record = module.build_model_record("AF-TEST-3", config, manifest_rows, entry_map, model_metadata)
    chain_records = module.build_chain_records("AF-TEST-3", config, manifest_rows, entry_map, model_metadata)

    assert "complexName" not in model_record
    assert all("complexName" not in record for record in chain_records)
