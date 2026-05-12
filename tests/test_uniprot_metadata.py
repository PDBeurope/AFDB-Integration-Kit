from __future__ import annotations

from uniprot.scripts.extract_subset import parse_de_sections, parse_organism
from uniprot.scripts.export_model_metadata import (
    ManifestRow,
    ModelMetadataRow,
    build_record,
)


def test_parse_de_sections_preserves_square_brackets_in_names() -> None:
    lines = [
        "DE   RecName: Full=Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase;",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == [
        "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase"
    ]
    assert short_names == []


def test_parse_de_sections_preserves_semicolons_inside_names() -> None:
    lines = [
        "DE   RecName: Full=nitrite reductase (cytochrome; ammonia-forming) {ECO:0000256|ARBA:ARBA00011887};",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == ["nitrite reductase (cytochrome; ammonia-forming)"]
    assert short_names == []


def test_parse_de_sections_uses_submitted_names_as_full_name_fallback() -> None:
    lines = [
        "DE   SubName: Full=Methyl-coenzyme M reductase alpha {ECO:0000313|EMBL:ABW73421.1};",
        "DE   Flags: Fragment;",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == ["Methyl-coenzyme M reductase alpha"]
    assert short_names == []


def test_parse_de_sections_keeps_recommended_names_before_submitted_names() -> None:
    lines = [
        "DE   RecName: Full=Curated recommended name;",
        "DE   SubName: Full=Submitted fallback name {ECO:0000313|EMBL:ABC123.1};",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == ["Curated recommended name", "Submitted fallback name"]
    assert short_names == []


def test_parse_organism_preserves_parenthetical_text_in_scientific_name() -> None:
    organism, common_names, synonyms = parse_organism("Escherichia coli (strain K12).")

    assert organism == "Escherichia coli (strain K12)"
    assert common_names == ["strain K12"]
    assert synonyms == []


def test_build_record_adds_model_level_complex_fields() -> None:
    config = {
        "latestVersion": 1,
        "allVersions": [1],
        "providerId": "AF-TEST",
        "entityType": "protein",
    }
    manifest_rows = [
        ManifestRow(
            model_entity_id="AF-0000000000000001",
            entity_id="1",
            chain_id="A",
            uniprot_ac="P11111",
            sequence_start=None,
            sequence_end=None,
            is_fragment=False,
            is_isoform=False,
            entity_type="protein",
            average_plddt=90.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.0,
            fraction_plddt_confident=0.2,
            fraction_plddt_very_high=0.8,
        ),
        ManifestRow(
            model_entity_id="AF-0000000000000001",
            entity_id="2",
            chain_id="B",
            uniprot_ac="Q22222",
            sequence_start=None,
            sequence_end=None,
            is_fragment=False,
            is_isoform=False,
            entity_type="protein",
            average_plddt=88.0,
            fraction_plddt_very_low=0.0,
            fraction_plddt_low=0.1,
            fraction_plddt_confident=0.2,
            fraction_plddt_very_high=0.7,
        ),
    ]
    entry_map = {
        "P11111": {
            "primary_ac": "P11111",
            "entry_name": "PROT1_TEST",
            "protein_full_names": ["Protein one"],
            "protein_short_names": None,
            "gene_names": "GENE1",
            "gene_synonyms": None,
            "gene_ordered_locus_names": None,
            "gene_orf_names": None,
            "organism": "Organismus exampleus",
            "organism_common_names": None,
            "organism_synonyms": None,
            "taxid": 1111,
            "sequence_version_date": "2024-01-01",
            "sequence": "ACDE",
            "is_uniprot_reference_proteome": True,
            "reviewed": True,
        },
        "Q22222": {
            "primary_ac": "Q22222",
            "entry_name": "PROT2_TEST",
            "protein_full_names": ["Protein two"],
            "protein_short_names": None,
            "gene_names": "GENE2",
            "gene_synonyms": None,
            "gene_ordered_locus_names": None,
            "gene_orf_names": None,
            "organism": "Organismus exampleus",
            "organism_common_names": None,
            "organism_synonyms": None,
            "taxid": 1111,
            "sequence_version_date": "2024-01-01",
            "sequence": "FGHI",
            "is_uniprot_reference_proteome": True,
            "reviewed": True,
        },
    }
    model_metadata = {
        "AF-0000000000000001": ModelMetadataRow(
            iptm=0.91,
            average_plddt=89.3,
            complex_name=None,
            is_am_data=False,
        )
    }

    record = build_record(
        "AF-0000000000000001",
        config,
        manifest_rows,
        entry_map,
        model_metadata,
    )

    assert record["assemblyType"] == "Hetero"
    assert record["oligomericState"] == "dimer"
    assert record["oligomericStateDescription"] == "Heterodimer"
    assert record["complexComposition"] == ["P11111_1", "Q22222_1"]
