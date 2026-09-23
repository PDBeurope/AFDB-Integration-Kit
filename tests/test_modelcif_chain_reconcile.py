"""Tests for structure-based chain<->entity reconciliation in the ModelCIF exporter.

Covers the heterodimer chain-swap fix: chain->entity (and thus chain->accession
provenance) must follow the coordinates, not the manifest order.
"""
import pytest

from afdb_integration_kit.modelcif.generate import (
    _assign_chains_to_entities,
    _chain_seqs_from_atom_site,
    _entity_ref_seqs_from_input,
    assert_chain_entity_consistency,
    reconcile_entities_with_structure,
)

CAT_STRUCT_REF = "_struct_ref"
CAT_STRUCT_REF_SEQ = "_struct_ref_seq"
CAT_TARGET_ENTITY_INSTANCE = "_ma_target_entity_instance"
from afdb_integration_kit.utils.cifstorage import CifDataStorage
from afdb_integration_kit.utils.constant import (
    CAT_ATOM_SITE,
    CAT_STRUCT_ASYM,
    ITEM_LABEL_ASYM_ID,
    ITEM_LABEL_COMP_ID,
    ITEM_LABEL_SEQ_ID,
)


def test_assign_chains_distinct_length_swapped():
    # chain A is the short protein, chain B the long one; entity refs are reversed.
    chain_seqs = {"A": "A" * 94, "B": "G" * 117}
    entity_seqs = {"1": "G" * 117, "2": "A" * 94}  # manifest order: entity1=long
    assert _assign_chains_to_entities(chain_seqs, entity_seqs) == {"A": "2", "B": "1"}


def test_assign_chains_equal_length_by_sequence():
    chain_seqs = {"A": "A" * 100, "B": "G" * 100}
    entity_seqs = {"1": "G" * 100, "2": "A" * 100}
    assert _assign_chains_to_entities(chain_seqs, entity_seqs) == {"A": "2", "B": "1"}


def test_assign_chains_already_correct():
    chain_seqs = {"A": "A" * 94, "B": "G" * 117}
    entity_seqs = {"1": "A" * 94, "2": "G" * 117}
    assert _assign_chains_to_entities(chain_seqs, entity_seqs) == {"A": "1", "B": "2"}


def test_entity_ref_seqs_from_entity_poly_seq_fallback():
    md = {
        "categories": {
            "_entity_poly_seq": {
                "entity_id": ["1", "1", "1", "2", "2"],
                "num": ["1", "2", "3", "1", "2"],
                "mon_id": ["ALA", "GLY", "SER", "TRP", "LYS"],
            }
        }
    }
    assert _entity_ref_seqs_from_input(md) == {"1": "AGS", "2": "WK"}


def test_entity_ref_seqs_from_entity_poly_primary():
    # This is what export_modelcif_input.py actually writes (per-entity one-letter seq).
    md = {
        "categories": {
            "_entity_poly": {
                "entity_id": ["1", "2"],
                "pdbx_seq_one_letter_code": ["MASSG", "MEAT\nDIW"],  # may contain newlines
            }
        }
    }
    assert _entity_ref_seqs_from_input(md) == {"1": "MASSG", "2": "MEATDIW"}


def _atom_rows(chain, comp3, n):
    return (
        [chain] * n,
        [comp3] * n,
        [str(i) for i in range(1, n + 1)],
    )


def test_reconcile_fixes_swapped_chain_entity():
    # Structure: chain A = 3 ALA residues, chain B = 5 GLY residues.
    a_asym, a_comp, a_seq = _atom_rows("A", "ALA", 3)
    b_asym, b_comp, b_seq = _atom_rows("B", "GLY", 5)
    cif = CifDataStorage()
    cif.set_items(
        CAT_ATOM_SITE,
        {
            ITEM_LABEL_ASYM_ID: a_asym + b_asym,
            ITEM_LABEL_COMP_ID: a_comp + b_comp,
            ITEM_LABEL_SEQ_ID: a_seq + b_seq,
        },
    )
    # Manifest-order (WRONG) mapping: chain A -> entity 1, chain B -> entity 2,
    # but entity 1 is really the 5-residue protein and entity 2 the 3-residue one.
    cif.set_items(CAT_STRUCT_ASYM, {"id": ["A", "B"], "entity_id": ["1", "2"]})
    # Entity-keyed _struct_ref (correct regardless of chain order).
    cif.set_items(CAT_STRUCT_REF, {
        "id": ["1", "2"], "entity_id": ["1", "2"], "db_name": ["UNP", "UNP"],
        "pdbx_db_accession": ["P-LONG", "P-SHORT"],
    })
    # Chain-attributed categories in WRONG (manifest) order.
    cif.set_items(CAT_TARGET_ENTITY_INSTANCE, {
        "asym_id": ["A", "B"], "entity_id": ["1", "2"],
        "details": ["Chain A from UniProt P-LONG", "Chain B from UniProt P-SHORT"],
    })
    cif.set_items(CAT_STRUCT_REF_SEQ, {
        "align_id": ["1", "2"], "ref_id": ["1", "2"], "pdbx_PDB_id_code": ["m", "m"],
        "pdbx_strand_id": ["A", "B"],
        "seq_align_beg": ["1", "1"], "seq_align_end": ["5", "3"],
        "db_align_beg": ["1", "1"], "db_align_end": ["5", "3"],
    })
    entity_ref_seqs = {"1": "G" * 5, "2": "A" * 3}

    # Sanity: coordinate sequences extracted correctly.
    assert _chain_seqs_from_atom_site(cif.get_data()[CAT_ATOM_SITE]) == {"A": "AAA", "B": "GGGGG"}

    reconcile_entities_with_structure(cif, entity_ref_seqs)

    asym = cif.get_data()[CAT_STRUCT_ASYM]
    mapping = dict(zip(asym["id"], asym["entity_id"]))
    # chain A (3 res) must now map to entity 2 (the 3-res protein), chain B to entity 1.
    assert mapping == {"A": "2", "B": "1"}
    assert cif.get_data()[CAT_ATOM_SITE]["label_entity_id"] == ["2", "2", "2", "1", "1", "1", "1", "1"]

    # _ma_target_entity_instance re-attributed (asym fixed; entity + accession follow structure).
    tei = cif.get_data()[CAT_TARGET_ENTITY_INSTANCE]
    assert tei["asym_id"] == ["A", "B"]
    assert tei["entity_id"] == ["2", "1"]
    assert tei["details"] == ["Chain A from UniProt P-SHORT", "Chain B from UniProt P-LONG"]

    # _struct_ref_seq: each strand now points at the right ref_id, and the entity's
    # alignment ranges travelled with it (chain A -> 3 res, chain B -> 5 res).
    srs = cif.get_data()[CAT_STRUCT_REF_SEQ]
    assert srs["pdbx_strand_id"] == ["A", "B"]
    assert srs["ref_id"] == ["2", "1"]
    assert srs["seq_align_end"] == ["3", "5"]
    assert srs["db_align_end"] == ["3", "5"]

    # Whole cif is internally consistent after the fix.
    assert_chain_entity_consistency(cif)


def test_reconcile_remaps_categories_absent_gracefully():
    # Swapped struct_asym/atom_site but NO struct_ref family seeded -> still fixes the
    # core mapping and does not raise.
    a_asym, a_comp, a_seq = _atom_rows("A", "ALA", 3)
    b_asym, b_comp, b_seq = _atom_rows("B", "GLY", 5)
    cif = CifDataStorage()
    cif.set_items(CAT_ATOM_SITE, {
        ITEM_LABEL_ASYM_ID: a_asym + b_asym,
        ITEM_LABEL_COMP_ID: a_comp + b_comp,
        ITEM_LABEL_SEQ_ID: a_seq + b_seq,
    })
    cif.set_items(CAT_STRUCT_ASYM, {"id": ["A", "B"], "entity_id": ["1", "2"]})
    reconcile_entities_with_structure(cif, {"1": "G" * 5, "2": "A" * 3})
    asym = cif.get_data()[CAT_STRUCT_ASYM]
    assert dict(zip(asym["id"], asym["entity_id"])) == {"A": "2", "B": "1"}
    assert_chain_entity_consistency(cif)  # passes (absent categories skipped)


def test_consistency_guard_raises_on_inconsistent_cif():
    cif = CifDataStorage()
    a_asym, a_comp, a_seq = _atom_rows("A", "ALA", 3)
    b_asym, b_comp, b_seq = _atom_rows("B", "GLY", 5)
    cif.set_items(CAT_ATOM_SITE, {
        ITEM_LABEL_ASYM_ID: a_asym + b_asym,
        ITEM_LABEL_COMP_ID: a_comp + b_comp,
        ITEM_LABEL_SEQ_ID: a_seq + b_seq,
        "label_entity_id": ["2", "2", "2", "1", "1", "1", "1", "1"],
    })
    cif.set_items(CAT_STRUCT_ASYM, {"id": ["A", "B"], "entity_id": ["2", "1"]})
    # _ma_target_entity_instance left in the OLD (disagreeing) order.
    cif.set_items(CAT_TARGET_ENTITY_INSTANCE, {"asym_id": ["A", "B"], "entity_id": ["1", "2"]})
    with pytest.raises(ValueError, match="_ma_target_entity_instance"):
        assert_chain_entity_consistency(cif)


def test_consistency_guard_passes_when_consistent():
    cif = CifDataStorage()
    cif.set_items(CAT_ATOM_SITE, {
        ITEM_LABEL_ASYM_ID: ["A", "B"],
        "label_entity_id": ["2", "1"],
    })
    cif.set_items(CAT_STRUCT_ASYM, {"id": ["A", "B"], "entity_id": ["2", "1"]})
    cif.set_items(CAT_TARGET_ENTITY_INSTANCE, {"asym_id": ["A", "B"], "entity_id": ["2", "1"]})
    cif.set_items(CAT_STRUCT_REF, {"id": ["1", "2"], "entity_id": ["1", "2"]})
    cif.set_items(CAT_STRUCT_REF_SEQ, {"pdbx_strand_id": ["A", "B"], "ref_id": ["2", "1"]})
    assert_chain_entity_consistency(cif)  # no raise


def test_reconcile_noop_for_homomultimer():
    # Two chains, ONE entity (homodimer): reconcile must not touch it.
    cif = CifDataStorage()
    a_asym, a_comp, a_seq = _atom_rows("A", "ALA", 4)
    b_asym, b_comp, b_seq = _atom_rows("B", "ALA", 4)
    cif.set_items(
        CAT_ATOM_SITE,
        {
            ITEM_LABEL_ASYM_ID: a_asym + b_asym,
            ITEM_LABEL_COMP_ID: a_comp + b_comp,
            ITEM_LABEL_SEQ_ID: a_seq + b_seq,
        },
    )
    cif.set_items(CAT_STRUCT_ASYM, {"id": ["A", "B"], "entity_id": ["1", "1"]})
    reconcile_entities_with_structure(cif, {"1": "A" * 4})  # one entity != two chains
    asym = cif.get_data()[CAT_STRUCT_ASYM]
    assert dict(zip(asym["id"], asym["entity_id"])) == {"A": "1", "B": "1"}
