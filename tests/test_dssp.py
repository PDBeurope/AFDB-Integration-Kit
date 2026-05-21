from pathlib import Path

import gemmi
import numpy as np
import pytest

from afdb_integration_kit.dssp.dssp import (
    DEFAULT_ALGORITHM,
    _compute_sse_tmalign,
    run_dssp,
)


EXAMPLE_CIF = Path("examples/AF-0000000000000001-model_v1.cif")


def _category_rows(block, category):
    table = block.find_mmcif_category(category)
    return [
        [row[i] for i in range(table.width())]
        for row in table
    ]


def test_default_algorithm_matches_production_docs():
    assert DEFAULT_ALGORITHM == "pydssp"


def test_tmalign_classifies_secondary_structure_from_ca_distances():
    doc = gemmi.cif.read(str(EXAMPLE_CIF))
    structure = gemmi.make_structure_from_block(doc.sole_block())
    ca_coords = []
    for chain in structure[0]:
        for residue in chain.get_polymer():
            ca = residue.find_atom("CA", "\0")
            ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])

    sse = _compute_sse_tmalign(np.array(ca_coords, dtype=np.float32))

    assert 2 in sse
    assert 4 in sse


def test_tmalign_writes_struct_conf_and_preserves_existing_categories(tmp_path):
    output = tmp_path / "with_dssp.cif"

    assert run_dssp(EXAMPLE_CIF, output, algorithm="tmalign")

    block = gemmi.cif.read(str(output)).sole_block()
    conf_type_rows = _category_rows(block, "_struct_conf_type.")
    conf_rows = _category_rows(block, "_struct_conf.")

    assert ["HELX_P", "TM-align"] in conf_type_rows
    assert ["STRN", "TM-align"] in conf_type_rows
    assert conf_rows[0] == [
        "HELX_P1", "HELX_P", "ILE", "A", "3", "?",
        "THR", "A", "25", "?", "ILE", "A", "3", "THR", "A", "25",
    ]
    assert block.find_value("_entry.id") == "AF-0000000000000001-model_v1"
    assert block.find_mmcif_category("_ma_qa_metric_local.").width() > 0
    entity_rows = _category_rows(block, "_entity.")
    assert entity_rows[0][1] == "polymer"


def test_no_polymer_residues_fail_without_writing_output(tmp_path):
    input_cif = tmp_path / "nonpolymer.cif"
    output_cif = tmp_path / "out.cif"
    input_cif.write_text(
        """data_nonpolymer
_entry.id nonpolymer
#
_entity.id 1
_entity.type non-polymer
#
_struct_asym.id A
_struct_asym.entity_id 1
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
_atom_site.auth_comp_id
HETATM 1 O O . HOH A 1 . ? 0.0 0.0 0.0 1.0 10.0 ? 1 A 1 HOH
""",
        encoding="utf-8",
    )

    assert not run_dssp(input_cif, output_cif, algorithm="tmalign")
    assert not output_cif.exists()


def test_psea_writes_struct_conf_when_biotite_is_installed(tmp_path):
    pytest.importorskip("biotite")
    output = tmp_path / "psea.cif"

    assert run_dssp(EXAMPLE_CIF, output, algorithm="psea")

    block = gemmi.cif.read(str(output)).sole_block()
    assert _category_rows(block, "_struct_conf_type.")


def test_pydssp_writes_struct_conf_when_optional_package_is_installed(tmp_path):
    pytest.importorskip("pydssp")
    output = tmp_path / "pydssp.cif"

    assert run_dssp(EXAMPLE_CIF, output, algorithm="pydssp")

    block = gemmi.cif.read(str(output)).sole_block()
    assert ["HELX_P", "PyDSSP"] in _category_rows(block, "_struct_conf_type.")
