import datetime
import pytest

from afdb_integration_kit.utils.pdbeditor import PDBFileEditor


@pytest.fixture
def bare_pdb_content():
    """Fixture to provide sample ATOM data for testing."""
    return [
        "ATOM      1  N   ALA A   1      11.459   0.117   0.145  1.00 70.81           N  \n",
        "ATOM      2  CA  ALA A   1      12.227  -0.966  -0.536  1.00 72.93           C  \n",
        "ATOM      3  C   ALA A   1      13.593  -0.457  -0.970  1.00 74.34           C  \n",
        "ATOM      4  O   ALA A   1      14.390  -1.229  -1.503  1.00 74.22           O  \n",
        "ATOM      5  CB  ALA A   1      11.594  -1.921  -1.520  1.00 74.96           C  \n",
        "ATOM      6  N   LYS A   2      13.882   0.803  -0.784  1.00 75.31           N  \n",
        "ATOM      7  CA  LYS A   2      15.197   1.416  -1.127  1.00 77.06           C  \n",
        "ATOM      8  C   LYS A   2      15.701   0.887  -2.450  1.00 78.89           C  \n",
        "ATOM      9  O   LYS A   2      16.719   1.382  -2.853  1.00 79.79           O  \n",
        "ATOM     10  CB  LYS A   2      15.011   2.946  -0.925  1.00 78.50           C  \n",
    ]


@pytest.fixture
def create_temp_pdb_file(tmp_path, bare_pdb_content):
    """Fixture to create a temporary PDB file for testing file I/O."""
    temp_pdb = tmp_path / "temp.pdb"
    temp_pdb.write_text("".join(bare_pdb_content))
    return temp_pdb


def test_initialization_with_content(bare_pdb_content):
    """Test if PDBFileEditor initializes correctly with content."""
    editor = PDBFileEditor(pdb_content=bare_pdb_content)
    assert editor.lines == bare_pdb_content


def test_initialization_without_content():
    """Test if PDBFileEditor initializes with an empty list."""
    editor = PDBFileEditor()
    assert editor.lines == []


def test_load_pdb(create_temp_pdb_file, bare_pdb_content):
    """Test the load_pdb method."""
    editor = PDBFileEditor()
    editor.load_pdb(create_temp_pdb_file)
    assert editor.lines == bare_pdb_content


def test_write_pdb(bare_pdb_content, tmp_path):
    """Test the write_pdb method with a full set of headers."""
    output_file = tmp_path / "test_output.pdb"

    # Set up the editor with bare ATOM lines
    editor = PDBFileEditor(pdb_content=bare_pdb_content)
    
    # Add a comprehensive set of header records in a non-canonical order
    editor.add_model(model_number=1)
    editor.add_title("A TEST TITLE")
    editor.add_remark_blank_line(remark_number=1)
    editor.add_dbref(
        pdb_id="TEST", chain_id="A", seq_begin=1, seq_end=10,
        db="UNP", db_accession="A1B2C3", db_id="ID_TEST",
        db_seq_begin=1, db_seq_end=10
    )
    editor.add_header(pdb_id="TEST", date=datetime.date(2023, 1, 1))
    
    # Add a multi-line remark
    multi_line_remark_text = (
        "This is a very long remark that needs to be wrapped across multiple lines. "
        "It should be handled correctly by the editor's text wrapping logic."
    )
    editor.add_remark(remark_number=1, text=multi_line_remark_text, numbered=True)
    
    # Add crystal info
    editor.add_crystal_info(10.0, 20.0, 30.0, 90.0, 90.0, 90.0)

    editor.write_pdb(output_file)

    # Read the output file and check content and order
    with open(output_file, 'r') as f:
        lines = f.readlines()

    # Check for correct order
    header_index = next((i for i, line in enumerate(lines) if line.startswith("HEADER")), -1)
    title_index = next((i for i, line in enumerate(lines) if line.startswith("TITLE")), -1)
    dbref_index = next((i for i, line in enumerate(lines) if line.startswith("DBREF")), -1)
    remark_index = next((i for i, line in enumerate(lines) if line.startswith("REMARK")), -1)
    crystal_index = next((i for i, line in enumerate(lines) if line.startswith("CRYST1")), -1)
    model_index = next((i for i, line in enumerate(lines) if line.startswith("MODEL")), -1)
    atom_index = next((i for i, line in enumerate(lines) if line.startswith("ATOM")), -1)
    
    # Assert that records appear in the expected order
    assert header_index < title_index
    assert title_index < remark_index
    assert remark_index < dbref_index
    assert dbref_index < crystal_index
    assert crystal_index < model_index
    assert model_index < atom_index
    
    # Check that the number of ATOM lines is correct
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    assert len(atom_lines) == len(bare_pdb_content)
    
    # Check that the title and remarks are present and wrapped
    assert any("A TEST TITLE" in line for line in lines)
    assert any("This is a very long remark" in line for line in lines)


def test_add_title():
    """Test title wrapping and formatting."""
    editor = PDBFileEditor()
    title = "This is a very long title that needs to be wrapped across multiple lines to fit the PDB file specification."
    editor.add_title(title)
    
    # The title should be split into two lines
    assert len(editor._header_lines_to_insert) == 2
    assert "TITLE    This is a very long title that needs to be wrapped across multiple" in editor._header_lines_to_insert[0]
    assert "TITLE    2 lines to fit the PDB file specification." in editor._header_lines_to_insert[1]
    

def test_add_compnd():
    """Test COMPND record formatting and multi-line handling."""
    editor = PDBFileEditor()
    molecule_name = "A VERY LONG MOLECULE NAME FOR TESTING THE WRAPPING LOGIC OF THE COMPND RECORD IN THE PDB FILE."
    editor.add_compnd(molecule_id=1, molecule=molecule_name, chain="A")
    
    # Check the number of lines generated
    expected_lines = 1 + 2 + 1 # MOL_ID + wrapped molecule (2 lines) + CHAIN
    assert len(editor._header_lines_to_insert) == expected_lines

    # Check the content of the lines
    lines = [line.strip() for line in editor._header_lines_to_insert]
    assert "COMPND    MOL_ID: 1" in lines[0]
    assert "COMPND   2 MOLECULE: A VERY LONG MOLECULE NAME FOR TESTING THE WRAPPI;" in lines[1]
    assert "COMPND   3 NG LOGIC OF THE COMPND RECORD IN THE PDB FILE.;" in lines[2]
    assert "COMPND   4 CHAIN: A" in lines[3]


def test_add_source():
    """Test SOURCE record formatting."""
    editor = PDBFileEditor()
    editor.add_source(
        molecule_id=1,
        organism_scientific="TESTUS ORGANISMUS",
        organism_taxid=9999
    )
    
    assert len(editor._header_lines_to_insert) == 3
    lines = [line.strip() for line in editor._header_lines_to_insert]
    assert "SOURCE    MOL_ID: 1;" in lines[0]
    assert "SOURCE   2 ORGANISM_SCIENTIFIC: TESTUS ORGANISMUS;" in lines[1]
    assert "SOURCE   3 ORGANISM_TAXID: 9999" in lines[2]


def test_add_header():
    """Test HEADER record formatting."""
    editor = PDBFileEditor()
    test_date = datetime.date(2025, 9, 3)
    editor.add_header(pdb_id="TEST", date=test_date)
    
    line = editor._header_lines_to_insert[0]
    assert line.startswith("HEADER")
    assert "TEST" in line
    assert "03-SEP-25" in line.upper()


def test_add_remark_reference():
    """Test structured REMARK reference block."""
    editor = PDBFileEditor()
    editor.add_remark_reference(
        authors="First Author, Second Author, Third Author",
        title="A Title",
        journal="Journal",
        volume="1",
        page="100",
        year="2023",
        issn="1234-5678",
        pmid="12345678",
        doi="10.1234/test.doi"
    )
    
    # Just check the number of lines and a few key pieces of content
    assert len(editor._header_lines_to_insert) == 6
    lines = [line.strip() for line in editor._header_lines_to_insert]
    assert "REMARK   1 REFERENCE  1" in lines[0]
    assert "REMARK   1  AUTH First Author, Second Author, Third Author" in lines[1]
    assert "REMARK   1  TITL A Title" in lines[2]
    assert "REMARK   1  REF    Journal              V. 1      100 2023" in lines[3]
    assert "REMARK   1  PMID  12345678" in lines[4]
    assert "REMARK   1  DOI   10.1234/test.doi" in lines[5]
    

def test_add_dbref():
    """Test DBREF record formatting."""
    editor = PDBFileEditor()
    editor.add_dbref(
        pdb_id="TEST",
        chain_id="A",
        seq_begin=1,
        seq_end=10,
        db="UNP",
        db_accession="P12345",
        db_id="P12345_TEST",
        db_seq_begin=50,
        db_seq_end=59
    )
    
    line = editor._header_lines_to_insert[0]
    assert line.startswith("DBREF")
    assert "TEST   A    1    10   UNP    P12345  P12345_TEST      50      59" in line


def test_add_seqres():
    """Test SEQRES record formatting and wrapping."""
    editor = PDBFileEditor()
    test_sequence = ["ALA", "LYS", "VAL"] * 5 # 15 residues
    editor.add_seqres(chain_id="A", sequence_list=test_sequence)
    
    # Expect 2 lines since 15 residues > 13 per line limit
    assert len(editor._header_lines_to_insert) == 2
    
    # Check the first line content
    line1 = editor._header_lines_to_insert[0].strip()
    assert line1.startswith("SEQRES   1 A   15  ALA LYS VAL ALA LYS VAL ALA LYS VAL ALA LYS VAL")
    
    # Check the second line content
    line2 = editor._header_lines_to_insert[1].strip()
    assert line2.startswith("SEQRES   2 A   15  ALA LYS VAL")
    
def test_record_order_key():
    """Test the internal method for sorting records."""
    editor = PDBFileEditor()
    
    # Check the correct order for known records
    assert editor._get_record_order_key("HEADER    ") < editor._get_record_order_key("TITLE     ")
    assert editor._get_record_order_key("REMARK    ") < editor._get_record_order_key("DBREF     ")
    
    # Check that unknown records are placed at the end
    unknown_key = editor._get_record_order_key("UNKNOWN_RECORD")
    atom_key = editor._get_record_order_key("ATOM      ")
    assert unknown_key > atom_key