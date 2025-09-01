from afdb_integration_kit.modelcif.pdb import PDBFileEditor
import datetime

# Create a sample "bare" PDB file content with only ATOM data.
BARE_PDB_CONTENT = [
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

# Create an instance of the editor with the bare PDB content
# pdb_editor = PDBFileEditor(pdb_content=BARE_PDB_CONTENT)

pdb_editor = PDBFileEditor()
pdb_editor.load_pdb("../model_cif_conversion/external_files/alphafold-upload-v4/alphafold-upload/AF-0000000065710001-v1.pdb")


# 1. Add HEADER record
pdb_editor.add_header(pdb_id="XXXX", date=datetime.date(2022, 6, 1))

# 2. Add TITLE records
title = "ALPHAFOLD MONOMER V2.0 PREDICTION FOR PROBABLE DISEASE RESISTANCE PROTEIN AT1G58602 (Q8W3K0)"
pdb_editor.add_title(title)

# 3. Add COMPND records
pdb_editor.add_compnd(
    molecule_id=1,
    molecule="PROBABLE DISEASE RESISTANCE PROTEIN AT1G58602",
    chain="A"
)

# 4. Add SOURCE records
pdb_editor.add_source(
    molecule_id=1,
    organism_scientific="ARABIDOPSIS THALIANA",
    organism_taxid=3702
)

# 5. Use the new add_remark_reference method to add the citation
pdb_editor.add_remark_reference(
    authors="JOHN JUMPER, RICHARD EVANS, ALEXANDER PRITZEL, TIM GREEN, MICHAEL FIGURNOV, OLAF RONNEBERGER, KATHRYN TUNYASUVUNAKOOL, RUSS BATES, AUGUSTIN ZIDEK, ANNA POTAPENKO, ALEX BRIDGLAND, CLEMENS MEYER, SIMON A A KOHL, ANDREW J BALLARD, ANDREW COWIE, BERNARDINO ROMERA-PAREDES, STANISLAV NIKOLOV, RISHUB JAIN, JONAS ADLER, TREVOR BACK, STIG PETERSEN, DAVID REIMAN, ELLEN CLANCY, MICHAL ZIELINSKI, MARTIN STEINEGGER, MICHALINA PACHOLSKA, TAMAS BERGHAMMER, DAVID SILVER, ORIOL VINYALS, ANDREW W SENIOR, KORAY KAVUKCUOGLU, PUSHMEET KOHLI, DEMIS HASSABIS",
    title="HIGHLY ACCURATE PROTEIN STRUCTURE PREDICTION WITH ALPHAFOLD",
    journal="NATURE",
    volume="596",
    page="583",
    year="2021",
    issn="0028-0836",
    pmid="34265844",
    doi="10.1038/s41586-021-03819-2"
)

# 6. Add a blank line using the new method
pdb_editor.add_remark_blank_line()

# 7. Add other REMARK 1 records
pdb_editor.add_remark(
    remark_number=1,
    text=
    'DISCLAIMERS ALPHAFOLD DATA, COPYRIGHT (2021) DEEPMIND TECHNOLOGIES LIMITED. THE INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY AND CAUTION SHOULD BE EXERCISED IN ITS USE. IT IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, WHETHER EXPRESSED OR IMPLIED. NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION SHALL NOT INFRINGE THE RIGHTS OF ANY THIRD PARTY. THE INFORMATION IS NOT INTENDED TO BE A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE, DIAGNOSIS, OR TREATMENT, AND DOES NOT CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. IT IS AVAILABLE FOR ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE.',
    numbered=False
)

# 8. Add DBREF record
pdb_editor.add_dbref(
    pdb_id="XXXX",
    chain_id="A",
    seq_begin=1,
    seq_end=1138,
    db="UNP",
    db_accession="Q8W3K0",
    db_id="DRL9_ARATH",
    db_seq_begin=1,
    db_seq_end=1138
)

# 9. Add a SEQRES record with three-letter codes
three_letter_codes = ["ALA", "LYS", "VAL", "LEU", "GLY"] * 227 + ["ALA", "LYS", "VAL"]
pdb_editor.add_seqres(
    chain_id="A",
    sequence_list=three_letter_codes
)

# 10. Add CRYST1 and matrix records
pdb_editor.add_crystal_info(
    a=1.0,
    b=1.0,
    c=1.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0
)

# 11. Add a MODEL record just before the ATOM data
pdb_editor.add_model(model_number=1)

pdb_editor.write_pdb("modified_example.pdb")