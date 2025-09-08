import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict
import datetime

import gemmi
from afdb_integration_kit.utils.pdbeditor import PDBFileEditor
from afdb_integration_kit.utils.cifstorage import CifDataStorage

CAT_ATOM_SITE = "_atom_site."
CAT_ENTITY_POLY_SEQ = "_entity_poly_seq."
CAT_ENTITY_POLY_SEQ_SCHEME = "_pdbx_poly_seq_scheme."
CAT_STRUCT_ASYM = "_struct_asym."
CAT_TARGET_REF_DB = "_ma_target_ref_db_details."
ITEM_AUTH_ASYM_ID = "auth_asym_id"

# Configure logger
logger = logging.getLogger("afdb_integration_kit")

def load_cif_file(cif_path: str) -> CifDataStorage:
    """Loads a mmCIF file and returns a CifDataStorage instance."""
    logger.info(f"Reading mmCIF file: {cif_path}")
    try:
        doc = gemmi.cif.read(cif_path)
        cif_block = doc.sole_block()
        cif_data = CifDataStorage()
        cif_data.populate_from_cif_block(cif_block)
        return cif_data
    except Exception as e:
        logger.error(f"Could not read or parse mmCIF file '{cif_path}': {e}")
        sys.exit(1)


def load_json_file(path: str) -> Dict[str, Any]:
    """Loads and returns data from a JSON file."""
    logger.info(f"Loading JSON file: {path}")
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error reading or parsing JSON file '{path}': {e}")
        sys.exit(1)


def add_pdb_headers(pdb_editor: PDBFileEditor, cif_data: CifDataStorage, output_file_name: str, copyrights_str: str):
    """
    Enriches the PDB file with header information extracted from the mmCIF data.

    Args:
        pdb_editor: An instance of PDBFileEditor, pre-loaded with ATOM coordinates.
        cif_data: An instance of CifDataStorage containing all mmCIF data.
    """
    logger.info("Adding PDB header information from mmCIF data.")
    data = cif_data.get_data()

    # Get data from common categories
    citation = data.get("_citation.", {})
    citation_author = data.get("_citation_author.", {})
    entity_poly = data.get("_entity_poly.", {})
    entity = data.get("_entity.", {})
    target_ref_db = data.get("_ma_target_ref_db_details.", {})
    data_usage = data.get("_pdbx_data_usage.", {})
    seq_scheme = data.get("_pdbx_poly_seq_scheme.", {})
    db_status = data.get("_pdbx_database_status.", {})

    # 1. HEADER and AUDIT
    # The date is the most recent revision date from recvd_initial_deposition_date
    # The PDB ID is the entry ID from _pdbx_database_status
    date_str = db_status.get("recvd_initial_deposition_date", ["?"])[-1]
    # pdb_id = db_status.get("entry_id", ["?"])[-1].split("-")[-1]
    dbref_id = output_file_name

    try:
        if date_str != "?":
            date_obj = datetime.datetime.strptime(date_str, "%Y-%m-%d").date()
            pdb_editor.add_header(pdb_id=dbref_id.upper(), date=date_obj)
        else:
            logger.warning(
                "Cannot determine PDB header date. The recvd_initial_deposition_date "
                "category is missing revision_date."
            )
            pdb_editor.add_header(pdb_id=dbref_id, date=datetime.datetime.now().date())
    except (ValueError, IndexError):
        logger.warning("Invalid date format in recvd_initial_deposition_date. Skipping date.")
        pdb_editor.add_header(pdb_id=dbref_id)

    # 2. TITLE and COMPND
    # Use _entity.pdbx_description and _entity_poly.pdbx_seq_one_letter_code to build the title
    if entity.get("pdbx_description") and entity_poly.get("pdbx_seq_one_letter_code"):
        description = entity["pdbx_description"][0]
        seq = entity_poly["pdbx_seq_one_letter_code"][0]
        title_text = f"{description} ({target_ref_db.get('db_accession', ['?'])[0]})"
        pdb_editor.add_title(title_text.upper())

        # COMPND information
        # Get chain_id from _struct_ref_seq to handle multiple chains correctly
        chains_present = set(cif_data.data.get(CAT_ATOM_SITE, {}).get(ITEM_AUTH_ASYM_ID, []))
        for chain_id in sorted(chains_present):
            entity_id = cif_data.data.get(CAT_STRUCT_ASYM, {}).get("entity_id", ["?"])[0] # Assuming one entity for simplicity based on prompt
            pdb_editor.add_compnd(
                molecule_id=entity_id,
                molecule=description.upper(),
                chain=chain_id
            )
    else:
        logger.warning(
            "Missing _entity or _entity_poly data. Cannot generate TITLE and COMPND records."
        )

    # 3. SOURCE
    # Using the first entry from _ma_target_ref_db_details which contains UNIPROT info
    if target_ref_db.get("organism_scientific") and target_ref_db.get("ncbi_taxonomy_id"):
        scientific_name = target_ref_db["organism_scientific"][0]
        taxonomy_id = target_ref_db["ncbi_taxonomy_id"][0]
        pdb_editor.add_source(
            molecule_id=target_ref_db.get("target_entity_id", ["?"])[0],
            organism_scientific=scientific_name.upper(),
            organism_taxid=int(taxonomy_id) if taxonomy_id != "?" else None
        )
    else:
        logger.warning(
            "Missing _ma_target_ref_db_details for SOURCE record. UniProt data might be missing."
        )

    # 4. REMARK section
    # Find all non primary citations (primary citations go in JRNL)
    non_primary_citation = None
    if "_citation." in data and "_citation_author." in data:
        citation = data["_citation."]
        reference_number = 1
        for citation_id, doi in zip(citation.get("id", []), citation.get("pdbx_database_id_DOI", [])):
            if citation_id != "primary":
                try:
                    idx = citation.get("id", []).index(citation_id)
                except ValueError:
                    continue  # skip if not found

            def safe_get(key: str) -> str:
                values = citation.get(key, [])
                if values[idx] is None or idx >= len(values):
                    return ""                    
                return values[idx]

            non_primary_citation_data = {
                "title": safe_get("title"),
                "journal": safe_get("journal_full"),
                "volume": safe_get("journal_volume"),
                "page": safe_get("page_first"),
                "year": safe_get("year"),
                "doi": doi or "",
                "pmid": safe_get("pdbx_database_id_PubMed"),
                "issn": safe_get("journal_id_ISSN"),
            }
                
            # Get the authors for this specific citation
            citation_authors = [
                name for name, cit_id in zip(
                    citation_author.get("name", []),
                    citation_author.get("citation_id", [])
                ) if str(cit_id) == str(citation_id)
            ]
            # remove commas from name in citation author list
            citation_authors = [author.replace(",", "") for author in citation_authors]
            non_primary_citation_data["authors"] = ",".join(citation_authors)
            non_primary_citation = non_primary_citation_data

            if non_primary_citation:
                pdb_editor.add_remark_reference(
                    authors=non_primary_citation.get("authors", "?").upper(),
                    title=non_primary_citation.get("title", "?").upper(),
                    journal=non_primary_citation.get("journal", "?").upper(),
                    volume=non_primary_citation.get("volume", "?").upper(),
                    page=non_primary_citation.get("page", "?").upper(),
                    year=non_primary_citation.get("year", "?").upper(),
                    doi=non_primary_citation.get("doi", "?").upper(),
                    pmid=non_primary_citation.get("pmid", "?").upper(),
                    issn=non_primary_citation.get("issn", "?").upper(),
                    reference_number=reference_number
                )
                reference_number += 1
            pdb_editor.add_remark_blank_line()
 
    if data_usage.get("details"):
        # disclaimer = data_usage.get("details", ["?"])[0]
        cleaned_copyrights_str = copyrights_str.encode('ascii', errors='ignore').decode('ascii')
        disclaimer = (
            f'{cleaned_copyrights_str.upper()} THE INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY '
            'AND CAUTION SHOULD BE EXERCISED IN ITS USE. IT IS PROVIDED "AS-IS" '
            'WITHOUT ANY WARRANTY OF ANY KIND, WHETHER EXPRESSED OR IMPLIED. '
            'NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION SHALL NOT INFRINGE THE RIGHTS OF ANY THIRD PARTY. '
            'THE INFORMATION IS NOT INTENDED TO BE A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE '
            'DIAGNOSIS, OR TREATMENT, AND DOES NOT CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. '
            'IT IS AVAILABLE FOR ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE. '

        )
        pdb_editor.add_remark(remark_number=1, text=disclaimer, numbered=False)

    # 5. DBREF records
    if target_ref_db:
        # 1. Build a mapping from entity_id to UniProt details
        uniprot_by_entity = {}
        if "target_entity_id" in target_ref_db:
            for i, entity_id in enumerate(target_ref_db["target_entity_id"]):
                uniprot_by_entity[entity_id] = {
                    "db_accession": target_ref_db.get("db_accession", ["?"])[i],
                    "db_name": target_ref_db.get("db_name", ["?"])[i],
                    "db_code": target_ref_db.get("db_code", ["?"])[i],
                    "seq_db_align_begin": target_ref_db.get("seq_db_align_begin", ["?"])[i],
                    "seq_db_align_end": target_ref_db.get("seq_db_align_end", ["?"])[i],
                }
        
        # 2. Get the chain-to-entity mapping from _struct_asym and
        #    the min/max seq_id from _pdbx_poly_seq_scheme for the PDB sequence range.
        struct_asym_data = data.get(CAT_STRUCT_ASYM, {})
        chain_to_entity_map = {}
        if "id" in struct_asym_data and "entity_id" in struct_asym_data:
            for chain_id, entity_id in zip(struct_asym_data["id"], struct_asym_data["entity_id"]):
                chain_to_entity_map[chain_id] = entity_id

        pdb_seq_range_by_chain = defaultdict(list)
        seq_scheme = data.get(CAT_ENTITY_POLY_SEQ_SCHEME, {})
        if "asym_id" in seq_scheme and "seq_id" in seq_scheme:
            for asym_id, seq_id in zip(seq_scheme["asym_id"], seq_scheme["seq_id"]):
                pdb_seq_range_by_chain[asym_id].append(int(seq_id))

        # 3. Iterate through chains and create DBREF records
        dbref_id = "XXXX"
        chains_in_pdb = sorted(set(cif_data.data.get(CAT_ATOM_SITE, {}).get(ITEM_AUTH_ASYM_ID, [])))
        
        for pdb_chain_id in chains_in_pdb:
            entity_id = chain_to_entity_map.get(pdb_chain_id)
            if entity_id and entity_id in uniprot_by_entity:
                uniprot_info = uniprot_by_entity[entity_id]
                pdb_seq_ids = pdb_seq_range_by_chain.get(pdb_chain_id, [])
                
                if not pdb_seq_ids:
                    logger.warning(f"No sequence data found for chain '{pdb_chain_id}'. Skipping DBREF.")
                    continue
                
                try:
                    pdb_editor.add_dbref(
                        pdb_id=dbref_id,
                        chain_id=pdb_chain_id,
                        seq_begin=min(pdb_seq_ids),
                        seq_end=max(pdb_seq_ids),
                        db=uniprot_info["db_name"],
                        db_accession=uniprot_info["db_accession"] if len(uniprot_info["db_accession"]) < 8 else "XXXX",
                        db_id=uniprot_info["db_code"] if len(uniprot_info["db_code"]) < 12 else "XXXX",
                        db_seq_begin=uniprot_info["seq_db_align_begin"],
                        db_seq_end=uniprot_info["seq_db_align_end"]
                    )
                except (KeyError, IndexError):
                    logger.warning(
                        f"Incomplete UniProt data for entity '{entity_id}'. "
                        "Cannot add DBREF record."
                    )
            else:
                logger.warning(
                    f"No UniProt data found for chain '{pdb_chain_id}' (entity '{entity_id}'). "
                    "Cannot add DBREF record."
                )
    else:
        logger.warning("Missing _ma_target_ref_db_details for DBREF records. Skipping.")



    # 6. SEQRES records
    # Get sequence data from _entity_poly_seq and group by chain
    if entity_poly and entity_poly.get("pdbx_seq_one_letter_code"):
        seq_scheme = data.get(CAT_ENTITY_POLY_SEQ_SCHEME, {})
        if seq_scheme and 'entity_id' in seq_scheme:
            seq_by_chain = defaultdict(list)
            for entity_id, asym_id, mon_id in zip(
                seq_scheme.get('entity_id', []),
                seq_scheme.get('asym_id', []),
                seq_scheme.get('mon_id', [])
            ):
                seq_by_chain[asym_id].append(mon_id)
            
            for chain_id, sequence in seq_by_chain.items():
                pdb_editor.add_seqres(chain_id=chain_id, sequence_list=sequence)
        else:
            logger.warning("Missing _pdbx_poly_seq_scheme for SEQRES records. Skipping.")
    else:
        logger.warning("Missing _entity_poly data. Cannot add SEQRES records.")

    # 7. CRYST1 record
    pdb_editor.add_crystal_info(
    a=1.0,
    b=1.0,
    c=1.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0
)

    # 8. Add MODEL record before ATOMs
    pdb_editor.add_model(model_number=1)
    

def generate_pdb_headers(
    cif_file: str, pdb_file: str, output_pdb_file: str, provider_json_file: str
):
    """
    Main function to orchestrate the generation of a PDB file with headers from a mmCIF file.
    """
    cif_data = load_cif_file(cif_file)
    provider_json = load_json_file(provider_json_file)
    copyrights = provider_json.get("copyrights", [])
    copyrights_str = ". ".join(copyrights)
    if copyrights_str[-1] != ".":
        copyrights_str += "."

    pdb_editor = PDBFileEditor()
    pdb_editor.load_pdb(pdb_file)

    add_pdb_headers(pdb_editor, cif_data, Path(output_pdb_file).stem, copyrights_str)
    pdb_editor.validate_ter_records()
    pdb_editor.write_pdb(output_pdb_file)
    logger.info(f"PDB file with headers written to: {output_pdb_file}")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(
            "Usage: python script_name.py <mmcif_file> <input_pdb_file> <output_pdb_file> <provider_json_file>"
        )
        sys.exit(1)

    mmcif_file = sys.argv[1]
    input_pdb_file = sys.argv[2]
    output_pdb_file = sys.argv[3]
    provider_json_file = sys.argv[4]
    generate_pdb_headers(mmcif_file, input_pdb_file, output_pdb_file, provider_json_file)
