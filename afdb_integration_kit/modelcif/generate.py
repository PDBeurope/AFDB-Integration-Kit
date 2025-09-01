import hashlib
import json
import logging
import shutil
import subprocess
import sys
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict, List, Optional, TypedDict

import gemmi
import jsonschema
import numpy as np
import requests
from afdb_integration_kit.modelcif.pdb import PDBFileEditor
import datetime

# --- Configuration & Constants ---

# Configure logger
logger = logging.getLogger("modelcif")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)s] %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
JSON_SCHEMA_PATH = str(
    files("afdb_integration_kit.modelcif.resources").joinpath("schema.json")
)

# MMCIF Category and Item Constants
CAT_ATOM_SITE = "_atom_site."
CAT_CHEM_COMP = "_chem_comp"
CAT_SOFTWARE = "_software"
CAT_MODEL_LIST = "_ma_model_list"
CAT_TARGET_REF_DB = "_ma_target_ref_db_details"
CAT_GLOBAL_QA = "_ma_qa_metric_global"
CAT_LOCAL_QA = "_ma_qa_metric_local"
CAT_ENTITY_POLY_SEQ = "_entity_poly_seq"
CAT_ENTITY_POLY_SEQ_SCHEME = "_pdbx_poly_seq_scheme"
CAT_STRUCT_ASYM = "_struct_asym."

ITEM_B_FACTOR = "B_iso_or_equiv"
ITEM_LABEL_ASYM_ID = "label_asym_id"
ITEM_AUTH_ASYM_ID = "auth_asym_id"
ITEM_LABEL_COMP_ID = "label_comp_id"
ITEM_AUTH_COMP_ID = "auth_comp_id"
ITEM_LABEL_SEQ_ID = "label_seq_id"
ITEM_AUTH_SEQ_ID = "auth_seq_id"
ITEM_PDB_INS_CODE = "pdbx_PDB_ins_code"


UNIPROT_API_BASE_URL = "https://rest.uniprot.org/uniprotkb"

# --- Data Structures (Type Definitions) ---


class ChainMetadata(TypedDict):
    uniprot_accession: str
    chain_id: str


class InputMetadata(TypedDict):
    metadata: Dict[str, Any]
    categories: Dict[str, Any]
    chains: List[ChainMetadata]


# --- Data Storage Class ---


class CifDataStorage:
    """A container for holding and writing mmCIF data."""

    def __init__(self):
        self.data: Dict[str, Dict[str, List[Any]]] = {}

    def set_items(self, category_name: str, items_dict: Dict[str, List[Any]]):
        if category_name not in self.data:
            self.data[category_name] = {}
        for item, values in items_dict.items():
            self.data[category_name][item] = values

    def set_item(self, category_name: str, item_name: str, item_value: Any):
        if category_name not in self.data:
            self.data[category_name] = {}
        self.data[category_name][item_name] = item_value

    def get_data(self) -> Dict[str, Dict[str, List[Any]]]:
        return self.data

    def populate_from_cif_block(self, cif_block: gemmi.cif.Block):
        """Initializes the storage from a gemmi cif.Block."""
        for category in cif_block.get_mmcif_category_names():
            items = cif_block.get_mmcif_category(category)
            self.set_items(category, items)

        # Perform initial data mappings required for consistency
        self.data[CAT_ATOM_SITE][ITEM_LABEL_ASYM_ID] = self.data[CAT_ATOM_SITE][
            ITEM_AUTH_ASYM_ID
        ]
        self.data[CAT_ATOM_SITE][ITEM_AUTH_COMP_ID] = self.data[CAT_ATOM_SITE][
            ITEM_LABEL_COMP_ID
        ]
        self.data[CAT_ATOM_SITE][ITEM_LABEL_SEQ_ID] = self.data[CAT_ATOM_SITE][
            ITEM_AUTH_SEQ_ID
        ]

    def write_to_cif(self, output_file: str, block_name: str = "model"):
        """Writes the stored data to an mmCIF file."""
        logger.info("Writing CIF file...")
        doc = gemmi.cif.Document()
        block = doc.add_new_block(block_name)

        for category, items in self.data.items():
            block.set_mmcif_category(category, items)

        write_options = gemmi.cif.WriteOptions()
        write_options.align_loops = 50
        write_options.align_pairs = 50
        write_options.prefer_pairs = True
        doc.write_file(output_file, write_options)
        logger.info(f"mmCIF file written to: {output_file}")


# --- External Service Client ---


class UniprotAPIClient:
    """Client for fetching data from the UniProt API."""

    def __init__(self, session: requests.Session, base_url: str = UNIPROT_API_BASE_URL):
        self.session = session
        self.base_url = base_url

    def fetch_metadata(self, uniprot_accession: str) -> Dict[str, Any]:
        """Fetches metadata for a given UniProt ID."""
        url = f"{self.base_url}/{uniprot_accession}.json"
        logger.info(f"Fetching UniProt metadata for ID: {uniprot_accession}")
        try:
            response = self.session.get(url, timeout=15)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            logger.warning(
                f"Failed to fetch UniProt ID {uniprot_accession}: {e}. "
                "Please check the ID and your network connection."
            )
            return {}


# --- Core Processing Functions (Pure Functions) ---


def process_uniprot_response(data: Dict[str, Any]) -> Dict[str, str]:
    """Processes the JSON response from UniProt into a flat dictionary."""
    if not data:
        return {
            key: "?"
            for key in [
                "db_accession",
                "db_code",
                "db_name",
                "gene_name",
                "ncbi_taxonomy_id",
                "organism_scientific",
                "seq_db_align_begin",
                "seq_db_align_end",
                "seq_db_isoform",
                "seq_db_sequence_checksum",
                "seq_db_sequence_version_date",
            ]
        }

    seq = data.get("sequence", {}).get("value", "")
    crc_64 = data.get("sequence", {}).get("crc64", None)

    return {
        "db_accession": data.get("primaryAccession", None),
        "db_code": data.get("uniProtkbId", None),
        "db_name": "UNP",
        "gene_name": data.get("genes", [{}])[0].get("geneName", {}).get("value", None),
        "ncbi_taxonomy_id": str(data.get("organism", {}).get("taxonId", None)),
        "organism_scientific": data.get("organism", {}).get("scientificName", None),
        "seq_db_align_begin": "1",
        "seq_db_align_end": str(len(seq)) if seq else None,
        "seq_db_isoform": None,
        "seq_db_sequence_checksum": crc_64,
        "seq_db_sequence_version_date": data.get("entryAudit", {}).get(
            "lastSequenceUpdateDate", None
        ),
    }


def compute_global_plddt(b_factors: List[str]) -> float:
    """Computes the average pLDDT from a list of B-factor strings."""
    logger.info("Computing global pLDDT...")
    try:
        b_floats = [float(v) for v in b_factors if v not in ("?", ".", "")]
        return float(np.mean(b_floats)) if b_floats else -1.0
    except (ValueError, TypeError) as e:
        logger.error(f"Error parsing B-factors for global pLDDT: {e}")
        return -1.0


def compute_local_plddt_metrics(
    asym_ids: List[str], comp_ids: List[str], seq_ids: List[str], b_factors: List[str]
) -> Dict[str, List[Any]]:
    """Computes local pLDDT per residue and returns a dictionary for the local
    QA category."""
    logger.info("Computing local pLDDT...")
    residue_b_factors = defaultdict(list)
    residue_info = {}

    for asym_id, comp_id, seq_id, b_factor in zip(
        asym_ids, comp_ids, seq_ids, b_factors
    ):
        if b_factor not in ("?", ".", ""):
            try:
                residue_key = (asym_id, comp_id, seq_id)
                residue_b_factors[residue_key].append(float(b_factor))
                residue_info[residue_key] = {
                    "label_asym_id": asym_id,
                    "label_comp_id": comp_id,
                    "label_seq_id": seq_id,
                }
            except (ValueError, TypeError):
                continue

    local_metrics = defaultdict(list)
    for i, (residue_key, b_list) in enumerate(residue_b_factors.items()):
        mean_plddt = np.mean(b_list)
        info = residue_info[residue_key]
        local_metrics["label_asym_id"].append(info["label_asym_id"])
        local_metrics["label_comp_id"].append(info["label_comp_id"])
        local_metrics["label_seq_id"].append(info["label_seq_id"])
        local_metrics["metric_id"].append("2")
        local_metrics["metric_value"].append(f"{mean_plddt:.2f}")
        local_metrics["model_id"].append("1")
        local_metrics["ordinal_id"].append(str(i + 1))

    logger.info(
        f"Computed local pLDDT for {len(local_metrics['ordinal_id'])} residues."
    )
    return dict(local_metrics)


def create_polymer_sequence_categories(
    atom_site_data: Dict[str, List[Any]],
) -> Dict[str, Dict[str, List[Any]]]:
    """
    Creates both _entity_poly_seq and _ma_entity_poly_seq_scheme categories
    by processing unique residues from the _atom_site data.
    """
    logger.info(
        "Creating polymer sequence categories "
        "(_entity_poly_seq and _ma_entity_poly_seq_scheme)..."
    )
    if not all(
        k in atom_site_data
        for k in [ITEM_LABEL_ASYM_ID, ITEM_LABEL_SEQ_ID, "label_entity_id"]
    ):
        logger.error(
            "Cannot create polymer sequence categories: "
            "required columns are missing from _atom_site."
        )
        return {}

    unique_residues = {}
    # Use .get() to safely access columns that might be missing from some files
    ins_codes = atom_site_data.get(
        ITEM_PDB_INS_CODE, ["?"] * len(atom_site_data[ITEM_LABEL_ASYM_ID])
    )

    for i, asym_id in enumerate(atom_site_data[ITEM_LABEL_ASYM_ID]):
        seq_id = atom_site_data[ITEM_LABEL_SEQ_ID][i]
        residue_key = (asym_id, seq_id)
        if residue_key not in unique_residues:
            unique_residues[residue_key] = {
                "asym_id": asym_id,
                "entity_id": atom_site_data["label_entity_id"][i],
                "seq_id": seq_id,
                "mon_id": atom_site_data[ITEM_LABEL_COMP_ID][i],
                "auth_seq_num": atom_site_data[ITEM_AUTH_SEQ_ID][i],
                "pdb_ins_code": ins_codes[i],
                "pdb_mon_id": atom_site_data[ITEM_AUTH_COMP_ID][i],
                "pdb_strand_id": atom_site_data[ITEM_AUTH_ASYM_ID][i],
            }

    # Sort residues by entity, chain, and sequence number for consistent ordering
    sorted_residues = sorted(
        unique_residues.values(),
        key=lambda r: (r["entity_id"], r["asym_id"], int(r["seq_id"])),
    )

    # --- Build _ma_entity_poly_seq_scheme ---
    scheme_data = defaultdict(list)
    for r in sorted_residues:
        scheme_data["asym_id"].append(r["asym_id"])
        scheme_data["auth_seq_num"].append(r["auth_seq_num"])
        scheme_data["entity_id"].append(r["entity_id"])
        scheme_data["hetero"].append("n")
        scheme_data["mon_id"].append(r["mon_id"])
        scheme_data["pdb_ins_code"].append(r["pdb_ins_code"])
        scheme_data["pdb_mon_id"].append(r["pdb_mon_id"])
        scheme_data["pdb_seq_num"].append(r["auth_seq_num"])  # Same as auth_seq_num
        scheme_data["pdb_strand_id"].append(r["pdb_strand_id"])
        scheme_data["seq_id"].append(r["seq_id"])

    # --- Build _entity_poly_seq ---
    poly_seq_data = defaultdict(list)
    entity_sequences = defaultdict(dict)
    for r in sorted_residues:
        entity_sequences[r["entity_id"]][int(r["seq_id"])] = r["mon_id"]

    for entity_id, residues in sorted(entity_sequences.items()):
        for seq_num, mon_id in sorted(residues.items()):
            poly_seq_data["entity_id"].append(entity_id)
            poly_seq_data["num"].append(str(seq_num))
            poly_seq_data["mon_id"].append(mon_id)
            poly_seq_data["hetero"].append("n")

    return {
        CAT_ENTITY_POLY_SEQ: dict(poly_seq_data),
        CAT_ENTITY_POLY_SEQ_SCHEME: dict(scheme_data),
    }


def map_entities_and_chains(
    cif_data: CifDataStorage, json_chains_info: Optional[List[ChainMetadata]]
):
    """
    Assigns entity IDs based on chain information, handling single vs. multiple chains.
    Updates cif_data in place.
    """
    atom_site_data = cif_data.get_data().get(CAT_ATOM_SITE, {})
    if not atom_site_data:
        logger.error("Cannot map entities: _atom_site data is missing.")
        return

    pdb_asym_ids = sorted(set(atom_site_data.get(ITEM_LABEL_ASYM_ID, [])))
    num_pdb_chains = len(pdb_asym_ids)
    entity_id_map: Dict[str, str] = {}

    if num_pdb_chains == 0:
        logger.warning("No chains found in the PDB data to map.")
        return

    if num_pdb_chains == 1:
        logger.info("Single chain detected. Assigning entity_id '1'.")
        entity_id_map = {pdb_asym_ids[0]: "1"}
    else:  # Multiple chains
        logger.info(
            f"Multiple chains ({num_pdb_chains}) detected. Using JSON for"
            "entity mapping."
        )

        if not json_chains_info:
            logger.error(
                "CRITICAL: Multiple chains found in PDB, but the 'chains' section"
                " is missing in the JSON metadata. Cannot map chains to entities. "
                "Please provide this mapping."
            )
            sys.exit(1)

        if len(json_chains_info) != num_pdb_chains:
            logger.warning(
                f"Mismatch: PDB file has {num_pdb_chains} chains,"
                f"but JSON provides info for {len(json_chains_info)}."
            )

        json_chain_ids = {chain.get("chain_id") for chain in json_chains_info}
        if set(pdb_asym_ids) != json_chain_ids:
            pdb_set = set(pdb_asym_ids)
            json_set = json_chain_ids
            logger.warning(
                f"Chain ID mismatch! "
                f"Chains in PDB but not JSON: {pdb_set - json_set}. "
                f"Chains in JSON but not PDB: {json_set - pdb_set}."
            )

        for chain in json_chains_info:
            chain_id, entity_id = chain.get("chain_id"), chain.get("entity_id")
            if not chain_id or not entity_id:
                logger.warning(f"Skipping incomplete chain entry in JSON: {chain}")
                continue
            entity_id_map[chain_id] = entity_id

    if not entity_id_map:
        logger.error(
            "CRITICAL: Failed to create a valid chain_id to entity_id map. "
            "Please check the 'chains' section in your JSON file."
        )
        sys.exit(1)

    logger.info(f"Applying entity map: {entity_id_map}")

    # 1. Update _atom_site.label_entity_id
    all_asym_ids_in_order = atom_site_data.get(ITEM_LABEL_ASYM_ID, [])
    cif_data.data[CAT_ATOM_SITE]["label_entity_id"] = [
        entity_id_map.get(asym, "?") for asym in all_asym_ids_in_order
    ]

    # 2. Create/Update _struct_asym
    sorted_asym_ids = sorted(entity_id_map.keys())
    cif_data.set_items(
        CAT_STRUCT_ASYM,
        {
            "id": sorted_asym_ids,
            "entity_id": [entity_id_map[asym_id] for asym_id in sorted_asym_ids],
        },
    )


def add_standard_chem_comp_data(cif_data: CifDataStorage):
    """Adds the hardcoded _chem_comp category for the 20 standard amino acids."""
    logger.info("Adding hardcoded _chem_comp category for standard amino acids.")
    chem_comp_data = {
        "formula": [
            "C3 H7 N O2",
            "C6 H15 N4 O2",
            "C4 H8 N2 O3",
            "C4 H7 N O4",
            "C3 H7 N O2 S",
            "C5 H10 N2 O3",
            "C5 H9 N O4",
            "C2 H5 N O2",
            "C6 H10 N3 O2",
            "C6 H13 N O2",
            "C6 H13 N O2",
            "C6 H15 N2 O2",
            "C5 H11 N O2 S",
            "C9 H11 N O2",
            "C5 H9 N O2",
            "C3 H7 N O3",
            "C4 H9 N O3",
            "C11 H12 N2 O2",
            "C9 H11 N O3",
            "C5 H11 N O2",
        ],
        "formula_weight": [
            "89.093",
            "175.209",
            "132.118",
            "133.103",
            "121.158",
            "146.144",
            "147.129",
            "75.067",
            "156.162",
            "131.173",
            "131.173",
            "147.195",
            "149.211",
            "165.189",
            "115.130",
            "105.093",
            "119.119",
            "204.225",
            "181.189",
            "117.146",
        ],
        "id": [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ],
        "mon_nstd_flag": ["y"] * 20,
        "name": [
            "ALANINE",
            "ARGININE",
            "ASPARAGINE",
            "ASPARTIC ACID",
            "CYSTEINE",
            "GLUTAMINE",
            "GLUTAMIC ACID",
            "GLYCINE",
            "HISTIDINE",
            "ISOLEUCINE",
            "LEUCINE",
            "LYSINE",
            "METHIONINE",
            "PHENYLALANINE",
            "PROLINE",
            "SERINE",
            "THREONINE",
            "TRYPTOPHAN",
            "TYROSINE",
            "VALINE",
        ],
        "pdbx_synonyms": [None] * 20,
        "type": [
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
            "L-PEPTIDE LINKING",
        ],
    }
    cif_data.set_items(CAT_CHEM_COMP, chem_comp_data)


def load_json_file(path: str) -> InputMetadata:
    logger.info(f"Loading metadata from: {path}")
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error reading or parsing JSON file '{path}': {e}")
        sys.exit(1)


def pdb_to_cif_block(pdb_path: str) -> gemmi.cif.Block:
    logger.info(f"Reading PDB file and converting to CIF block: {pdb_path}")
    try:
        structure = gemmi.read_structure(pdb_path)
        return structure.make_mmcif_block()
    except Exception as e:
        logger.error(f"Could not read or process PDB file '{pdb_path}': {e}")
        sys.exit(1)


def validate_json_with_schema(data: Dict[str, Any], schema_path: str):
    """Validates the given data against a JSON schema."""
    logger.info(f"Validating metadata against schema: {schema_path}")

    try:
        with open(schema_path, "r") as f:
            schema = json.load(f)
    except FileNotFoundError:
        logger.error(f"CRITICAL: JSON schema file not found at '{schema_path}'.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        logger.error(f"CRITICAL: Could not parse JSON schema file '{schema_path}': {e}")
        sys.exit(1)

    try:
        jsonschema.validate(instance=data, schema=schema)
        logger.info("Metadata JSON validation successful.")
    except jsonschema.ValidationError as e:
        logger.error("CRITICAL: Metadata JSON validation failed.")
        logger.error(f"Validation Error: {e.message}")
        logger.error(f"Location: {' -> '.join(map(str, e.path))}")
        sys.exit(1)


def validate_with_gemmi(cif_path: str, dict_path: str):
    """Validates the CIF file using the external 'gemmi validate' command."""
    if not shutil.which("gemmi"):
        logger.error(
            "Validation skipped: The 'gemmi' command-line program is not in your PATH. "
            "Please run 'pip install gemmi-program'."
        )
        return
    if not Path(dict_path).is_file():
        logger.error(f"Validation skipped: Dictionary not found at '{dict_path}'.")
        return

    logger.info(f"Validating '{cif_path}' against dictionary '{dict_path}'...")
    command = ["gemmi", "validate", "-p", "-d", dict_path, cif_path]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            logger.error(f"Gemmi validator exited with error code {result.returncode}.")
            print(f"\n--- STDERR ---\n{result.stderr.strip()}\n--------------")

        if result.stdout.strip():
            logger.warning("Validation found the following issues:")
            print(
                f"\n--- GEMMI VALIDATION REPORT ---\n{result.stdout.strip()}"
                "\n-----------------------------"
            )
        else:
            logger.info("Validation successful: No issues found by gemmi.")

    except subprocess.TimeoutExpired:
        logger.error("Validation command timed out after 60 seconds.")
    except Exception as e:
        logger.error(f"An unexpected error occurred during validation: {e}")

def add_pdb_headers(pdb_editor: PDBFileEditor, cif_data: CifDataStorage, output_file_name: str):
    """
    Enriches the PDB file with header information extracted from the mmCIF data.

    Args:
        pdb_editor: An instance of PDBFileEditor, pre-loaded with ATOM coordinates.
        cif_data: An instance of CifDataStorage containing all mmCIF data.
    """
    logger.info("Adding PDB header information from mmCIF data.")
    data = cif_data.get_data()

    # Get data from common categories
    audit_author = data.get("_audit_author", {})
    citation = data.get("_citation", {})
    citation_author = data.get("_citation_author", {})
    entity_poly = data.get("_entity_poly", {})
    entity = data.get("_entity", {})
    target_ref_db = data.get("_ma_target_ref_db_details", {})
    data_usage = data.get("_pdbx_data_usage", {})
    seq_scheme = data.get("_pdbx_poly_seq_scheme", {})
    db_status = data.get("_pdbx_database_status", {})
    software = data.get("_software", {})

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
            organism_scientific=scientific_name,
            organism_taxid=int(taxonomy_id) if taxonomy_id != "?" else None
        )
    else:
        logger.warning(
            "Missing _ma_target_ref_db_details for SOURCE record. UniProt data might be missing."
        )

    # 4. REMARK section
    # Find all non primary citations (primary citations go in JRNL)
    non_primary_citation = None
    if "_citation" in data and "_citation_author" in data:
        citation = data["_citation"]
        reference_number = 1
        for citation_id, doi in zip(citation.get("id", []), citation.get("pdbx_database_id_DOI", [])):
            if citation_id != "primary":
                try:
                    idx = citation.get("id", []).index(citation_id)
                except ValueError:
                    continue  # skip if not found

            def safe_get(key: str) -> str:
                values = citation.get(key, [])
                return values[idx] if idx < len(values) else ""

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
        disclaimer = 'THE INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY \
            AND CAUTION SHOULD BE EXERCISED IN ITS USE. IT IS PROVIDED "AS-IS"\
            WITHOUT ANY WARRANTY OF ANY KIND, WHETHER EXPRESSED OR IMPLIED.\
            NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION SHALL NOT \
            INFRINGE THE RIGHTS OF ANY THIRD PARTY. THE INFORMATION IS NOT\
            INTENDED TO BE A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE, \
            DIAGNOSIS, OR TREATMENT, AND DOES NOT CONSTITUTE MEDICAL OR OTHER \
            PROFESSIONAL ADVICE. IT IS AVAILABLE FOR ACADEMIC AND COMMERCIAL \
            PURPOSES, UNDER CC-BY 4.0 LICENCE.'
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
    



# --- Main Orchestration ---
def generate(
    pdb_file: str, metadata_file: str, output_file: str, validate_dict_path: str
):
    """Main function to orchestrate the PDB to mmCIF conversion and enrichment."""
    # 1. Load initial data
    input_metadata = load_json_file(metadata_file)
    validate_json_with_schema(input_metadata, JSON_SCHEMA_PATH)
    cif_block = pdb_to_cif_block(pdb_file)

    # 2. Initialize data storage and populate from PDB
    cif_data = CifDataStorage()
    cif_data.populate_from_cif_block(cif_block)
    map_entities_and_chains(cif_data, input_metadata.get("chains"))
    add_standard_chem_comp_data(cif_data)

    # 3. Add metadata from JSON file
    model_meta = input_metadata.get("metadata", {})
    cif_data.set_item(
        CAT_MODEL_LIST,
        "model_group_name",
        [f"AlphaFold {model_meta.get(
            'model_type', 'Monomer')} v{model_meta.get('version', '3.1')} model"],
    )
    for category, items in input_metadata.get("categories", {}).items():
        if items:
            cif_data.set_items(category, items)
    # check if _software is alphafold. If yes get version from model_meta. If softwares are multiple log that versions should be correctly provided
    software_category = input_metadata.get("categories", {}).get(CAT_SOFTWARE, {})
    versions = []
    if software_category.get("version", None) is None:
        software_names = software_category.get("name", None)
        for software_name in software_names:
            if software_name == "AlphaFold":
                versions.append(f"v{model_meta.get('version', '3.1')}")
            else:
                versions.append(None)
                logger.warning("Different softwares found in JSON file. Consider adding versions for all softwares you have mentioned.")
        cif_data.set_item(CAT_SOFTWARE, "version", versions)

    # 4. Fetch and process external data (UniProt)
    uniprot_details = defaultdict(list)
    with requests.Session() as session:
        api_client = UniprotAPIClient(session)
        for i, chain_info in enumerate(input_metadata.get("chains", [])):
            if uniprot_accession := chain_info.get("uniprot_accession"):
                response_data = api_client.fetch_metadata(uniprot_accession)
                processed_data = process_uniprot_response(response_data)
                for key, value in processed_data.items():
                    uniprot_details[key].append(value)
                # for multiple chains entity_id will always be present.
                # For single chain it is assumed 1.
                entity_id = chain_info.get("entity_id", i + 1)
                uniprot_details["target_entity_id"].append(str(entity_id))

    if uniprot_details:
        cif_data.set_items(CAT_TARGET_REF_DB, dict(uniprot_details))

    # 5. Compute metrics from atomic data
    atom_site_data = cif_data.get_data().get(CAT_ATOM_SITE, {})

    global_plddt = compute_global_plddt(atom_site_data.get(ITEM_B_FACTOR, []))
    if global_plddt >= 0:
        cif_data.set_item(CAT_GLOBAL_QA, "metric_value", [f"{global_plddt:.2f}"])

    local_plddt_metrics = compute_local_plddt_metrics(
        atom_site_data.get(ITEM_LABEL_ASYM_ID, []),
        atom_site_data.get(ITEM_LABEL_COMP_ID, []),
        atom_site_data.get(ITEM_LABEL_SEQ_ID, []),
        atom_site_data.get(ITEM_B_FACTOR, []),
    )
    if local_plddt_metrics:
        cif_data.set_items(CAT_LOCAL_QA, local_plddt_metrics)

    polymer_seq_cats = create_polymer_sequence_categories(atom_site_data)

    if polymer_seq_cats:
        for cat_name, cat_data in polymer_seq_cats.items():
            if cat_data:
                cif_data.set_items(cat_name, cat_data)

    # 6. Write the final mmCIF file
    block_name = Path(output_file).stem
    mmcif_output_file = Path(output_file)
    cif_data.write_to_cif(str(mmcif_output_file), block_name=block_name)

    # 7. Optionally validate the output file
    if validate_dict_path:
        validate_with_gemmi(str(mmcif_output_file), validate_dict_path)

    # 8. Create enriched PDB file
    pdb_output_file = mmcif_output_file.with_suffix(".pdb")

    pdb_editor = PDBFileEditor()
    pdb_editor.load_pdb(pdb_file)
    add_pdb_headers(pdb_editor, cif_data,mmcif_output_file.stem)
    pdb_editor.write_pdb(str(pdb_output_file))

    
