import json
import logging
import shutil
import subprocess
import sys
from collections import defaultdict
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict, List, Optional, TypedDict
import datetime
import gemmi
import jsonschema
import numpy as np
import requests
from afdb_integration_kit.utils.pdbeditor import PDBFileEditor
from afdb_integration_kit.utils.cifstorage import CifDataStorage
from afdb_integration_kit.utils.uniprot import UniprotAPIClient
from afdb_integration_kit.utils.constant import (
    CAT_ATOM_SITE,
    CAT_CHEM_COMP,
    CAT_SOFTWARE,
    CAT_MODEL_LIST,
    CAT_TARGET_REF_DB,
    CAT_GLOBAL_QA,
    CAT_LOCAL_QA,
    CAT_ENTITY_POLY_SEQ,
    CAT_ENTITY_POLY_SEQ_SCHEME,
    CAT_STRUCT_ASYM,
    ITEM_B_FACTOR,
    ITEM_LABEL_ASYM_ID,
    ITEM_AUTH_ASYM_ID,
    ITEM_LABEL_COMP_ID,
    ITEM_AUTH_COMP_ID,
    ITEM_LABEL_SEQ_ID,
    ITEM_AUTH_SEQ_ID,
    ITEM_PDB_INS_CODE,
)

# --- Configuration & Constants ---

# Configure logger
logger = logging.getLogger("afdb_integration_kit")

JSON_SCHEMA_PATH = str(
    files("afdb_integration_kit.modelcif.resources").joinpath("schema.json")
)

class ChainMetadata(TypedDict):
    uniprot_accession: str
    chain_id: str

class InputMetadata(TypedDict):
    metadata: Dict[str, Any]
    categories: Dict[str, Any]
    chains: List[ChainMetadata]

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

    
