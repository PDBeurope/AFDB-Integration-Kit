#!/usr/bin/env python3
"""
Set up mock test data for local Nextflow workflow testing.

This creates:
- A mock DuckDB database with UniProt entries
- Sample config files (dataset_config.json, provider.json, etc.)
- Sample PDB files in POSIX layout
- All required validation files per model:
  - *-model_v1.pdb (structure)
  - *-meta_v1.json (ColabFold metadata)
  - *-confidence_v1.json (pLDDT scores)
  - *-predicted_aligned_error_v1.json (PAE matrix)
  - *-msa_v1.a3m (multiple sequence alignment)
- sequences.fasta at dataset root
- Manifest files

Usage:
    python tests/fixtures/setup_mock_data.py --output-dir /tmp/mock_afdb_data
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False


# ============================================================================
# Constants
# ============================================================================

SAMPLE_MODEL_IDS = [
    "AF-0000000000000001",
    "AF-0000000000000002",
    "AF-0000000000000003",
    "AF-0000000000000004",
    "AF-0000000000000005",
]

SAMPLE_UNIPROT_ENTRIES = [
    {
        "primary_ac": "P12345",
        "entry_name": "TEST1_HUMAN",
        "organism": "Homo sapiens",
        "taxid": 9606,
        "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALVHQQIFEVTEALVKYGNPGGIPIYTVMNGKEDLSSLIIRSRSSTVFQQIFEQLQPILTHLNQQQQQQQQQQQQLQQQQILQVQDQQAQQLQVQQEQAQQVQQVQQQPLQVQVQQQQQQQQQQQQQQQLQVQ",
        "protein_full_names": ["Test protein 1"],
        "protein_short_names": ["TP1"],
        "gene_names": ["TEST1"],
        "gene_synonyms": ["TP1", "TESTP1"],
        "organism_common_names": ["Human"],
        "organism_synonyms": [],
        "reviewed": True,
        "is_uniprot_reference_proteome": True,
        "sequence_version_date": "2024-01-01",
    },
    {
        "primary_ac": "P67890",
        "entry_name": "TEST2_MOUSE",
        "organism": "Mus musculus",
        "taxid": 10090,
        "sequence": "MKFLILLFNILCLFPVLAADNHGVGPQGASNAFSNFISGVKGIDFKDGTLTLKDLKGKF" * 2,
        "protein_full_names": ["Test protein 2"],
        "protein_short_names": ["TP2"],
        "gene_names": ["Test2"],
        "gene_synonyms": [],
        "organism_common_names": ["Mouse"],
        "organism_synonyms": ["House mouse"],
        "reviewed": True,
        "is_uniprot_reference_proteome": True,
        "sequence_version_date": "2024-02-01",
    },
    {
        "primary_ac": "Q11111",
        "entry_name": "TEST3_YEAST",
        "organism": "Saccharomyces cerevisiae",
        "taxid": 559292,
        "sequence": "MSSQVPVTVSPLRASTHGLHHHHHHHSPALTPSATSTGSDSAASGSG" * 3,
        "protein_full_names": ["Test protein 3"],
        "protein_short_names": [],
        "gene_names": ["TST3"],
        "gene_synonyms": ["YAL001C"],
        "organism_common_names": ["Baker's yeast"],
        "organism_synonyms": [],
        "reviewed": False,
        "is_uniprot_reference_proteome": False,
        "sequence_version_date": "2024-03-01",
    },
]


# ============================================================================
# Data Generators
# ============================================================================

def generate_plddt_scores(seq_length: int) -> List[float]:
    """Generate realistic pLDDT scores."""
    return [round(random.uniform(40, 95), 2) for _ in range(seq_length)]


def generate_pae_matrix(seq_length: int) -> List[List[float]]:
    """Generate a realistic PAE matrix."""
    matrix = []
    for i in range(seq_length):
        row = []
        for j in range(seq_length):
            # Diagonal is always ~0, distance from diagonal correlates with PAE
            dist = abs(i - j)
            base = min(dist * 0.5, 25)
            noise = random.uniform(-2, 2)
            value = max(0.25, min(31.75, base + noise))
            row.append(round(value, 2))
        matrix.append(row)
    return matrix


def generate_pdb_file(model_id: str, seq_length: int) -> str:
    """Generate a minimal valid PDB file."""
    lines = [
        f"HEADER    PREDICTION                              {datetime.now().strftime('%d-%b-%y').upper()}   {model_id}",
        f"TITLE     ALPHAFOLD PREDICTION FOR {model_id}",
    ]

    atom_num = 1
    for res_num in range(1, seq_length + 1):
        # Generate coordinates along a helix-like path
        import math
        theta = res_num * 0.3
        x = math.cos(theta) * 5 + res_num * 0.1
        y = math.sin(theta) * 5
        z = res_num * 1.5

        # CA atom
        line = f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{random.uniform(40, 95):6.2f}           C"
        lines.append(line)
        atom_num += 1

    lines.append("END")
    return "\n".join(lines)


def generate_meta_json(model_id: str, seq_length: int, uniprot_ac: str) -> Dict[str, Any]:
    """Generate the ColabFold-style meta JSON.

    ColabFold converter expects keys: plddt, pae, max_pae
    """
    plddt = generate_plddt_scores(seq_length)
    pae = generate_pae_matrix(seq_length)
    max_pae = max(max(row) for row in pae)

    return {
        "plddt": plddt,
        "pae": pae,  # ColabFold uses 'pae', not 'predicted_aligned_error'
        "max_pae": max_pae,
        "ptm": round(random.uniform(0.5, 0.95), 3),
        "iptm": round(random.uniform(0.4, 0.9), 3),  # For multimers
    }


def generate_confidence_json(model_id: str, seq_length: int) -> Dict[str, Any]:
    """Generate AFDB confidence JSON format.

    Category codes per AFDB spec:
    - V (Very High): >90
    - H (High): 70-90
    - M (Medium): 50-70
    - L (Low): 30-50
    - D (Disordered): <30
    """
    plddt = generate_plddt_scores(seq_length)

    def score_to_category(score: float) -> str:
        if score >= 90:
            return "V"
        elif score >= 70:
            return "H"
        elif score >= 50:
            return "M"
        elif score >= 30:
            return "L"
        else:
            return "D"

    return {
        "residueNumber": list(range(1, seq_length + 1)),
        "confidenceScore": plddt,
        "confidenceCategory": [score_to_category(s) for s in plddt],
    }


def generate_pae_json(seq_length: int) -> List[Dict[str, Any]]:
    """Generate AFDB PAE JSON format."""
    pae = generate_pae_matrix(seq_length)
    return [{
        "predicted_aligned_error": pae,
        "max_predicted_aligned_error": max(max(row) for row in pae),
    }]


def generate_msa_a3m(model_id: str, sequence: str) -> str:
    """Generate a minimal MSA in A3M format.

    A3M is a compressed multiple sequence alignment format.
    We generate a simple alignment with the query sequence and a few homologs.
    """
    lines = [f">{model_id}"]
    lines.append(sequence)

    # Add a few "homologous" sequences with some mutations
    for i in range(3):
        lines.append(f">homolog_{i+1}")
        # Create a mutated version of the sequence
        mutated = list(sequence)
        for j in range(0, len(mutated), 10 + i):
            if j < len(mutated):
                mutated[j] = random.choice("ACDEFGHIKLMNPQRSTVWY")
        lines.append("".join(mutated))

    return "\n".join(lines)


# ============================================================================
# DuckDB Setup
# ============================================================================

def create_mock_duckdb(db_path: Path, entries: List[Dict[str, Any]]) -> None:
    """Create a mock DuckDB database with UniProt entries."""
    if not HAS_DUCKDB:
        raise ImportError("duckdb is required. Install with: pip install duckdb")

    con = duckdb.connect(str(db_path))

    # Create the entry table matching the real schema
    con.execute("""
        CREATE TABLE entry (
            primary_ac VARCHAR PRIMARY KEY,
            entry_name VARCHAR,
            organism VARCHAR,
            taxid INTEGER,
            sequence VARCHAR,
            protein_full_names VARCHAR[],
            protein_short_names VARCHAR[],
            gene_names VARCHAR[],
            gene_synonyms VARCHAR[],
            gene_ordered_locus_names VARCHAR[],
            gene_orf_names VARCHAR[],
            organism_common_names VARCHAR[],
            organism_synonyms VARCHAR[],
            reviewed BOOLEAN,
            is_uniprot_reference_proteome BOOLEAN,
            sequence_version_date DATE
        )
    """)

    # Insert entries
    for entry in entries:
        con.execute("""
            INSERT INTO entry VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, [
            entry["primary_ac"],
            entry["entry_name"],
            entry["organism"],
            entry["taxid"],
            entry["sequence"],
            entry.get("protein_full_names", []),
            entry.get("protein_short_names", []),
            entry.get("gene_names", []),
            entry.get("gene_synonyms", []),
            entry.get("gene_ordered_locus_names", []),
            entry.get("gene_orf_names", []),
            entry.get("organism_common_names", []),
            entry.get("organism_synonyms", []),
            entry.get("reviewed", False),
            entry.get("is_uniprot_reference_proteome", False),
            entry.get("sequence_version_date"),
        ])

    con.close()
    print(f"Created DuckDB with {len(entries)} entries at {db_path}")


# ============================================================================
# File Structure Setup
# ============================================================================

def setup_posix_layout(
    output_dir: Path,
    model_ids: List[str],
    entries: List[Dict[str, Any]],
) -> List[tuple]:
    """Set up the POSIX shard layout with model files.

    Creates all required files for AFDB validation:
    - PDB model file
    - Meta JSON (ColabFold format)
    - pLDDT/confidence JSON
    - PAE JSON
    - MSA A3M file

    Returns list of (model_id, sequence) tuples for sequences.fasta generation.
    """

    input_dir = output_dir / "input"
    sequences_for_fasta = []

    for i, model_id in enumerate(model_ids):
        entry = entries[i % len(entries)]
        sequence = entry["sequence"]
        seq_length = len(sequence)

        # Create shard path (AF-XXXX/XXXX/XXXX/XXXX)
        digits = model_id.replace("AF-", "")[:16].ljust(16, "0")
        shard_parts = [digits[j:j+4] for j in range(0, 16, 4)]
        shard_path = input_dir.joinpath(*shard_parts)
        shard_path.mkdir(parents=True, exist_ok=True)

        # Create meta JSON (ColabFold format with plddt/pae embedded)
        meta = generate_meta_json(model_id, seq_length, entry["primary_ac"])
        meta_path = shard_path / f"{model_id}-meta_v1.json"
        meta_path.write_text(json.dumps(meta, indent=2))

        # Create PDB file
        pdb_content = generate_pdb_file(model_id, seq_length)
        pdb_path = shard_path / f"{model_id}-model_v1.pdb"
        pdb_path.write_text(pdb_content)

        # Create pLDDT/confidence JSON (AFDB format)
        confidence = generate_confidence_json(model_id, seq_length)
        confidence_path = shard_path / f"{model_id}-confidence_v1.json"
        confidence_path.write_text(json.dumps(confidence, indent=2))

        # Create PAE JSON (AFDB format)
        pae = generate_pae_json(seq_length)
        pae_path = shard_path / f"{model_id}-predicted_aligned_error_v1.json"
        pae_path.write_text(json.dumps(pae, indent=2))

        # Create MSA A3M file
        msa_content = generate_msa_a3m(model_id, sequence)
        msa_path = shard_path / f"{model_id}-msa_v1.a3m"
        msa_path.write_text(msa_content)

        # Store for sequences.fasta
        sequences_for_fasta.append((model_id, sequence))

        print(f"Created files for {model_id} at {shard_path}")

    return sequences_for_fasta


def create_sequences_fasta(output_dir: Path, sequences: List[tuple]) -> None:
    """Create sequences.fasta file at dataset root.

    Args:
        output_dir: Dataset root directory
        sequences: List of (model_id, sequence) tuples
    """
    fasta_path = output_dir / "sequences.fasta"
    lines = []
    for model_id, sequence in sequences:
        lines.append(f">{model_id}")
        # Wrap sequence at 80 characters
        for i in range(0, len(sequence), 80):
            lines.append(sequence[i:i+80])

    fasta_path.write_text("\n".join(lines) + "\n")
    print(f"Created sequences.fasta with {len(sequences)} sequences at {fasta_path}")


def create_config_files(output_dir: Path, model_ids: List[str], entries: List[Dict[str, Any]]) -> None:
    """Create configuration files for the workflow."""
    config_dir = output_dir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    # Dataset config
    dataset_config = {
        "providerId": "test-provider",
        "toolUsed": "AlphaFold",
        "latestVersion": 1,
        "allVersions": [1],
        "entityType": "protein",
        "modelCreatedDate": "2024-01-01T00:00:00Z",
        "uniqueIdTemplate": "{model_entity_id}",
        "versionTag": "v1",
    }
    (config_dir / "dataset_config.json").write_text(json.dumps(dataset_config, indent=2))

    # Provider JSON
    provider = {
        "providerId": "test-provider",
        "providerName": "Test Provider",
        "providerUrl": "https://example.com",
        "copyrights": ["Copyright 2024 Test Provider. All rights reserved."],
    }
    (config_dir / "provider.json").write_text(json.dumps(provider, indent=2))

    # AF mapping file (TSV) - single column, NO header (workflow reads all rows as data)
    with (config_dir / "af_mapping.tsv").open("w") as f:
        for model_id in model_ids:
            f.write(f"{model_id}\n")

    # ColabFold manifest (CSV)
    colabfold_rows = []
    for i, model_id in enumerate(model_ids):
        entry = entries[i % len(entries)]
        colabfold_rows.append({
            "model_entity_id": model_id,
            "uniprot_ac": entry["primary_ac"],
            "chain_id": "A",
        })

    with (config_dir / "colabfold_manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["model_entity_id", "uniprot_ac", "chain_id"])
        writer.writeheader()
        writer.writerows(colabfold_rows)

    # ModelCIF template
    modelcif_template = {
        "metadata": {
            "title": "Test Model",
        },
        "categories": {
            "_entry": {"id": ["PLACEHOLDER"]},
            "_database_2": {"database_code": ["PLACEHOLDER"]},
            "_pdbx_database_status": {"entry_id": ["PLACEHOLDER"]},
        },
    }
    (config_dir / "modelcif_template.json").write_text(json.dumps(modelcif_template, indent=2))

    print(f"Created config files in {config_dir}")


def create_nextflow_config(output_dir: Path) -> None:
    """Create a minimal Nextflow config for local testing."""
    config_content = """// Local test configuration for AFDB Integration Kit

process {
    executor = 'local'
    cpus = 2
    memory = '4 GB'
}

docker.enabled = false
singularity.enabled = false
"""
    (output_dir / "config" / "nextflow_local.config").write_text(config_content)
    print(f"Created Nextflow config at {output_dir / 'config' / 'nextflow_local.config'}")


# ============================================================================
# Main
# ============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Set up mock test data for Nextflow workflow")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/mock_afdb_data"),
        help="Output directory for mock data",
    )
    parser.add_argument(
        "--num-models",
        type=int,
        default=5,
        help="Number of test models to create",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing output directory before creating",
    )
    args = parser.parse_args()

    output_dir = args.output_dir

    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
        print(f"Cleaned {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate model IDs
    model_ids = [f"AF-{str(i).zfill(16)}" for i in range(1, args.num_models + 1)]

    # Create DuckDB
    db_path = output_dir / "uniprot_test.duckdb"
    create_mock_duckdb(db_path, SAMPLE_UNIPROT_ENTRIES)

    # Create POSIX layout with all required files (PDB, meta, pLDDT, PAE, MSA)
    sequences = setup_posix_layout(output_dir, model_ids, SAMPLE_UNIPROT_ENTRIES)

    # Create sequences.fasta at dataset root
    create_sequences_fasta(output_dir, sequences)

    # Create config files
    create_config_files(output_dir, model_ids, SAMPLE_UNIPROT_ENTRIES)

    # Create Nextflow config
    create_nextflow_config(output_dir)

    print("\n" + "=" * 60)
    print("Mock data setup complete!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print(f"DuckDB database:  {db_path}")
    print(f"Model count:      {len(model_ids)}")
    print("\nTo run the workflow:")
    print(f"""
nextflow run workflow/end_to_end_with_validation_multibatch.nf \\
  -c {output_dir}/config/nextflow_local.config \\
  --posix_root {output_dir}/input \\
  --mapping_file {output_dir}/config/af_mapping.tsv \\
  --colabfold_manifest {output_dir}/config/colabfold_manifest.csv \\
  --dataset_dir {output_dir} \\
  --dataset_config {output_dir}/config/dataset_config.json \\
  --provider_json {output_dir}/config/provider.json \\
  --archive_root {output_dir}/input \\
  --converter_duckdb {output_dir}/uniprot_test.duckdb \\
  --uniprot_db {output_dir}/uniprot_test.duckdb \\
  --toolkit_cmd 'python main.py' \\
  --sample_size 0 \\
  -work-dir {output_dir}/work \\
  -with-report {output_dir}/report.html \\
  -with-timeline {output_dir}/timeline.html
""")


if __name__ == "__main__":
    main()
