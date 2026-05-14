"""UniProt REST API client for fetching entries and building DuckDB databases.

Provides batch fetching of UniProt entries via the REST API, parsing into the
standard entry schema used by the pipeline, and helpers to write Parquet and
build indexed DuckDB databases.  Intended as a dev/test convenience when a
pre-built DuckDB is not available.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
DEFAULT_BATCH_SIZE = 25
DEFAULT_RELEASE = "2025_01"

ENTRY_SCHEMA = pa.schema([
    ("primary_ac", pa.string()),
    ("entry_name", pa.string()),
    ("reviewed", pa.bool_()),
    ("protein_full_names", pa.list_(pa.string())),
    ("protein_short_names", pa.list_(pa.string())),
    ("gene_names", pa.string()),
    ("gene_synonyms", pa.list_(pa.string())),
    ("gene_ordered_locus_names", pa.list_(pa.string())),
    ("gene_orf_names", pa.list_(pa.string())),
    ("organism", pa.string()),
    ("organism_common_names", pa.list_(pa.string())),
    ("organism_synonyms", pa.list_(pa.string())),
    ("taxid", pa.int64()),
    ("length", pa.int32()),
    ("sequence_version_date", pa.string()),
    ("is_uniprot_reference_proteome", pa.bool_()),
    ("md5", pa.string()),
    ("sequence", pa.string()),
    ("release", pa.string()),
    ("is_isoform", pa.bool_()),
])


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------

def fetch_entries(
    accessions: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    release: str = DEFAULT_RELEASE,
) -> List[Dict[str, Any]]:
    """Fetch UniProt entries via REST API in batches.

    Args:
        accessions: UniProt accession IDs to fetch.
        batch_size: Accessions per API request (max ~100).
        release: Release tag to store with each entry.

    Returns:
        List of parsed entry dictionaries matching ``ENTRY_SCHEMA``.
    """
    entries: List[Dict[str, Any]] = []

    for i in tqdm(range(0, len(accessions), batch_size), desc="Fetching UniProt"):
        batch = accessions[i : i + batch_size]
        query = " OR ".join(f"accession:{ac}" for ac in batch)
        params = {
            "query": query,
            "format": "json",
            "fields": (
                "accession,id,protein_name,gene_names,organism_name,"
                "organism_id,length,sequence,reviewed,xref_proteomes"
            ),
            "size": len(batch),
        }
        try:
            response = requests.get(
                f"{UNIPROT_API_BASE}/search", params=params, timeout=30,
            )
            response.raise_for_status()
            data = response.json()
            for result in data.get("results", []):
                entry = parse_entry(result, release=release)
                if entry is not None:
                    entries.append(entry)
        except requests.RequestException as exc:
            logger.warning("Failed to fetch batch starting at %d: %s", i, exc)

    return entries


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_entry(
    result: Dict[str, Any],
    release: str = DEFAULT_RELEASE,
) -> Optional[Dict[str, Any]]:
    """Parse a single UniProt API JSON result into the DuckDB entry schema."""
    try:
        primary_ac = result.get("primaryAccession", "")
        entry_name = result.get("uniProtkbId", "")
        reviewed = result.get("entryType", "").startswith("UniProtKB reviewed")

        # -- Protein names --
        protein_desc = result.get("proteinDescription", {})
        recommended = protein_desc.get("recommendedName", {})
        protein_full_names: List[str] = []
        protein_short_names: List[str] = []

        if recommended:
            full = recommended.get("fullName", {}).get("value")
            if full:
                protein_full_names.append(full)
            for s in recommended.get("shortNames", []):
                val = s.get("value")
                if val:
                    protein_short_names.append(val)

        # Submission names (unreviewed entries without recommended names)
        if not protein_full_names:
            for sub in protein_desc.get("submissionNames", []):
                full = sub.get("fullName", {}).get("value")
                if full:
                    protein_full_names.append(full)

        for alt in protein_desc.get("alternativeNames", []):
            full = alt.get("fullName", {}).get("value")
            if full:
                protein_full_names.append(full)
            for s in alt.get("shortNames", []):
                val = s.get("value")
                if val:
                    protein_short_names.append(val)

        # -- Gene names --
        genes = result.get("genes", [])
        gene_names = ""
        gene_synonyms: List[str] = []
        gene_ordered_locus: List[str] = []
        gene_orf: List[str] = []
        if genes:
            g = genes[0]
            if "geneName" in g:
                gene_names = g["geneName"].get("value", "")
            gene_synonyms = [x.get("value", "") for x in g.get("synonyms", [])]
            gene_ordered_locus = [
                x.get("value", "") for x in g.get("orderedLocusNames", [])
            ]
            gene_orf = [x.get("value", "") for x in g.get("orfNames", [])]

        # -- Organism --
        org = result.get("organism", {})
        organism = org.get("scientificName", "")
        organism_common = [org["commonName"]] if org.get("commonName") else []
        organism_syn = org.get("synonyms", [])
        taxid = org.get("taxonId", 0)

        # -- Sequence --
        seq_data = result.get("sequence", {})
        sequence = seq_data.get("value", "")
        length = seq_data.get("length", len(sequence))
        md5 = hashlib.md5(sequence.encode()).hexdigest() if sequence else ""

        # -- Reference proteome --
        xrefs = result.get("uniProtKBCrossReferences", [])
        is_ref_proteome = any(
            x.get("database") == "Proteomes" and "UP" in x.get("id", "")
            for x in xrefs
        )

        return {
            "primary_ac": primary_ac,
            "entry_name": entry_name,
            "reviewed": reviewed,
            "protein_full_names": protein_full_names or [""],
            "protein_short_names": protein_short_names or [""],
            "gene_names": gene_names,
            "gene_synonyms": gene_synonyms or [""],
            "gene_ordered_locus_names": gene_ordered_locus or [""],
            "gene_orf_names": gene_orf or [""],
            "organism": organism,
            "organism_common_names": organism_common or [""],
            "organism_synonyms": organism_syn or [""],
            "taxid": taxid,
            "length": length,
            "sequence_version_date": "",
            "is_uniprot_reference_proteome": is_ref_proteome,
            "md5": md5,
            "sequence": sequence,
            "release": release,
            "is_isoform": "-" in primary_ac,
        }
    except Exception as exc:
        logger.warning("Failed to parse entry: %s", exc)
        return None


# ---------------------------------------------------------------------------
# Parquet / DuckDB
# ---------------------------------------------------------------------------

def entries_to_parquet(
    entries: List[Dict[str, Any]],
    output_path: Path,
) -> None:
    """Write UniProt entries to a Parquet file matching ``ENTRY_SCHEMA``."""
    columns = {field.name: [] for field in ENTRY_SCHEMA}
    for entry in entries:
        for field in ENTRY_SCHEMA:
            columns[field.name].append(entry.get(field.name))

    arrays = [pa.array(columns[f.name], type=f.type) for f in ENTRY_SCHEMA]
    table = pa.Table.from_arrays(arrays, schema=ENTRY_SCHEMA)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, output_path)
    logger.info("Wrote %d entries to %s", len(entries), output_path)


def build_duckdb(parquet_path: Path, db_path: Path) -> None:
    """Build an indexed DuckDB database from an ``entry.parquet`` file."""
    if db_path.exists():
        db_path.unlink()

    con = duckdb.connect(str(db_path))
    try:
        con.execute(
            f"CREATE TABLE entry AS SELECT * FROM read_parquet('{parquet_path}')"
        )
        con.execute("CREATE INDEX idx_entry_primary ON entry(primary_ac)")
        con.execute("CREATE INDEX idx_entry_release ON entry(release)")
        con.commit()
        count = con.execute("SELECT COUNT(*) FROM entry").fetchone()[0]
        logger.info("Built DuckDB with %d entries at %s", count, db_path)
    finally:
        con.close()
