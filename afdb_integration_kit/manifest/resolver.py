"""Resolve model IDs to per-chain UniProt accessions.

Supports two model ID patterns without any PDB file I/O:

- **Homodimers / homomultimers**: ``AF_XXXX`` (single AF ID).  The AFCDB
  manifest already stores the chain-to-accession mapping for these.
- **Heterodimers**: ``AF_XXXX_AF_YYYY`` (two AF IDs joined by underscore).
  Each component AF ID is looked up separately and assigned to chains A / B
  in order.

The AFCDB manifest (``afdb_toolkit_manifest_file.csv``) has columns::

    model_entity_id,entity_id,chain_id,uniprot_ac

and uses hyphens (``AF-XXXX``).  This module handles the underscore/hyphen
conversion transparently.
"""

from __future__ import annotations

import csv
import orjson
import logging
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple

import duckdb

try:
    import cudf

    _HAS_CUDF = True
except ImportError:
    _HAS_CUDF = False

logger = logging.getLogger(__name__)

# Matches heterodimer composite IDs like AF_0000_AF_1111
_HETERODIMER_PATTERN = re.compile(r"^(AF_\d+)_(AF_\d+)$")

# Matches a single AF ID: AF_0000 (underscore) or AF-0000 (canonical hyphen)
_SINGLE_AF_PATTERN = re.compile(r"^(AF[_-]\d+)$")


# ---------------------------------------------------------------------------
# Model ID classification
# ---------------------------------------------------------------------------

def classify_model_ids(
    model_ids: List[str],
) -> Tuple[Dict[str, List[str]], List[str]]:
    """Classify model IDs into homodimer (single AF) and heterodimer (composite).

    Returns:
        (classified, unrecognised) where *classified* maps each model_id to a
        list of its component AF IDs (length 1 for homodimers, 2 for
        heterodimers), and *unrecognised* lists IDs matching neither pattern.
    """
    classified: Dict[str, List[str]] = {}
    unrecognised: List[str] = []

    for mid in model_ids:
        het_match = _HETERODIMER_PATTERN.match(mid)
        if het_match:
            classified[mid] = [het_match.group(1), het_match.group(2)]
            continue

        single_match = _SINGLE_AF_PATTERN.match(mid)
        if single_match:
            classified[mid] = [single_match.group(1)]
            continue

        unrecognised.append(mid)

    return classified, unrecognised


def collect_unique_af_ids(
    classified: Dict[str, List[str]],
) -> Set[str]:
    """Collect all unique AF IDs across all classified model IDs."""
    af_ids: Set[str] = set()
    for components in classified.values():
        af_ids.update(components)
    return af_ids


# ---------------------------------------------------------------------------
# AFCDB manifest streaming
# ---------------------------------------------------------------------------

def _af_underscore_to_hyphen(af_id: str) -> str:
    """``AF_XXXX`` -> ``AF-XXXX`` (manifest uses hyphens)."""
    return af_id.replace("_", "-", 1)


def _resolve_via_cudf(
    manifest_path: Path,
    hyphen_to_underscore: Dict[str, str],
    target_ids: list,
) -> Dict[str, Set[str]]:
    """GPU-accelerated manifest lookup using cuDF."""
    df = cudf.read_csv(
        str(manifest_path),
        usecols=["model_entity_id", "uniprot_ac"],
    )
    target_series = cudf.Series(target_ids)
    matched = df[df["model_entity_id"].isin(target_series)]
    matched_pd = matched.to_pandas()

    mapping: Dict[str, Set[str]] = {}
    for mid, ac in zip(matched_pd["model_entity_id"], matched_pd["uniprot_ac"]):
        mapping.setdefault(hyphen_to_underscore[mid], set()).add(ac)
    return mapping


def _resolve_via_duckdb(
    manifest_path: Path,
    hyphen_to_underscore: Dict[str, str],
    target_ids: list,
) -> Dict[str, Set[str]]:
    """Fast CPU manifest lookup using DuckDB's parallel CSV reader."""
    con = duckdb.connect()
    con.execute("PRAGMA memory_limit = '8GB'")
    result = con.execute(
        """
        SELECT DISTINCT model_entity_id, uniprot_ac
        FROM read_csv(?, delim=',', header=true, parallel=true)
        WHERE model_entity_id IN (SELECT unnest(?::VARCHAR[]))
        """,
        [str(manifest_path), target_ids],
    ).fetchall()
    con.close()

    mapping: Dict[str, Set[str]] = {}
    for mid, ac in result:
        mapping.setdefault(hyphen_to_underscore[mid], set()).add(ac)
    return mapping


def _resolve_via_csv(
    manifest_path: Path,
    hyphen_to_underscore: Dict[str, str],
) -> Dict[str, Set[str]]:
    """Legacy single-threaded csv.DictReader streaming (no extra deps)."""
    target_hyphens = set(hyphen_to_underscore.keys())
    mapping: Dict[str, Set[str]] = {}

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            mid = row["model_entity_id"]
            if mid not in target_hyphens:
                continue
            mapping.setdefault(
                hyphen_to_underscore[mid], set()
            ).add(row["uniprot_ac"])

    return mapping


def _batch_query_sequence_lengths(
    db_path: Path,
    accessions: Set[str],
) -> Dict[str, int]:
    """Query DuckDB for sequence lengths of all accessions in one batch."""
    if not accessions:
        return {}
    con = duckdb.connect(str(db_path), read_only=True)
    con.execute("PRAGMA memory_limit = '512MB'")
    rows = con.execute(
        """
        SELECT primary_ac, LENGTH(sequence)
        FROM entry
        WHERE primary_ac IN (SELECT unnest(?::VARCHAR[]))
        """,
        [list(accessions)],
    ).fetchall()
    con.close()
    return {ac: length for ac, length in rows}


def _read_plddt_length(meta_path: Path) -> int | None:
    """Read pLDDT array length from a ColabFold meta JSON using orjson."""
    with open(meta_path, "rb") as fh:
        meta = orjson.loads(fh.read())
    plddt = meta.get("plddt")
    if plddt is None:
        return None
    return len(plddt)


def _deduplicate_accessions(
    raw_mapping: Dict[str, Set[str]],
    uniprot_db_path: Path | None = None,
    input_dir: Path | None = None,
) -> Tuple[Dict[str, str], List[str]]:
    """Resolve models with multiple candidate UniProt accessions.

    For models with a single accession, passes through unchanged.
    For models with multiple accessions, validates against DuckDB sequence
    lengths and the ColabFold pLDDT to pick the correct one.

    Returns:
        ``(clean_mapping, failed_ids)`` where *clean_mapping* maps each AF ID
        to exactly one accession, and *failed_ids* lists AF IDs that could
        not be disambiguated.
    """
    clean: Dict[str, str] = {}
    failed: List[str] = []

    ambiguous: Dict[str, Set[str]] = {}
    for af_id, accs in raw_mapping.items():
        if len(accs) == 1:
            clean[af_id] = next(iter(accs))
        else:
            ambiguous[af_id] = accs

    if not ambiguous:
        return clean, failed

    logger.warning(
        "%d AF IDs have multiple accessions in manifest; "
        "attempting deduplication",
        len(ambiguous),
    )

    if uniprot_db_path is None or input_dir is None:
        logger.warning(
            "No uniprot_db_path or input_dir provided; using first accession "
            "alphabetically for %d ambiguous models",
            len(ambiguous),
        )
        for af_id, accs in ambiguous.items():
            picked = sorted(accs)[0]
            logger.warning(
                "  %s: picked %s from %s (alphabetical fallback)",
                _af_underscore_to_hyphen(af_id), picked, sorted(accs),
            )
            clean[af_id] = picked
        return clean, failed

    all_candidates: Set[str] = set()
    for accs in ambiguous.values():
        all_candidates.update(accs)
    ac_lengths = _batch_query_sequence_lengths(uniprot_db_path, all_candidates)

    resolved_count = 0
    input_path = Path(input_dir) if not isinstance(input_dir, Path) else input_dir

    for af_id, accs in sorted(ambiguous.items()):
        model_id = _af_underscore_to_hyphen(af_id)
        meta_path = input_path / f"{model_id}-meta_v1.json"

        if not meta_path.exists():
            logger.warning(
                "%s: meta file not found at %s; failing model",
                model_id, meta_path,
            )
            failed.append(af_id)
            continue

        plddt_length = _read_plddt_length(meta_path)
        if plddt_length is None:
            logger.warning(
                "%s: no pLDDT in meta file %s; failing model",
                model_id, meta_path,
            )
            failed.append(af_id)
            continue
        monomer_length = plddt_length // 2

        matching = [
            ac for ac in accs if ac_lengths.get(ac) == monomer_length
        ]

        if len(matching) == 1:
            picked = matching[0]
            rejected = sorted(accs - {picked})
            logger.warning(
                "%s: picked %s (len=%d) over %s "
                "(pLDDT=%d, monomer=%d)",
                model_id, picked, ac_lengths[picked],
                rejected, plddt_length, monomer_length,
            )
            clean[af_id] = picked
            resolved_count += 1
        elif len(matching) > 1:
            logger.warning(
                "%s: AMBIGUOUS - multiple accessions match monomer "
                "length %d: %s; failing model",
                model_id, monomer_length, sorted(matching),
            )
            failed.append(af_id)
        else:
            logger.warning(
                "%s: NO accession matches monomer length %d "
                "(candidates: %s); failing model",
                model_id, monomer_length,
                {ac: ac_lengths.get(ac) for ac in sorted(accs)},
            )
            failed.append(af_id)

    logger.info(
        "Deduplication complete: %d resolved, %d failed "
        "out of %d ambiguous",
        resolved_count, len(failed), len(ambiguous),
    )
    return clean, failed


def stream_afcdb_manifest(
    manifest_path: Path,
    af_ids: Set[str],
    *,
    backend: str = "auto",
    uniprot_db_path: Path | None = None,
    input_dir: Path | None = None,
) -> Tuple[Dict[str, str], List[str]]:
    """Look up AF IDs in the AFCDB manifest CSV.

    Backend selection (``backend`` kwarg):

    * ``"auto"`` (default): cuDF if available, else DuckDB.
    * ``"cudf"``: Force GPU via cuDF (requires ``cudf`` package).
    * ``"duckdb"``: Force DuckDB's parallel CSV reader (CPU).
    * ``"csv"``: Legacy single-threaded ``csv.DictReader`` streaming.

    When *uniprot_db_path* and *input_dir* are provided, models that map to
    multiple UniProt accessions in the manifest are disambiguated by
    comparing DuckDB sequence lengths against the ColabFold pLDDT.

    Returns:
        ``({af_id: uniprot_ac}, [failed_af_ids])``
    """
    hyphen_to_underscore = {_af_underscore_to_hyphen(af): af for af in af_ids}
    for af in af_ids:
        hyphen_to_underscore[af] = af
    target_ids = list(hyphen_to_underscore.keys())

    if backend == "auto":
        backend = "cudf" if _HAS_CUDF else "duckdb"

    logger.info(
        "Resolving %d AF IDs from %s (backend: %s)...",
        len(af_ids), manifest_path, backend,
    )

    if backend == "cudf":
        if not _HAS_CUDF:
            raise ImportError("cudf is not installed; use backend='duckdb' or 'csv'")
        raw_mapping = _resolve_via_cudf(manifest_path, hyphen_to_underscore, target_ids)
    elif backend == "duckdb":
        raw_mapping = _resolve_via_duckdb(manifest_path, hyphen_to_underscore, target_ids)
    elif backend == "csv":
        raw_mapping = _resolve_via_csv(manifest_path, hyphen_to_underscore)
    else:
        raise ValueError(f"Unknown backend: {backend!r}")

    missing = af_ids - set(raw_mapping.keys())
    if missing:
        logger.warning(
            "%d AF IDs not found in manifest: %s",
            len(missing),
            sorted(_af_underscore_to_hyphen(af) for af in missing)[:10],
        )

    mapping, dedup_failed = _deduplicate_accessions(
        raw_mapping,
        uniprot_db_path=uniprot_db_path,
        input_dir=input_dir,
    )

    logger.info(
        "Resolved %d / %d AF IDs to UniProt accessions",
        len(mapping), len(af_ids),
    )
    return mapping, dedup_failed


# ---------------------------------------------------------------------------
# ColabFold manifest building
# ---------------------------------------------------------------------------

def build_colabfold_manifest(
    model_ids: List[str],
    af_to_uniprot: Dict[str, str],
    classified: Dict[str, List[str]],
) -> Tuple[List[Dict[str, str]], List[str]]:
    """Build ColabFold manifest rows from classified model IDs.

    For homodimers (single AF ID), both chains get the same accession.
    For heterodimers (two AF IDs), chain A gets the first, chain B the second.

    Returns:
        ``(manifest_rows, skipped_model_ids)``
    """
    chain_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    rows: List[Dict[str, str]] = []
    skipped: List[str] = []

    for mid in model_ids:
        components = classified.get(mid)
        if components is None:
            skipped.append(mid)
            continue

        accessions = [af_to_uniprot.get(af) for af in components]
        if any(ac is None for ac in accessions):
            skipped.append(mid)
            continue

        if len(components) == 1:
            canonical_mid = _af_underscore_to_hyphen(mid)
        else:
            canonical_mid = "_".join(_af_underscore_to_hyphen(c) for c in components)
        if len(components) == 1:
            # Homodimer: replicate accession for 2 chains
            for chain_idx in range(2):
                rows.append({
                    "model_entity_id": canonical_mid,
                    "chain_id": chain_labels[chain_idx],
                    "uniprot_ac": accessions[0],
                })
        else:
            # Heterodimer: one chain per component
            for chain_idx, ac in enumerate(accessions):
                if chain_idx >= len(chain_labels):
                    logger.warning("Model %s has more chains than labels", mid)
                    break
                rows.append({
                    "model_entity_id": canonical_mid,
                    "chain_id": chain_labels[chain_idx],
                    "uniprot_ac": ac,
                })

    return rows, skipped


def resolve_and_build_manifest(
    model_ids: List[str],
    afcdb_manifest_path: Path,
    *,
    uniprot_db_path: Path | None = None,
    input_dir: Path | None = None,
) -> Tuple[List[Dict[str, str]], List[str], Dict[str, str]]:
    """End-to-end: classify IDs, stream AFCDB manifest, build ColabFold manifest.

    When *uniprot_db_path* and *input_dir* are provided, models that map to
    multiple UniProt accessions in the manifest are disambiguated by
    comparing DuckDB sequence lengths against the ColabFold pLDDT.

    Returns:
        ``(manifest_rows, skipped_model_ids, af_to_uniprot_mapping)``
    """
    classified, unrecognised = classify_model_ids(model_ids)
    if unrecognised:
        logger.warning(
            "%d model IDs do not match any known pattern: %s",
            len(unrecognised), unrecognised[:5],
        )

    af_ids = collect_unique_af_ids(classified)
    af_to_uniprot, dedup_failed = stream_afcdb_manifest(
        afcdb_manifest_path,
        af_ids,
        uniprot_db_path=uniprot_db_path,
        input_dir=input_dir,
    )

    rows, skipped = build_colabfold_manifest(model_ids, af_to_uniprot, classified)
    all_skipped = unrecognised + dedup_failed + skipped

    return rows, all_skipped, af_to_uniprot
