from __future__ import annotations

import argparse
import csv
import logging
import os
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, TypedDict

import gemmi
import numpy as np
import orjson
import threading

from afdb_integration_kit.uniprot.naming import protein_description
from afdb_integration_kit.utils.rounding import round_float

logger = logging.getLogger(__name__)

# Module-level manifest cache to avoid re-parsing the same CSV file
_MANIFEST_CACHE: dict[str, list[dict[str, str]]] = {}
_MANIFEST_CACHE_LOCK = threading.Lock()

# Thread-local storage for DuckDB connections (DuckDB connections are not thread-safe)
_THREAD_LOCAL = threading.local()

# Module-level cache for prefetched DuckDB metadata (thread-safe read after initial population)
_DUCKDB_METADATA_CACHE: dict[str, dict[str, dict[str, Any]]] = {}
_DUCKDB_METADATA_LOCK = threading.Lock()


def _query_duckdb_entries(
    con: "duckdb.DuckDBPyConnection",
    accessions: Sequence[str],
) -> dict[str, dict[str, Any]]:
    """Fetch entry-table metadata for the requested accessions."""
    unique_accs = list(dict.fromkeys(accessions))
    if not unique_accs:
        return {}

    placeholders = ",".join("?" for _ in unique_accs)
    available_columns = {
        str(row[0]) for row in con.execute("DESCRIBE entry").fetchall()
    }
    selected_columns = [
        name
        for name in (
            "primary_ac",
            "protein_full_names",
            "protein_short_names",
            "entry_name",
            "sequence",
        )
        if name in available_columns
    ]
    query = "SELECT {} FROM entry WHERE primary_ac IN ({})".format(
        ", ".join(selected_columns),
        placeholders,
    )
    rows_rel = con.execute(query, unique_accs)
    rows = rows_rel.fetchall()
    col_index = {
        name: idx
        for idx, name in enumerate([col[0] for col in (rows_rel.description or [])])
    }

    entry_lookup: dict[str, dict[str, Any]] = {}
    for row in rows:
        entry_lookup[str(row[col_index["primary_ac"]])] = {
            name: row[col_index[name]] for name in col_index
        }
    return entry_lookup


def prefetch_duckdb_metadata(db_path: str, accessions: list[str]) -> None:
    """
    Prefetch metadata for all accessions from DuckDB into module-level cache.
    Call this ONCE before batch processing to avoid per-model queries.
    Thread-safe for initial population; subsequent reads don't need locks.
    """
    if not accessions:
        return

    cache_key = str(Path(db_path).resolve())

    with _DUCKDB_METADATA_LOCK:
        if cache_key in _DUCKDB_METADATA_CACHE:
            return  # Already populated

        try:
            import duckdb
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "duckdb Python package is required to read DuckDB manifests. "
                "Install with `pip install duckdb`."
            ) from exc

        con = duckdb.connect(db_path, read_only=True)
        duckdb_mem = os.environ.get("DUCKDB_MEMORY_LIMIT", "512MB")
        con.execute(f"SET memory_limit = '{duckdb_mem}'")
        try:
            unique_accs = list(set(accessions))
            placeholders = ",".join("?" for _ in unique_accs)
            query = (
                "SELECT primary_ac, protein_full_names, sequence "
                "FROM entry WHERE primary_ac IN ({})"
            ).format(placeholders)
            rows_rel = con.execute(query, unique_accs)
            rows = rows_rel.fetchall()
            col_index = {name: idx for idx, name in enumerate([col[0] for col in (rows_rel.description or [])])}

            entry_lookup: dict[str, dict[str, Any]] = {}
            for row in rows:
                entry_lookup[str(row[col_index["primary_ac"]])] = {
                    name: row[col_index[name]] for name in col_index
                }

            _DUCKDB_METADATA_CACHE[cache_key] = entry_lookup
            logger.info("Prefetched DuckDB metadata for %d/%d accessions", len(entry_lookup), len(unique_accs))
        finally:
            con.close()


def _get_prefetched_metadata(db_path: str) -> dict[str, dict[str, Any]] | None:
    """Get prefetched metadata if available, otherwise return None."""
    cache_key = str(Path(db_path).resolve())
    return _DUCKDB_METADATA_CACHE.get(cache_key)


def _get_duckdb_connection(db_path: str) -> "duckdb.DuckDBPyConnection":
    """Get or create a thread-local DuckDB connection."""
    try:
        import duckdb
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "duckdb Python package is required to read DuckDB manifests. "
            "Install with `pip install duckdb`."
        ) from exc

    if not hasattr(_THREAD_LOCAL, "duckdb_conns"):
        _THREAD_LOCAL.duckdb_conns = {}

    if db_path not in _THREAD_LOCAL.duckdb_conns:
        conn = duckdb.connect(db_path, read_only=True)
        duckdb_mem = os.environ.get("DUCKDB_MEMORY_LIMIT", "512MB")
        conn.execute(f"SET memory_limit = '{duckdb_mem}'")
        _THREAD_LOCAL.duckdb_conns[db_path] = conn

    return _THREAD_LOCAL.duckdb_conns[db_path]


def cleanup_caches() -> None:
    """Close all cached DuckDB connections and clear caches. Call after batch processing."""
    # Close thread-local DuckDB connections if they exist
    if hasattr(_THREAD_LOCAL, "duckdb_conns"):
        for conn in _THREAD_LOCAL.duckdb_conns.values():
            try:
                conn.close()
            except Exception:
                pass
        _THREAD_LOCAL.duckdb_conns.clear()
    with _MANIFEST_CACHE_LOCK:
        _MANIFEST_CACHE.clear()
    with _DUCKDB_METADATA_LOCK:
        _DUCKDB_METADATA_CACHE.clear()


class ChainMetadata(TypedDict):
    name: str
    label_asym_id: str
    sequenceStart: int
    sequenceEnd: int


class PAEItem(TypedDict):
    predicted_aligned_error: list[list[float]]
    max_predicted_aligned_error: float
    chains: list[ChainMetadata]


def _categorise_confidence(score: float) -> str:
    """
    Map a single pLDDT score to category:
        "V" : > 90
        "H" : 70–90  (inclusive of 70 and 90)
        "M" : 50–70
        "L" : 30–50
        "D" : < 30
    """
    if score > 90.0:
        return "V"
    if 70.0 <= score <= 90.0:
        return "H"
    if 50.0 <= score < 70.0:
        return "M"
    if 30.0 <= score < 50.0:
        return "L"
    return "D"


def _iterate_pdb_residues(pdb_path: Path) -> Iterable[tuple[str, int, str]]:
    """
    Yield (chain_id, resseq, insertion_code) for each residue in the first MODEL.
    Only ATOM records are considered to align with sequence-derived matrices.
    """
    in_first_model = True
    with pdb_path.open("r") as handle:
        for line in handle:
            if len(line) < 26:
                continue
            record = line[:6].strip()
            if record == "MODEL":
                try:
                    model_idx = int(line.split()[1])
                except (IndexError, ValueError):
                    model_idx = 1
                in_first_model = model_idx == 1
                continue
            if record == "ENDMDL":
                if in_first_model:
                    break
                in_first_model = False
                continue
            if record != "ATOM" or not in_first_model:
                continue

            chain_id = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            insertion_code = line[26].strip()
            yield chain_id, resseq, insertion_code


def _chain_spans_from_pdb(
    pdb_path: Path,
    chain_display_names: Dict[str, str] | None = None,
) -> tuple[list[ChainMetadata], int]:
    """
    Parse the PDB to derive chain metadata using gemmi for speed.
    Falls back to line-by-line parsing if gemmi fails.
    """
    try:
        return _chain_spans_from_pdb_gemmi(pdb_path, chain_display_names)
    except Exception:
        logger.warning("gemmi parsing failed, falling back to line parser")
        return _chain_spans_from_pdb_legacy(pdb_path, chain_display_names)


def _chain_spans_from_pdb_gemmi(
    pdb_path: Path,
    chain_display_names: Dict[str, str] | None = None,
) -> tuple[list[ChainMetadata], int]:
    """Fast PDB parsing using gemmi."""
    structure = gemmi.read_structure(str(pdb_path))
    if not structure:
        raise ValueError(f"No structure found in {pdb_path}")
    
    model = structure[0]
    chains: list[ChainMetadata] = []
    total_residues = 0
    
    for idx, chain in enumerate(model, start=1):
        # Count unique residues (excluding water)
        seen_residues: set[tuple[int, str]] = set()
        for residue in chain:
            if residue.name == "HOH":
                continue
            seen_residues.add((residue.seqid.num, residue.seqid.icode))
        
        if not seen_residues:
            continue
        
        chain_id = chain.name if chain.name.strip() else f"Chain{idx}"
        label = chain_id if chain_id != "_" else f"Chain{idx}"
        display_name = chain_display_names.get(chain_id, label) if chain_display_names else label
        
        # Expose per-chain local residue ranges in JSON metadata.
        start = 1
        end = len(seen_residues)
        chains.append({
            "name": display_name,
            "label_asym_id": label,
            "sequenceStart": start,
            "sequenceEnd": end,
        })
        total_residues += len(seen_residues)
    
    if not chains:
        raise ValueError(f"No chains found in {pdb_path}")
    
    return chains, total_residues


def _get_pdb_chain_ids(pdb_path: Path) -> list[str]:
    """
    Extract ordered list of chain IDs from PDB file.
    Uses gemmi for fast parsing, falls back to line-by-line if needed.
    Returns chain IDs in the order they appear in the PDB.
    """
    try:
        structure = gemmi.read_structure(str(pdb_path))
        if not structure:
            raise ValueError(f"No structure found in {pdb_path}")
        model = structure[0]
        chain_ids = []
        for chain in model:
            # Skip chains with no residues (excluding water)
            has_residues = any(r.name != "HOH" for r in chain)
            if has_residues:
                chain_id = chain.name if chain.name.strip() else "_"
                chain_ids.append(chain_id)
        return chain_ids
    except Exception:
        # Fallback to line-by-line parsing
        seen_chains: list[str] = []
        for chain_id, _, _ in _iterate_pdb_residues(pdb_path):
            if chain_id not in seen_chains:
                seen_chains.append(chain_id)
        return seen_chains


def _chain_spans_from_pdb_legacy(
    pdb_path: Path,
    chain_display_names: Dict[str, str] | None = None,
) -> tuple[list[ChainMetadata], int]:
    """
    Parse the PDB to derive chain metadata aligned to pLDDT/PAE indices.
    Returns (chains, total_residues). Legacy line-by-line fallback.
    """
    residues: OrderedDict[str, list[tuple[int, str]]] = OrderedDict()
    for chain_id, resseq, insertion_code in _iterate_pdb_residues(pdb_path):
        chain_residues = residues.setdefault(chain_id, [])
        resid = (resseq, insertion_code)
        if resid not in chain_residues:
            chain_residues.append(resid)

    if not residues:
        raise ValueError(f"No ATOM records found in {pdb_path} to derive chains.")

    chains: list[ChainMetadata] = []
    total_residues = 0
    for idx, (chain_id, res_list) in enumerate(residues.items(), start=1):
        label = chain_id if chain_id != "_" else f"Chain{idx}"
        display_name = chain_display_names.get(chain_id) if chain_display_names else label
        # Expose per-chain local residue ranges in JSON metadata.
        start = 1
        end = len(res_list)
        chains.append(
            {
                "name": display_name,
                "label_asym_id": label,
                "sequenceStart": start,
                "sequenceEnd": end,
            }
        )
        total_residues += len(res_list)

    return chains, total_residues


def _first_value(row: dict[str, Any], keys: Sequence[str]) -> str | None:
    """Return the first non-empty value for the provided keys in a manifest row."""
    for key in keys:
        if key in row and row[key] not in (None, ""):
            value = str(row[key]).strip()
            if value:
                return value
    return None


def _parse_int(value: str | None) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"Invalid integer value in manifest: {value!r}") from exc


def _load_manifest_chains(
    manifest_path: Path,
    model_entity_id: str | None = None,
) -> tuple[str, list[dict[str, str]]]:
    """
    Load chain_id + uniprot_ac pairs for the requested model_entity_id.
    Returns the resolved model_entity_id and ordered chain rows.
    Uses module-level cache to avoid re-parsing the same manifest file.
    """
    cache_key = str(manifest_path.resolve())

    with _MANIFEST_CACHE_LOCK:
        if cache_key not in _MANIFEST_CACHE:
            with manifest_path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                _MANIFEST_CACHE[cache_key] = list(reader)
        rows = _MANIFEST_CACHE[cache_key]

    if not rows:
        raise ValueError(f"Manifest {manifest_path} is empty.")

    model_ids = {row.get("model_entity_id") for row in rows}
    if model_entity_id is None:
        if len(model_ids) == 1:
            model_entity_id = next(iter(model_ids))
            logger.info("Using sole model_entity_id from manifest: %s", model_entity_id)
        else:
            raise ValueError(
                "manifest contains multiple model_entity_id values; "
                "provide --model-entity-id to disambiguate."
            )

    filtered = [row for row in rows if row.get("model_entity_id") == model_entity_id]
    if not filtered:
        raise ValueError(
            f"manifest {manifest_path} has no rows for model_entity_id={model_entity_id}."
        )

    chain_rows: list[dict[str, str]] = []
    for idx, row in enumerate(filtered, start=1):
        chain_id = (row.get("chain_id") or "").strip() or f"Chain{idx}"
        uniprot_ac = _first_value(row, ["uniprot_ac", "uniprotAccession"])
        if not uniprot_ac:
            raise ValueError(f"Manifest row for chain {chain_id} is missing uniprot_ac.")
        entity_id = (row.get("entity_id") or "").strip()
        chain_row = {
            "chain_id": chain_id,
            "uniprot_ac": uniprot_ac,
            "entity_id": entity_id,
        }
        for optional_key in (
            "sequence_start",
            "sequence_end",
            "protein_name",
            "is_fragment",
            "is_isoform",
            "entity_type",
        ):
            if optional_key == "protein_name" and optional_key in row:
                chain_row[optional_key] = str(
                    row.get(optional_key) or ""
                ).strip()
                continue
            value = row.get(optional_key)
            if value is None:
                continue
            text = str(value).strip()
            if text:
                chain_row[optional_key] = text
        chain_rows.append(chain_row)
    # Assign missing entity IDs by component identity.  Different fragments of
    # one UniProt accession are different components, while repeated copies of
    # the same fragment remain one entity.  Explicit manifest IDs are retained.
    component_to_entity: Dict[tuple[str, int | None, int | None], str] = {}
    used_ids: set[str] = set()
    for row in chain_rows:
        key = _manifest_component_key(row)
        eid = (row.get("entity_id") or "").strip()
        if eid:
            existing = component_to_entity.get(key)
            if existing is None:
                component_to_entity[key] = eid
                used_ids.add(eid)
            elif existing != eid:
                # Swapped duplicate rows can list the same component under two
                # explicit entity IDs; keep the first one so dedup downstream
                # still collapses them. Genuine conflicts are logged, not fatal.
                logger.warning(
                    "Manifest gives one component inconsistent entity IDs: "
                    "%r and %r; keeping the first (%r).",
                    existing, eid, existing,
                )
    next_id = 1
    for row in chain_rows:
        if row.get("entity_id"):
            continue
        key = _manifest_component_key(row)
        if key not in component_to_entity:
            while str(next_id) in used_ids:
                next_id += 1
            component_to_entity[key] = str(next_id)
            used_ids.add(str(next_id))
            next_id += 1
        row["entity_id"] = component_to_entity[key]
    return model_entity_id, chain_rows


def _manifest_component_key(
    row: dict[str, str],
) -> tuple[str, int | None, int | None]:
    """Return the identity used for entity assignment in chain manifests."""
    accession = row["uniprot_ac"]
    is_fragment = str(row.get("is_fragment") or "").strip().lower()
    start = _parse_int(row.get("sequence_start"))
    end = _parse_int(row.get("sequence_end"))
    if is_fragment in {"true", "1", "yes"} or start is not None or end is not None:
        return accession, start, end
    return accession, None, None


def _as_string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        items: list[str] = []
        for item in value:
            if item is None:
                continue
            text = str(item).strip()
            if text:
                items.append(text)
        return items
    if isinstance(value, str):
        text = value.strip()
        # Try to parse JSON array strings as emitted by some DuckDB exports.
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = orjson.loads(text)
                return _as_string_list(parsed)
            except orjson.JSONDecodeError:
                pass
        cleaned = value.strip()
        return [cleaned] if cleaned else []
    return []


# Standard + common modified residue 3-letter -> 1-letter codes, used as a
# fallback when gemmi cannot tabulate a residue while deriving chain sequences.
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def _structure_chain_seqs(pdb_path: Path) -> list[dict[str, Any]]:
    """Return the model's chains in coordinate (file) order.

    Each entry is ``{"label": str, "length": int, "seq": str}`` where ``seq`` is
    the per-chain one-letter sequence (empty if it cannot be derived). This order
    and these per-chain lengths mirror the pLDDT/PAE arrays, which follow the
    predicted structure — NOT the manifest. It is the authoritative ordering for
    assigning chain names/boundaries.
    """
    try:
        structure = gemmi.read_structure(str(pdb_path))
        if not structure:
            raise ValueError(f"No structure found in {pdb_path}")
        model = structure[0]
        out: list[dict[str, Any]] = []
        for idx, chain in enumerate(model, start=1):
            seen: set[tuple[int, str]] = set()
            seq_chars: list[str] = []
            for residue in chain:
                if residue.name == "HOH":
                    continue
                key = (residue.seqid.num, residue.seqid.icode)
                if key in seen:
                    continue
                seen.add(key)
                info = gemmi.find_tabulated_residue(residue.name)
                code = (
                    info.one_letter_code.upper()
                    if info and info.one_letter_code
                    else _AA3TO1.get(residue.name, "X")
                )
                seq_chars.append(code if code.isalpha() else "X")
            if not seen:
                continue
            label = chain.name if chain.name.strip() else f"Chain{idx}"
            out.append({"label": label, "length": len(seen), "seq": "".join(seq_chars)})
        if out:
            return out
    except Exception:
        logger.warning("gemmi parsing failed for %s; using line parser for chain order", pdb_path)

    # Fallback: line parser yields chain order + lengths but no residue names.
    residues: "OrderedDict[str, list[tuple[int, str]]]" = OrderedDict()
    for chain_id, resseq, icode in _iterate_pdb_residues(pdb_path):
        res_list = residues.setdefault(chain_id, [])
        rid = (resseq, icode)
        if rid not in res_list:
            res_list.append(rid)
    fallback: list[dict[str, Any]] = []
    for idx, (chain_id, res_list) in enumerate(residues.items(), start=1):
        label = chain_id if chain_id != "_" else f"Chain{idx}"
        fallback.append({"label": label, "length": len(res_list), "seq": ""})
    return fallback


def _resolve_chain_name(entry: dict[str, Any], acc: str) -> str:
    """Best chain display name for an accession from a DuckDB entry row."""
    names = _as_string_list(entry.get("protein_full_names"))
    if names:
        return names[0]
    if entry.get("entry_name"):
        return str(entry["entry_name"])
    if entry.get("gene_names"):
        return str(entry["gene_names"])
    return acc  # Ultimate fallback to accession itself


def _bind_chains_to_structure(
    physical: list[dict[str, Any]],
    candidates: list[dict[str, Any]],
) -> tuple[list[ChainMetadata], list[int], list[dict[str, str]]]:
    """Assign each physical chain (in structure order) to a candidate accession.

    Matching priority: exact sequence, then sequence length, then positional
    (last resort, logged). Names + boundaries are emitted in PHYSICAL order so
    they always line up with the coordinates and the pLDDT/PAE arrays. This is
    what prevents chain swaps when the manifest order disagrees with the order
    ColabFold actually emitted the chains.

    ``physical``   : [{"label", "length", "seq"}] in coordinate order.
    ``candidates`` : [{"uniprot_ac", "entity_id", "name", "seq"}] (manifest order).
    """
    remaining: list[dict[str, Any] | None] = list(candidates)
    assigned: list[dict[str, Any] | None] = [None] * len(physical)

    def _take(pred) -> dict[str, Any] | None:
        for j, cand in enumerate(remaining):
            if cand is not None and pred(cand):
                remaining[j] = None
                return cand
        return None

    # Pass 1: exact sequence match (disambiguates equal-length pairs).
    for i, pc in enumerate(physical):
        if pc["seq"]:
            cand = _take(lambda c, s=pc["seq"]: c["seq"] == s)
            if cand is not None:
                assigned[i] = cand
    # Pass 1b: ColabFold includes unknown (``X``) residues in pLDDT/PAE, but
    # those residues have no atoms and can therefore be absent from the PDB.
    # Match the remaining coordinate sequence against the UniProt sequence
    # with unknown residues removed so heteromers still bind deterministically.
    for i, pc in enumerate(physical):
        if assigned[i] is None and pc["seq"]:
            cand = _take(lambda c, s=pc["seq"]: c["seq"].replace("X", "") == s)
            if cand is not None:
                assigned[i] = cand
    # Pass 2: sequence-length match.
    for i, pc in enumerate(physical):
        if assigned[i] is None:
            cand = _take(lambda c, n=pc["length"]: len(c["seq"]) == n)
            if cand is not None:
                assigned[i] = cand
    # Pass 3: positional fallback for anything still unmatched.
    leftovers = [c for c in remaining if c is not None]
    li = 0
    for i in range(len(physical)):
        if assigned[i] is None:
            assigned[i] = leftovers[li]
            li += 1
            logger.warning(
                "Chain %s (len %s) not matched to an accession by sequence/length; "
                "assigned positionally to %s.",
                physical[i]["label"], physical[i]["length"], assigned[i]["uniprot_ac"],
            )

    chains: list[ChainMetadata] = []
    residue_numbers: list[int] = []
    effective: list[dict[str, str]] = []
    for pc, cand in zip(physical, assigned):
        # Preserve the full scored sequence when the only coordinate omissions
        # are unresolved UniProt ``X`` residues.  The confidence arrays include
        # these positions even though the PDB cannot provide atoms for them.
        candidate_without_unknowns = cand["seq"].replace("X", "")
        length = (
            len(cand["seq"])
            if pc["seq"] and candidate_without_unknowns == pc["seq"]
            else pc["length"]
        )
        chains.append(
            {
                "name": cand["name"],
                "label_asym_id": pc["label"],
                "sequenceStart": 1,
                "sequenceEnd": length,
            }
        )
        effective_row = {
            "chain_id": pc["label"],
            "uniprot_ac": cand["uniprot_ac"],
            "entity_id": cand.get("entity_id", ""),
            "is_fragment": cand.get("is_fragment", ""),
            "is_isoform": cand.get("is_isoform", ""),
            "entity_type": cand.get("entity_type", "protein"),
            "sequence_start": cand.get("sequence_start", ""),
            "sequence_end": cand.get("sequence_end", ""),
        }
        if cand.get("protein_name"):
            effective_row["protein_name"] = cand["protein_name"]
        effective.append(effective_row)
        residue_numbers.extend(range(len(residue_numbers) + 1, len(residue_numbers) + length + 1))

    # Tripwire: a boundary may exceed the coordinate count only by UniProt
    # unknown residues that are absent from the PDB.
    for chain, pc, cand in zip(chains, physical, assigned):
        emitted = chain["sequenceEnd"] - chain["sequenceStart"] + 1
        allowed = {pc["length"]}
        if pc["seq"] and cand["seq"].replace("X", "") == pc["seq"]:
            allowed.add(len(cand["seq"]))
        if emitted not in allowed:
            raise AssertionError("Chain boundary/order does not match the structure after binding.")
    return chains, residue_numbers, effective


def _dedup_manifest_chains(chains: list[dict[str, str]]) -> list[dict[str, str]]:
    """Drop duplicate chain rows keyed by chain_id (label_asym_id), keeping the first.

    A renamed heterodimer can receive duplicate chain rows in the work manifest
    (two compound IDs -- forward ``AF_A_AF_B`` and swapped ``AF_B_AF_A``, or two
    component pairs resolving to the same accession pair -- renaming to the same
    single AF-ID), which inflates the residue count in stage_03 by exactly 2x.

    The key is ``chain_id`` (the physical chain / ``label_asym_id``), not
    ``entity_id``: ``_load_manifest_chains`` reassigns ``entity_id`` per accession,
    so a swapped duplicate gets a scrambled ``entity_id`` and a composite key would
    not collapse it. Keeping the first occurrence per ``chain_id`` preserves the
    correct accession set; ``_bind_chains_to_structure`` re-derives the per-chain
    assignment from the physical structure, so first-vs-swapped order is irrelevant.
    """
    seen: set[str] = set()
    deduped: list[dict[str, str]] = []
    for row in chains:
        key = row.get("chain_id", "")
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    if len(deduped) != len(chains):
        logger.warning(
            "Dropped %d duplicate manifest chain row(s); kept %d unique chain_id "
            "row(s): %s",
            len(chains) - len(deduped),
            len(deduped),
            [r.get("chain_id", "") for r in deduped],
        )
    return deduped


def _load_chain_metadata_from_duckdb(
    db_path: Path,
    manifest_chains: list[dict[str, str]],
    pdb_path: Path | None = None,
) -> tuple[list[ChainMetadata], list[int], list[dict[str, str]]]:
    """
    Resolve chain names and residue ranges from DuckDB using accessions from the CSV manifest.
    - Uses uniprot_ac to find matching rows in the entry table.
    - Uses the first protein_full_names entry as the chain name (required).
    - Derives sequenceStart/sequenceEnd as per-chain local ranges 1..len(sequence).

    If pdb_path is provided, detects actual chains in the PDB and auto-expands
    the manifest for homomultimers (multiple chains with the same UniProt accession).

    Returns:
        (chains, residue_numbers, effective_manifest_chains)
        - effective_manifest_chains may be expanded for homomultimers

    Uses prefetched metadata cache if available (populated by prefetch_duckdb_metadata).
    """
    # Determine the physical chain order/lengths/sequences from the predicted
    # structure. This is the authoritative ordering for the pLDDT/PAE arrays and
    # for assigning chain names/boundaries.
    physical: list[dict[str, Any]] | None = None
    effective_chains = _dedup_manifest_chains(manifest_chains)
    if pdb_path is not None:
        physical = _structure_chain_seqs(pdb_path)
        pdb_chain_ids = [p["label"] for p in physical]
        if len(pdb_chain_ids) > len(manifest_chains):
            # Check if this is a homomultimer (single accession for multiple chains)
            unique_accs = {c["uniprot_ac"] for c in manifest_chains}
            if len(unique_accs) == 1:
                # Homodimer/homomultimer: replicate the single accession for each PDB chain
                base_acc = list(unique_accs)[0]
                base_entity_id = manifest_chains[0].get("entity_id", "1")
                base_chain = manifest_chains[0]
                effective_chains = []
                for chain_id in pdb_chain_ids:
                    expanded_chain = dict(base_chain)
                    expanded_chain.update({
                        "chain_id": chain_id,
                        "uniprot_ac": base_acc,
                        "entity_id": base_entity_id,
                    })
                    effective_chains.append(expanded_chain)
                logger.info(
                    "Auto-expanded manifest for homomultimer: %d chains -> %d chains (accession: %s)",
                    len(manifest_chains), len(effective_chains), base_acc
                )
            else:
                # Heteromultimer with incomplete manifest - warn but proceed
                logger.warning(
                    "PDB has %d chains but manifest has %d entries with %d unique accessions. "
                    "Cannot auto-expand heteromultimer manifests.",
                    len(pdb_chain_ids), len(manifest_chains), len(unique_accs)
                )

    cache_key = str(db_path.resolve())
    accs = [c["uniprot_ac"] for c in effective_chains]

    # Try to use prefetched metadata first (fast path)
    prefetched = _get_prefetched_metadata(cache_key)
    if prefetched is not None:
        entry_lookup = prefetched
    else:
        # Fallback to per-model query (slow path)
        con = _get_duckdb_connection(cache_key)
        unique_accs_list = list(set(accs))
        placeholders = ",".join("?" for _ in unique_accs_list)
        query = (
            "SELECT primary_ac, protein_full_names, sequence "
            "FROM entry WHERE primary_ac IN ({})"
        ).format(placeholders)
        rows_rel = con.execute(query, unique_accs_list)
        rows = rows_rel.fetchall()
        if not rows:
            raise ValueError(f"No matching accessions found in DuckDB entry table for {unique_accs_list}.")
        col_index = {name: idx for idx, name in enumerate([col[0] for col in (rows_rel.description or [])])}

        entry_lookup: Dict[str, Dict[str, Any]] = {}
        for row in rows:
            entry_lookup[str(row[col_index["primary_ac"]])] = {name: row[col_index[name]] for name in col_index}

    # Verify all unique accessions are in the lookup
    unique_accs_set = set(accs)
    missing = [acc for acc in unique_accs_set if acc not in entry_lookup]
    if missing:
        raise ValueError(f"No matching accessions found in DuckDB entry table for {missing}.")

    # Build candidate accession metadata (name + sequence), preserving manifest order.
    candidates: list[dict[str, Any]] = []
    for chain in effective_chains:
        acc = chain["uniprot_ac"]
        entry = entry_lookup.get(acc)
        if entry is None:
            raise ValueError(f"Accession {acc} not found in DuckDB entry table.")
        desc = protein_description(chain.get("protein_name"), entry, acc)
        seq = entry.get("sequence") or ""
        if len(seq) == 0:
            raise ValueError(f"No sequence found in DuckDB entry table for accession {acc}.")
        candidate = {
            "uniprot_ac": acc,
            "entity_id": chain.get("entity_id", ""),
            "name": desc,
            "seq": seq,
            "is_fragment": chain.get("is_fragment", ""),
            "is_isoform": chain.get("is_isoform", ""),
            "entity_type": chain.get("entity_type", "protein"),
            "sequence_start": chain.get("sequence_start", ""),
            "sequence_end": chain.get("sequence_end", ""),
        }
        if chain.get("protein_name"):
            candidate["protein_name"] = chain["protein_name"]
        candidates.append(candidate)

    # Preferred path: bind names/boundaries to the physical structure order so
    # they always match the coordinates and PAE/pLDDT arrays. This fixes chain
    # swaps caused by the manifest order disagreeing with ColabFold's output.
    if physical is not None and len(physical) == len(candidates):
        return _bind_chains_to_structure(physical, candidates)

    # Fallback (no structure, or chain-count mismatch): legacy manifest-order build.
    if physical is not None:
        logger.warning(
            "Structure chain count (%d) != accession count (%d); using manifest order for %s.",
            len(physical), len(candidates), accs,
        )
    chains: list[ChainMetadata] = []
    residue_numbers: list[int] = []
    for chain, cand in zip(effective_chains, candidates):
        seqlen = len(cand["seq"])
        chains.append(
            {
                "name": cand["name"],
                "label_asym_id": chain["chain_id"],
                "sequenceStart": 1,
                "sequenceEnd": seqlen,
            }
        )
        residue_numbers.extend(range(len(residue_numbers) + 1, len(residue_numbers) + seqlen + 1))

    return chains, residue_numbers, effective_chains


def _compute_plddt_metrics(
    plddt: Sequence[float],
    chains: list[ChainMetadata],
    manifest_chains: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], float]:
    """
    Compute per-chain average pLDDT and model-level average.
    Returns (chain_rows, model_avg_plddt).
    """
    row_lookup = {c["chain_id"]: c for c in manifest_chains}
    all_scores = np.asarray(plddt, dtype=np.float64)
    chain_rows: list[dict[str, Any]] = []
    total_sum = 0.0
    total_count = 0
    running_offset = 0
    for chain in chains:
        start = chain["sequenceStart"]
        end = chain["sequenceEnd"]
        length = end - start + 1
        values = all_scores[running_offset : running_offset + length]
        if len(values) == 0:
            raise ValueError(f"No pLDDT values found for chain {chain['label_asym_id']}.")
        running_offset += length
        n = len(values)
        chain_sum = float(np.sum(values))
        avg = round_float(chain_sum / n, 2)
        total_sum += chain_sum
        total_count += n
        fraction_plddt_very_low = round_float(int(np.sum(values < 50.0)) / n, 3)
        fraction_plddt_low = round_float(int(np.sum((values >= 50.0) & (values < 70.0))) / n, 3)
        fraction_plddt_confident = round_float(int(np.sum((values >= 70.0) & (values <= 90.0))) / n, 3)
        fraction_plddt_very_high = round_float(int(np.sum(values > 90.0)) / n, 3)
        manifest_row = row_lookup.get(chain["label_asym_id"], {})
        parsed_start = _parse_int(manifest_row.get("sequence_start"))
        sequence_start = parsed_start if parsed_start is not None else 1
        parsed_end = _parse_int(manifest_row.get("sequence_end"))
        sequence_end = parsed_end if parsed_end is not None else length
        chain_row = {
            "model_entity_id": None,  # filled later
            "entity_id": manifest_row.get("entity_id", ""),
            "chain_id": chain["label_asym_id"],
            "uniprot_ac": manifest_row.get("uniprot_ac", ""),
            "is_fragment": manifest_row.get("is_fragment", ""),
            "is_isoform": manifest_row.get("is_isoform", ""),
            "entity_type": manifest_row.get("entity_type", "protein"),
            "sequence_start": sequence_start,
            "sequence_end": sequence_end,
            "average_plddt": avg,
            "fraction_plddt_very_low": fraction_plddt_very_low,
            "fraction_plddt_low": fraction_plddt_low,
            "fraction_plddt_confident": fraction_plddt_confident,
            "fraction_plddt_very_high": fraction_plddt_very_high,
        }
        if "protein_name" in manifest_row:
            chain_row["protein_name"] = manifest_row["protein_name"]
        chain_rows.append(chain_row)
    model_avg = round_float(total_sum / total_count, 2) if total_count else 0.0
    return chain_rows, model_avg


def _write_manifest_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _merge_manifest_csv(
    path: Path,
    base_fieldnames: list[str],
    model_entity_id: str,
    new_rows: list[dict[str, Any]],
) -> None:
    """
    Merge manifest rows by replacing any existing rows for model_entity_id with new_rows.
    Preserves other models and unions any extra columns present in the existing file.
    """
    existing_rows: list[dict[str, Any]] = []
    existing_fields: list[str] = []
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            existing_fields = reader.fieldnames or []
            for row in reader:
                if (row.get("model_entity_id") or "") != model_entity_id:
                    existing_rows.append(row)

    merged_fields_set = set(base_fieldnames) | set(existing_fields)
    merged_fieldnames = [fn for fn in base_fieldnames if fn in merged_fields_set] + [
        fn for fn in existing_fields if fn not in base_fieldnames
    ]

    def _normalize(row: dict[str, Any]) -> dict[str, Any]:
        return {fn: row.get(fn, "") for fn in merged_fieldnames}

    normalized_rows = [_normalize(r) for r in existing_rows] + [_normalize(r) for r in new_rows]
    _write_manifest_csv(path, merged_fieldnames, normalized_rows)


def plddt_to_ingest(
    plddt: Sequence[float],
    chains: list[ChainMetadata],
    residue_numbers: Sequence[int] | None = None,
) -> Dict[str, Any]:
    """Build the AFDB pLDDT payload."""
    residue_numbers = (
        list(residue_numbers)
        if residue_numbers is not None
        else list(range(1, len(plddt) + 1))
    )  # 1-based indexing
    if len(residue_numbers) != len(plddt):
        raise ValueError(
            f"Residue numbers length ({len(residue_numbers)}) does not match pLDDT length ({len(plddt)})."
        )
    scores = np.asarray(plddt, dtype=np.float64)
    rounded_scores = np.round(scores, 2)

    categories = np.full(len(rounded_scores), "D", dtype="U1")
    categories[rounded_scores >= 30.0] = "L"
    categories[rounded_scores >= 50.0] = "M"
    categories[rounded_scores >= 70.0] = "H"
    categories[rounded_scores > 90.0] = "V"

    return {
        "residueNumber": residue_numbers,
        "confidenceScore": rounded_scores,
        "confidenceCategory": categories.tolist(),
        "chains": chains,
    }


def pae_to_ingest(pae: Sequence[Sequence[float]], max_pae: float, chains: list[ChainMetadata]) -> list[PAEItem]:
    """Build the AFDB PAE payload with light validation."""
    pae_arr = np.asarray(pae, dtype=np.float64)
    if pae_arr.ndim != 2 or pae_arr.shape[0] != pae_arr.shape[1] or pae_arr.shape[0] == 0:
        raise ValueError("PAE must be a non-empty square matrix (NxN).")
    rounded_pae = np.round(pae_arr, 2)
    return [
        {
            "predicted_aligned_error": rounded_pae,
            "max_predicted_aligned_error": round(max_pae, 2),
            "chains": chains,
        }
    ]


def convert_file(
    scores_json_path: str,
    pdb_path: str,
    out_plddt_path: str | None = None,
    out_pae_path: str | None = None,
    outdir: str | None = None,
    manifest_path: str | None = None,
    model_entity_id: str | None = None,
    duckdb_path: str | None = None,
    out_chain_manifest: str | None = None,
    out_model_manifest: str | None = None,
    chain_manifest_dir: str | None = None,
    model_manifest_dir: str | None = None,
) -> Dict[str, str]:
    """
    Convert ColabFold score JSON + PDB into AFDB-format JSONs.
    Returns written paths.

    If a CSV manifest and DuckDB path are provided, chain names (uniprotDescription)
    and residue ranges are resolved by:
      1) reading chain_id/uniprot_ac from the CSV manifest (filtered by model_entity_id)
      2) looking up those accessions in the DuckDB entry table to fetch uniprotDescription
         and sequence length (start=1, end=len(sequence))
    If DuckDB lookup is unavailable, the converter falls back to PDB parsing (using
    manifest-provided chain names when possible).
    """
    scores_path = Path(scores_json_path)
    pdb = Path(pdb_path)

    data = orjson.loads(scores_path.read_bytes())

    try:
        plddt = data["plddt"]
        pae = data["pae"]
        max_pae = data["max_pae"]
    except KeyError as e:
        raise KeyError(f"Input JSON is missing required key: {e}")

    # pLDDT and the PAE matrix must describe the same residues, in the same order.
    n_res = len(plddt)
    if not hasattr(pae, "__len__") or len(pae) != n_res:
        raise ValueError(
            f"PAE matrix size ({len(pae) if hasattr(pae, '__len__') else 'NA'}) "
            f"does not match pLDDT length ({n_res})."
        )

    resolved_model_id = model_entity_id
    manifest_chains: list[dict[str, str]] = []
    if manifest_path:
        resolved_model_id, manifest_chains = _load_manifest_chains(Path(manifest_path), model_entity_id)
    elif duckdb_path:
        raise ValueError("CSV manifest is required when using --duckdb to map chains to accessions.")

    chains: list[ChainMetadata]
    residue_numbers: list[int] | None = None
    effective_manifest_chains: list[dict[str, str]] = manifest_chains

    if manifest_chains and duckdb_path:
        chains, residue_numbers, effective_manifest_chains = _load_chain_metadata_from_duckdb(
            Path(duckdb_path),
            manifest_chains=manifest_chains,
            pdb_path=pdb,
        )
        if len(residue_numbers) != len(plddt):
            raise ValueError(
                f"DuckDB residue count ({len(residue_numbers)}) does not match pLDDT length ({len(plddt)})."
            )
    else:
        display_names = (
            {
                c["chain_id"]: c.get("protein_name") or c["uniprot_ac"]
                for c in manifest_chains
            }
            if manifest_chains
            else None
        )
        chains, pdb_residue_total = _chain_spans_from_pdb(pdb, display_names)
        if pdb_residue_total != len(plddt):
            raise ValueError(
                f"PDB residue count ({pdb_residue_total}) does not match pLDDT length ({len(plddt)})."
            )
        residue_numbers = None

    # Order-sensitive guard: the chain annotation must cover exactly the residues
    # present in the pLDDT/PAE arrays. Catches any future regression that lets the
    # chain boundaries drift out of sync with the structure (the original bug).
    chain_total = sum(c["sequenceEnd"] - c["sequenceStart"] + 1 for c in chains)
    if chain_total != n_res:
        raise ValueError(
            f"Chain boundaries cover {chain_total} residues but pLDDT/PAE have {n_res}."
        )

    plddt_payload = plddt_to_ingest(plddt, chains, residue_numbers)
    pae_payload = pae_to_ingest(pae, max_pae, chains)

    chain_manifest_rows: list[dict[str, Any]] | None = None
    model_avg_plddt: float | None = None
    if out_chain_manifest or out_model_manifest or chain_manifest_dir or model_manifest_dir:
        if not effective_manifest_chains:
            raise ValueError("--manifest is required when writing pLDDT manifests.")
        chain_manifest_rows, model_avg_plddt = _compute_plddt_metrics(plddt, chains, effective_manifest_chains)
        if model_entity_id:
            for row in chain_manifest_rows:
                row["model_entity_id"] = model_entity_id
                row.setdefault("is_fragment", "")
                row.setdefault("is_isoform", "")
                row.setdefault("entity_type", "protein")
                row.setdefault("protein_name", "")
        else:
            raise ValueError("--model-entity-id is required when writing pLDDT manifests.")

    base_name = model_entity_id or scores_path.stem
    default_plddt_name = f"{base_name}-confidence_v1.json"
    default_pae_name = f"{base_name}-predicted_aligned_error_v1.json"

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out_plddt_path = out_plddt_path or os.path.join(outdir, default_plddt_name)
        out_pae_path = out_pae_path or os.path.join(outdir, default_pae_name)
    else:
        out_plddt_path = out_plddt_path or default_plddt_name
        out_pae_path = out_pae_path or default_pae_name

    def _dump(obj: Any, path: str) -> None:
        with open(path, "wb") as f:
            f.write(orjson.dumps(obj, option=orjson.OPT_SERIALIZE_NUMPY))

    _dump(plddt_payload, out_plddt_path)
    _dump(pae_payload, out_pae_path)

    if chain_manifest_rows is not None:
        chain_manifest_fieldnames = [
            "model_entity_id",
            "entity_id",
            "chain_id",
            "uniprot_ac",
            "protein_name",
        ]
        chain_manifest_fieldnames.extend(
            [
                "is_fragment",
                "is_isoform",
                "entity_type",
                "sequence_start",
                "sequence_end",
                "average_plddt",
                "fraction_plddt_very_low",
                "fraction_plddt_low",
                "fraction_plddt_confident",
                "fraction_plddt_very_high",
            ]
        )
        chain_manifest_path: Path | None = None
        if chain_manifest_dir:
            chain_manifest_path = Path(chain_manifest_dir) / f"{base_name}_afid_mapping.csv"
            _write_manifest_csv(
                chain_manifest_path,
                chain_manifest_fieldnames,
                chain_manifest_rows,
            )
        elif out_chain_manifest:
            _merge_manifest_csv(
                Path(out_chain_manifest),
                chain_manifest_fieldnames,
                model_entity_id,
                chain_manifest_rows,
            )

    if model_avg_plddt is not None:
        model_manifest_path: Path | None = None
        model_rows = [
            {
                "model_entity_id": model_entity_id,
                "ipTM": data.get("iptm", data.get("ipTM", "")),
                "average_plddt": model_avg_plddt,
                "complexName": "",
                "isAMdata": "",
            }
        ]
        if model_manifest_dir:
            model_manifest_path = Path(model_manifest_dir) / f"{base_name}_model_metadata.csv"
            _write_manifest_csv(
                model_manifest_path,
                [
                    "model_entity_id",
                    "ipTM",
                    "average_plddt",
                    "complexName",
                    "isAMdata",
                ],
                model_rows,
            )
        elif out_model_manifest:
            _merge_manifest_csv(
                Path(out_model_manifest),
                [
                    "model_entity_id",
                    "ipTM",
                    "average_plddt",
                    "complexName",
                    "isAMdata",
                ],
                model_entity_id,
                model_rows,
            )

    return {"plddt": out_plddt_path, "pae": out_pae_path}


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Convert ColabFold score JSON and accompanying PDB to AFDB ingestion format "
            "(outputs: *-confidence_v1.json and *-predicted_aligned_error_v1.json)."
        )
    )
    p.add_argument("input", help="Path to ColabFold score JSON with keys: plddt, pae, max_pae")
    p.add_argument("pdb", help="Path to the corresponding PDB (for chain metadata)")
    p.add_argument("--outdir", help="Directory to write outputs (defaults use AFDB names)")
    p.add_argument("--plddt", help="Explicit output path for pLDDT JSON")
    p.add_argument("--pae", help="Explicit output path for PAE JSON")
    p.add_argument(
        "--manifest",
        help=(
            "Optional CSV manifest mapping model_entity_id,chain_id -> chain metadata. "
            "If sequence_start/sequence_end columns are present, they are used instead of parsing the PDB."
        ),
    )
    p.add_argument(
        "--duckdb",
        help=(
            "Optional DuckDB file containing the 'entry' table. "
            "When supplied, chain names (uniprotDescription) and sequence lengths are read from 'entry'."
        ),
    )
    p.add_argument("--model-entity-id", help="Model entity ID to select rows from manifest when provided")
    p.add_argument(
        "--chain-manifest-out",
        help="Optional output CSV for chain-level metrics (adds average_plddt). Requires --manifest and --model-entity-id.",
    )
    p.add_argument(
        "--model-manifest-out",
        help="Optional output CSV for model-level metrics (average_plddt). Requires --manifest and --model-entity-id.",
    )
    p.add_argument(
        "--chain-manifest-dir",
        help="Optional directory to write per-model chain manifests (<model_entity_id>_afid_mapping.csv).",
    )
    p.add_argument(
        "--model-manifest-dir",
        help="Optional directory to write per-model model manifests (<model_entity_id>_model_metadata.csv).",
    )
    return p


def main(argv: List[str] | None = None) -> None:
    args = _build_arg_parser().parse_args(argv)
    paths = convert_file(
        args.input,
        args.pdb,
        out_plddt_path=args.plddt,
        out_pae_path=args.pae,
        outdir=args.outdir,
        manifest_path=args.manifest,
        model_entity_id=args.model_entity_id,
        duckdb_path=args.duckdb,
        out_chain_manifest=args.chain_manifest_out,
        out_model_manifest=args.model_manifest_out,
        chain_manifest_dir=args.chain_manifest_dir,
        model_manifest_dir=args.model_manifest_dir,
    )
    print(orjson.dumps(paths, option=orjson.OPT_INDENT_2).decode())


if __name__ == "__main__":
    main()
