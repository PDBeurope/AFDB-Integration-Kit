from __future__ import annotations

import argparse
import csv
import logging
import os
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, TypedDict

import orjson
logger = logging.getLogger(__name__)


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
    Parse the PDB to derive chain metadata aligned to pLDDT/PAE indices.
    Returns (chains, total_residues).
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
    running_total = 1
    for idx, (chain_id, res_list) in enumerate(residues.items(), start=1):
        label = chain_id if chain_id != "_" else f"Chain{idx}"
        display_name = chain_display_names.get(chain_id) if chain_display_names else label
        start = running_total
        end = running_total + len(res_list) - 1
        chains.append(
            {
                "name": display_name,
                "label_asym_id": label,
                "sequenceStart": start,
                "sequenceEnd": end,
            }
        )
        running_total = end + 1

    total_residues = running_total - 1
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
    """
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

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
        chain_rows.append(
            {
                "chain_id": chain_id,
                "uniprot_ac": uniprot_ac,
            }
        )
    return model_entity_id, chain_rows


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


def _load_chain_metadata_from_duckdb(
    db_path: Path,
    manifest_chains: list[dict[str, str]],
) -> tuple[list[ChainMetadata], list[int]]:
    """
    Resolve chain names and residue ranges from DuckDB using accessions from the CSV manifest.
    - Uses uniprot_ac to find matching rows in the entry table.
    - Uses the first protein_full_names entry as the chain name (required).
    - Derives sequenceStart/sequenceEnd as 1..len(sequence) (required).
    """
    try:
        import duckdb  # type: ignore
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "duckdb Python package is required to read DuckDB manifests. "
            "Install with `pip install duckdb`."
        ) from exc

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        accs = [c["uniprot_ac"] for c in manifest_chains]
        placeholders = ",".join("?" for _ in accs)
        query = (
            "SELECT primary_ac, protein_full_names, sequence "
            "FROM entry WHERE primary_ac IN ({})"
        ).format(placeholders)
        rows_rel = con.execute(query, accs)
        rows = rows_rel.fetchall()
        if not rows:
            raise ValueError(f"No matching accessions found in DuckDB entry table for {accs}.")
        col_index = {name: idx for idx, name in enumerate([col[0] for col in (rows_rel.description or [])])}

        entry_lookup: Dict[str, Dict[str, Any]] = {}
        for row in rows:
            entry_lookup[str(row[col_index["primary_ac"]])] = {name: row[col_index[name]] for name in col_index}

        chains: list[ChainMetadata] = []
        residue_numbers: list[int] = []
        for chain in manifest_chains:
            acc = chain["uniprot_ac"]
            chain_id = chain["chain_id"]
            entry = entry_lookup.get(acc)
            if entry is None:
                raise ValueError(f"Accession {acc} not found in DuckDB entry table.")
            names = _as_string_list(entry.get("protein_full_names"))
            if not names:
                raise ValueError(f"No protein_full_names found in DuckDB entry table for accession {acc}.")
            desc = names[0]
            seq = entry.get("sequence") or ""
            seqlen = len(seq)
            if seqlen == 0:
                raise ValueError(f"No sequence found in DuckDB entry table for accession {acc}.")
            seq_start = 1
            seq_end = seqlen
            chains.append(
                {
                    "name": desc,
                    "label_asym_id": chain_id,
                    "sequenceStart": seq_start,
                    "sequenceEnd": seq_end,
                }
            )
            residue_numbers.extend(range(seq_start, seq_end + 1))

        return chains, residue_numbers
    finally:
        con.close()


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
    categories = [_categorise_confidence(x) for x in plddt]
    return {
        "residueNumber": residue_numbers,
        "confidenceScore": list(plddt),
        "confidenceCategory": categories,
        "chains": chains,
    }


def pae_to_ingest(pae: Sequence[Sequence[float]], max_pae: float, chains: list[ChainMetadata]) -> list[PAEItem]:
    """Build the AFDB PAE payload with light validation."""
    if not pae or any(len(row) != len(pae) for row in pae):
        raise ValueError("PAE must be a non-empty square matrix (NxN).")
    return [
        {
            "predicted_aligned_error": pae,
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

    resolved_model_id = model_entity_id
    manifest_chains: list[dict[str, str]] = []
    if manifest_path:
        resolved_model_id, manifest_chains = _load_manifest_chains(Path(manifest_path), model_entity_id)
    elif duckdb_path:
        raise ValueError("CSV manifest is required when using --duckdb to map chains to accessions.")

    chains: list[ChainMetadata]
    residue_numbers: list[int] | None = None

    if manifest_chains and duckdb_path:
        chains, residue_numbers = _load_chain_metadata_from_duckdb(
            Path(duckdb_path),
            manifest_chains=manifest_chains,
        )
        if len(residue_numbers) != len(plddt):
            raise ValueError(
                f"DuckDB residue count ({len(residue_numbers)}) does not match pLDDT length ({len(plddt)})."
            )
    else:
        display_names = {c["chain_id"]: c["uniprot_ac"] for c in manifest_chains} if manifest_chains else None
        chains, pdb_residue_total = _chain_spans_from_pdb(pdb, display_names)
        if pdb_residue_total != len(plddt):
            raise ValueError(
                f"PDB residue count ({pdb_residue_total}) does not match pLDDT length ({len(plddt)})."
            )
        residue_numbers = None

    plddt_payload = plddt_to_ingest(plddt, chains, residue_numbers)
    pae_payload = pae_to_ingest(pae, max_pae, chains)

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
            f.write(orjson.dumps(obj))

    _dump(plddt_payload, out_plddt_path)
    _dump(pae_payload, out_pae_path)

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
    )
    print(orjson.dumps(paths, option=orjson.OPT_INDENT_2).decode())


if __name__ == "__main__":
    main()
