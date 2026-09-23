#!/usr/bin/env python3
"""Create deterministic synthetic PDB/score fixtures for fragment E2E tests.

The target identities and residue ranges come from the fragment E2E manifest.
Coordinates and confidence values are mechanically remapped from two explicit
donor models.  The outputs are suitable only for testing software behaviour.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi


ALGORITHM = "deterministic-backbone-crop-tile-v1"
WARNING = "synthetic_for_software_testing_only"
BACKBONE_ATOMS = ("N", "CA", "C", "O")
MONOMER_DONOR_ID = "AF-0000000300000001"
DIMER_DONOR_ID = "AF-0000000065760046"
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DONOR_DIR = Path("/mnt/disks/toolkit-data/viruses/sample_data")
DEFAULT_EXAMPLE_DIR = REPO_ROOT / "examples" / "fragment_metadata_e2e"
DEFAULT_OUTPUT_DIR = Path(
    "/mnt/disks/toolkit-data/viruses/fragment_metadata_synthetic_e2e"
    "/generated/input"
)


@dataclass(frozen=True)
class AtomData:
    name: str
    x: float
    y: float
    z: float
    element: str


@dataclass(frozen=True)
class ResidueData:
    atoms: tuple[AtomData, ...]


@dataclass(frozen=True)
class Donor:
    donor_id: str
    pdb_path: Path
    metadata_path: Path
    chains: tuple[tuple[ResidueData, ...], ...]
    chain_offsets: tuple[int, ...]
    plddt: tuple[float, ...]
    pae: tuple[tuple[float, ...], ...]
    max_pae: float
    ptm: float
    iptm: float | None


@dataclass(frozen=True)
class TargetChain:
    chain_id: str
    uniprot_ac: str
    is_fragment: bool
    sequence_start: int
    sequence_end: int
    length: int
    entity_id: str
    protein_name: str


@dataclass(frozen=True)
class TargetModel:
    model_id: str
    chains: tuple[TargetChain, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric.")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result


def _load_score_payload(path: Path) -> tuple[
    tuple[float, ...], tuple[tuple[float, ...], ...], float, float,
    float | None,
]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read donor score JSON {path}: {exc}")
    if not isinstance(payload, dict):
        raise ValueError(f"Donor score JSON must be an object: {path}")
    raw_plddt = payload.get("plddt")
    raw_pae = payload.get("pae", payload.get("predicted_aligned_error"))
    if not isinstance(raw_plddt, list) or not raw_plddt:
        raise ValueError(f"Donor {path} has no non-empty plddt vector.")
    if not isinstance(raw_pae, list) or not raw_pae:
        raise ValueError(f"Donor {path} has no non-empty PAE matrix.")
    plddt = tuple(
        _finite_number(value, f"plddt[{index}]")
        for index, value in enumerate(raw_plddt)
    )
    pae_rows: list[tuple[float, ...]] = []
    for row_index, raw_row in enumerate(raw_pae):
        if not isinstance(raw_row, list) or len(raw_row) != len(plddt):
            raise ValueError(
                f"Donor PAE row {row_index} does not have {len(plddt)} "
                "columns."
            )
        pae_rows.append(tuple(
            _finite_number(value, f"pae[{row_index}][{column_index}]")
            for column_index, value in enumerate(raw_row)
        ))
    if len(pae_rows) != len(plddt):
        raise ValueError("Donor PAE row count does not match plddt length.")
    max_pae = _finite_number(payload.get("max_pae"), "max_pae")
    observed_max = max(max(row) for row in pae_rows)
    # ColabFold JSONs commonly serialize ``max_pae`` and matrix entries at
    # slightly different precision.  Accept only that small rounding gap and
    # normalize upward so the generated invariant remains exact.
    if max_pae + 0.01 < observed_max:
        raise ValueError(
            f"Donor max_pae {max_pae} is below matrix maximum "
            f"{observed_max}."
        )
    ptm = _finite_number(payload.get("ptm"), "ptm")
    iptm_raw = payload.get("iptm")
    iptm = None if iptm_raw is None else _finite_number(iptm_raw, "iptm")
    return plddt, tuple(pae_rows), max(max_pae, observed_max), ptm, iptm


def _structure_chains(path: Path) -> tuple[tuple[ResidueData, ...], ...]:
    try:
        structure = gemmi.read_structure(str(path))
    except (RuntimeError, OSError) as exc:
        raise ValueError(f"Cannot read donor PDB {path}: {exc}")
    if len(structure) != 1:
        raise ValueError(f"Donor PDB must contain exactly one model: {path}")
    chains: list[tuple[ResidueData, ...]] = []
    for chain in structure[0]:
        residues: list[ResidueData] = []
        for residue in chain:
            if residue.het_flag != "A":
                continue
            selected: list[AtomData] = []
            for atom_name in BACKBONE_ATOMS:
                matches = [
                    atom for atom in residue
                    if atom.name.strip() == atom_name
                    and atom.altloc in ("\x00", "A")
                ]
                if len(matches) != 1:
                    raise ValueError(
                        f"Donor residue {chain.name}:{residue.seqid.num} "
                        f"must have one {atom_name} atom."
                    )
                atom = matches[0]
                coords = (atom.pos.x, atom.pos.y, atom.pos.z)
                if not all(math.isfinite(value) for value in coords):
                    raise ValueError("Donor coordinates must be finite.")
                selected.append(AtomData(
                    atom_name, *coords, atom.element.name or atom_name[0]
                ))
            residues.append(ResidueData(tuple(selected)))
        if residues:
            chains.append(tuple(residues))
    if not chains:
        raise ValueError(f"Donor PDB has no polymer chains: {path}")
    return tuple(chains)


def load_donor(donor_dir: Path, donor_id: str, chain_count: int) -> Donor:
    """Load one explicitly named and internally consistent donor pair."""
    pdb_path = donor_dir / f"{donor_id}-model_v1.pdb"
    metadata_path = donor_dir / f"{donor_id}-meta_v1.json"
    if not pdb_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(
            f"Required donor pair is missing for {donor_id}: "
            f"{pdb_path}, {metadata_path}"
        )
    chains = _structure_chains(pdb_path)
    if len(chains) != chain_count:
        raise ValueError(
            f"Donor {donor_id} has {len(chains)} chains; expected "
            f"{chain_count}."
        )
    plddt, pae, max_pae, ptm, iptm = _load_score_payload(metadata_path)
    lengths = [len(chain) for chain in chains]
    if sum(lengths) != len(plddt):
        raise ValueError(
            f"Donor {donor_id} has {sum(lengths)} PDB residues but "
            f"{len(plddt)} plddt values."
        )
    offsets: list[int] = []
    offset = 0
    for length in lengths:
        offsets.append(offset)
        offset += length
    if chain_count > 1 and iptm is None:
        raise ValueError(f"Dimer donor {donor_id} lacks iptm.")
    return Donor(
        donor_id, pdb_path, metadata_path, chains, tuple(offsets), plddt,
        pae, max_pae, ptm, iptm,
    )


def load_targets(manifest_path: Path, seed_path: Path) -> tuple[TargetModel, ...]:
    """Resolve target chain lengths from canonical fragment semantics."""
    seed_payload = json.loads(seed_path.read_text(encoding="utf-8"))
    entries = seed_payload.get("entries")
    if not isinstance(entries, list):
        raise ValueError("Mock seed must contain an entries list.")
    lengths = {
        str(item["primary_ac"]): int(item["sequence_length"])
        for item in entries
    }
    by_model: dict[str, list[TargetChain]] = {}
    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "model_entity_id", "entity_id", "chain_id", "uniprot_ac",
            "is_fragment", "sequence_start", "sequence_end",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("Canonical manifest lacks required columns.")
        for row_number, row in enumerate(reader, 2):
            model_id = row["model_entity_id"].strip()
            chain_id = row["chain_id"].strip()
            accession = row["uniprot_ac"].strip()
            if chain_id not in ("A", "B"):
                raise ValueError(
                    f"Row {row_number}: only A/B chains are supported."
                )
            if accession not in lengths:
                raise ValueError(
                    f"Row {row_number}: unknown seed accession {accession}."
                )
            raw_fragment = row["is_fragment"].strip().lower()
            if raw_fragment not in ("true", "false"):
                raise ValueError(
                    f"Row {row_number}: is_fragment must be true or false."
                )
            is_fragment = raw_fragment == "true"
            if is_fragment:
                try:
                    start = int(row["sequence_start"])
                    end = int(row["sequence_end"])
                except ValueError as exc:
                    raise ValueError(
                        f"Row {row_number}: fragment range is invalid."
                    ) from exc
                if start < 1 or end < start or end > lengths[accession]:
                    raise ValueError(
                        f"Row {row_number}: fragment range is out of bounds."
                    )
                length = end - start + 1
            else:
                start = 1
                end = lengths[accession]
                length = lengths[accession]
            by_model.setdefault(model_id, []).append(TargetChain(
                chain_id, accession, is_fragment, start, end, length,
                row["entity_id"].strip(), row.get("protein_name", "").strip(),
            ))
    targets: list[TargetModel] = []
    for model_id, chains in by_model.items():
        if not model_id or not chains or len(chains) > 2:
            raise ValueError(f"Invalid target topology for {model_id!r}.")
        expected_ids = [chr(ord("A") + index) for index in range(len(chains))]
        actual_ids = [chain.chain_id for chain in chains]
        if actual_ids != expected_ids:
            raise ValueError(
                f"Target {model_id} chains must be ordered {expected_ids}."
            )
        targets.append(TargetModel(model_id, tuple(chains)))
    if not targets:
        raise ValueError("Canonical manifest contains no target models.")
    return tuple(targets)


def _translation_vector(chain: tuple[ResidueData, ...]) -> tuple[float, ...]:
    positions = [
        (atom.x, atom.y, atom.z)
        for residue in chain for atom in residue.atoms
    ]
    spans = tuple(
        max(position[axis] for position in positions)
        - min(position[axis] for position in positions)
        for axis in range(3)
    )
    # A positive shift in every axis makes repeated copies unambiguous even
    # when a donor happens to be flat along one dimension.
    return tuple(max(span + 10.0, 10.0) for span in spans)


def synthesize_model(
    target: TargetModel, donor: Donor
) -> tuple[list[str], dict[str, Any], dict[str, Any]]:
    """Return PDB lines, raw score JSON, and model provenance."""
    if len(target.chains) != len(donor.chains):
        raise ValueError("Target and donor topology do not match.")
    pdb_lines: list[str] = []
    global_source_indices: list[int] = []
    atom_serial = 1
    chain_provenance: list[dict[str, Any]] = []
    for chain_index, target_chain in enumerate(target.chains):
        donor_chain = donor.chains[chain_index]
        donor_length = len(donor_chain)
        translation = _translation_vector(donor_chain)
        source_indices: list[int] = []
        for target_index in range(target_chain.length):
            source_index = target_index % donor_length
            tile_index = target_index // donor_length
            source_indices.append(source_index)
            global_source_index = (
                donor.chain_offsets[chain_index] + source_index
            )
            global_source_indices.append(global_source_index)
            score = donor.plddt[global_source_index]
            residue = donor_chain[source_index]
            for atom in residue.atoms:
                coords = (
                    atom.x + tile_index * translation[0],
                    atom.y + tile_index * translation[1],
                    atom.z + tile_index * translation[2],
                )
                if any(abs(value) >= 9999.999 for value in coords):
                    raise ValueError("Synthetic coordinates exceed PDB limits.")
                pdb_lines.append(
                    f"ATOM  {atom_serial:5d} {atom.name:^4s} ALA "
                    f"{target_chain.chain_id:1s}{target_index + 1:4d}    "
                    f"{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}"
                    f"  1.00{score:6.2f}          {atom.element:>2s}  "
                )
                atom_serial += 1
        pdb_lines.append(
            f"TER   {atom_serial:5d}      ALA "
            f"{target_chain.chain_id:1s}{target_chain.length:4d}"
        )
        atom_serial += 1
        chain_provenance.append({
            "chain_id": target_chain.chain_id,
            "entity_id": target_chain.entity_id,
            "uniprot_ac": target_chain.uniprot_ac,
            "is_fragment": target_chain.is_fragment,
            "sequence_start": target_chain.sequence_start,
            "sequence_end": target_chain.sequence_end,
            "target_length": target_chain.length,
            "donor_chain_index": chain_index,
            "donor_chain_length": donor_length,
            "mapping": "target_index_modulo_donor_chain_length",
            "tile_count": math.ceil(target_chain.length / donor_length),
            "cropped": target_chain.length < donor_length,
            "tiled": target_chain.length > donor_length,
            "translation_per_tile": list(translation),
        })
    pdb_lines.append("END")
    target_plddt = [donor.plddt[index] for index in global_source_indices]
    target_pae = [
        [donor.pae[row_index][column_index]
         for column_index in global_source_indices]
        for row_index in global_source_indices
    ]
    scores: dict[str, Any] = {
        "plddt": target_plddt,
        "max_pae": donor.max_pae,
        "pae": target_pae,
        "ptm": donor.ptm,
    }
    if len(target.chains) > 1:
        scores["iptm"] = donor.iptm
    provenance = {
        "model_id": target.model_id,
        "donor_id": donor.donor_id,
        "total_residue_count": len(global_source_indices),
        "chains": chain_provenance,
    }
    return pdb_lines, scores, provenance


def validate_asset_pair(
    pdb_path: Path, metadata_path: Path, target: TargetModel
) -> None:
    """Strongly validate a generated target pair."""
    structure = gemmi.read_structure(str(pdb_path))
    if len(structure) != 1:
        raise AssertionError(f"{pdb_path} must contain one model.")
    observed_chains = list(structure[0])
    if [chain.name for chain in observed_chains] != [
        chain.chain_id for chain in target.chains
    ]:
        raise AssertionError(f"Unexpected chain IDs in {pdb_path}.")
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    plddt = payload.get("plddt")
    pae = payload.get("pae")
    expected_total = sum(chain.length for chain in target.chains)
    if not isinstance(plddt, list) or len(plddt) != expected_total:
        raise AssertionError(f"pLDDT length mismatch for {target.model_id}.")
    if (
        not isinstance(pae, list) or len(pae) != expected_total
        or any(not isinstance(row, list) or len(row) != expected_total
               for row in pae)
    ):
        raise AssertionError(f"PAE dimensions mismatch for {target.model_id}.")
    observed_scores: list[float] = []
    for expected_chain, observed_chain in zip(
        target.chains, observed_chains
    ):
        residues = [
            residue for residue in observed_chain
            if residue.het_flag == "A"
        ]
        if len(residues) != expected_chain.length:
            raise AssertionError(
                f"Residue count mismatch in {target.model_id} "
                f"chain {expected_chain.chain_id}."
            )
        if [residue.seqid.num for residue in residues] != list(
            range(1, expected_chain.length + 1)
        ):
            raise AssertionError("Synthetic residue numbering is not 1..N.")
        for residue in residues:
            if residue.name != "ALA":
                raise AssertionError("Synthetic residues must all be ALA.")
            atom_names = [atom.name.strip() for atom in residue]
            if atom_names != list(BACKBONE_ATOMS):
                raise AssertionError(
                    f"Backbone atoms are incomplete or reordered: {atom_names}"
                )
            values = [
                coordinate
                for atom in residue
                for coordinate in (atom.pos.x, atom.pos.y, atom.pos.z)
            ]
            if not all(math.isfinite(value) for value in values):
                raise AssertionError("Synthetic coordinates are not finite.")
            if any(abs(atom.b_iso - residue[0].b_iso) > 0.01
                   for atom in residue):
                raise AssertionError("Residue B-factors are inconsistent.")
            observed_scores.append(residue[0].b_iso)
    for index, (observed, expected) in enumerate(
        zip(observed_scores, plddt)
    ):
        if not math.isfinite(float(expected)):
            raise AssertionError("Synthetic pLDDT contains non-finite values.")
        if abs(observed - float(expected)) > 0.011:
            raise AssertionError(
                f"B-factor/pLDDT mismatch at residue {index + 1}."
            )
    all_pae = [float(value) for row in pae for value in row]
    if not all(math.isfinite(value) for value in all_pae):
        raise AssertionError("Synthetic PAE contains non-finite values.")
    if float(payload["max_pae"]) + 1e-6 < max(all_pae):
        raise AssertionError("Synthetic max_pae does not cover the matrix.")


def write_synthetic_assets(
    donor_dir: Path,
    manifest_path: Path,
    seed_path: Path,
    output_dir: Path,
    provenance_path: Path,
) -> dict[str, Any]:
    donors = {
        1: load_donor(donor_dir, MONOMER_DONOR_ID, 1),
        2: load_donor(donor_dir, DIMER_DONOR_ID, 2),
    }
    targets = load_targets(manifest_path, seed_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError(f"Synthetic asset output directory is not empty: {output_dir}")
    model_records: list[dict[str, Any]] = []
    for target in targets:
        donor = donors[len(target.chains)]
        pdb_lines, scores, record = synthesize_model(target, donor)
        pdb_path = output_dir / f"{target.model_id}-model_v1.pdb"
        metadata_path = output_dir / f"{target.model_id}-meta_v1.json"
        pdb_path.write_text("\n".join(pdb_lines) + "\n", encoding="ascii")
        metadata_path.write_text(
            json.dumps(scores, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        validate_asset_pair(pdb_path, metadata_path, target)
        record["outputs"] = {
            "pdb": str(pdb_path),
            "pdb_sha256": sha256(pdb_path),
            "metadata": str(metadata_path),
            "metadata_sha256": sha256(metadata_path),
        }
        model_records.append(record)
    provenance = {
        "warning": WARNING,
        "algorithm": ALGORITHM,
        "source_manifest": str(manifest_path.resolve()),
        "source_manifest_sha256": sha256(manifest_path),
        "source_seed": str(seed_path.resolve()),
        "source_seed_sha256": sha256(seed_path),
        "donors": {
            donor.donor_id: {
                "pdb": str(donor.pdb_path.resolve()),
                "pdb_sha256": sha256(donor.pdb_path),
                "metadata": str(donor.metadata_path.resolve()),
                "metadata_sha256": sha256(donor.metadata_path),
                "chain_lengths": [len(chain) for chain in donor.chains],
            }
            for donor in donors.values()
        },
        "models": model_records,
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return provenance


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-dir", type=Path, default=DEFAULT_DONOR_DIR)
    parser.add_argument(
        "--manifest", type=Path,
        default=DEFAULT_EXAMPLE_DIR / "config/canonical_input_manifest.csv",
    )
    parser.add_argument(
        "--seed", type=Path,
        default=DEFAULT_EXAMPLE_DIR / "config/mock_uniprot_seed.json",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--provenance", type=Path,
        default=DEFAULT_OUTPUT_DIR.parent / "reports/synthetic_provenance.json",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        provenance = write_synthetic_assets(
            args.donor_dir.resolve(), args.manifest.resolve(),
            args.seed.resolve(), args.output_dir.resolve(),
            args.provenance.resolve(),
        )
    except (OSError, ValueError, AssertionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(
        f"Generated {len(provenance['models'])} deterministic synthetic "
        f"asset pairs in {args.output_dir.resolve()}"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
