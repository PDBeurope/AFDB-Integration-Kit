# Copyright 2026 Maciej Majewski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# SPDX-License-Identifier: Apache-2.0

"""
PDB parser using fastpdb with multi-CPU parallelization.

Parses PDB files into Protein objects with heavy atoms only (no hydrogens).
Supports parallel parsing of multiple files using ProcessPoolExecutor.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import fastpdb
from biotite.structure.atoms import AtomArray

try:
    from .protein import (
        Protein,
        empty_protein,
        MAX_HEAVY_ATOMS,
        AA_ORDER_HEAVY,
        AA_ALIASES,
        _CANON_AA,
        _CANON_OR_ALIAS,
        _ATOM_SLOT_LOOKUP,
    )
except ImportError:
    from protein import (
        Protein,
        empty_protein,
        MAX_HEAVY_ATOMS,
        AA_ORDER_HEAVY,
        AA_ALIASES,
        _CANON_AA,
        _CANON_OR_ALIAS,
        _ATOM_SLOT_LOOKUP,
    )


def _read_structure(path: str | Path, model: Optional[int] = 1) -> AtomArray:
    """Read PDB file using fastpdb."""
    path = Path(path)
    pdb = fastpdb.PDBFile.read(str(path))
    if model is None:
        return pdb.get_structure()
    else:
        return pdb.get_structure(model=model)


def parse_protein(
    path: str | Path,
    model: Optional[int] = 1,
) -> Protein:
    """
    Parse a PDB file into a Protein object.

    Args:
        path: Path to PDB file
        model: Model number (1-indexed) for multi-model files.
               Use None to load all models concatenated.

    Returns:
        Protein object with heavy atom coordinates
    """
    path = Path(path)

    try:
        arr = _read_structure(path, model=model)
    except Exception as e:
        print(f"Warning: Could not parse {path}: {e}")
        return empty_protein(str(path))

    if arr.array_length() == 0:
        return empty_protein(str(path))

    return _parse_atomarray(arr, path)


def _parse_atomarray(arr: AtomArray, path: Path) -> Protein:
    """Convert biotite AtomArray to Protein object."""
    cats = set(arr.get_annotation_categories())
    n_atoms_total = arr.array_length()

    # Get raw arrays
    res_name_raw = arr.res_name
    atom_name_all = arr.atom_name
    chain_all = arr.chain_id
    res_id_all = arr.res_id.astype(np.int32, copy=False)
    coord_all = arr.coord
    ins_all = arr.ins_code if "ins_code" in cats else np.full(n_atoms_total, "", dtype="<U1")

    # Filter to protein residues only
    keep_mask = np.isin(res_name_raw, list(_CANON_OR_ALIAS))
    if not keep_mask.any():
        return empty_protein(str(path))

    arr = arr[keep_mask]
    n_atoms_total = arr.array_length()
    res_name_raw = arr.res_name
    atom_name_all = arr.atom_name
    chain_all = arr.chain_id
    res_id_all = arr.res_id.astype(np.int32, copy=False)
    coord_all = arr.coord
    ins_all = arr.ins_code if "ins_code" in cats else np.full(n_atoms_total, "", dtype="<U1")

    # Handle altloc - keep only first conformer
    if "altloc" in cats:
        alt = arr.altloc
        keep = (alt == "") | (alt == ".") | (alt == "A")
        arr = arr[keep]
        n_atoms_total = arr.array_length()
        res_name_raw = arr.res_name
        atom_name_all = arr.atom_name
        chain_all = arr.chain_id
        res_id_all = arr.res_id.astype(np.int32, copy=False)
        coord_all = arr.coord
        ins_all = arr.ins_code if "ins_code" in cats else np.full(n_atoms_total, "", dtype="<U1")

    if n_atoms_total == 0:
        return empty_protein(str(path))

    # Filter out hydrogens
    if "element" in cats:
        is_hydrogen = arr.element == "H"
    else:
        is_hydrogen = np.char.startswith(atom_name_all.astype(str), "H")

    keep_heavy = ~is_hydrogen
    res_name_raw = res_name_raw[keep_heavy]
    atom_name_all = atom_name_all[keep_heavy]
    chain_all = chain_all[keep_heavy]
    res_id_all = res_id_all[keep_heavy]
    coord_all = coord_all[keep_heavy]
    ins_all = ins_all[keep_heavy]
    n_atoms_total = len(res_name_raw)

    if n_atoms_total == 0:
        return empty_protein(str(path))

    # Map aliases to canonical names
    res_name_mapped = res_name_raw.copy()
    for alias, canonical in AA_ALIASES.items():
        mask_alias = res_name_raw == alias
        if mask_alias.any():
            res_name_mapped[mask_alias] = canonical

    # Find residue boundaries
    key_change = np.ones(n_atoms_total, dtype=bool)
    key_change[1:] = (
        (chain_all[1:] != chain_all[:-1]) |
        (res_id_all[1:] != res_id_all[:-1]) |
        (ins_all[1:] != ins_all[:-1])
    )
    starts = np.flatnonzero(key_change)
    ends = np.concatenate([starts[1:], [n_atoms_total]])
    n_res = len(starts)

    # Allocate output arrays
    coords = np.zeros((n_res, MAX_HEAVY_ATOMS, 3), dtype=np.float32)
    mask = np.zeros((n_res, MAX_HEAVY_ATOMS), dtype=bool)
    atom_names = np.full((n_res, MAX_HEAVY_ATOMS), "", dtype="<U4")
    res_names = np.empty(n_res, dtype="<U3")
    chain_ids = np.empty(n_res, dtype="<U1")
    res_ids = np.empty(n_res, dtype=np.int32)

    # Fill residue data
    for i, (s, e) in enumerate(zip(starts, ends)):
        res_name = res_name_mapped[s]
        res_names[i] = res_name
        chain_ids[i] = chain_all[s]
        res_ids[i] = res_id_all[s]

        if res_name not in _ATOM_SLOT_LOOKUP:
            continue

        slot_lookup = _ATOM_SLOT_LOOKUP[res_name]

        for j in range(s, e):
            atom_name = atom_name_all[j]
            if atom_name in slot_lookup:
                slot = slot_lookup[atom_name]
                coords[i, slot] = coord_all[j]
                mask[i, slot] = True
                atom_names[i, slot] = atom_name

    return Protein(
        coords=coords,
        mask=mask,
        atom_names=atom_names,
        res_names=res_names,
        chain_ids=chain_ids,
        res_ids=res_ids,
        path=str(path),
    )


def parse_proteins(
    paths: Sequence[str | Path],
    model: Optional[int] = 1,
    n_jobs: int = 1,
) -> List[Protein]:
    """
    Parse multiple PDB files, optionally in parallel.

    Args:
        paths: List of paths to PDB files
        model: Model number for multi-model files
        n_jobs: Number of parallel workers. Use -1 for all CPUs.

    Returns:
        List of Protein objects (same order as input paths)
    """
    paths = [Path(p) for p in paths]

    if n_jobs == 1:
        # Sequential processing
        return [parse_protein(p, model=model) for p in paths]

    # Parallel processing
    if n_jobs == -1:
        import os
        n_jobs = os.cpu_count() or 1

    # Use ProcessPoolExecutor for parallel parsing
    proteins = [None] * len(paths)
    parse_fn = partial(parse_protein, model=model)

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Submit all jobs
        future_to_idx = {
            executor.submit(parse_fn, path): i
            for i, path in enumerate(paths)
        }

        # Collect results
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                proteins[idx] = future.result()
            except Exception as e:
                print(f"Warning: Failed to parse {paths[idx]}: {e}")
                proteins[idx] = empty_protein(str(paths[idx]))

    return proteins
