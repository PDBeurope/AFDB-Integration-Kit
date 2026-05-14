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
GPU batching utilities for protein structure analysis.

Creates GPU-resident tensors with coordinates and metadata for efficient
batched computation of clashes, interfaces, and other metrics.

Optimized for high-throughput processing of millions of structures.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Iterator

import numpy as np
import torch

try:
    from .protein import Protein, MAX_HEAVY_ATOMS
except ImportError:
    from protein import Protein, MAX_HEAVY_ATOMS


@dataclass
class ProteinBatch:
    """
    GPU-resident batch of protein structures with metadata.
    """
    coords: torch.Tensor          # (B, max_res, max_atoms, 3)
    mask: torch.Tensor            # (B, max_res, max_atoms)
    ca_coords: torch.Tensor       # (B, max_res, 3)
    ca_mask: torch.Tensor         # (B, max_res)
    chain_ids: torch.Tensor       # (B, max_res)
    res_ids: torch.Tensor         # (B, max_res)
    vdw_radii: torch.Tensor       # (B, max_res, max_atoms)
    batch_indices: List[int]
    proteins: List[Protein]
    device: str

    @property
    def batch_size(self) -> int:
        return self.coords.shape[0]

    @property
    def max_residues(self) -> int:
        return self.coords.shape[1]

    @property
    def max_atoms(self) -> int:
        return self.coords.shape[2]


# Van der Waals radii
VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "default": 1.70,
}

# Sidechain heteroatom positions per residue type
# Format: {res_name: [(slot_index, vdw_radius), ...]}
_SIDECHAIN_HETEROATOMS = {
    "ARG": [(7, 1.55), (9, 1.55), (10, 1.55)],   # NE, NH1, NH2
    "ASN": [(6, 1.52), (7, 1.55)],                # OD1, ND2
    "ASP": [(6, 1.52), (7, 1.52)],                # OD1, OD2
    "CYS": [(5, 1.80)],                           # SG
    "GLN": [(7, 1.52), (8, 1.55)],                # OE1, NE2
    "GLU": [(7, 1.52), (8, 1.52)],                # OE1, OE2
    "HIS": [(6, 1.55), (9, 1.55)],                # ND1, NE2
    "LYS": [(8, 1.55)],                           # NZ
    "MET": [(6, 1.80)],                           # SD
    "SER": [(5, 1.52)],                           # OG
    "THR": [(5, 1.52)],                           # OG1
    "TRP": [(8, 1.55)],                           # NE1
    "TYR": [(11, 1.52)],                          # OH
}


def create_batch(
    proteins: Sequence[Protein],
    indices: List[int],
    device: str = "cuda",
) -> ProteinBatch:
    """
    Create a GPU-resident batch from proteins.

    Optimized using array slicing instead of loops.
    """
    B = len(proteins)
    if B == 0:
        raise ValueError("Cannot create batch from empty protein list")

    # Get dimensions
    max_res = max(p.n_residues for p in proteins)
    max_atoms = max(p.coords.shape[1] for p in proteins)

    if max_res == 0:
        raise ValueError("All proteins have 0 residues")

    # Allocate arrays
    coords_np = np.zeros((B, max_res, max_atoms, 3), dtype=np.float32)
    mask_np = np.zeros((B, max_res, max_atoms), dtype=bool)
    res_ids_np = np.zeros((B, max_res), dtype=np.int32)
    chain_ids_np = np.zeros((B, max_res), dtype=np.int32)
    res_names_np = np.empty((B, max_res), dtype='<U3')
    res_names_np[:] = ''

    # Pre-fill VDW with carbon (most common)
    vdw_np = np.full((B, max_res, max_atoms), 1.70, dtype=np.float32)

    # Set backbone heteroatoms (always at fixed positions)
    vdw_np[:, :, 0] = 1.55  # N at slot 0
    vdw_np[:, :, 3] = 1.52  # O at slot 3

    # Build chain ID mapping
    all_chains = sorted(set(c for p in proteins for c in np.unique(p.chain_ids)))
    chain_to_int = {c: i for i, c in enumerate(all_chains)}

    # Fill per-protein data (single loop, array ops inside)
    for i, p in enumerate(proteins):
        nr = p.n_residues
        na = min(p.coords.shape[1], max_atoms)

        # Array slicing - fast
        coords_np[i, :nr, :na] = p.coords[:, :na]
        mask_np[i, :nr, :na] = p.mask[:, :na]
        res_ids_np[i, :nr] = p.res_ids
        res_names_np[i, :nr] = p.res_names

        # Chain IDs - vectorized with np.vectorize
        chain_ids_np[i, :nr] = np.array([chain_to_int[c] for c in p.chain_ids], dtype=np.int32)

    # Set sidechain heteroatoms using residue type masks
    for res_type, heteroatoms in _SIDECHAIN_HETEROATOMS.items():
        # Find all (protein, residue) positions with this residue type
        res_mask = res_names_np == res_type
        if not res_mask.any():
            continue

        # Set each heteroatom position
        for slot, vdw_val in heteroatoms:
            if slot < max_atoms:
                vdw_np[res_mask, slot] = vdw_val

    # Extract CA coords/mask (slot 1)
    ca_coords_np = coords_np[:, :, 1, :]
    ca_mask_np = mask_np[:, :, 1]

    # Transfer to GPU
    if device == "cuda" and torch.cuda.is_available():
        coords = torch.from_numpy(coords_np).pin_memory().to(device, non_blocking=True)
        mask = torch.from_numpy(mask_np).pin_memory().to(device, non_blocking=True)
        ca_coords = torch.from_numpy(ca_coords_np.copy()).pin_memory().to(device, non_blocking=True)
        ca_mask = torch.from_numpy(ca_mask_np.copy()).pin_memory().to(device, non_blocking=True)
        chain_ids = torch.from_numpy(chain_ids_np).pin_memory().to(device, non_blocking=True)
        res_ids = torch.from_numpy(res_ids_np).pin_memory().to(device, non_blocking=True)
        vdw_radii = torch.from_numpy(vdw_np).pin_memory().to(device, non_blocking=True)
    else:
        coords = torch.from_numpy(coords_np).to(device)
        mask = torch.from_numpy(mask_np).to(device)
        ca_coords = torch.from_numpy(ca_coords_np.copy()).to(device)
        ca_mask = torch.from_numpy(ca_mask_np.copy()).to(device)
        chain_ids = torch.from_numpy(chain_ids_np).to(device)
        res_ids = torch.from_numpy(res_ids_np).to(device)
        vdw_radii = torch.from_numpy(vdw_np).to(device)

    return ProteinBatch(
        coords=coords,
        mask=mask,
        ca_coords=ca_coords,
        ca_mask=ca_mask,
        chain_ids=chain_ids,
        res_ids=res_ids,
        vdw_radii=vdw_radii,
        batch_indices=list(indices),
        proteins=list(proteins),
        device=device,
    )


def iter_batches(
    proteins: Sequence[Protein],
    batch_size: int = 32,
    device: str = "cuda",
) -> Iterator[ProteinBatch]:
    """Iterate over proteins in batches."""
    n = len(proteins)
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch_proteins = proteins[start:end]
        indices = list(range(start, end))
        yield create_batch(batch_proteins, indices, device)
