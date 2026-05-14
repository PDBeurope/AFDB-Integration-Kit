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
GPU-accelerated atomic clash detection.

Uses torch_cluster.radius_graph when available (CUDA); otherwise falls back
to a pure-PyTorch implementation so clash detection works on systems where
the torch_cluster wheel cannot be loaded (e.g. older glibc).

Handles disulfide bonds: S-S pairs within ~2.5 Å are excluded from
clash detection since they represent covalent bonds, not steric clashes.
"""
from __future__ import annotations

from typing import List, Literal, Sequence

import numpy as np
import torch

try:
    from torch_cluster import radius_graph as _radius_graph_ext
    _RADIUS_GRAPH_USE_FALLBACK = False
except (ImportError, OSError):
    _radius_graph_ext = None
    _RADIUS_GRAPH_USE_FALLBACK = True

try:
    from .protein import Protein, MAX_HEAVY_ATOMS
except ImportError:
    from protein import Protein, MAX_HEAVY_ATOMS

# Chunk size for fallback radius graph to limit cdist memory
_RADIUS_GRAPH_CHUNK_SIZE = 2000


def _radius_graph_pytorch(
    x: torch.Tensor,
    r: float,
    batch: torch.Tensor | None = None,
    loop: bool = False,
    max_num_neighbors: int = 128,
) -> torch.Tensor:
    """Pure-PyTorch fallback for radius_graph when torch_cluster is unavailable."""
    device = x.device
    n = x.size(0)
    if batch is None:
        batch = torch.zeros(n, dtype=torch.long, device=device)
    batch_ids = batch.unique(sorted=True)
    src_list: List[torch.Tensor] = []
    dst_list: List[torch.Tensor] = []
    for b in batch_ids:
        mask_b = batch == b
        idx_b = mask_b.nonzero(as_tuple=True)[0]
        x_b = x[idx_b]
        n_b = x_b.size(0)
        if n_b < 2:
            continue
        chunk_size = min(_RADIUS_GRAPH_CHUNK_SIZE, n_b)
        for start in range(0, n_b, chunk_size):
            end = min(start + chunk_size, n_b)
            x_chunk = x_b[start:end]
            d = torch.cdist(x_chunk, x_b)
            if not loop:
                diag_idx = torch.arange(end - start, device=device) + start
                d[torch.arange(end - start, device=device), diag_idx] = float("inf")
            in_radius = d <= r
            for i in range(end - start):
                j_candidates = in_radius[i].nonzero(as_tuple=True)[0]
                if j_candidates.numel() == 0:
                    continue
                dists = d[i, j_candidates]
                _, order = dists.sort()
                j_keep = j_candidates[order[:max_num_neighbors]]
                global_i = idx_b[start + i]
                global_j = idx_b[j_keep]
                src_list.append(global_i.expand_as(global_j))
                dst_list.append(global_j)
    if not src_list:
        return torch.empty(2, 0, dtype=torch.long, device=device)
    return torch.stack([torch.cat(src_list), torch.cat(dst_list)], dim=0)


def radius_graph(
    x: torch.Tensor,
    r: float,
    batch: torch.Tensor | None = None,
    loop: bool = False,
    max_num_neighbors: int = 128,
) -> torch.Tensor:
    """Build edge index for pairs of points within radius r. Uses torch_cluster if available, else PyTorch fallback."""
    if _RADIUS_GRAPH_USE_FALLBACK or _radius_graph_ext is None:
        return _radius_graph_pytorch(x, r, batch=batch, loop=loop, max_num_neighbors=max_num_neighbors)
    return _radius_graph_ext(x, r, batch, loop, max_num_neighbors)


# Van der Waals radii for common protein atoms (in Angstroms)
VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "default": 1.70,
}

# Disulfide bond parameters
SULFUR_VDW = 1.80  # Used to identify sulfur atoms by VDW radius
DISULFIDE_MAX = 2.5  # Max S-S distance to consider a disulfide bond (Å)


@torch.no_grad()
def count_clashes(
    proteins: Sequence[Protein],
    selection: Literal["heavy", "backbone"] = "heavy",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    device: str = "cpu",
    batch_size: int = 32,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
) -> List[int]:
    """
    Count atomic clashes for a list of proteins.

    A clash occurs when two atoms are closer than the sum of their VDW radii
    minus a tolerance (defined by clash_cutoff).

    Args:
        proteins: List of Protein objects
        selection: Atom selection
            - "heavy": All heavy atoms
            - "backbone": Only backbone atoms (N, CA, C, O)
        clash_cutoff: Fraction of VDW overlap to consider a clash (default 0.4 = 40%)
        min_seq_sep: Minimum sequence separation between residues (default 3)
        device: Computation device ("cpu" or "cuda")
        batch_size: Number of proteins per batch
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, S-S pairs within 2.5 Å are not counted as clashes

    Returns:
        List of clash counts, one per protein
    """
    n_proteins = len(proteins)

    if n_proteins == 0:
        return []

    clash_counts = []

    # Process in batches
    for batch_start in range(0, n_proteins, batch_size):
        batch_end = min(batch_start + batch_size, n_proteins)
        batch_proteins = proteins[batch_start:batch_end]

        batch_counts = _compute_clashes_batch(
            batch_proteins,
            selection=selection,
            clash_cutoff=clash_cutoff,
            min_seq_sep=min_seq_sep,
            device=device,
            max_num_neighbors=max_num_neighbors,
            exclude_disulfides=exclude_disulfides,
        )

        clash_counts.extend(batch_counts)

    return clash_counts


@torch.no_grad()
def _compute_clashes_batch(
    proteins: Sequence[Protein],
    selection: Literal["heavy", "backbone"],
    clash_cutoff: float,
    min_seq_sep: int,
    device: str,
    max_num_neighbors: int,
    exclude_disulfides: bool = True,
) -> List[int]:
    """Compute clashes for a batch of proteins."""
    B = len(proteins)

    if B == 0:
        return []

    # Determine atom slots to use
    if selection == "backbone":
        n_slots = 4  # N, CA, C, O
    else:
        n_slots = MAX_HEAVY_ATOMS

    # Find max residue count for padding
    max_res = max(p.n_residues for p in proteins)

    if max_res == 0:
        return [0] * B

    # Pad and stack coordinates and masks
    padded_coords = np.zeros((B, max_res, n_slots, 3), dtype=np.float32)
    padded_mask = np.zeros((B, max_res, n_slots), dtype=bool)
    padded_vdw = np.full((B, max_res, n_slots), VDW_RADII["default"], dtype=np.float32)

    for i, p in enumerate(proteins):
        n_res = p.n_residues
        if n_res == 0:
            continue

        if selection == "backbone":
            padded_coords[i, :n_res] = p.coords[:, :4]
            padded_mask[i, :n_res] = p.mask[:, :4]
            # Set VDW radii based on atom type (N, CA, C, O)
            padded_vdw[i, :n_res, 0] = VDW_RADII["N"]  # N
            padded_vdw[i, :n_res, 1] = VDW_RADII["C"]  # CA
            padded_vdw[i, :n_res, 2] = VDW_RADII["C"]  # C
            padded_vdw[i, :n_res, 3] = VDW_RADII["O"]  # O
        else:
            padded_coords[i, :n_res] = p.coords
            padded_mask[i, :n_res] = p.mask
            # Set VDW radii based on atom names
            for r in range(n_res):
                for a in range(n_slots):
                    if p.mask[r, a]:
                        atom_name = p.atom_names[r, a]
                        element = atom_name[0] if atom_name else "C"
                        padded_vdw[i, r, a] = VDW_RADII.get(element, VDW_RADII["default"])

    # Move to device
    coords = torch.tensor(padded_coords, dtype=torch.float32, device=device)
    mask = torch.tensor(padded_mask, dtype=torch.bool, device=device)
    vdw = torch.tensor(padded_vdw, dtype=torch.float32, device=device)

    N = max_res
    A = n_slots

    # Flatten all data
    coords_flat = coords.reshape(B * N * A, 3)
    mask_flat = mask.reshape(B * N * A)
    vdw_flat = vdw.reshape(B * N * A)

    # Create protein, residue, and atom indices
    protein_idx = torch.arange(B, device=device).repeat_interleave(N * A)
    res_idx = torch.arange(N, device=device).repeat_interleave(A).repeat(B)
    atom_idx = torch.arange(A, device=device).repeat(B * N)

    # Filter to valid atoms
    valid = mask_flat
    coords_valid = coords_flat[valid]
    vdw_valid = vdw_flat[valid]
    prot_valid = protein_idx[valid]
    res_valid = res_idx[valid]
    atom_valid = atom_idx[valid]

    if coords_valid.shape[0] < 2:
        return [0] * B

    # Sort by protein index for radius_graph (required for correct batch handling)
    sort_idx = prot_valid.argsort(stable=True)
    coords_sorted = coords_valid[sort_idx]
    vdw_sorted = vdw_valid[sort_idx]
    prot_sorted = prot_valid[sort_idx]
    res_sorted = res_valid[sort_idx]
    atom_sorted = atom_valid[sort_idx]

    # Use radius_graph with batch indices
    max_vdw_sum = 2 * VDW_RADII["default"] * 2  # Conservative upper bound

    edge_index = radius_graph(
        coords_sorted,
        r=max_vdw_sum,
        batch=prot_sorted,
        loop=False,
        max_num_neighbors=max_num_neighbors,
    )

    if edge_index.shape[1] == 0:
        return [0] * B

    src, dst = edge_index[0], edge_index[1]

    # Get attributes for edge endpoints
    prot_src = prot_sorted[src]
    res_src = res_sorted[src]
    res_dst = res_sorted[dst]
    atom_src = atom_sorted[src]
    atom_dst = atom_sorted[dst]
    vdw_src = vdw_sorted[src]
    vdw_dst = vdw_sorted[dst]

    # Compute distances
    d = torch.linalg.norm(coords_sorted[src] - coords_sorted[dst], dim=-1)

    # Compute VDW sum and threshold
    vdw_sum = vdw_src + vdw_dst
    clash_threshold = vdw_sum * (1.0 - clash_cutoff)

    # Filter criteria:
    # 1. Avoid double counting (only keep src < dst by some ordering)
    # 2. Minimum sequence separation
    # 3. Distance less than clash threshold
    seq_sep = torch.abs(res_dst - res_src)
    order_mask = (res_src < res_dst) | ((res_src == res_dst) & (atom_src < atom_dst))
    sep_mask = seq_sep >= min_seq_sep
    clash_mask = d < clash_threshold

    keep = order_mask & sep_mask & clash_mask

    # Exclude disulfide bonds (S-S pairs within bond distance)
    # Sulfur identified by exact VDW radius match (faster than isclose)
    if exclude_disulfides:
        is_sulfur_src = vdw_src == SULFUR_VDW
        is_sulfur_dst = vdw_dst == SULFUR_VDW
        is_disulfide = is_sulfur_src & is_sulfur_dst & (d <= DISULFIDE_MAX)
        keep = keep & ~is_disulfide

    # Count clashes per protein
    prot_clash = prot_src[keep]
    clash_counts = torch.bincount(prot_clash, minlength=B).tolist()

    return clash_counts


def get_clash_pairs(
    protein: Protein,
    selection: Literal["heavy", "backbone"] = "heavy",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    device: str = "cpu",
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
) -> List[tuple]:
    """
    Get detailed clash information for a single protein.

    Args:
        protein: Protein object
        selection: "heavy" or "backbone"
        clash_cutoff: Fraction of VDW overlap
        min_seq_sep: Minimum sequence separation
        device: Computation device
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, S-S pairs within 2.5 Å are not counted as clashes

    Returns:
        List of tuples: (res_i, atom_name_i, res_j, atom_name_j, distance, overlap)
    """
    if protein.n_residues == 0:
        return []

    # Determine atom slots
    if selection == "backbone":
        n_slots = 4
        coords = protein.coords[:, :4]
        mask = protein.mask[:, :4]
        atom_names = protein.atom_names[:, :4]
    else:
        n_slots = MAX_HEAVY_ATOMS
        coords = protein.coords
        mask = protein.mask
        atom_names = protein.atom_names

    N = protein.n_residues

    # Flatten
    coords_flat = coords.reshape(-1, 3)
    mask_flat = mask.reshape(-1)

    # Create VDW radii
    vdw_flat = np.full(N * n_slots, VDW_RADII["default"], dtype=np.float32)
    for r in range(N):
        for a in range(n_slots):
            if mask[r, a]:
                atom_name = atom_names[r, a]
                element = atom_name[0] if atom_name else "C"
                vdw_flat[r * n_slots + a] = VDW_RADII.get(element, VDW_RADII["default"])

    # Filter to valid atoms
    valid_idx = np.where(mask_flat)[0]
    coords_valid = coords_flat[valid_idx]
    vdw_valid = vdw_flat[valid_idx]

    if len(coords_valid) < 2:
        return []

    # Move to device
    coords_t = torch.tensor(coords_valid, dtype=torch.float32, device=device)
    vdw_t = torch.tensor(vdw_valid, dtype=torch.float32, device=device)

    # Create index mapping
    res_idx = np.arange(N).repeat(n_slots)[valid_idx]
    atom_idx = np.tile(np.arange(n_slots), N)[valid_idx]

    # radius_graph
    max_vdw_sum = 2 * VDW_RADII["default"] * 2

    edge_index = radius_graph(
        coords_t,
        r=max_vdw_sum,
        loop=False,
        max_num_neighbors=max_num_neighbors,
    )

    if edge_index.shape[1] == 0:
        return []

    src, dst = edge_index[0].cpu().numpy(), edge_index[1].cpu().numpy()

    # Compute distances
    d = np.linalg.norm(coords_valid[src] - coords_valid[dst], axis=-1)

    # Get residue/atom indices
    ri = res_idx[src]
    rj = res_idx[dst]
    ai = atom_idx[src]
    aj = atom_idx[dst]

    # VDW sum and threshold
    vdw_sum = vdw_valid[src] + vdw_valid[dst]
    clash_threshold = vdw_sum * (1.0 - clash_cutoff)

    # Filter
    seq_sep = np.abs(rj - ri)
    order_mask = (ri < rj) | ((ri == rj) & (ai < aj))
    sep_mask = seq_sep >= min_seq_sep
    clash_mask = d < clash_threshold

    keep = order_mask & sep_mask & clash_mask

    # Exclude disulfide bonds
    if exclude_disulfides:
        is_sulfur_src = np.isclose(vdw_valid[src], SULFUR_VDW)
        is_sulfur_dst = np.isclose(vdw_valid[dst], SULFUR_VDW)
        is_disulfide = is_sulfur_src & is_sulfur_dst & (d <= DISULFIDE_MAX)
        keep = keep & ~is_disulfide

    # Build result
    clashes = []
    for idx in np.where(keep)[0]:
        res_i = int(ri[idx])
        res_j = int(rj[idx])
        atom_i = int(ai[idx])
        atom_j = int(aj[idx])
        distance = float(d[idx])
        overlap = float(vdw_sum[idx] - d[idx])

        name_i = atom_names[res_i, atom_i] if atom_names[res_i, atom_i] else f"slot{atom_i}"
        name_j = atom_names[res_j, atom_j] if atom_names[res_j, atom_j] else f"slot{atom_j}"

        clashes.append((
            protein.res_ids[res_i],
            name_i,
            protein.res_ids[res_j],
            name_j,
            distance,
            overlap,
        ))

    return clashes


@torch.no_grad()
def compute_clashes_from_batch(
    coords: torch.Tensor,          # (B, N, A, 3)
    mask: torch.Tensor,            # (B, N, A)
    vdw_radii: torch.Tensor,       # (B, N, A)
    selection: Literal["heavy", "backbone"] = "heavy",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Compute clashes from pre-batched GPU tensors.

    This is the core clash detection kernel optimized for GPU.

    Args:
        coords: (B, N, A, 3) coordinates on device
        mask: (B, N, A) atom validity mask
        vdw_radii: (B, N, A) VDW radii per atom
        selection: "heavy" or "backbone"
        clash_cutoff: Fraction of VDW overlap
        min_seq_sep: Minimum sequence separation
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, exclude S-S pairs within DISULFIDE_MAX distance

    Returns:
        protein_idx: (E,) protein index for each clash
        res_i: (E,) residue index i (within protein)
        res_j: (E,) residue index j (within protein)
        atom_i: (E,) atom slot i
        atom_j: (E,) atom slot j
        distances: (E,) inter-atom distances in Angstroms
        overlaps: (E,) VDW overlap (vdw_sum - distance) in Angstroms
    """
    device = coords.device
    B, N, A, _ = coords.shape

    # Apply selection
    if selection == "backbone":
        coords = coords[:, :, :4, :]
        mask = mask[:, :, :4]
        vdw_radii = vdw_radii[:, :, :4]
        A = 4

    # Flatten
    coords_flat = coords.reshape(B * N * A, 3)
    mask_flat = mask.reshape(B * N * A)
    vdw_flat = vdw_radii.reshape(B * N * A)

    protein_idx = torch.arange(B, device=device).repeat_interleave(N * A)
    res_idx = torch.arange(N, device=device).repeat_interleave(A).repeat(B)
    atom_idx = torch.arange(A, device=device).repeat(B * N)

    # Filter valid
    valid = mask_flat
    coords_valid = coords_flat[valid]
    vdw_valid = vdw_flat[valid]
    prot_valid = protein_idx[valid]
    res_valid = res_idx[valid]
    atom_valid = atom_idx[valid]

    if coords_valid.shape[0] < 2:
        empty = torch.empty((0,), dtype=torch.long, device=device)
        empty_f = torch.empty((0,), dtype=torch.float32, device=device)
        return empty, empty, empty, empty, empty, empty_f, empty_f

    # Sort by protein
    sort_idx = prot_valid.argsort(stable=True)
    coords_sorted = coords_valid[sort_idx]
    vdw_sorted = vdw_valid[sort_idx]
    prot_sorted = prot_valid[sort_idx]
    res_sorted = res_valid[sort_idx]
    atom_sorted = atom_valid[sort_idx]

    # radius_graph
    max_vdw_sum = 2 * 1.80 * 2  # S-S worst case

    edge_index = radius_graph(
        coords_sorted,
        r=max_vdw_sum,
        batch=prot_sorted,
        loop=False,
        max_num_neighbors=max_num_neighbors,
    )

    if edge_index.shape[1] == 0:
        empty = torch.empty((0,), dtype=torch.long, device=device)
        empty_f = torch.empty((0,), dtype=torch.float32, device=device)
        return empty, empty, empty, empty, empty, empty_f, empty_f

    src, dst = edge_index[0], edge_index[1]

    # Get attributes
    prot_src = prot_sorted[src]
    res_src = res_sorted[src]
    res_dst = res_sorted[dst]
    atom_src = atom_sorted[src]
    atom_dst = atom_sorted[dst]
    vdw_src = vdw_sorted[src]
    vdw_dst = vdw_sorted[dst]

    # Distances
    d = torch.linalg.norm(coords_sorted[src] - coords_sorted[dst], dim=-1)

    # Threshold
    vdw_sum = vdw_src + vdw_dst
    clash_threshold = vdw_sum * (1.0 - clash_cutoff)

    # Filters
    seq_sep = torch.abs(res_dst - res_src)
    order_mask = (res_src < res_dst) | ((res_src == res_dst) & (atom_src < atom_dst))
    sep_mask = seq_sep >= min_seq_sep
    clash_mask = d < clash_threshold

    keep = order_mask & sep_mask & clash_mask

    # Exclude disulfide bonds (S-S pairs within bond distance)
    # Sulfur identified by exact VDW radius match (faster than isclose)
    if exclude_disulfides:
        is_sulfur_src = vdw_src == SULFUR_VDW
        is_sulfur_dst = vdw_dst == SULFUR_VDW
        is_disulfide = is_sulfur_src & is_sulfur_dst & (d <= DISULFIDE_MAX)
        keep = keep & ~is_disulfide

    return (
        prot_src[keep].long(),
        res_src[keep].long(),
        res_dst[keep].long(),
        atom_src[keep].long(),
        atom_dst[keep].long(),
        d[keep],
        (vdw_sum - d)[keep],
    )


@torch.no_grad()
def compute_clashing_residues_from_batch(
    coords: torch.Tensor,          # (B, N, A, 3)
    mask: torch.Tensor,            # (B, N, A)
    vdw_radii: torch.Tensor,       # (B, N, A)
    selection: Literal["heavy", "backbone"] = "heavy",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Compute clash counts and clashing residue masks.

    Args:
        coords: (B, N, A, 3) coordinates on device
        mask: (B, N, A) atom validity mask
        vdw_radii: (B, N, A) VDW radii per atom
        selection: "heavy" or "backbone"
        clash_cutoff: Fraction of VDW overlap
        min_seq_sep: Minimum sequence separation
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, exclude S-S pairs within DISULFIDE_MAX distance

    Returns:
        clash_counts: (B,) number of clashing atom pairs per protein
        clashing_residue_mask: (B, N) bool mask of residues involved in clashes
    """
    device = coords.device
    B, N, _, _ = coords.shape

    prot_idx, res_i, res_j, _, _, _, _ = compute_clashes_from_batch(
        coords, mask, vdw_radii,
        selection=selection,
        clash_cutoff=clash_cutoff,
        min_seq_sep=min_seq_sep,
        max_num_neighbors=max_num_neighbors,
        exclude_disulfides=exclude_disulfides,
    )

    # Count clashes per protein
    if prot_idx.numel() == 0:
        clash_counts = torch.zeros(B, dtype=torch.long, device=device)
        clashing_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
        return clash_counts, clashing_mask

    clash_counts = torch.bincount(prot_idx, minlength=B)

    # Mark clashing residues
    clashing_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
    clashing_mask[prot_idx, res_i] = True
    clashing_mask[prot_idx, res_j] = True

    return clash_counts, clashing_mask
