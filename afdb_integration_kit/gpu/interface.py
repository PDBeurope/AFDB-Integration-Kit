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
Interface residue detection using CA distances.

A residue is an interface residue if it has CA within contact_cutoff of
a CA from a different chain.
"""
from __future__ import annotations

from typing import List, Tuple

from ._runtime import require_torch
from .clashes import radius_graph

torch = require_torch("GPU clash/interface analysis")

try:
    from .batch import ProteinBatch
except ImportError:
    from batch import ProteinBatch


@torch.no_grad()
def compute_interface_residues(
    batch: ProteinBatch,
    contact_cutoff: float = 8.0,
    max_num_neighbors: int = 64,
) -> List[List[Tuple[int, str]]]:
    """
    Compute interface residues for a batch of proteins.

    Interface residues are those with CA within contact_cutoff of a CA
    from a different chain.

    Args:
        batch: ProteinBatch with CA coordinates and chain IDs on GPU
        contact_cutoff: Maximum CA-CA distance to consider a contact (default 8.0 Å)
        max_num_neighbors: Maximum neighbors per atom for radius_graph

    Returns:
        List of lists, one per protein in batch. Each inner list contains
        (res_id, chain_id) tuples for interface residues.
    """
    B = batch.batch_size
    device = batch.device

    # Flatten CA coordinates
    ca_flat = batch.ca_coords.reshape(B * batch.max_residues, 3)
    ca_mask_flat = batch.ca_mask.reshape(B * batch.max_residues)
    chain_flat = batch.chain_ids.reshape(B * batch.max_residues)

    # Protein index for each residue
    protein_idx = torch.arange(B, device=device).repeat_interleave(batch.max_residues)
    # Residue index within protein
    res_idx = torch.arange(batch.max_residues, device=device).repeat(B)

    # Filter to valid CA atoms
    valid = ca_mask_flat
    ca_valid = ca_flat[valid]
    chain_valid = chain_flat[valid]
    protein_valid = protein_idx[valid]
    res_idx_valid = res_idx[valid]

    if ca_valid.shape[0] == 0:
        return [[] for _ in range(B)]

    # Find neighbors within cutoff
    edge_index = radius_graph(
        ca_valid,
        r=contact_cutoff,
        batch=protein_valid,
        loop=False,
        max_num_neighbors=max_num_neighbors,
    )

    if edge_index.shape[1] == 0:
        return [[] for _ in range(B)]

    src, dst = edge_index[0], edge_index[1]

    # Filter to inter-chain contacts (different chain IDs)
    chain_src = chain_valid[src]
    chain_dst = chain_valid[dst]
    inter_chain = chain_src != chain_dst

    src = src[inter_chain]
    dst = dst[inter_chain]

    if src.shape[0] == 0:
        return [[] for _ in range(B)]

    # Get protein indices for contacts
    protein_src = protein_valid[src]

    # Both src and dst residues are interface residues
    # Combine them and get unique per protein
    interface_protein = torch.cat([protein_src, protein_valid[dst]])
    interface_res_idx = torch.cat([res_idx_valid[src], res_idx_valid[dst]])

    # Move to CPU for output construction
    interface_protein = interface_protein.cpu().numpy()
    interface_res_idx = interface_res_idx.cpu().numpy()

    # Build output: unique (res_id, chain_id) per protein
    results: List[List[Tuple[int, str]]] = [[] for _ in range(B)]

    # Group by protein
    for i in range(B):
        protein_mask = interface_protein == i
        if not protein_mask.any():
            continue

        # Get unique residue indices for this protein
        res_indices = set(interface_res_idx[protein_mask])

        # Map back to original res_id and chain_id
        p = batch.proteins[i]
        for r in res_indices:
            if r < p.n_residues:
                res_id = int(p.res_ids[r])
                chain_id = str(p.chain_ids[r])
                results[i].append((res_id, chain_id))

        # Sort by chain, then res_id
        results[i].sort(key=lambda x: (x[1], x[0]))

    return results


@torch.no_grad()
def compute_interface_residues_flat(
    ca_coords: torch.Tensor,       # (B, N, 3)
    ca_mask: torch.Tensor,         # (B, N)
    chain_ids: torch.Tensor,       # (B, N)
    contact_cutoff: float = 8.0,
    max_num_neighbors: int = 64,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Compute interface residue masks (pure tensor version).

    Returns boolean mask of interface residues for each protein.

    Args:
        ca_coords: (B, N, 3) CA coordinates
        ca_mask: (B, N) validity mask
        chain_ids: (B, N) chain IDs as integers
        contact_cutoff: Maximum CA-CA distance
        max_num_neighbors: Maximum neighbors for radius_graph

    Returns:
        interface_mask: (B, N) bool mask of interface residues
        contact_src: (E,) source residue flat indices for all contacts
        contact_dst: (E,) destination residue flat indices for all contacts
        contact_distances: (E,) CA-CA distances in Angstroms for all contacts
    """
    B, N, _ = ca_coords.shape
    device = ca_coords.device

    # Flatten
    ca_flat = ca_coords.reshape(B * N, 3)
    ca_mask_flat = ca_mask.reshape(B * N)
    chain_flat = chain_ids.reshape(B * N)

    protein_idx = torch.arange(B, device=device).repeat_interleave(N)

    # Filter valid
    valid = ca_mask_flat
    valid_indices = torch.where(valid)[0]
    ca_valid = ca_flat[valid]
    chain_valid = chain_flat[valid]
    protein_valid = protein_idx[valid]

    # Initialize output mask
    interface_mask = torch.zeros((B, N), dtype=torch.bool, device=device)

    if ca_valid.shape[0] == 0:
        empty = torch.empty((0,), dtype=torch.long, device=device)
        empty_dist = torch.empty((0,), dtype=torch.float32, device=device)
        return interface_mask, empty, empty, empty_dist

    # Find neighbors
    edge_index = radius_graph(
        ca_valid,
        r=contact_cutoff,
        batch=protein_valid,
        loop=False,
        max_num_neighbors=max_num_neighbors,
    )

    if edge_index.shape[1] == 0:
        empty = torch.empty((0,), dtype=torch.long, device=device)
        empty_dist = torch.empty((0,), dtype=torch.float32, device=device)
        return interface_mask, empty, empty, empty_dist

    src, dst = edge_index[0], edge_index[1]

    # Filter inter-chain
    inter_chain = chain_valid[src] != chain_valid[dst]
    src = src[inter_chain]
    dst = dst[inter_chain]

    if src.shape[0] == 0:
        empty = torch.empty((0,), dtype=torch.long, device=device)
        empty_dist = torch.empty((0,), dtype=torch.float32, device=device)
        return interface_mask, empty, empty, empty_dist

    # Compute CA-CA distances for inter-chain contacts
    distances = torch.linalg.norm(ca_valid[src] - ca_valid[dst], dim=-1)

    # Map back to original flat indices
    src_orig = valid_indices[src]
    dst_orig = valid_indices[dst]

    # Mark interface residues
    src_protein = src_orig // N
    src_res = src_orig % N
    dst_protein = dst_orig // N
    dst_res = dst_orig % N

    interface_mask[src_protein, src_res] = True
    interface_mask[dst_protein, dst_res] = True

    return interface_mask, src_orig, dst_orig, distances
