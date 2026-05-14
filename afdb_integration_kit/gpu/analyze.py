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
Master analysis script for protein structure quality metrics.

Computes clashes (backbone and all atoms) and interface residues
in a batched, GPU-accelerated pipeline.

Designed for high-throughput processing of millions of structures.

Pipeline architecture:
  [Parser Pool] --> [Parse Queue] --> [GPU] --> [Write Queue] --> [Writer Pool]

All three stages run concurrently for maximum throughput.
"""
from __future__ import annotations

import json
import queue
import threading
from concurrent.futures import ThreadPoolExecutor, Future
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import List, Optional, Sequence, Set, Union

import torch
from tqdm import tqdm

try:
    from .protein import Protein
    from .parse import parse_protein, parse_proteins
    from .batch import ProteinBatch, create_batch, iter_batches
    from .clashes import compute_clashes_from_batch
    from .interface import compute_interface_residues_flat
    from .schema import result_to_interface_schema, result_to_clash_schema
except ImportError:
    from protein import Protein
    from parse import parse_protein, parse_proteins
    from batch import ProteinBatch, create_batch, iter_batches
    from clashes import compute_clashes_from_batch
    from interface import compute_interface_residues_flat
    from schema import result_to_interface_schema, result_to_clash_schema


@dataclass
class InterfaceContact:
    """A single inter-chain CA-CA contact."""
    chain_1: str
    res_1_id: int
    aa_1_type: str
    chain_2: str
    res_2_id: int
    aa_2_type: str
    distance: float


@dataclass
class ClashContact:
    """A single atomic clash between two residues."""
    chain_1: str
    res_1_id: int
    aa_1_type: str
    atom_1_name: str
    chain_2: str
    res_2_id: int
    aa_2_type: str
    atom_2_name: str
    distance: float
    overlap: float  # vdw_sum - distance (Angstroms)


@dataclass
class ProteinAnalysisResult:
    """Analysis results for a single protein."""
    # Identification
    path: str
    n_residues: int
    n_atoms: int

    # Clash counts
    n_backbone_clashes: int
    n_heavy_clashes: int

    # Clashing residues (list of {"res_id": int, "chain_id": str, "aa_type": str})
    backbone_clashing_residues: List[dict]
    heavy_clashing_residues: List[dict]

    # Interface residues (list of {"res_id": int, "chain_id": str, "aa_type": str})
    interface_residues: List[dict]

    # Inter-chain contacts with distances
    interface_contacts: List[InterfaceContact] = field(default_factory=list)

    # Clash contacts with distances and overlaps
    backbone_clash_contacts: List[ClashContact] = field(default_factory=list)
    heavy_clash_contacts: List[ClashContact] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _materialize_contacts(result: ProteinAnalysisResult) -> None:
    """Build contact lists from raw arrays attached by analyze_batch.

    Called in the writer thread so that dataclass construction does not
    block the GPU thread between batches.
    """
    if not hasattr(result, '_protein'):
        return  # already materialized or built directly

    p = result._protein
    n_res = p.n_residues

    # --- Interface contacts ---
    src, dst, dist = result._iface_src, result._iface_dst, result._iface_dist
    if len(src) > 0:
        c1, r1, a1 = p.chain_ids[src], p.res_ids[src], p.res_names[src]
        c2, r2, a2 = p.chain_ids[dst], p.res_ids[dst], p.res_names[dst]
        result.interface_contacts = [
            InterfaceContact(
                chain_1=str(c1[k]), res_1_id=int(r1[k]), aa_1_type=str(a1[k]),
                chain_2=str(c2[k]), res_2_id=int(r2[k]), aa_2_type=str(a2[k]),
                distance=float(dist[k]),
            )
            for k in range(len(src))
        ]

    # --- Clash contacts ---
    result.backbone_clash_contacts = _clash_list_from_arrays(p, result._bb_clash)
    result.heavy_clash_contacts = _clash_list_from_arrays(p, result._hv_clash)

    # Free raw data (no longer needed)
    del result._protein
    del result._iface_src, result._iface_dst, result._iface_dist
    del result._bb_clash, result._hv_clash


def _clash_list_from_arrays(
    p: Protein, arrays: tuple,
) -> List[ClashContact]:
    """Build ClashContact list from per-protein clash arrays."""
    ri, rj, ai, aj, dists, overlaps = arrays
    if len(ri) == 0:
        return []

    valid = (ri < p.n_residues) & (rj < p.n_residues)
    if not valid.all():
        ri, rj, ai, aj = ri[valid], rj[valid], ai[valid], aj[valid]
        dists, overlaps = dists[valid], overlaps[valid]

    if len(ri) == 0:
        return []

    chains_i, res_ids_i, aa_i = p.chain_ids[ri], p.res_ids[ri], p.res_names[ri]
    anames_i = p.atom_names[ri, ai]
    chains_j, res_ids_j, aa_j = p.chain_ids[rj], p.res_ids[rj], p.res_names[rj]
    anames_j = p.atom_names[rj, aj]

    return [
        ClashContact(
            chain_1=str(chains_i[k]),
            res_1_id=int(res_ids_i[k]),
            aa_1_type=str(aa_i[k]),
            atom_1_name=str(anames_i[k]) if anames_i[k] else f"slot{ai[k]}",
            chain_2=str(chains_j[k]),
            res_2_id=int(res_ids_j[k]),
            aa_2_type=str(aa_j[k]),
            atom_2_name=str(anames_j[k]) if anames_j[k] else f"slot{aj[k]}",
            distance=float(dists[k]),
            overlap=float(overlaps[k]),
        )
        for k in range(len(ri))
    ]


def _write_result_json(
    result: ProteinAnalysisResult,
    output_dir: Path,
    interface_cutoff: float = 8.0,
) -> List[str]:
    """Write interface and clash schema JSONs as separate files.

    Produces up to two files per structure:
      - {stem}_interface.json  (if interface contacts exist)
      - {stem}_clashes.json    (if clashes exist)

    Returns list of output paths written.
    """
    # Materialize contact lists from raw arrays (deferred from GPU thread)
    _materialize_contacts(result)

    stem = Path(result.path).stem
    paths: List[str] = []

    # Interface file
    iface = result_to_interface_schema(result, interface_cutoff=interface_cutoff)
    if iface is not None:
        p = output_dir / f"{stem}_interface.json"
        with open(p, "w") as f:
            json.dump(iface, f)
        paths.append(str(p))

    # Clashes file (always written — includes both sites with n_clashes=0 when clean)
    clash = result_to_clash_schema(result)
    p = output_dir / f"{stem}_clashes.json"
    with open(p, "w") as f:
        json.dump(clash, f)
    paths.append(str(p))

    return paths


def _write_batch_results(
    results: List[ProteinAnalysisResult],
    output_dir: Path,
    interface_cutoff: float = 8.0,
) -> List[str]:
    """Write a batch of results to individual JSON files."""
    paths = []
    for result in results:
        paths.extend(_write_result_json(result, output_dir, interface_cutoff=interface_cutoff))
    return paths


ALL_ANALYSES = {"interface", "backbone_clashes", "heavy_atom_clashes"}


@torch.no_grad()
def analyze_batch(
    batch: ProteinBatch,
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    interface_cutoff: float = 8.0,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
    analyses: Optional[Set[str]] = None,
) -> List[ProteinAnalysisResult]:
    """
    Analyze a batch of proteins on GPU.

    Computes (controlled by *analyses*):
    - Backbone clashes (N, CA, C, O atoms)
    - Heavy atom clashes (all non-hydrogen atoms, excluding disulfide S-S bonds)
    - Interface residues (CA-CA distance < interface_cutoff between chains)

    Args:
        analyses: Set of analyses to run. Valid values are
            ``"interface"``, ``"backbone_clashes"``, ``"heavy_atom_clashes"``.
            ``None`` (default) runs all three.
    """
    if analyses is None:
        analyses = ALL_ANALYSES

    B = batch.batch_size
    N = batch.max_residues
    device = batch.coords.device
    _empty = torch.empty(0, dtype=torch.long, device=device)
    _empty_f = torch.empty(0, dtype=torch.float, device=device)

    # ===== BACKBONE CLASHES =====
    if "backbone_clashes" in analyses:
        # Note: backbone doesn't include S, so exclude_disulfides has no effect
        bb_prot, bb_ri, bb_rj, bb_ai, bb_aj, bb_dists, bb_overlaps = compute_clashes_from_batch(
            batch.coords,
            batch.mask,
            batch.vdw_radii,
            selection="backbone",
            clash_cutoff=clash_cutoff,
            min_seq_sep=min_seq_sep,
            max_num_neighbors=max_num_neighbors,
            exclude_disulfides=False,  # No S in backbone
        )

        # Compute backbone counts and mask
        if bb_prot.numel() > 0:
            bb_counts = torch.bincount(bb_prot, minlength=B)
            bb_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
            bb_mask[bb_prot, bb_ri] = True
            bb_mask[bb_prot, bb_rj] = True
        else:
            bb_counts = torch.zeros(B, dtype=torch.long, device=device)
            bb_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
    else:
        bb_counts = torch.zeros(B, dtype=torch.long, device=device)
        bb_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
        bb_prot = bb_ri = bb_rj = bb_ai = bb_aj = _empty
        bb_dists = bb_overlaps = _empty_f

    # ===== HEAVY ATOM CLASHES =====
    if "heavy_atom_clashes" in analyses:
        hv_prot, hv_ri, hv_rj, hv_ai, hv_aj, hv_dists, hv_overlaps = compute_clashes_from_batch(
            batch.coords,
            batch.mask,
            batch.vdw_radii,
            selection="heavy",
            clash_cutoff=clash_cutoff,
            min_seq_sep=min_seq_sep,
            max_num_neighbors=max_num_neighbors,
            exclude_disulfides=exclude_disulfides,
        )

        # Compute heavy counts and mask
        if hv_prot.numel() > 0:
            heavy_counts = torch.bincount(hv_prot, minlength=B)
            heavy_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
            heavy_mask[hv_prot, hv_ri] = True
            heavy_mask[hv_prot, hv_rj] = True
        else:
            heavy_counts = torch.zeros(B, dtype=torch.long, device=device)
            heavy_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
    else:
        heavy_counts = torch.zeros(B, dtype=torch.long, device=device)
        heavy_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
        hv_prot = hv_ri = hv_rj = hv_ai = hv_aj = _empty
        hv_dists = hv_overlaps = _empty_f

    # ===== INTERFACE RESIDUES =====
    if "interface" in analyses:
        interface_mask, contact_src, contact_dst, contact_dist = compute_interface_residues_flat(
            batch.ca_coords,
            batch.ca_mask,
            batch.chain_ids,
            contact_cutoff=interface_cutoff,
            max_num_neighbors=max_num_neighbors,
        )
    else:
        interface_mask = torch.zeros((B, N), dtype=torch.bool, device=device)
        contact_src = contact_dst = _empty
        contact_dist = _empty_f

    # ===== MOVE TO CPU =====
    bb_counts_np = bb_counts.cpu().numpy()
    bb_mask_np = bb_mask.cpu().numpy()
    bb_prot_np = bb_prot.cpu().numpy()
    bb_ri_np = bb_ri.cpu().numpy()
    bb_rj_np = bb_rj.cpu().numpy()
    bb_ai_np = bb_ai.cpu().numpy()
    bb_aj_np = bb_aj.cpu().numpy()
    bb_dists_np = bb_dists.cpu().numpy()
    bb_overlaps_np = bb_overlaps.cpu().numpy()

    heavy_counts_np = heavy_counts.cpu().numpy()
    heavy_mask_np = heavy_mask.cpu().numpy()
    hv_prot_np = hv_prot.cpu().numpy()
    hv_ri_np = hv_ri.cpu().numpy()
    hv_rj_np = hv_rj.cpu().numpy()
    hv_ai_np = hv_ai.cpu().numpy()
    hv_aj_np = hv_aj.cpu().numpy()
    hv_dists_np = hv_dists.cpu().numpy()
    hv_overlaps_np = hv_overlaps.cpu().numpy()

    interface_mask_np = interface_mask.cpu().numpy()
    contact_src_np = contact_src.cpu().numpy()
    contact_dst_np = contact_dst.cpu().numpy()
    contact_dist_np = contact_dist.cpu().numpy()

    # ===== BUILD RESULTS =====
    results = []

    for i in range(B):
        p = batch.proteins[i]
        n_res = p.n_residues

        # Extract clashing residue indices via numpy (vectorized)
        bb_res_indices = bb_mask_np[i, :n_res].nonzero()[0]
        heavy_res_indices = heavy_mask_np[i, :n_res].nonzero()[0]
        interface_indices = interface_mask_np[i, :n_res].nonzero()[0]

        # Map to (res_id, chain_id, aa_type) format (vectorized lookups)
        def indices_to_residue_list(idx_arr):
            if len(idx_arr) == 0:
                return []
            ids = p.res_ids[idx_arr]
            chains = p.chain_ids[idx_arr]
            names = p.res_names[idx_arr]
            return [
                {"res_id": int(ids[k]), "chain_id": str(chains[k]), "aa_type": str(names[k])}
                for k in range(len(idx_arr))
            ]

        # --- Slice interface arrays (no contact construction yet) ---
        lo = i * N
        hi = (i + 1) * N
        edge_mask = (contact_src_np >= lo) & (contact_src_np < hi)
        p_src = contact_src_np[edge_mask] - lo
        p_dst = contact_dst_np[edge_mask] - lo
        p_dist = contact_dist_np[edge_mask]

        dedup = p_src < p_dst
        p_src, p_dst, p_dist = p_src[dedup], p_dst[dedup], p_dist[dedup]
        valid_if = (p_src < n_res) & (p_dst < n_res)
        if not valid_if.all():
            p_src, p_dst, p_dist = p_src[valid_if], p_dst[valid_if], p_dist[valid_if]

        # --- Slice clash arrays per-protein ---
        bb_sel = bb_prot_np == i
        hv_sel = hv_prot_np == i

        result = ProteinAnalysisResult(
            path=p.path,
            n_residues=n_res,
            n_atoms=int(p.mask.sum()),
            n_backbone_clashes=int(bb_counts_np[i]),
            n_heavy_clashes=int(heavy_counts_np[i]),
            backbone_clashing_residues=indices_to_residue_list(bb_res_indices),
            heavy_clashing_residues=indices_to_residue_list(heavy_res_indices),
            interface_residues=indices_to_residue_list(interface_indices),
            # Contact lists are materialized later by _materialize_contacts
            # (runs in writer thread, concurrent with next GPU batch)
        )

        # Attach raw arrays for deferred contact materialization
        result._protein = p
        result._iface_src = p_src
        result._iface_dst = p_dst
        result._iface_dist = p_dist
        result._bb_clash = (
            bb_ri_np[bb_sel], bb_rj_np[bb_sel], bb_ai_np[bb_sel],
            bb_aj_np[bb_sel], bb_dists_np[bb_sel], bb_overlaps_np[bb_sel],
        )
        result._hv_clash = (
            hv_ri_np[hv_sel], hv_rj_np[hv_sel], hv_ai_np[hv_sel],
            hv_aj_np[hv_sel], hv_dists_np[hv_sel], hv_overlaps_np[hv_sel],
        )
        results.append(result)

    return results


def _extract_clash_contacts(
    p: Protein,
    protein_idx: int,
    prot_np, ri_np, rj_np, ai_np, aj_np, dists_np, overlaps_np,
    selection: str,
) -> List[ClashContact]:
    """Extract ClashContact list for a single protein from batch clash arrays."""
    mask = prot_np == protein_idx
    if not mask.any():
        return []

    indices = mask.nonzero()[0]

    # Batch-extract all arrays at once (one fancy-index op each)
    ri = ri_np[indices]
    rj = rj_np[indices]
    ai = ai_np[indices]
    aj = aj_np[indices]
    dists = dists_np[indices]
    overlaps = overlaps_np[indices]

    # Filter out-of-bounds residues
    valid = (ri < p.n_residues) & (rj < p.n_residues)
    if not valid.all():
        ri, rj, ai, aj = ri[valid], rj[valid], ai[valid], aj[valid]
        dists, overlaps = dists[valid], overlaps[valid]

    if len(ri) == 0:
        return []

    # Vectorized lookups into protein arrays
    chains_i = p.chain_ids[ri]
    res_ids_i = p.res_ids[ri]
    aa_i = p.res_names[ri]
    anames_i = p.atom_names[ri, ai]

    chains_j = p.chain_ids[rj]
    res_ids_j = p.res_ids[rj]
    aa_j = p.res_names[rj]
    anames_j = p.atom_names[rj, aj]

    contacts = []
    for k in range(len(ri)):
        a1 = str(anames_i[k]) if anames_i[k] else f"slot{ai[k]}"
        a2 = str(anames_j[k]) if anames_j[k] else f"slot{aj[k]}"
        contacts.append(ClashContact(
            chain_1=str(chains_i[k]),
            res_1_id=int(res_ids_i[k]),
            aa_1_type=str(aa_i[k]),
            atom_1_name=a1,
            chain_2=str(chains_j[k]),
            res_2_id=int(res_ids_j[k]),
            aa_2_type=str(aa_j[k]),
            atom_2_name=a2,
            distance=float(dists[k]),
            overlap=float(overlaps[k]),
        ))

    return contacts


def analyze_proteins(
    proteins: Sequence[Protein],
    batch_size: int = 32,
    device: str = "cuda",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    interface_cutoff: float = 8.0,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
    progress: bool = True,
    output_dir: Optional[Path] = None,
    n_write_workers: int = 4,
    analyses: Optional[Set[str]] = None,
) -> List[ProteinAnalysisResult]:
    """
    Analyze a list of proteins with optional parallel output writing.

    Args:
        proteins: List of Protein objects
        batch_size: Number of proteins per GPU batch
        device: Computation device ("cuda" or "cpu")
        clash_cutoff: VDW overlap fraction for clashes
        min_seq_sep: Minimum sequence separation for clashes
        interface_cutoff: CA-CA distance cutoff for interface
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, S-S pairs within 2.5 Å are not counted as clashes
        progress: Show progress bar
        output_dir: If provided, write individual JSON files here
        n_write_workers: Number of parallel writers (overlaps I/O with compute)

    Returns:
        List of ProteinAnalysisResult
    """
    results = []
    n_proteins = len(proteins)

    iterator = range(0, n_proteins, batch_size)
    if progress:
        iterator = tqdm(iterator, desc="Analyzing", unit="batch")

    # Set up parallel writer if output_dir provided
    writer_pool = None
    pending_writes: List[Future] = []

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        writer_pool = ThreadPoolExecutor(max_workers=n_write_workers)

    try:
        for batch_start in iterator:
            batch_end = min(batch_start + batch_size, n_proteins)
            batch_proteins = proteins[batch_start:batch_end]
            indices = list(range(batch_start, batch_end))

            batch = create_batch(batch_proteins, indices, device)
            batch_results = analyze_batch(
                batch,
                clash_cutoff=clash_cutoff,
                min_seq_sep=min_seq_sep,
                interface_cutoff=interface_cutoff,
                max_num_neighbors=max_num_neighbors,
                exclude_disulfides=exclude_disulfides,
                analyses=analyses,
            )
            results.extend(batch_results)

            # Submit batch results to writer pool (overlaps with next batch compute)
            if writer_pool is not None:
                future = writer_pool.submit(
                    _write_batch_results, batch_results, output_dir,
                    interface_cutoff,
                )
                pending_writes.append(future)

            # Clear GPU cache periodically
            if device == "cuda" and batch_start % (batch_size * 10) == 0:
                torch.cuda.empty_cache()

        # Wait for all writes to complete
        if writer_pool is not None:
            for future in pending_writes:
                future.result()  # Raises any exceptions

    finally:
        if writer_pool is not None:
            writer_pool.shutdown(wait=True)

    return results


def analyze_pdb_files(
    paths: Sequence[Union[str, Path]],
    output_dir: Optional[Union[str, Path]] = None,
    batch_size: int = 32,
    n_workers: int = 4,
    device: str = "cuda",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    interface_cutoff: float = 8.0,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
    progress: bool = True,
    analyses: Optional[Set[str]] = None,
) -> List[ProteinAnalysisResult]:
    """
    Analyze PDB files and write individual JSON results.

    This is the main entry point for batch processing.

    Args:
        paths: List of PDB file paths
        output_dir: Directory for output JSON files (one per PDB)
        batch_size: Number of proteins per GPU batch
        n_workers: Number of CPUs for parallel PDB parsing and JSON writing
        device: Computation device
        clash_cutoff: VDW overlap fraction for clashes
        min_seq_sep: Minimum sequence separation for clashes
        interface_cutoff: CA-CA distance cutoff for interface
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, S-S pairs within 2.5 Å are not counted as clashes
        progress: Show progress bars

    Returns:
        List of ProteinAnalysisResult

    Example:
        >>> from clashes import analyze_pdb_files
        >>> results = analyze_pdb_files(
        ...     ["complex1.pdb", "complex2.pdb"],
        ...     output_dir="results/",
        ...     batch_size=64,
        ...     device="cuda",
        ... )
        # Creates: results/complex1.json, results/complex2.json
    """
    # Parse PDB files
    if progress:
        print(f"Parsing {len(paths)} PDB files...")
    proteins = parse_proteins(paths, n_jobs=n_workers)

    if progress:
        print(f"Parsed {len(proteins)} proteins successfully")

    # Convert output_dir to Path if provided
    output_path = Path(output_dir) if output_dir is not None else None

    # Analyze with parallel writing
    results = analyze_proteins(
        proteins,
        batch_size=batch_size,
        device=device,
        clash_cutoff=clash_cutoff,
        min_seq_sep=min_seq_sep,
        interface_cutoff=interface_cutoff,
        max_num_neighbors=max_num_neighbors,
        exclude_disulfides=exclude_disulfides,
        progress=progress,
        output_dir=output_path,
        n_write_workers=n_workers,  # Use same number of workers
        analyses=analyses,
    )

    if progress and output_dir is not None:
        print(f"Results written to {output_dir}/ ({len(results)} files)")

    return results


def analyze_pdb_files_pipelined(
    paths: Sequence[Union[str, Path]],
    output_dir: Optional[Union[str, Path]] = None,
    batch_size: int = 32,
    n_workers: int = 4,
    device: str = "cuda",
    clash_cutoff: float = 0.12,
    min_seq_sep: int = 3,
    interface_cutoff: float = 8.0,
    max_num_neighbors: int = 128,
    exclude_disulfides: bool = True,
    progress: bool = True,
    prefetch_batches: int = 2,
    analyses: Optional[Set[str]] = None,
) -> List[ProteinAnalysisResult]:
    """
    Analyze PDB files with pipelined parsing, computation, and writing.

    Pipeline stages overlap:
      - Batch N-1: being written to disk (ThreadPoolExecutor)
      - Batch N: being computed on GPU
      - Batch N+1: being parsed from disk (multiprocessing.Pool)

    Performance notes:
        - For small datasets (<200 proteins), use analyze_pdb_files instead
          (multiprocessing overhead dominates)
        - For large datasets (>500 proteins), this is ~1.4x faster
        - At 1000 proteins: ~490 proteins/s vs ~350 proteins/s sequential
        - For 10M proteins: ~5.7 hours vs ~7.8 hours sequential

    Args:
        paths: List of PDB file paths
        output_dir: Directory for output JSON files (one per PDB)
        batch_size: Number of proteins per GPU batch
        n_workers: Number of CPUs for parsing and writing
        device: Computation device
        clash_cutoff: VDW overlap fraction for clashes
        min_seq_sep: Minimum sequence separation for clashes
        interface_cutoff: CA-CA distance cutoff for interface
        max_num_neighbors: Max neighbors for radius_graph
        exclude_disulfides: If True, S-S pairs within 2.5 Å are not counted as clashes
        progress: Show progress bars
        prefetch_batches: Number of batches to parse ahead (currently unused)

    Returns:
        List of ProteinAnalysisResult
    """
    import multiprocessing as mp

    n_paths = len(paths)
    if n_paths == 0:
        return []

    # Set up output directory
    output_path = None
    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

    all_results: List[ProteinAnalysisResult] = []
    pending_writes: List[Future] = []

    # Set up writer pool for async writing
    writer_pool = None
    if output_path is not None:
        writer_pool = ThreadPoolExecutor(max_workers=n_workers)

    try:
        # Use multiprocessing.Pool with imap_unordered for streaming parsing
        # chunksize controls how many files each worker processes at once
        chunksize = max(1, batch_size // n_workers)

        with mp.Pool(processes=n_workers) as parse_pool:
            # Start parsing all files asynchronously
            protein_iter = parse_pool.imap_unordered(
                parse_protein,
                paths,
                chunksize=chunksize,
            )

            pbar = None
            if progress:
                pbar = tqdm(total=n_paths, desc="Processing", unit="protein")

            current_batch: List[Protein] = []

            for protein in protein_iter:
                if protein is None:
                    if pbar:
                        pbar.update(1)
                    continue

                current_batch.append(protein)

                # Process batch when full
                if len(current_batch) >= batch_size:
                    batch_results = _process_batch(
                        current_batch, device,
                        clash_cutoff, min_seq_sep, interface_cutoff,
                        max_num_neighbors, exclude_disulfides,
                        analyses=analyses,
                    )
                    all_results.extend(batch_results)

                    # Async write (overlaps with next batch parse/compute)
                    if writer_pool is not None and output_path is not None:
                        future = writer_pool.submit(
                            _write_batch_results, batch_results, output_path,
                            interface_cutoff,
                        )
                        pending_writes.append(future)

                    if pbar:
                        pbar.update(len(current_batch))

                    current_batch = []

            # Process remaining proteins
            if current_batch:
                batch_results = _process_batch(
                    current_batch, device,
                    clash_cutoff, min_seq_sep, interface_cutoff,
                    max_num_neighbors, exclude_disulfides,
                    analyses=analyses,
                )
                all_results.extend(batch_results)

                if writer_pool is not None and output_path is not None:
                    future = writer_pool.submit(
                        _write_batch_results, batch_results, output_path,
                        interface_cutoff,
                    )
                    pending_writes.append(future)

                if pbar:
                    pbar.update(len(current_batch))

            if pbar:
                pbar.close()

        # Wait for all writes to complete
        for future in pending_writes:
            future.result()

        if progress and output_path is not None:
            print(f"Results written to {output_dir}/ ({len(all_results)} files)")

        return all_results

    finally:
        if writer_pool is not None:
            writer_pool.shutdown(wait=True)


def _process_batch(
    proteins: List[Protein],
    device: str,
    clash_cutoff: float,
    min_seq_sep: int,
    interface_cutoff: float,
    max_num_neighbors: int,
    exclude_disulfides: bool,
    analyses: Optional[Set[str]] = None,
) -> List[ProteinAnalysisResult]:
    """Create batch, move to GPU, compute metrics."""
    indices = list(range(len(proteins)))
    batch = create_batch(proteins, indices, device)
    return analyze_batch(
        batch,
        clash_cutoff=clash_cutoff,
        min_seq_sep=min_seq_sep,
        interface_cutoff=interface_cutoff,
        max_num_neighbors=max_num_neighbors,
        exclude_disulfides=exclude_disulfides,
        analyses=analyses,
    )


def main():
    """Command-line interface."""
    import argparse
    import glob

    parser = argparse.ArgumentParser(
        description="Analyze protein structures for clashes and interfaces"
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="PDB files or glob patterns",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for JSON files (one per PDB)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Proteins per GPU batch (default: 32)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="CPUs for PDB parsing and JSON writing (default: 4)",
    )
    parser.add_argument(
        "--device",
        default="cuda",
        choices=["cuda", "cpu"],
        help="Computation device (default: cuda)",
    )
    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=0.12,
        help="VDW overlap fraction for clashes (default: 0.12, ~0.4 A overlap)",
    )
    parser.add_argument(
        "--interface-cutoff",
        type=float,
        default=8.0,
        help="CA-CA distance for interface (default: 8.0 Å)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress output",
    )
    parser.add_argument(
        "--no-disulfide-exclusion",
        action="store_true",
        help="Count S-S bonds as clashes (default: exclude disulfide bonds)",
    )
    parser.add_argument(
        "--analyses",
        nargs="+",
        choices=["interface", "backbone_clashes", "heavy_atom_clashes"],
        default=["interface", "backbone_clashes", "heavy_atom_clashes"],
        help="Which analyses to run (default: all three)",
    )
    parser.add_argument(
        "--no-pipeline",
        action="store_true",
        help="Disable pipelined I/O (parse all first, then compute)",
    )
    parser.add_argument(
        "--prefetch",
        type=int,
        default=2,
        help="Number of batches to prefetch in pipeline mode (default: 2)",
    )

    args = parser.parse_args()

    # Expand glob patterns
    paths = []
    for pattern in args.inputs:
        expanded = glob.glob(pattern)
        if expanded:
            paths.extend(expanded)
        else:
            paths.append(pattern)  # Might be a literal path

    if not paths:
        parser.error("No input files found")

    # Run analysis (pipelined by default)
    analyses_set = set(args.analyses)
    if args.no_pipeline:
        analyze_pdb_files(
            paths,
            output_dir=args.output_dir,
            batch_size=args.batch_size,
            n_workers=args.workers,
            device=args.device,
            clash_cutoff=args.clash_cutoff,
            interface_cutoff=args.interface_cutoff,
            exclude_disulfides=not args.no_disulfide_exclusion,
            progress=not args.quiet,
            analyses=analyses_set,
        )
    else:
        analyze_pdb_files_pipelined(
            paths,
            output_dir=args.output_dir,
            batch_size=args.batch_size,
            n_workers=args.workers,
            device=args.device,
            clash_cutoff=args.clash_cutoff,
            interface_cutoff=args.interface_cutoff,
            exclude_disulfides=not args.no_disulfide_exclusion,
            progress=not args.quiet,
            prefetch_batches=args.prefetch,
            analyses=analyses_set,
        )


if __name__ == "__main__":
    main()
