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
Clashes - GPU-accelerated protein structure analysis.

Computes atomic clashes and interface residues using efficient
batched processing on GPU with torch_cluster.

Example (small datasets):
    >>> from clashes import analyze_pdb_files
    >>> results = analyze_pdb_files(
    ...     ["complex1.pdb", "complex2.pdb"],
    ...     output_dir="results/",
    ...     device="cuda",
    ... )

Example (large datasets, >500 proteins):
    >>> from clashes import analyze_pdb_files_pipelined
    >>> results = analyze_pdb_files_pipelined(
    ...     pdb_paths,  # List of 10M paths
    ...     output_dir="results/",
    ...     device="cuda",
    ... )
    # ~1.4x faster due to pipelined parsing/compute/writing
"""

# Core data structures
from .protein import Protein, empty_protein, MAX_HEAVY_ATOMS, AA_ORDER_HEAVY

# Parsing
from .parse import parse_protein, parse_proteins

# Batching
from .batch import ProteinBatch, create_batch, iter_batches, VDW_RADII

# Clash detection
from .clashes import (
    count_clashes,
    get_clash_pairs,
    compute_clashes_from_batch,
    compute_clashing_residues_from_batch,
    SULFUR_VDW,
    DISULFIDE_MAX,
)

# Interface detection
from .interface import (
    compute_interface_residues,
    compute_interface_residues_flat,
)

# High-level analysis
from .analyze import (
    ALL_ANALYSES,
    ClashContact,
    InterfaceContact,
    ProteinAnalysisResult,
    analyze_batch,
    analyze_proteins,
    analyze_pdb_files,
    analyze_pdb_files_pipelined,
)

# Schema conversion
from .schema import result_to_interface_schema, result_to_clash_schema

__all__ = [
    # Protein
    "Protein",
    "empty_protein",
    "MAX_HEAVY_ATOMS",
    "AA_ORDER_HEAVY",
    # Parsing
    "parse_protein",
    "parse_proteins",
    # Batching
    "ProteinBatch",
    "create_batch",
    "iter_batches",
    "VDW_RADII",
    # Clashes
    "count_clashes",
    "get_clash_pairs",
    "compute_clashes_from_batch",
    "compute_clashing_residues_from_batch",
    "SULFUR_VDW",
    "DISULFIDE_MAX",
    # Interface
    "compute_interface_residues",
    "compute_interface_residues_flat",
    # Analysis
    "ALL_ANALYSES",
    "ClashContact",
    "InterfaceContact",
    "ProteinAnalysisResult",
    "analyze_batch",
    "analyze_proteins",
    "analyze_pdb_files",
    "analyze_pdb_files_pipelined",
    # Schema
    "result_to_interface_schema",
    "result_to_clash_schema",
]
