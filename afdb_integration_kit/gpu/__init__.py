# Copyright 2026 Maciej Majewski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# SPDX-License-Identifier: Apache-2.0

"""GPU clash/interface analysis with lazy optional imports."""
from __future__ import annotations

from importlib import import_module

from .protein import AA_ORDER_HEAVY, MAX_HEAVY_ATOMS, Protein, empty_protein
from .schema import result_to_clash_schema, result_to_interface_schema

_LAZY_EXPORTS = {
    "parse_protein": ".parse",
    "parse_proteins": ".parse",
    "ProteinBatch": ".batch",
    "create_batch": ".batch",
    "iter_batches": ".batch",
    "VDW_RADII": ".batch",
    "count_clashes": ".clashes",
    "get_clash_pairs": ".clashes",
    "compute_clashes_from_batch": ".clashes",
    "compute_clashing_residues_from_batch": ".clashes",
    "SULFUR_VDW": ".clashes",
    "DISULFIDE_MAX": ".clashes",
    "compute_interface_residues": ".interface",
    "compute_interface_residues_flat": ".interface",
    "ALL_ANALYSES": ".analyze",
    "ClashContact": ".analyze",
    "InterfaceContact": ".analyze",
    "ProteinAnalysisResult": ".analyze",
    "analyze_batch": ".analyze",
    "analyze_proteins": ".analyze",
    "analyze_pdb_files": ".analyze",
    "analyze_pdb_files_pipelined": ".analyze",
}


def __getattr__(name: str):
    """Load production-analysis modules only when requested."""
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))


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
