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
Protein data structure for molecular analysis.

Simple dataclass to hold protein heavy atom coordinates and metadata.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np


# Maximum heavy atoms per residue
MAX_HEAVY_ATOMS = 14

# Amino acid heavy atom ordering (canonical)
AA_ORDER_HEAVY: Dict[str, Tuple[str, ...]] = {
    "ALA": ("N", "CA", "C", "O", "CB"),
    "ARG": ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "ASN": ("N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"),
    "ASP": ("N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"),
    "CYS": ("N", "CA", "C", "O", "CB", "SG"),
    "GLN": ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"),
    "GLU": ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"),
    "GLY": ("N", "CA", "C", "O"),
    "HIS": ("N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"),
    "ILE": ("N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"),
    "LEU": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"),
    "LYS": ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"),
    "MET": ("N", "CA", "C", "O", "CB", "CG", "SD", "CE"),
    "PHE": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "PRO": ("N", "CA", "C", "O", "CB", "CG", "CD"),
    "SER": ("N", "CA", "C", "O", "CB", "OG"),
    "THR": ("N", "CA", "C", "O", "CB", "OG1", "CG2"),
    "TRP": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "TYR": ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "VAL": ("N", "CA", "C", "O", "CB", "CG1", "CG2"),
}

# Amino acid name aliases (non-standard to standard)
AA_ALIASES: Dict[str, str] = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "CYX": "CYS", "CYM": "CYS",
    "ASH": "ASP", "GLH": "GLU",
    "LYN": "LYS", "ARGN": "ARG",
    "MSE": "MET", "SEC": "CYS",
}

_CANON_AA = frozenset(AA_ORDER_HEAVY.keys())
_CANON_OR_ALIAS = _CANON_AA | frozenset(AA_ALIASES.keys())

# Pre-build slot lookup: res_type -> {atom_name: slot_index}
_ATOM_SLOT_LOOKUP: Dict[str, Dict[str, int]] = {
    res: {atom: i for i, atom in enumerate(atoms)}
    for res, atoms in AA_ORDER_HEAVY.items()
}


@dataclass(slots=True)
class Protein:
    """
    Protein structure with heavy atom coordinates.

    Attributes:
        coords: (N_res, MAX_HEAVY_ATOMS, 3) float32 coordinates in Angstroms
        mask: (N_res, MAX_HEAVY_ATOMS) bool mask for valid atoms
        atom_names: (N_res, MAX_HEAVY_ATOMS) str atom names
        res_names: (N_res,) str residue names (3-letter code)
        chain_ids: (N_res,) str chain identifiers
        res_ids: (N_res,) int residue numbers
        path: str source file path
    """
    coords: np.ndarray          # (N_res, MAX_HEAVY_ATOMS, 3) float32
    mask: np.ndarray            # (N_res, MAX_HEAVY_ATOMS) bool
    atom_names: np.ndarray      # (N_res, MAX_HEAVY_ATOMS) <U4
    res_names: np.ndarray       # (N_res,) <U3
    chain_ids: np.ndarray       # (N_res,) <U1
    res_ids: np.ndarray         # (N_res,) int32
    path: str = ""

    @property
    def n_residues(self) -> int:
        """Number of residues."""
        return self.coords.shape[0]

    @property
    def n_atoms(self) -> int:
        """Number of valid atoms."""
        return int(self.mask.sum())

    def get_backbone_coords(self) -> np.ndarray:
        """Get backbone coordinates (N, CA, C, O)."""
        return self.coords[:, :4, :]

    def get_backbone_mask(self) -> np.ndarray:
        """Get backbone atom mask."""
        return self.mask[:, :4]

    def get_ca_coords(self) -> np.ndarray:
        """Get CA coordinates."""
        return self.coords[:, 1, :]

    def __repr__(self) -> str:
        return f"Protein(n_residues={self.n_residues}, n_atoms={self.n_atoms}, path='{self.path}')"


def empty_protein(path: str = "") -> Protein:
    """Create an empty protein."""
    return Protein(
        coords=np.zeros((0, MAX_HEAVY_ATOMS, 3), dtype=np.float32),
        mask=np.zeros((0, MAX_HEAVY_ATOMS), dtype=bool),
        atom_names=np.zeros((0, MAX_HEAVY_ATOMS), dtype="<U4"),
        res_names=np.zeros((0,), dtype="<U3"),
        chain_ids=np.zeros((0,), dtype="<U1"),
        res_ids=np.zeros((0,), dtype=np.int32),
        path=path,
    )
