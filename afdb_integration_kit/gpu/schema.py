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
Convert ProteinAnalysisResult to FunPDBe-compatible schema dicts.

Produces two separate outputs:
  - Interface schema (one per chain pair)
  - Clash schema (backbone + side-chain sites)

Both conform to interface_annotations/interface_schema.json.
"""
from __future__ import annotations

import re
from collections import defaultdict
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from .analyze import ClashContact, ProteinAnalysisResult

# Regex to extract AF-XXXXXXXXXXXXXXXX from filenames
_AF_ID_RE = re.compile(r"(AF-\d{16})")


def _extract_af_id(path: str) -> Optional[str]:
    """Extract AlphaFold DB identifier from a file path."""
    m = _AF_ID_RE.search(path)
    return m.group(1) if m else None


def _get_af_id_or_stem(path: str) -> str:
    """Return the AF ID if present, otherwise the filename stem."""
    af_id = _extract_af_id(path)
    if af_id is not None:
        return af_id
    from pathlib import Path
    return Path(path).stem


def _canon_chain_pair(c1: str, c2: str) -> Tuple[str, str]:
    """Return chain pair in canonical (alphabetical) order."""
    return (c1, c2) if c1 <= c2 else (c2, c1)


_ECO = [{"eco_term": "computational evidence", "eco_code": "ECO_0000246"}]


# ---------------------------------------------------------------------------
# Interface schema
# ---------------------------------------------------------------------------

def _precompute_closest_partners(contacts) -> Dict[Tuple[str, int], dict]:
    """Build (chain, res_id) -> closest-partner dict in one pass over contacts."""
    best: Dict[Tuple[str, int], Tuple[float, dict]] = {}
    for c in contacts:
        # Side 1 -> partner is side 2
        key1 = (c.chain_1, c.res_1_id)
        if key1 not in best or c.distance < best[key1][0]:
            best[key1] = (c.distance, {
                "partner_chain": c.chain_2,
                "partner_res_label": str(c.res_2_id),
                "partner_aa_type": c.aa_2_type,
                "distance_angstrom": round(c.distance, 2),
                "method": "ca_ca_distance",
            })
        # Side 2 -> partner is side 1
        key2 = (c.chain_2, c.res_2_id)
        if key2 not in best or c.distance < best[key2][0]:
            best[key2] = (c.distance, {
                "partner_chain": c.chain_1,
                "partner_res_label": str(c.res_1_id),
                "partner_aa_type": c.aa_1_type,
                "distance_angstrom": round(c.distance, 2),
                "method": "ca_ca_distance",
            })
    return {k: v[1] for k, v in best.items()}


def result_to_interface_schema(
    result: "ProteinAnalysisResult",
    interface_cutoff: float = 8.0,
    data_resource: str = "AFDB-Integration-Kit",
) -> Optional[dict]:
    """
    Build the interface-only FunPDBe schema dict.

    Returns None if there are no interface contacts.
    """
    # Group contacts by canonical chain pair
    pair_contacts: Dict[Tuple[str, str], list] = defaultdict(list)
    for c in result.interface_contacts:
        pair = _canon_chain_pair(c.chain_1, c.chain_2)
        pair_contacts[pair].append(c)

    if not pair_contacts:
        return None

    sites: List[dict] = []
    site_id_counter = 1

    residue_site_data: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    residue_aa: Dict[Tuple[str, str], str] = {}
    residue_annotations: Dict[Tuple[str, str], dict] = {}

    for pair in sorted(pair_contacts.keys()):
        contacts = pair_contacts[pair]
        site_id = site_id_counter
        site_id_counter += 1

        # Build interactions list
        interactions = []
        for c in contacts:
            if _canon_chain_pair(c.chain_1, c.chain_2) == (c.chain_1, c.chain_2):
                interaction = {
                    "chain_1": c.chain_1,
                    "res_1_label": str(c.res_1_id),
                    "aa_1_type": c.aa_1_type,
                    "chain_2": c.chain_2,
                    "res_2_label": str(c.res_2_id),
                    "aa_2_type": c.aa_2_type,
                    "distance_angstrom": round(c.distance, 2),
                }
            else:
                interaction = {
                    "chain_1": c.chain_2,
                    "res_1_label": str(c.res_2_id),
                    "aa_1_type": c.aa_2_type,
                    "chain_2": c.chain_1,
                    "res_2_label": str(c.res_1_id),
                    "aa_2_type": c.aa_1_type,
                    "distance_angstrom": round(c.distance, 2),
                }
            interactions.append(interaction)

        label = f"{pair[0]}{pair[1]}_interface"
        sites.append({
            "site_id": site_id,
            "label": label,
            "additional_site_annotations": {
                "contact_cutoff_angstrom": interface_cutoff,
                "description": f"{pair[0]}{pair[1]} interface",
                "interactions": interactions,
            },
        })

        # Pre-compute closest partner per residue in one pass (O(C))
        closest = _precompute_closest_partners(contacts)

        # Collect interface residues with partner annotations
        seen: set = set()
        for c in contacts:
            for chain, res_id, aa in [
                (c.chain_1, c.res_1_id, c.aa_1_type),
                (c.chain_2, c.res_2_id, c.aa_2_type),
            ]:
                key = (chain, str(res_id))
                residue_aa[key] = aa
                if key not in seen:
                    seen.add(key)
                    residue_site_data[key].append({
                        "site_id_ref": site_id,
                        "confidence_classification": "null",
                    })
                    partner = closest.get((chain, res_id))
                    if partner and key not in residue_annotations:
                        residue_annotations[key] = partner

    # Build chains array
    chain_residues: Dict[str, List[dict]] = defaultdict(list)
    for (chain, res_id_str), sd_list in sorted(residue_site_data.items()):
        rd: dict = {
            "residue_number": res_id_str,
            "aa_type": residue_aa.get((chain, res_id_str), "UNK"),
            "site_data": sd_list,
        }
        ann = residue_annotations.get((chain, res_id_str))
        if ann:
            rd["additional_residue_annotations"] = ann
        chain_residues[chain].append(rd)

    chains = []
    for cl in sorted(chain_residues.keys()):
        residues = chain_residues[cl]
        residues.sort(key=lambda r: int(r["residue_number"]))
        chains.append({"chain_label": cl, "residues": residues})

    return {
        "data_resource": data_resource,
        "af_id": _get_af_id_or_stem(result.path),
        "chains": chains,
        "sites": sites,
        "evidence_code_ontology": _ECO,
    }


# ---------------------------------------------------------------------------
# Clash schema
# ---------------------------------------------------------------------------

def _clash_severity(overlap: float) -> str:
    """Classify clash severity based on VDW overlap (Angstroms).

    Overlap = vdw_sum - distance.  With clash_cutoff=0.12 the minimum
    detected overlap is ~0.4 A, giving two tiers:

        >= 0.8 A  -> "high"   (severe steric clash)
        <  0.8 A  -> "medium" (moderate clash)
    """
    if overlap >= 0.8:
        return "high"
    return "medium"


def _precompute_worst_overlaps(
    clash_contacts: List["ClashContact"],
) -> Dict[Tuple[str, int], float]:
    """Build (chain, res_id) -> max overlap dict in one pass over contacts."""
    worst: Dict[Tuple[str, int], float] = {}
    for c in clash_contacts:
        key1 = (c.chain_1, c.res_1_id)
        if c.overlap > worst.get(key1, 0.0):
            worst[key1] = c.overlap
        key2 = (c.chain_2, c.res_2_id)
        if c.overlap > worst.get(key2, 0.0):
            worst[key2] = c.overlap
    return worst


def _build_clash_site(
    site_id: int,
    label: str,
    description: str,
    n_clashes: int,
    clashing_residues: List[dict],
    clash_contacts: List["ClashContact"],
    residue_site_data: Dict[Tuple[str, str], List[dict]],
    residue_aa: Dict[Tuple[str, str], str],
) -> dict:
    """Build a single clash site and populate per-residue site_data."""
    # Pre-compute max overlap per residue in one pass (O(C))
    worst_overlaps = _precompute_worst_overlaps(clash_contacts)

    for res in clashing_residues:
        key = (res["chain_id"], str(res["res_id"]))
        if key not in residue_aa:
            residue_aa[key] = res.get("aa_type", "UNK")

        worst = worst_overlaps.get((res["chain_id"], res["res_id"]), 0.0)
        residue_site_data[key].append({
            "site_id_ref": site_id,
            "raw_score": round(worst, 2),
            "raw_score_unit": "vdw_overlap_angstrom",
            "confidence_classification": _clash_severity(worst),
        })

    return {
        "site_id": site_id,
        "label": label,
        "additional_site_annotations": {
            "description": description,
            "n_clashes": n_clashes,
        },
    }


def result_to_clash_schema(
    result: "ProteinAnalysisResult",
    data_resource: str = "AFDB-Integration-Kit",
) -> dict:
    """
    Build the clash-only FunPDBe schema dict.

    Always includes both backbone_clashes and heavy_atom_clashes sites (with
    n_clashes=0 when none found) so that counts are always recorded in the JSON
    and downstream enrichment can distinguish "analyzed, 0 clashes" from "not analyzed".
    """
    sites: List[dict] = []
    site_id_counter = 1

    residue_site_data: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    residue_aa: Dict[Tuple[str, str], str] = {}

    site = _build_clash_site(
        site_id=site_id_counter,
        label="backbone_clashes",
        description="Residues with backbone atom clashes (N, CA, C, O)",
        n_clashes=result.n_backbone_clashes,
        clashing_residues=result.backbone_clashing_residues,
        clash_contacts=result.backbone_clash_contacts,
        residue_site_data=residue_site_data,
        residue_aa=residue_aa,
    )
    sites.append(site)
    site_id_counter += 1

    site = _build_clash_site(
        site_id=site_id_counter,
        label="heavy_atom_clashes",
        description="Residues with heavy atom clashes (all non-hydrogen atoms)",
        n_clashes=result.n_heavy_clashes,
        clashing_residues=result.heavy_clashing_residues,
        clash_contacts=result.heavy_clash_contacts,
        residue_site_data=residue_site_data,
        residue_aa=residue_aa,
    )
    sites.append(site)
    site_id_counter += 1

    # Build chains array
    chain_residues: Dict[str, List[dict]] = defaultdict(list)
    for (chain, res_id_str), sd_list in sorted(residue_site_data.items()):
        chain_residues[chain].append({
            "residue_number": res_id_str,
            "aa_type": residue_aa.get((chain, res_id_str), "UNK"),
            "site_data": sd_list,
        })

    chains = []
    for cl in sorted(chain_residues.keys()):
        residues = chain_residues[cl]
        residues.sort(key=lambda r: int(r["residue_number"]))
        chains.append({"chain_label": cl, "residues": residues})

    return {
        "data_resource": data_resource,
        "af_id": _get_af_id_or_stem(result.path),
        "chains": chains,
        "sites": sites,
        "evidence_code_ontology": _ECO,
    }
