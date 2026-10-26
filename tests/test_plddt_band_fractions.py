"""Regression tests for per-chain pLDDT band fractions (AFPDBE-1511).

AFDB bands: very low < 50; low 50 to < 70; confident 70 to <= 90; very high > 90.
All scores below are SYNTHETIC; expected values are derived by hand.
"""

import pytest

from afdb_integration_kit.colabfold.converter import (
    _categorise_confidence,
    _compute_plddt_metrics,
)

FRACTION_KEYS = (
    "fraction_plddt_very_low",
    "fraction_plddt_low",
    "fraction_plddt_confident",
    "fraction_plddt_very_high",
)


def _band_fractions(scores_by_chain: dict[str, list[float]]) -> dict[str, tuple[float, ...]]:
    flat: list[float] = []
    chains = []
    manifest = []
    for chain_id, scores in scores_by_chain.items():
        start = len(flat) + 1
        flat.extend(scores)
        chains.append(
            {"label_asym_id": chain_id, "sequenceStart": start, "sequenceEnd": len(flat)}
        )
        manifest.append(
            {
                "chain_id": chain_id,
                "uniprot_ac": "",
                "entity_id": "1",
                "sequence_start": "1",
                "sequence_end": str(len(scores)),
            }
        )
    rows, _ = _compute_plddt_metrics(flat, chains, manifest)
    return {row["chain_id"]: tuple(row[key] for key in FRACTION_KEYS) for row in rows}


def test_boundary_scores_fall_in_the_documented_band() -> None:
    scores = [49.99, 50.0, 69.99, 70.0, 89.99, 90.0, 90.01, 30.0, 95.0, 60.0]
    # very low {49.99, 30.0}; low {50.0, 69.99, 60.0};
    # confident {70.0, 89.99, 90.0}; very high {90.01, 95.0}
    assert _band_fractions({"A": scores})["A"] == (0.2, 0.3, 0.3, 0.2)


def test_fractions_are_not_derived_from_five_letter_categories() -> None:
    # Same band composition as the VADB chain investigated in AFPDBE-1511:
    # D=41 (<30), L=141 (30-<50), M=72 (50-<70), H=30 (70-90), V=0.
    scores = [25.0] * 41 + [40.0] * 141 + [60.0] * 72 + [80.0] * 30
    assert [_categorise_confidence(s) for s in (25.0, 40.0, 60.0, 80.0)] == ["D", "L", "M", "H"]

    # 182/284, 72/284, 30/284, 0/284 rounded to 3 d.p.
    assert _band_fractions({"A": scores})["A"] == (0.641, 0.254, 0.106, 0.0)
    # Regrouping letters as D / L / M+H / V gives 0.144 / 0.496 / 0.359 / 0.0.
    assert _band_fractions({"A": scores})["A"] != (0.144, 0.496, 0.359, 0.0)


def test_dimer_chains_get_their_own_fractions_not_pooled_model_values() -> None:
    fractions = _band_fractions({"A": [40.0, 60.0, 80.0, 95.0], "B": [95.0, 95.0, 95.0, 95.0]})

    assert fractions["A"] == (0.25, 0.25, 0.25, 0.25)
    assert fractions["B"] == (0.0, 0.0, 0.0, 1.0)
    # Pooling both chains would give 0.125 / 0.125 / 0.125 / 0.625 for each chain.


@pytest.mark.parametrize(
    ("score", "category"),
    [(29.99, "D"), (30.0, "L"), (49.99, "L"), (50.0, "M"), (69.99, "M"),
     (70.0, "H"), (90.0, "H"), (90.01, "V")],
)
def test_category_letter_boundaries(score: float, category: str) -> None:
    assert _categorise_confidence(score) == category
