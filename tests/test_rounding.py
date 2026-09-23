"""Tests for the deterministic fixed-decimal rounding helpers.

These guard against the ModelCIF/BCIF payload-parity drift where a recomputed float
mean sitting on a half-cent boundary (e.g. 81.585) rounded to 81.58 on one run and
81.59 on another. The helpers must be deterministic and independent of summation
order, while remaining byte-identical to the previous round-half-even behaviour for
every non-boundary value.
"""
import random

import pytest

from afdb_integration_kit.utils.rounding import (
    mean_2dp_str,
    quantize_half_even,
    round_float,
)


# --- mean_2dp_str: the site that actually drifted (CIF global/local pLDDT) ---------

def test_mean_half_cent_boundary_is_half_even():
    # 163.17 / 2 = 81.585 exactly -> half-to-even keeps the even hundredth (81.58).
    assert mean_2dp_str(["81.58", "81.59"]) == "81.58"
    # 163.19 / 2 = 81.595 exactly -> half-to-even rounds up to the even hundredth (81.60).
    assert mean_2dp_str(["81.59", "81.60"]) == "81.60"


def test_mean_is_order_independent_on_boundary():
    # The core regression: the boundary result must not depend on summation order.
    assert mean_2dp_str(["81.58", "81.59"]) == mean_2dp_str(["81.59", "81.58"])

    base = ["80.12", "81.58", "81.59", "79.44", "85.00", "62.33", "90.01"]
    expected = mean_2dp_str(base)
    rng = random.Random(0)
    for _ in range(50):
        shuffled = base[:]
        rng.shuffle(shuffled)
        assert mean_2dp_str(shuffled) == expected


def test_mean_matches_legacy_for_non_boundary_values():
    # Away from a half-cent boundary, the helper is byte-identical to the old f"{:.2f}".
    for values in (["70.00", "80.00", "90.00"], ["12.34", "56.78"], ["88.88"]):
        legacy = f"{sum(map(float, values)) / len(values):.2f}"
        assert mean_2dp_str(values) == legacy


def test_mean_filters_sentinels_and_handles_empty():
    assert mean_2dp_str(["80.00", "?", "82.00", ".", ""]) == "81.00"
    assert mean_2dp_str([]) is None
    assert mean_2dp_str(["?", ".", ""]) is None
    assert mean_2dp_str(["73.21"]) == "73.21"


def test_mean_reproduces_854_855_case_deterministically():
    # An AF-0000000210947144-style population whose exact mean lands on a half-cent
    # boundary: many atoms at 81.58 and 81.59 averaging to 81.585. Must be stable and
    # order-independent regardless of how many atoms or in what order.
    values = ["81.58"] * 250 + ["81.59"] * 250
    result = mean_2dp_str(values)
    assert result == "81.58"  # 81.585 -> half-even
    rng = random.Random(42)
    shuffled = values[:]
    rng.shuffle(shuffled)
    assert mean_2dp_str(shuffled) == result


# --- round_float / quantize_half_even: the recommended JSON/CSV sites ---------------

def test_round_float_half_even_on_boundary():
    assert round_float(0.125, 2) == 0.12   # exact 0.125 -> even hundredth
    assert round_float(0.135, 2) == 0.14   # -> even hundredth
    assert round_float(2.5, 0) == 2.0
    assert round_float(3.5, 0) == 4.0


def test_round_float_matches_builtin_for_non_boundary():
    for value in (81.58, 0.0107492, 7.43432, 0.6, 0.23, 99.994, 12.341):
        assert round_float(value, 2) == round(value, 2)
    for value in (0.0107492, 0.106888, 0.769506):
        assert round_float(value, 3) == round(value, 3)


def test_round_float_accepts_int_and_returns_float():
    assert round_float(5, 2) == 5.0
    assert isinstance(round_float(5, 2), float)


def test_quantize_supports_three_decimals():
    assert str(quantize_half_even(0.7695, 3)) == "0.770"   # half-even -> even last digit (0)
    assert str(quantize_half_even("0.7685", 3)) == "0.768"  # half-even -> even last digit (8)
