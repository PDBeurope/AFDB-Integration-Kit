"""Deterministic fixed-decimal rounding helpers.

Every numeric value written into a deposited payload (mmCIF/BCIF QA metrics, the
metadata JSONs, the chain/model metric manifests, the clash/interface JSONs) is a
float that gets recomputed (a mean, a fraction, a distance) and then rounded to a
fixed number of decimals. Plain ``round(x, n)`` / ``f"{x:.nf}"`` / ``np.round`` make
that decision on the *binary* float, so a value sitting exactly on a half-cent
boundary (e.g. a mean of ``81.585``) can round to ``81.58`` on one run and ``81.59``
on another after a tiny change in summation order or code path. That produces
non-deterministic bytes for the same model and breaks payload-parity checks.

These helpers make the rounding a function of the exact decimal value instead of the
fragile float, using ``Decimal`` with ``ROUND_HALF_EVEN`` (the same half-to-even rule
Python's ``round``/``:.nf``/``np.round`` already use), so:

* every non-boundary value is byte-identical to the previous behaviour, and
* exact half-cent boundaries become deterministic and order-independent.
"""

from __future__ import annotations

from decimal import Decimal, ROUND_HALF_EVEN
from typing import Iterable, Optional, Union

Number = Union[int, float, str, Decimal]


def quantize_half_even(value: Number, decimals: int) -> Decimal:
    """Round ``value`` to ``decimals`` places via exact Decimal half-to-even.

    ``value`` is parsed via ``str()`` first so floats are taken at their shortest
    decimal repr (matching what ``round``/``:.nf`` display), not their full binary
    expansion.
    """
    quantum = Decimal(1).scaleb(-decimals)  # 0.01 for decimals=2, 0.001 for 3, ...
    return Decimal(str(value)).quantize(quantum, rounding=ROUND_HALF_EVEN)


def round_float(value: Number, decimals: int) -> float:
    """Deterministic drop-in for ``round(float, n)`` returning a float.

    Use for JSON/CSV numeric fields so their type is unchanged but the rounding
    decision is made on the exact decimal, not the binary float.
    """
    return float(quantize_half_even(value, decimals))


def mean_2dp_str(string_values: Iterable[str]) -> Optional[str]:
    """Exact 2-decimal mean of textual numeric values, as a formatted string.

    Built for averaging PDB B-factors (textual, already 2-dp) into a ModelCIF
    ``metric_value``. Empty/sentinel-only inputs (``"?"``, ``"."``, ``""``) yield
    ``None`` so callers can skip writing the metric. Because the sum is taken over
    exact ``Decimal`` values, the result is independent of summation order.
    """
    vals = [Decimal(v) for v in string_values if v not in ("?", ".", "")]
    if not vals:
        return None
    mean = sum(vals) / Decimal(len(vals))
    return str(mean.quantize(Decimal("0.01"), rounding=ROUND_HALF_EVEN))
