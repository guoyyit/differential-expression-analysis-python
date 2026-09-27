"""Multiple-testing correction over a vector of p-values."""

from __future__ import annotations

import numpy as np

from app.models import CorrectionMethod


def adjust_pvalues(p_values: np.ndarray, method: CorrectionMethod) -> np.ndarray:
    """Return adjusted p-values aligned with the input order.

    BH: q_i = min over ranks j>=i of (n/j * p_(j)), clipped to [0, 1].
    Bonferroni: min(1, n * p_i). NONE: identity. All deterministic.
    """
    p = np.asarray(p_values, dtype=np.float64)
    n = p.size
    if n == 0:
        return p.copy()
    if method is CorrectionMethod.NONE:
        return p.copy()
    if method is CorrectionMethod.BONFERRONI:
        return np.minimum(1.0, p * n)
    if method is CorrectionMethod.BH:
        order = np.argsort(p, kind="stable")
        ranked = p[order]
        ranks = np.arange(1, n + 1, dtype=np.float64)
        adjusted_sorted = ranked * n / ranks
        # enforce monotonicity from the largest p-value downwards
        adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
        adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
        out = np.empty(n, dtype=np.float64)
        out[order] = adjusted_sorted
        return out
    raise ValueError(f"unknown correction method: {method!r}")
