"""Per-gene two-group differential tests. Deterministic, no randomness."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm, t as t_dist

from app.models import TestMethod


def run_test(a: np.ndarray, b: np.ndarray, method: TestMethod) -> tuple[float, float]:
    """Return (statistic, p_value) for group vectors ``a`` and ``b``."""
    if method is TestMethod.WELCH:
        return welch_t(a, b)
    if method is TestMethod.MANN_WHITNEY:
        return mann_whitney(a, b)
    raise ValueError(f"unknown test method: {method!r}")


def welch_t(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Welch's t-test: returns (t statistic, two-sided p value).

    Sign convention: positive t means group B > group A. Degenerate
    zero-variance cases return defined values instead of NaN so the
    pipeline output stays total and reproducible.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    na, nb = a.size, b.size
    mean_a, mean_b = a.mean(), b.mean()
    var_a = a.var(ddof=1) if na > 1 else 0.0
    var_b = b.var(ddof=1) if nb > 1 else 0.0
    se2 = var_a / na + var_b / nb
    diff = mean_b - mean_a
    if se2 <= 0.0:
        if diff == 0.0:
            return 0.0, 1.0
        return math.copysign(math.inf, diff), 0.0
    t_stat = diff / math.sqrt(se2)
    # Welch-Satterthwaite degrees of freedom
    num = se2 * se2
    den = 0.0
    if na > 1:
        den += (var_a / na) ** 2 / (na - 1)
    if nb > 1:
        den += (var_b / nb) ** 2 / (nb - 1)
    df = num / den if den > 0.0 else float(na + nb - 2)
    p = 2.0 * float(t_dist.sf(abs(t_stat), df))
    return float(t_stat), min(1.0, max(0.0, p))


def mann_whitney(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Rank-sum test: (z statistic, two-sided p), normal approximation.

    Uses average ranks for ties, a tie-corrected variance and a continuity
    correction. Positive z means group B tends to exceed group A.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    na, nb = a.size, b.size
    n = na + nb
    values = np.concatenate([a, b])
    ranks = _average_ranks(values)
    rank_b = ranks[na:].sum()
    u_b = rank_b - nb * (nb + 1) / 2.0
    mean_u = na * nb / 2.0
    # tie correction: sum(t^3 - t) over tied groups
    _, counts = np.unique(values, return_counts=True)
    tie_term = float((counts**3 - counts).sum())
    var_u = na * nb / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    if var_u <= 0.0:
        return 0.0, 1.0
    diff = u_b - mean_u
    # continuity correction towards the mean
    z = (diff - math.copysign(0.5, diff)) / math.sqrt(var_u) if diff != 0.0 else 0.0
    p = 2.0 * float(norm.sf(abs(z)))
    return float(z), min(1.0, max(0.0, p))


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Ranks from 1..n with ties assigned their average rank. Stable."""
    order = np.argsort(values, kind="stable")
    sorted_vals = values[order]
    ranks = np.empty(values.size, dtype=np.float64)
    i = 0
    n = values.size
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        ranks[order[i : j + 1]] = avg
        i = j + 1
    return ranks
