"""Effect sizes on normalized expression values."""

from __future__ import annotations

import math

import numpy as np

_PSEUDOCOUNT = 1.0


def log2_fold_change(mean_a: float, mean_b: float) -> float:
    """log2(B/A) with a pseudocount so zero means stay finite.

    Computed as log2(mean_b + 1) - log2(mean_a + 1): symmetric, defined at
    zero, and monotone in the mean difference.
    """
    return math.log2(mean_b + _PSEUDOCOUNT) - math.log2(mean_a + _PSEUDOCOUNT)


def cohens_d(a: np.ndarray, b: np.ndarray) -> float | None:
    """Pooled-SD standardized mean difference (b - a) / s_pooled.

    Returns None when the pooled standard deviation is zero or undefined
    (e.g. both groups constant, or a group of size < 2).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return None
    var_a = a.var(ddof=1)
    var_b = b.var(ddof=1)
    pooled = ((na - 1) * var_a + (nb - 1) * var_b) / (na + nb - 2)
    if pooled <= 0.0:
        return None
    return float((b.mean() - a.mean()) / math.sqrt(pooled))
