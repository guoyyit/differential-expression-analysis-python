"""Multiple testing correction procedures.

Each function takes p-values in gene order and returns adjusted p-values in
the same order. Sorting is deterministic: ties are broken by original index.
"""
from __future__ import annotations

import math

from .models import CorrectionMethod


def _ranked_indices(pvals: list[float]) -> list[int]:
    return sorted(range(len(pvals)), key=lambda i: (pvals[i], i))


def bonferroni(pvals: list[float]) -> list[float]:
    m = len(pvals)
    return [min(1.0, p * m) for p in pvals]


def holm(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = _ranked_indices(pvals)
    adj = [0.0] * m
    running = 0.0
    # smallest p first; adjusted = max over prior ranks of (m-k)*p, enforced monotone
    for k, idx in enumerate(order):
        val = (m - k) * pvals[idx]
        running = max(running, val)
        adj[idx] = min(1.0, running)
    return adj


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = _ranked_indices(pvals)
    adj = [0.0] * m
    running = 1.0
    # walk from largest p downwards, enforcing monotonicity
    for k in range(m - 1, -1, -1):
        idx = order[k]
        val = m / (k + 1) * pvals[idx]
        running = min(running, val)
        adj[idx] = min(1.0, running)
    return adj


def benjamini_yekutieli(pvals: list[float]) -> list[float]:
    m = len(pvals)
    c_m = math.fsum(1.0 / k for k in range(1, m + 1)) if m else 1.0
    order = _ranked_indices(pvals)
    adj = [0.0] * m
    running = 1.0
    for k in range(m - 1, -1, -1):
        idx = order[k]
        val = c_m * m / (k + 1) * pvals[idx]
        running = min(running, val)
        adj[idx] = min(1.0, running)
    return adj


_CORRECTIONS = {
    CorrectionMethod.BONFERRONI: bonferroni,
    CorrectionMethod.HOLM: holm,
    CorrectionMethod.BH: benjamini_hochberg,
    CorrectionMethod.BY: benjamini_yekutieli,
}


def correct(method: CorrectionMethod, pvals: list[float]) -> list[float]:
    return _CORRECTIONS[method](pvals)
