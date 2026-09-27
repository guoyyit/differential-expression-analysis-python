"""Per-gene two-group differential tests. Deterministic, no randomness."""

from __future__ import annotations

import numpy as np

from app.models import TestMethod


def run_test(a: np.ndarray, b: np.ndarray, method: TestMethod) -> tuple[float, float]:
    """Return (statistic, p_value) for group vectors ``a`` and ``b``."""
    raise NotImplementedError


def welch_t(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Welch's t-test: returns (t statistic, two-sided p value)."""
    raise NotImplementedError


def mann_whitney(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Rank-sum test, normal approximation with tie correction."""
    raise NotImplementedError
