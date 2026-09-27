"""Effect sizes on normalized expression values."""

from __future__ import annotations

import numpy as np


def log2_fold_change(mean_a: float, mean_b: float) -> float:
    """log2(B/A) with a pseudocount so zero means stay finite."""
    raise NotImplementedError


def cohens_d(a: np.ndarray, b: np.ndarray) -> float | None:
    """Pooled-SD standardized mean difference; None when undefined."""
    raise NotImplementedError
