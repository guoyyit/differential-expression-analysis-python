"""Library-size normalization strategies. Pure functions, deterministic."""

from __future__ import annotations

import numpy as np

from app.models import NormalizationMethod


def normalize(counts: np.ndarray, method: NormalizationMethod) -> np.ndarray:
    """Return a normalized, log-scale matrix with the same shape as ``counts``.

    Rows are genes, columns are samples. Output ordering matches input.
    """
    raise NotImplementedError


def size_factors_median_ratio(counts: np.ndarray) -> np.ndarray:
    """DESeq-style median-of-ratios size factors, one per sample."""
    raise NotImplementedError
