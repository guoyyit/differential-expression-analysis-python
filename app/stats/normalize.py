"""Library-size normalization strategies. Pure functions, deterministic."""

from __future__ import annotations

import numpy as np

from app.models import NormalizationMethod


def normalize(counts: np.ndarray, method: NormalizationMethod) -> np.ndarray:
    """Return a normalized matrix with the same shape as ``counts``.

    Rows are genes, columns are samples. CPM/LOG_CPM scale by library size;
    MEDIAN_RATIO uses DESeq-style size factors and returns log2(x + 1).
    """
    x = np.asarray(counts, dtype=np.float64)
    if method is NormalizationMethod.NONE:
        return x
    if method is NormalizationMethod.CPM:
        return _cpm(x)
    if method is NormalizationMethod.LOG_CPM:
        return np.log2(_cpm(x) + 1.0)
    if method is NormalizationMethod.MEDIAN_RATIO:
        sf = size_factors_median_ratio(x)
        return np.log2(x / sf[np.newaxis, :] + 1.0)
    raise ValueError(f"unknown normalization method: {method!r}")


def _cpm(x: np.ndarray) -> np.ndarray:
    lib = x.sum(axis=0)
    # A zero-library sample would divide by zero; keep its column at zero.
    safe = np.where(lib > 0.0, lib, 1.0)
    out = x / safe[np.newaxis, :] * 1.0e6
    out[:, lib <= 0.0] = 0.0
    return out


def size_factors_median_ratio(counts: np.ndarray) -> np.ndarray:
    """DESeq-style median-of-ratios size factors, one per sample.

    Uses genes with strictly positive counts in every sample. If no such
    gene exists, falls back to library-size factors so the pipeline still
    produces a defined, reproducible result.
    """
    x = np.asarray(counts, dtype=np.float64)
    positive = np.all(x > 0.0, axis=1)
    if np.any(positive):
        sub = x[positive]
        log_geo = np.log(sub).mean(axis=1)  # log geometric mean per gene
        ratios = sub / np.exp(log_geo)[:, np.newaxis]
        sf = np.median(ratios, axis=0)
    else:
        lib = x.sum(axis=0)
        mean_lib = lib.mean() if lib.size else 0.0
        sf = lib / mean_lib if mean_lib > 0.0 else np.ones(x.shape[1])
    sf = np.where(sf > 0.0, sf, 1.0)
    return sf
