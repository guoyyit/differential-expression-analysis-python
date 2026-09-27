"""Deterministic, dependency-free normalization strategies.

A normalizer receives the per-gene values split into control / treatment and
returns the same structure after scaling. Column order inside each group is
preserved; all computations are plain Python ``float`` so results are bit
identical for identical input regardless of thread scheduling.
"""
from __future__ import annotations

import math
from typing import Callable

from .models import AnalysisConfig, NormalizationMethod
from .validation import ValidatedInput

Matrix = tuple[tuple[tuple[float, ...], tuple[float, ...]], ...]


def normalize(
    method: NormalizationMethod,
    data: ValidatedInput,
    config: AnalysisConfig,
) -> Matrix:
    fn = _NORMALIZERS[method]
    return fn(data, config)


def _none(data: ValidatedInput, config: AnalysisConfig) -> Matrix:
    return data.values


def _columns(data: ValidatedInput) -> tuple[list[float], int]:
    """Build full sample-major matrix: columns ordered control then treatment."""
    nc = len(data.control_samples)
    nt = len(data.treatment_samples)
    cols: list[list[float]] = [
        [data.values[g][0][j] for g in range(len(data.gene_ids))]
        for j in range(nc)
    ] + [
        [data.values[g][1][j] for g in range(len(data.gene_ids))]
        for j in range(nt)
    ]
    return cols, nc


def _split_back(
    cols: list[list[float]], nc: int, n_genes: int
) -> Matrix:
    return tuple(
        (
            tuple(cols[j][g] for j in range(nc)),
            tuple(cols[j][g] for j in range(nc, len(cols))),
        )
        for g in range(n_genes)
    )


def _cpm(data: ValidatedInput, config: AnalysisConfig) -> Matrix:
    cols, nc = _columns(data)
    scaled: list[list[float]] = []
    for col in cols:
        lib = math.fsum(col)
        if lib <= 0:
            scaled.append([0.0] * len(col))
        else:
            factor = 1_000_000.0 / lib
            scaled.append([v * factor for v in col])
    return _split_back(scaled, nc, len(data.gene_ids))


def _tpm(data: ValidatedInput, config: AnalysisConfig) -> Matrix:
    assert config.gene_lengths is not None
    lengths = [config.gene_lengths[g] for g in data.gene_ids]
    cols, nc = _columns(data)
    scaled: list[list[float]] = []
    for col in cols:
        rpk = [v / lengths[i] * 1000.0 for i, v in enumerate(col)]
        total = math.fsum(rpk)
        if total <= 0:
            scaled.append([0.0] * len(col))
        else:
            factor = 1_000_000.0 / total
            scaled.append([r * factor for r in rpk])
    return _split_back(scaled, nc, len(data.gene_ids))


def _quantile(data: ValidatedInput, config: AnalysisConfig) -> Matrix:
    """Classical quantile normalization (mean of per-rank values).

    Ties inside a column receive the average of the rank means they occupy,
    so equal input values map to equal output values.
    """
    cols, nc = _columns(data)
    n_genes = len(data.gene_ids)
    # Rank means: sort every column (stable), average across columns per rank.
    sorted_cols = [sorted(col) for col in cols]
    rank_means = [
        math.fsum(sc[r] for sc in sorted_cols) / len(sorted_cols)
        for r in range(n_genes)
    ]
    normalized: list[list[float]] = []
    for col in cols:
        order = sorted(range(n_genes), key=lambda i: col[i])  # stable on ties
        out = [0.0] * n_genes
        k = 0
        while k < n_genes:
            end = k + 1
            while end < n_genes and col[order[end]] == col[order[k]]:
                end += 1
            if end - k == 1:
                out[order[k]] = rank_means[k]
            else:
                tie_mean = math.fsum(rank_means[k:end]) / (end - k)
                for r in range(k, end):
                    out[order[r]] = tie_mean
            k = end
        normalized.append(out)
    return _split_back(normalized, nc, n_genes)


def _upper_quartile(col: list[float]) -> float:
    s = sorted(col)
    n = len(s)
    if n < 4:
        return math.fsum(s) / n
    # median of the upper half (inclusive-of-middle convention, deterministic)
    half = s[n // 2:]
    m = len(half)
    return half[m // 2] if m % 2 else 0.5 * (half[m // 2 - 1] + half[m // 2])


def _tmm(data: ValidatedInput, config: AnalysisConfig) -> Matrix:
    """edgeR-style trimmed mean of M values followed by library scaling.

    Reference sample = upper-quartile closest to the mean upper quartile
    (ties broken by sample order). M trimmed 30%, A trimmed 5%; weights use
    binomial-type asymptotic variances; factors are rescaled so their
    geometric mean is 1.
    """
    cols, nc = _columns(data)
    n_samples = len(cols)
    n_genes = len(data.gene_ids)
    libs = [math.fsum(c) for c in cols]

    q3 = [_upper_quartile(c) for c in cols]
    target = math.fsum(q3) / n_samples
    ref_idx = min(range(n_samples), key=lambda i: (abs(q3[i] - target), i))
    ref = cols[ref_idx]
    ref_lib = libs[ref_idx]

    log2 = math.log2
    factors: list[float] = []
    for j in range(n_samples):
        if j == ref_idx or libs[j] <= 0:
            factors.append(1.0)
            continue
        obs: list[tuple[float, float, float]] = []
        for g in range(n_genes):
            y, yr = cols[j][g], ref[g]
            if y <= 0 or yr <= 0:
                continue
            yg = y / libs[j]
            yrg = yr / ref_lib
            m = log2(yg / yrg)
            a = 0.5 * log2(yg * yrg)
            # var ~ 1/y - 1/lib + 1/yr - 1/reflib (asymptotic binomial);
            # single-gene libraries give zero variance -> skip weighting
            w = (1.0 / y - 1.0 / libs[j]) + (1.0 / yr - 1.0 / ref_lib)
            if not (w > 0.0):
                continue
            w = 1.0 / w
            obs.append((m, a, w))
        if not obs:
            factors.append(1.0)
            continue
        ms = sorted(m for m, _, _ in obs)
        als = sorted(a for _, a, _ in obs)
        n = len(obs)
        lo_m = ms[int(0.30 * (n - 1))]
        hi_m = ms[n - 1 - int(0.30 * (n - 1))]
        lo_a = als[int(0.05 * (n - 1))]
        hi_a = als[n - 1 - int(0.05 * (n - 1))]
        num = 0.0
        den = 0.0
        for m, a, w in obs:
            if lo_m <= m <= hi_m and lo_a <= a <= hi_a:
                num += w * m
                den += w
        factors.append(2.0 ** (num / den) if den > 0 else 1.0)

    # rescale factors so their geometric mean is exactly 1
    log_mean = math.fsum(math.log(f) for f in factors) / n_samples
    factors = [math.exp(math.log(f) - log_mean) for f in factors]

    scaled = []
    for j in range(n_samples):
        if libs[j] <= 0:
            scaled.append([0.0] * n_genes)
            continue
        scaled.append([
            cols[j][g] / (libs[j] * factors[j]) * 1_000_000.0
            for g in range(n_genes)
        ])
    return _split_back(scaled, nc, n_genes)


_NORMALIZERS: dict[NormalizationMethod, Callable[[ValidatedInput, AnalysisConfig], Matrix]] = {
    NormalizationMethod.NONE: _none,
    NormalizationMethod.CPM: _cpm,
    NormalizationMethod.TPM: _tpm,
    NormalizationMethod.QUANTILE: _quantile,
    NormalizationMethod.TMM: _tmm,
}
