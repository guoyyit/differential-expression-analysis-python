"""Deterministic differential-expression tests and effect sizes.

All math is pure-Python (``math.fsum`` / continued fractions) so that the
numeric output is identical on every machine and under any thread schedule.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .models import TestMethod

PSEUDOCOUNT = 1e-300  # floor used when a group mean is exactly zero so that
# log2 fold change stays finite; documented in README


@dataclass(frozen=True)
class GeneStat:
    gene_id: str
    mean_control: float
    mean_treatment: float
    log2fc: float
    statistic: float
    p_value: float


# ---------------------------------------------------------------------------
# Statistical building blocks
# ---------------------------------------------------------------------------


def _mean(xs: tuple[float, ...]) -> float:
    return math.fsum(xs) / len(xs)


def _sample_variance(xs: tuple[float, ...], mean: float) -> float:
    if len(xs) < 2:
        return 0.0
    return math.fsum((x - mean) ** 2 for x in xs) / (len(xs) - 1)


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta (Lentz, Numerical Recipes)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t: float, df: float) -> float:
    """Two-sided survival p-value of Student's t: 2 * P(T >= |t|)."""
    if df <= 0 or not math.isfinite(t):
        return 0.0 if math.isinf(t) else 1.0
    x = df / (df + t * t)
    return _betai(0.5 * df, 0.5, x)


def _normal_sf_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


PSEUDOCOUNT = 1e-300  # floor used when a group mean is exactly zero so that
# log2 fold change stays finite; documented in README


def _effect(mc: float, mt: float) -> float:
    """log2((treatment+eps)/(control+eps)).

    Inputs are normalized expression values (or raw values when
    normalization is "none"). The tiny floor only matters for all-zero
    groups; for ordinary counts the fold change is the plain ratio.
    """
    return math.log2((mt + PSEUDOCOUNT) / (mc + PSEUDOCOUNT))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _welch(c: tuple[float, ...], t: tuple[float, ...]) -> tuple[float, float]:
    mc, mt = _mean(c), _mean(t)
    vc, vt = _sample_variance(c, mc), _sample_variance(t, mt)
    se2 = vc / len(c) + vt / len(t)
    diff = mt - mc
    if se2 <= 0.0:
        return _zero_scale_statistic(diff)
    stat = diff / math.sqrt(se2)
    df_num = se2 * se2
    df_den = (
        (vc / len(c)) ** 2 / (len(c) - 1) + (vt / len(t)) ** 2 / (len(t) - 1)
    )
    if df_den <= 0.0:
        return _zero_scale_statistic(diff)
    df = df_num / df_den
    return stat, t_sf_two_sided(stat, df)


def _zero_scale_statistic(diff: float) -> tuple[float, float]:
    """No within-group variance: exact equality -> p=1, else t=inf, p=0."""
    if diff == 0.0:
        return 0.0, 1.0
    return math.inf if diff > 0 else -math.inf, 0.0


def _student(c: tuple[float, ...], t: tuple[float, ...]) -> tuple[float, float]:
    nc, nt = len(c), len(t)
    mc, mt = _mean(c), _mean(t)
    vc, vt = _sample_variance(c, mc), _sample_variance(t, mt)
    sp2 = ((nc - 1) * vc + (nt - 1) * vt) / (nc + nt - 2)
    se2 = sp2 * (1.0 / nc + 1.0 / nt)
    diff = mt - mc
    if se2 <= 0.0:
        return _zero_scale_statistic(diff)
    stat = diff / math.sqrt(se2)
    return stat, t_sf_two_sided(stat, nc + nt - 2)


def _paired(c: tuple[float, ...], t: tuple[float, ...]) -> tuple[float, float]:
    diffs = tuple(t[i] - c[i] for i in range(len(c)))
    md = _mean(diffs)
    sd = math.sqrt(_sample_variance(diffs, md))
    if sd == 0.0:
        return _zero_scale_statistic(md)
    stat = md / (sd / math.sqrt(len(diffs)))
    return stat, t_sf_two_sided(stat, len(diffs) - 1)


def _mann_whitney(
    c: tuple[float, ...], t: tuple[float, ...]
) -> tuple[float, float]:
    """Wilcoxon rank-sum with tie correction + normal approximation.

    Statistic reported is the standardized z (signed, control -> treatment).
    """
    n1, n2 = len(c), len(t)
    pooled = list(c) + list(t)
    order = sorted(range(n1 + n2), key=lambda i: pooled[i])
    ranks = [0.0] * (n1 + n2)
    k = 0
    tie_sum_sq = 0.0
    while k < len(pooled):
        end = k + 1
        while end < len(pooled) and pooled[order[end]] == pooled[order[k]]:
            end += 1
        avg_rank = (k + 1 + end) / 2.0  # ranks are 1-based: (k+1 .. end) mean
        for r in range(k, end):
            ranks[order[r]] = avg_rank
        tie_n = end - k
        tie_sum_sq += tie_n ** 3 - tie_n
        k = end
    r1 = math.fsum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    n = n1 + n2
    var_u = n1 * n2 / 12.0 * ((n + 1) - tie_sum_sq / (n * (n - 1)))
    if var_u <= 0:
        return 0.0, 1.0
    # continuity correction
    if u1 > mean_u:
        z = (u1 - mean_u - 0.5) / math.sqrt(var_u)
    elif u1 < mean_u:
        z = (u1 - mean_u + 0.5) / math.sqrt(var_u)
    else:
        z = 0.0
    # return sign-positive statistic magnitude? Keep signed z for traceability
    return z, _normal_sf_two_sided(z)


_TESTS = {
    TestMethod.WELCH: _welch,
    TestMethod.STUDENT: _student,
    TestMethod.PAIRED: _paired,
    TestMethod.MANN_WHITNEY: _mann_whitney,
}


def gene_stat(
    method: TestMethod, gene_id: str, control: tuple[float, ...], treatment: tuple[float, ...]
) -> GeneStat:
    mc, mt = _mean(control), _mean(treatment)
    stat, p = _TESTS[method](control, treatment)
    p = min(1.0, max(0.0, p))
    return GeneStat(
        gene_id=gene_id,
        mean_control=mc,
        mean_treatment=mt,
        log2fc=_effect(mc, mt),
        statistic=stat,
        p_value=p,
    )
