"""End-to-end analysis: validate -> normalize -> test -> correct -> rank.

Determinism contract: per-gene statistics are computed independently in
input gene order; correction and ranking are total and stable; no random
number generation is used anywhere, so identical inputs always produce
byte-identical outputs regardless of worker scheduling.
"""

from __future__ import annotations

import threading
import time
from collections.abc import Callable

import numpy as np

from app.models import AnalysisConfig, GeneResult
from app.stats.correction import adjust_pvalues
from app.stats.effect import cohens_d, log2_fold_change
from app.stats.normalize import normalize
from app.stats.tests import run_test


class Cancelled(Exception):
    """Raised inside the pipeline when the cancel token is set."""


def check_cancel(token: threading.Event) -> None:
    if token.is_set():
        raise Cancelled()


def run_analysis(
    genes: list[str],
    matrix: list[list[float]],
    group_a_idx: list[int],
    group_b_idx: list[int],
    config: AnalysisConfig,
    cancel_token: threading.Event,
    progress: Callable[[int, int], None] | None = None,
    sleep_per_gene: float = 0.0,
) -> list[GeneResult]:
    """Run the full pipeline and return ranked results (rank field set)."""
    check_cancel(cancel_token)
    counts = np.asarray(matrix, dtype=np.float64)
    normalized = normalize(counts, config.normalization)
    check_cancel(cancel_token)

    totals = counts.sum(axis=1)
    keep = np.flatnonzero(totals >= config.min_total_count)
    n_kept = int(keep.size)

    stats = np.empty(n_kept, dtype=np.float64)
    pvals = np.empty(n_kept, dtype=np.float64)
    effects: list[float] = []
    cohens: list[float | None] = []
    means_a: list[float] = []
    means_b: list[float] = []

    step = max(1, n_kept // 100) if n_kept else 1
    for row, gene_idx in enumerate(keep):
        if row % step == 0:
            check_cancel(cancel_token)
            if progress is not None:
                progress(row, n_kept)
        if sleep_per_gene > 0.0:
            time.sleep(sleep_per_gene)
        vec = normalized[gene_idx]
        a = vec[group_a_idx]
        b = vec[group_b_idx]
        stat, p = run_test(a, b, config.test)
        stats[row] = stat
        pvals[row] = p
        ma = float(a.mean())
        mb = float(b.mean())
        means_a.append(ma)
        means_b.append(mb)
        effects.append(log2_fold_change(ma, mb))
        cohens.append(cohens_d(a, b))

    check_cancel(cancel_token)
    adjusted = adjust_pvalues(pvals, config.correction)

    # Total, stable ordering: significance first, then raw p, then magnitude
    # of effect, then gene name — so ties can never reorder output.
    order = sorted(
        range(n_kept),
        key=lambda r: (
            adjusted[r],
            pvals[r],
            -abs(effects[r]),
            genes[int(keep[r])],
        ),
    )

    results: list[GeneResult] = []
    for rank, r in enumerate(order, start=1):
        gene_idx = int(keep[r])
        results.append(
            GeneResult(
                gene=genes[gene_idx],
                statistic=float(stats[r]),
                p_value=float(pvals[r]),
                adjusted_p_value=float(adjusted[r]),
                effect_size=float(effects[r]),
                cohens_d=cohens[r],
                mean_group_a=means_a[r],
                mean_group_b=means_b[r],
                significant=bool(adjusted[r] <= config.alpha),
                rank=rank,
            )
        )
    if progress is not None:
        progress(n_kept, n_kept)
    return results
