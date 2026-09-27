"""Pure computation pipeline: normalization -> per-gene tests -> correction.

Output ordering is fully deterministic:
rank key = (adjusted p asc, |log2fc| desc, |statistic| desc, gene_id asc).
"""
from __future__ import annotations

from typing import Callable

from . import __version__, multipletest, normalization, statistics
from .models import AnalysisConfig, ResultResponse, ResultRow, StrategyInfo, TaskStatus
from .validation import ValidatedInput

CancellationChecker = Callable[[], bool]


class TaskCancelled(Exception):
    """Raised internally when the cancellation/timeout checker trips."""


def run_analysis(
    task_id: str,
    data: ValidatedInput,
    config: AnalysisConfig,
    is_cancelled: CancellationChecker | None = None,
) -> ResultResponse:
    matrix = normalization.normalize(config.normalization, data, config)

    gene_stats: list[statistics.GeneStat] = []
    for i, gene_id in enumerate(data.gene_ids):
        if is_cancelled is not None and is_cancelled():
            raise TaskCancelled()
        control, treatment = matrix[i]
        gene_stats.append(statistics.gene_stat(config.test, gene_id, control, treatment))

    pvals = [gs.p_value for gs in gene_stats]
    adjusted = multipletest.correct(config.correction, pvals)

    rows = [
        ResultRow(
            gene_id=gs.gene_id,
            mean_control=gs.mean_control,
            mean_treatment=gs.mean_treatment,
            log2fc=gs.log2fc,
            statistic=gs.statistic,
            p_value=gs.p_value,
            p_value_adjusted=adjusted[i],
            rank=0,
        )
        for i, gs in enumerate(gene_stats)
    ]
    rows.sort(
        key=lambda r: (
            r.p_value_adjusted,
            -abs(r.log2fc),
            -abs(r.statistic),
            r.gene_id,
        )
    )
    for i, row in enumerate(rows, start=1):
        row.rank = i

    strategy = StrategyInfo(
        normalization=config.normalization.value,
        test=config.test.value,
        correction=config.correction.value,
        effect=config.effect,
        seed=config.seed,
        parameters={
            "pseudocount": statistics.PSEUDOCOUNT,
            "control_group": data.control_group,
            "treatment_group": data.treatment_group,
            "n_control": len(data.control_samples),
            "n_treatment": len(data.treatment_samples),
            "engine_version": __version__,
        },
    )
    return ResultResponse(
        task_id=task_id,
        status=TaskStatus.SUCCEEDED,
        strategy=strategy,
        n_genes=len(rows),
        rows=rows,
    )
