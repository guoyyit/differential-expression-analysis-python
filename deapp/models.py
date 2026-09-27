"""Pydantic schemas for API requests, responses and persisted records."""
from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field, field_serializer

# ---------------------------------------------------------------------------
# Enums – strategy choices are part of the public, traceable configuration
# ---------------------------------------------------------------------------


class NormalizationMethod(str, Enum):
    NONE = "none"          # use raw counts / values as given
    CPM = "cpm"            # counts per million (library size scaling)
    TPM = "tpm"            # transcripts per million (needs gene_lengths)
    QUANTILE = "quantile"  # quantile normalization across samples
    TMM = "tmm"            # trimmed mean of M values (edgeR-style)


class TestMethod(str, Enum):
    WELCH = "welch-ttest"                 # Welch's unequal variance t-test
    STUDENT = "student-ttest"             # pooled variance t-test
    PAIRED = "paired-ttest"               # paired t-test (equal n, matched ids)
    MANN_WHITNEY = "mann-whitney"         # Wilcoxon rank-sum (no permutation)


class CorrectionMethod(str, Enum):
    BONFERRONI = "bonferroni"
    BH = "benjamini-hochberg"             # Benjamini-Hochberg FDR
    BY = "benjamini-yekutieli"            # Benjamini-Yekutieli FDR
    HOLM = "holm"                         # Holm step-down FWER


class TaskStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"


# ---------------------------------------------------------------------------
# Request payloads
# ---------------------------------------------------------------------------


class ExpressionMatrix(BaseModel):
    """Row = gene, column = sample."""

    gene_ids: list[str] = Field(min_length=1)
    sample_ids: list[str] = Field(min_length=1)
    values: list[list[float]]  # values[i][j] for gene_ids[i], sample_ids[j]


class AnalysisConfig(BaseModel):
    normalization: NormalizationMethod = NormalizationMethod.NONE
    test: TestMethod = TestMethod.WELCH
    correction: CorrectionMethod = CorrectionMethod.BH
    effect: str = "log2fc"  # currently the single supported effect size
    seed: int = 0
    timeout_seconds: float | None = Field(default=None, gt=0)
    gene_lengths: dict[str, float] | None = None  # required for TPM


class CreateTaskRequest(BaseModel):
    name: str | None = None
    # canonical design: group name -> ordered list of sample ids
    design: dict[str, list[str]]
    # which group is treated vs control; inferred alphabetically if omitted
    control: str | None = None
    treatment: str | None = None
    matrix: ExpressionMatrix
    config: AnalysisConfig = Field(default_factory=AnalysisConfig)


# ---------------------------------------------------------------------------
# Result / task records
# ---------------------------------------------------------------------------


class ResultRow(BaseModel):
    gene_id: str
    mean_control: float
    mean_treatment: float
    log2fc: float
    statistic: float
    p_value: float
    p_value_adjusted: float
    rank: int

    @field_serializer(
        "mean_control",
        "mean_treatment",
        "log2fc",
        "statistic",
        "p_value",
        "p_value_adjusted",
        when_used="json",
    )
    def _serialize_finite(self, v: float) -> float:
        """JSON has no infinity/NaN token.

        +/- infinity (only produced by zero-variance, clearly-separated
        groups) is encoded with a signed 1e308 sentinel; any defensive NaN
        becomes 0.0. Semantics: |statistic| == 1e308 with p == 0 means an
        infinitely separated, perfectly constant contrast.
        """
        import math

        if math.isnan(v):
            return 0.0
        if math.isinf(v):
            return 1e308 if v > 0 else -1e308
        return v


class StrategyInfo(BaseModel):
    normalization: str
    test: str
    correction: str
    effect: str
    seed: int
    parameters: dict[str, Any] = Field(default_factory=dict)


class TaskRecord(BaseModel):
    id: str
    name: str | None = None
    status: TaskStatus
    created_at: str
    started_at: str | None = None
    finished_at: str | None = None
    control: str | None = None
    treatment: str | None = None
    strategy: StrategyInfo | None = None
    error: dict[str, Any] | None = None
    n_genes: int | None = None
    n_samples: int | None = None


class TaskDetail(TaskRecord):
    design: dict[str, list[str]] | None = None


class ResultResponse(BaseModel):
    task_id: str
    status: TaskStatus
    strategy: StrategyInfo
    n_genes: int
    rows: list[ResultRow]
