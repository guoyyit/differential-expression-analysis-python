"""API and domain models. Strategy enums make every analysis choice traceable."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, Field


class NormalizationMethod(str, Enum):
    NONE = "none"
    CPM = "cpm"  # counts per million
    LOG_CPM = "log_cpm"  # log2(CPM + 1)
    MEDIAN_RATIO = "median_ratio"  # DESeq-style size factors, then log2(x + 1)


class TestMethod(str, Enum):
    WELCH = "welch"  # Welch two-sample t-test (unequal variances)
    MANN_WHITNEY = "mann_whitney"  # rank-sum, normal approximation with ties


class CorrectionMethod(str, Enum):
    BH = "bh"  # Benjamini-Hochberg FDR
    BONFERRONI = "bonferroni"
    NONE = "none"


class TaskStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    TIMED_OUT = "timed_out"
    INTERRUPTED = "interrupted"  # service stopped mid-task and policy=fail

    @property
    def terminal(self) -> bool:
        return self in _TERMINAL


_TERMINAL = frozenset(
    {
        TaskStatus.COMPLETED,
        TaskStatus.FAILED,
        TaskStatus.CANCELLED,
        TaskStatus.TIMED_OUT,
        TaskStatus.INTERRUPTED,
    }
)


class DesignGroup(BaseModel):
    """One condition: a name plus the sample columns that belong to it."""

    name: str = Field(min_length=1, max_length=128)
    samples: list[str] = Field(min_length=1)


class AnalysisConfig(BaseModel):
    """Analysis strategy. Echoed back in every task snapshot for traceability."""

    normalization: NormalizationMethod = NormalizationMethod.LOG_CPM
    test: TestMethod = TestMethod.WELCH
    correction: CorrectionMethod = CorrectionMethod.BH
    alpha: float = Field(default=0.05, gt=0.0, lt=1.0)
    min_total_count: float = Field(
        default=1.0, ge=0.0,
        description="Prefilter: drop genes whose total raw count is below this.",
    )
    timeout_seconds: float | None = Field(default=None, gt=0.0)
    seed: int = Field(default=0, description="Reserved; pipeline is deterministic.")


class AnalysisRequest(BaseModel):
    """Submission payload: expression matrix plus group design."""

    genes: list[str] = Field(min_length=1)
    samples: list[str] = Field(min_length=2)
    matrix: list[list[float]]  # shape: len(genes) x len(samples)
    groups: list[DesignGroup] = Field(min_length=2, max_length=2)
    config: AnalysisConfig = Field(default_factory=AnalysisConfig)
    task_id: str | None = Field(
        default=None,
        min_length=1,
        max_length=128,
        description="Optional client idempotency key. Re-submitting the same "
        "task_id returns the existing task instead of duplicating work.",
    )


class GeneResult(BaseModel):
    """One row of the differential expression table."""

    gene: str
    statistic: float
    p_value: float
    adjusted_p_value: float
    effect_size: float  # log2 fold change of normalized means, B over A
    cohens_d: float | None  # standardized effect; None when undefined
    mean_group_a: float
    mean_group_b: float
    significant: bool
    rank: int  # 1-based position in the significance/effect ordering


class AnalysisSummary(BaseModel):
    """Aggregate facts about a finished analysis."""

    n_genes_input: int
    n_genes_tested: int
    n_significant: int
    group_a: str
    group_b: str


class TaskSnapshot(BaseModel):
    """Persisted, client-visible state of a task. Never exposes partial results."""

    task_id: str
    status: TaskStatus
    created_at: str
    started_at: str | None = None
    finished_at: str | None = None
    config: AnalysisConfig
    error: dict[str, Any] | None = None
    summary: AnalysisSummary | None = None


class ResultPage(BaseModel):
    """A filtered, sorted slice of a completed task's results."""

    task_id: str
    total: int  # rows matching the filter, before pagination
    offset: int
    limit: int
    sort: str
    order: str
    results: list[GeneResult]


class CompareResponse(BaseModel):
    """Cross-task concordance over genes tested in both tasks."""

    task_a: str
    task_b: str
    n_common_genes: int
    n_both_significant: int
    n_concordant_sign: int  # both significant, same effect direction
    effect_size_correlation: float | None  # Pearson r over common genes
    common_genes: list[str]
