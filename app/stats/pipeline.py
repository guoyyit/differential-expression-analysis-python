"""End-to-end analysis: validate -> normalize -> test -> correct -> rank."""

from __future__ import annotations

import threading
from collections.abc import Callable

from app.models import AnalysisConfig, GeneResult


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
    raise NotImplementedError
