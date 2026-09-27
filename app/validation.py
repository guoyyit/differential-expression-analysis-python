"""Input validation. Raises ValidationError; never yields partial results."""

from __future__ import annotations

from app.config import Settings
from app.models import AnalysisRequest


def validate_request(req: AnalysisRequest, settings: Settings) -> None:
    """Check sample counts, group design, matrix shape and missing values."""
    raise NotImplementedError


def group_column_indices(req: AnalysisRequest) -> tuple[list[int], list[int]]:
    """Map each group's sample names to matrix column indices."""
    raise NotImplementedError
