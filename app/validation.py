"""Input validation. Raises ValidationError; never yields partial results.

All checks run before a task is persisted or enqueued, so an invalid
submission cannot leave any trace in the store.
"""

from __future__ import annotations

import math

from app.config import Settings
from app.errors import ValidationError
from app.models import AnalysisRequest


def validate_request(req: AnalysisRequest, settings: Settings) -> None:
    """Check sample counts, group design, matrix shape and missing values."""
    if len(req.genes) > settings.max_genes:
        raise ValidationError(
            "too many genes",
            {"field": "genes", "limit": settings.max_genes, "got": len(req.genes)},
        )
    if len(req.samples) > settings.max_samples:
        raise ValidationError(
            "too many samples",
            {"field": "samples", "limit": settings.max_samples, "got": len(req.samples)},
        )
    _require_unique(req.genes, "genes")
    _require_unique(req.samples, "samples")
    for g in req.genes:
        if not g.strip():
            raise ValidationError("gene names must be non-empty", {"field": "genes"})
    for s in req.samples:
        if not s.strip():
            raise ValidationError("sample names must be non-empty", {"field": "samples"})

    if len(req.matrix) != len(req.genes):
        raise ValidationError(
            "matrix row count must equal number of genes",
            {"field": "matrix", "expected_rows": len(req.genes), "got_rows": len(req.matrix)},
        )
    for i, row in enumerate(req.matrix):
        if len(row) != len(req.samples):
            raise ValidationError(
                "matrix row length must equal number of samples",
                {"field": "matrix", "row": i, "expected_cols": len(req.samples),
                 "got_cols": len(row)},
            )
        for j, v in enumerate(row):
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                raise ValidationError(
                    "matrix entries must be numbers",
                    {"field": "matrix", "row": i, "col": j},
                )
            if math.isnan(v) or math.isinf(v):
                raise ValidationError(
                    "matrix contains missing or non-finite values",
                    {"field": "matrix", "row": i, "col": j, "value": str(v)},
                )
            if v < 0:
                raise ValidationError(
                    "expression counts must be non-negative",
                    {"field": "matrix", "row": i, "col": j, "value": v},
                )

    names = [g.name for g in req.groups]
    _require_unique(names, "groups")
    sample_set = set(req.samples)
    seen: set[str] = set()
    for g in req.groups:
        for s in g.samples:
            if s not in sample_set:
                raise ValidationError(
                    "group references unknown sample",
                    {"field": "groups", "group": g.name, "sample": s},
                )
            if s in seen:
                raise ValidationError(
                    "sample assigned to more than one group",
                    {"field": "groups", "sample": s},
                )
            seen.add(s)

    if (
        req.config.timeout_seconds is not None
        and req.config.timeout_seconds > settings.max_timeout_seconds
    ):
        raise ValidationError(
            "timeout_seconds exceeds server limit",
            {"field": "config.timeout_seconds", "limit": settings.max_timeout_seconds},
        )


def group_column_indices(req: AnalysisRequest) -> tuple[list[int], list[int]]:
    """Map each group's sample names to matrix column indices."""
    pos = {name: i for i, name in enumerate(req.samples)}
    a = [pos[s] for s in req.groups[0].samples]
    b = [pos[s] for s in req.groups[1].samples]
    return a, b


def _require_unique(items: list[str], field: str) -> None:
    if len(set(items)) != len(items):
        raise ValidationError("duplicate entries", {"field": field})
