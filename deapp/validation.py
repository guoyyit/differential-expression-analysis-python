"""Strict semantic validation of submitted matrices and designs.

Validation happens *before* a task is created, so rejected input never
produces a task record or partial results.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .errors import ValidationError
from .models import CreateTaskRequest, NormalizationMethod, TestMethod


@dataclass(frozen=True)
class ValidatedInput:
    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    control_group: str
    treatment_group: str
    control_samples: tuple[str, ...]
    treatment_samples: tuple[str, ...]
    # values[gene_index] -> (control values, treatment values) in sample order
    values: tuple[tuple[tuple[float, ...], tuple[float, ...]], ...]


def validate_request(req: CreateTaskRequest) -> ValidatedInput:
    design = req.design
    if not design:
        raise ValidationError(
            "design must contain at least two non-empty groups",
            details={"reason": "empty_design"},
        )

    # --- group categories -------------------------------------------------
    groups = list(design.keys())
    if len(groups) != 2:
        raise ValidationError(
            "exactly two group categories are required",
            details={"reason": "invalid_group_count", "groups": groups},
        )
    for g, samples in design.items():
        if not samples:
            raise ValidationError(
                f"group {g!r} has no samples",
                details={"reason": "empty_group", "group": g},
            )

    control, treatment = _resolve_comparison(groups, req.control, req.treatment)

    # --- sample membership -------------------------------------------------
    matrix_samples = req.matrix.sample_ids
    sample_set = set(matrix_samples)
    if len(sample_set) != len(matrix_samples):
        raise ValidationError(
            "duplicate sample ids in matrix header",
            details={"reason": "duplicate_sample"},
        )

    seen: set[str] = set()
    ordered: dict[str, list[str]] = {}
    for group in (control, treatment):
        members = design[group]
        if len(set(members)) != len(members):
            raise ValidationError(
                f"duplicate sample inside group {group!r}",
                details={"reason": "duplicate_sample_in_group", "group": group},
            )
        for s in members:
            if s in seen:
                raise ValidationError(
                    f"sample {s!r} appears in more than one group",
                    details={"reason": "sample_in_multiple_groups", "sample": s},
                )
            if s not in sample_set:
                raise ValidationError(
                    f"sample {s!r} from design is missing in matrix",
                    details={"reason": "unknown_sample", "sample": s},
                )
            seen.add(s)
        ordered[group] = list(members)

    # paired test requires equal, positive sizes
    if req.config.test == TestMethod.PAIRED:
        if len(ordered[control]) != len(ordered[treatment]):
            raise ValidationError(
                "paired test requires equal sample counts in both groups",
                details={
                    "reason": "unequal_paired_groups",
                    "n_control": len(ordered[control]),
                    "n_treatment": len(ordered[treatment]),
                },
            )
    else:
        for group in (control, treatment):
            if len(ordered[group]) < 2:
                raise ValidationError(
                    f"group {group!r} needs at least 2 samples",
                    details={"reason": "not_enough_samples", "group": group},
                )

    # --- matrix dimensions -------------------------------------------------
    gene_ids = req.matrix.gene_ids
    if len(set(gene_ids)) != len(gene_ids):
        raise ValidationError(
            "duplicate gene ids",
            details={"reason": "duplicate_gene"},
        )
    n_genes, n_samples = len(gene_ids), len(matrix_samples)
    values = req.matrix.values
    if len(values) != n_genes:
        raise ValidationError(
            "matrix row count does not match gene id count",
            details={
                "reason": "dimension_mismatch",
                "expected_rows": n_genes,
                "actual_rows": len(values),
            },
        )
    for i, row in enumerate(values):
        if len(row) != n_samples:
            raise ValidationError(
                f"row {i} length does not match sample count",
                details={
                    "reason": "dimension_mismatch",
                    "row": i,
                    "expected_columns": n_samples,
                    "actual_columns": len(row),
                },
            )
        for j, v in enumerate(row):
            if not isinstance(v, (int, float)) or isinstance(v, bool) or math.isnan(float(v)):
                raise ValidationError(
                    f"missing/non-numeric value at gene {gene_ids[i]!r}, "
                    f"sample {matrix_samples[j]!r}",
                    details={
                        "reason": "missing_value",
                        "gene": gene_ids[i],
                        "sample": matrix_samples[j],
                    },
                )
            if not math.isfinite(float(v)):
                raise ValidationError(
                    f"non-finite value at gene {gene_ids[i]!r}, "
                    f"sample {matrix_samples[j]!r}",
                    details={
                        "reason": "non_finite_value",
                        "gene": gene_ids[i],
                        "sample": matrix_samples[j],
                    },
                )

    # --- strategy-specific requirements -----------------------------------
    if req.config.normalization == NormalizationMethod.TPM:
        lengths = req.config.gene_lengths
        if lengths is None:
            raise ValidationError(
                "TPM normalization requires gene_lengths",
                details={"reason": "missing_gene_lengths"},
            )
        missing = [g for g in gene_ids if g not in lengths]
        if missing:
            raise ValidationError(
                "gene_lengths missing entries for some genes",
                details={"reason": "incomplete_gene_lengths", "missing": missing[:10]},
            )
        bad = [g for g in gene_ids if not math.isfinite(lengths[g]) or lengths[g] <= 0]
        if bad:
            raise ValidationError(
                "gene_lengths must be positive finite numbers",
                details={"reason": "invalid_gene_length", "genes": bad[:10]},
            )

    # --- project values into control/treatment order ----------------------
    idx = {s: j for j, s in enumerate(matrix_samples)}
    ci = [idx[s] for s in ordered[control]]
    ti = [idx[s] for s in ordered[treatment]]
    projected = tuple(
        (
            tuple(float(values[g][j]) for j in ci),
            tuple(float(values[g][j]) for j in ti),
        )
        for g in range(n_genes)
    )

    return ValidatedInput(
        gene_ids=tuple(gene_ids),
        sample_ids=tuple(matrix_samples),
        control_group=control,
        treatment_group=treatment,
        control_samples=tuple(ordered[control]),
        treatment_samples=tuple(ordered[treatment]),
        values=projected,
    )


def _resolve_comparison(
    groups: list[str], control: str | None, treatment: str | None
) -> tuple[str, str]:
    if control is None and treatment is None:
        c, t = sorted(groups)
        return c, t
    if control is None or treatment is None:
        raise ValidationError(
            "control and treatment must be provided together",
            details={"reason": "partial_comparison"},
        )
    if control not in groups or treatment not in groups:
        raise ValidationError(
            "control/treatment must be group names in the design",
            details={"reason": "unknown_group", "groups": groups},
        )
    if control == treatment:
        raise ValidationError(
            "control and treatment must differ",
            details={"reason": "same_group"},
        )
    return control, treatment
