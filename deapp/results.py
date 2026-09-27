"""Result querying: deterministic filtering, ordering, comparison, export."""
from __future__ import annotations

import csv
import io
import json
from typing import Iterable

from .errors import InvalidFilterError, ValidationError
from .models import ResultResponse, ResultRow

_EXPORT_COLUMNS = [
    "rank",
    "gene_id",
    "mean_control",
    "mean_treatment",
    "log2fc",
    "statistic",
    "p_value",
    "p_value_adjusted",
]


def filter_rows(
    result: ResultResponse,
    *,
    gene: str | None = None,
    gene_contains: str | None = None,
    genes: list[str] | None = None,
    padj_max: float | None = None,
    pvalue_max: float | None = None,
    log2fc_min_abs: float | None = None,
    direction: str | None = None,  # "up" | "down"
    top: int | None = None,
) -> list[ResultRow]:
    """Rows are always returned in the stored rank order (significance)."""
    _check_range("padj_max", padj_max)
    _check_range("pvalue_max", pvalue_max)
    if log2fc_min_abs is not None and log2fc_min_abs < 0:
        raise InvalidFilterError("log2fc_min_abs must be >= 0")
    if direction not in (None, "up", "down"):
        raise InvalidFilterError("direction must be 'up' or 'down'")
    if top is not None and top < 0:
        raise InvalidFilterError("top must be >= 0")

    wanted = set(genes) if genes else None
    out: list[ResultRow] = []
    for row in result.rows:
        if gene is not None and row.gene_id != gene:
            continue
        if gene_contains is not None and gene_contains not in row.gene_id:
            continue
        if wanted is not None and row.gene_id not in wanted:
            continue
        if padj_max is not None and row.p_value_adjusted > padj_max:
            continue
        if pvalue_max is not None and row.p_value > pvalue_max:
            continue
        if log2fc_min_abs is not None and abs(row.log2fc) < log2fc_min_abs:
            continue
        if direction == "up" and not row.log2fc > 0:
            continue
        if direction == "down" and not row.log2fc < 0:
            continue
        out.append(row)
    return out[: top if top is not None else len(out)]


def _check_range(name: str, value: float | None) -> None:
    if value is None:
        return
    if not (0.0 <= value <= 1.0):
        raise InvalidFilterError(f"{name} must be within [0, 1]")


def compare(
    results: dict[str, ResultResponse],
    *,
    padj_max: float | None = None,
) -> list[dict]:
    """Join gene rows across tasks. Genes are emitted in a deterministic
    order (task-0 rank order, then gene_id), so concurrent store updates
    cannot reshuffle a single response.
    """
    if not results:
        raise ValidationError("at least one task is required for comparison")

    per_task_genes: dict[str, dict[str, ResultRow]] = {}
    for task_id, res in results.items():
        keep = {
            r.gene_id: r
            for r in res.rows
            if padj_max is None or r.p_value_adjusted <= padj_max
        }
        per_task_genes[task_id] = keep

    first_id = next(iter(results))
    base_order = [r.gene_id for r in results[first_id].rows if r.gene_id in per_task_genes[first_id]]
    extras: set[str] = set()
    for keep in per_task_genes.values():
        extras.update(keep.keys())
    extras.difference_update(base_order)
    ordered_genes = base_order + sorted(extras)

    out: list[dict] = []
    for gene in ordered_genes:
        entry: dict = {"gene_id": gene, "tasks": {}}
        for task_id, keep in per_task_genes.items():
            row = keep.get(gene)
            entry["tasks"][task_id] = (
                None
                if row is None
                else {
                    "rank": row.rank,
                    "log2fc": row.log2fc,
                    "statistic": row.statistic,
                    "p_value": row.p_value,
                    "p_value_adjusted": row.p_value_adjusted,
                }
            )
        out.append(entry)
    return out


def to_csv(rows: Iterable[ResultRow]) -> str:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=_EXPORT_COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: _format(getattr(row, k)) for k in _EXPORT_COLUMNS})
    return buf.getvalue()


def to_json(result: ResultResponse, rows: list[ResultRow]) -> str:
    payload = result.model_dump()
    payload["rows"] = [r.model_dump() for r in rows]
    payload["n_returned"] = len(rows)
    return json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True)


def _format(v: object) -> str:
    if isinstance(v, float):
        if v == float("inf"):
            return "inf"
        if v == float("-inf"):
            return "-inf"
        return repr(v)
    return str(v)
