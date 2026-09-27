"""FastAPI wiring: HTTP surface over the task manager."""

from __future__ import annotations

import csv
import io
import math
from contextlib import asynccontextmanager
from typing import Literal

from fastapi import FastAPI, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, Response

from app.config import Settings
from app.errors import ResultNotReadyError, ServiceError
from app.manager import TaskManager
from app.models import (
    AnalysisRequest,
    CompareResponse,
    GeneResult,
    ResultPage,
    TaskSnapshot,
    TaskStatus,
)
from app.store import TaskStore

_SORT_KEYS = {
    "rank": lambda r: (r.rank, r.gene),
    "adjusted_p_value": lambda r: (r.adjusted_p_value, r.p_value, r.gene),
    "p_value": lambda r: (r.p_value, r.adjusted_p_value, r.gene),
    "effect_size": lambda r: (r.effect_size, r.gene),
    "abs_effect_size": lambda r: (-abs(r.effect_size), r.gene),
    "statistic": lambda r: (r.statistic, r.gene),
    "gene": lambda r: (r.gene,),
}


def create_app(settings: Settings | None = None) -> FastAPI:
    """Build the app with its own manager/store (test-friendly factory)."""
    settings = settings or Settings()
    store = TaskStore(settings.db_path)
    manager = TaskManager(settings, store)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        manager.start()
        yield
        manager.shutdown(wait=False)
        store.close()

    app = FastAPI(title="Differential Expression Analysis Service", lifespan=lifespan)
    app.state.manager = manager
    app.state.store = store
    app.state.settings = settings

    @app.exception_handler(ServiceError)
    async def service_error_handler(_: Request, exc: ServiceError) -> JSONResponse:
        return JSONResponse(status_code=exc.http_status, content=exc.to_body())

    @app.exception_handler(RequestValidationError)
    async def request_validation_handler(
        _: Request, exc: RequestValidationError
    ) -> JSONResponse:
        body = {
            "error": {
                "code": "validation_error",
                "message": "request body failed schema validation",
                "details": {"issues": exc.errors()},
            }
        }
        return JSONResponse(status_code=422, content=body)

    # -- tasks ------------------------------------------------------------

    @app.get("/api/v1/health")
    def health() -> dict:
        return {"status": "ok"}

    @app.post("/api/v1/analyses", status_code=202)
    def submit(req: AnalysisRequest) -> TaskSnapshot:
        return manager.submit(req)

    @app.get("/api/v1/analyses")
    def list_analyses(status: TaskStatus | None = None) -> list[TaskSnapshot]:
        return manager.list(status)

    @app.get("/api/v1/analyses/{task_id}")
    def get_analysis(task_id: str) -> TaskSnapshot:
        return manager.get(task_id)

    @app.post("/api/v1/analyses/{task_id}/cancel")
    def cancel_analysis(task_id: str) -> TaskSnapshot:
        return manager.cancel(task_id)

    # -- results ------------------------------------------------------------

    def _completed_results(task_id: str) -> list[GeneResult]:
        snap = manager.get(task_id)
        if snap.status is not TaskStatus.COMPLETED:
            raise ResultNotReadyError(
                "results are only available for completed tasks",
                {"task_id": task_id, "status": snap.status.value},
            )
        return store.load_results(task_id)

    def _filter_sort(
        results: list[GeneResult],
        gene: list[str] | None,
        significant_only: bool,
        max_adjusted_p: float | None,
        min_abs_effect_size: float | None,
        sort: str,
        order: str,
    ) -> list[GeneResult]:
        out = results
        if gene:
            wanted = set(gene)
            out = [r for r in out if r.gene in wanted]
        if significant_only:
            out = [r for r in out if r.significant]
        if max_adjusted_p is not None:
            out = [r for r in out if r.adjusted_p_value <= max_adjusted_p]
        if min_abs_effect_size is not None:
            out = [r for r in out if abs(r.effect_size) >= min_abs_effect_size]
        # total order: every key ends with a unique tiebreak, so pagination
        # is stable even while other tasks are being updated concurrently
        out = sorted(out, key=_SORT_KEYS[sort], reverse=(order == "desc"))
        return out

    @app.get("/api/v1/analyses/{task_id}/results")
    def get_results(
        task_id: str,
        gene: list[str] | None = Query(default=None),
        significant_only: bool = False,
        max_adjusted_p: float | None = Query(default=None, ge=0.0, le=1.0),
        min_abs_effect_size: float | None = Query(default=None, ge=0.0),
        sort: Literal["rank", "adjusted_p_value", "p_value", "effect_size",
                      "abs_effect_size", "statistic", "gene"] = "rank",
        order: Literal["asc", "desc"] = "asc",
        offset: int = Query(default=0, ge=0),
        limit: int = Query(default=100, ge=1, le=10000),
    ) -> ResultPage:
        results = _completed_results(task_id)
        filtered = _filter_sort(results, gene, significant_only, max_adjusted_p,
                                min_abs_effect_size, sort, order)
        page = filtered[offset : offset + limit]
        return ResultPage(
            task_id=task_id, total=len(filtered), offset=offset, limit=limit,
            sort=sort, order=order, results=page,
        )

    @app.get("/api/v1/analyses/{task_id}/results/export")
    def export_results(
        task_id: str,
        format: Literal["csv", "json"] = "csv",
        gene: list[str] | None = Query(default=None),
        significant_only: bool = False,
        max_adjusted_p: float | None = Query(default=None, ge=0.0, le=1.0),
        min_abs_effect_size: float | None = Query(default=None, ge=0.0),
    ) -> Response:
        results = _completed_results(task_id)
        filtered = _filter_sort(results, gene, significant_only, max_adjusted_p,
                                min_abs_effect_size, "rank", "asc")
        if format == "json":
            payload = [r.model_dump() for r in filtered]
            return JSONResponse(
                content={"task_id": task_id, "results": payload},
                headers={"Content-Disposition":
                         f'attachment; filename="{task_id}_results.json"'},
            )
        buf = io.StringIO()
        writer = csv.writer(buf, lineterminator="\n")
        writer.writerow(["rank", "gene", "statistic", "p_value", "adjusted_p_value",
                         "effect_size", "cohens_d", "mean_group_a", "mean_group_b",
                         "significant"])
        for r in filtered:
            writer.writerow([r.rank, r.gene, r.statistic, r.p_value,
                             r.adjusted_p_value, r.effect_size,
                             "" if r.cohens_d is None else r.cohens_d,
                             r.mean_group_a, r.mean_group_b, r.significant])
        return Response(
            content=buf.getvalue(),
            media_type="text/csv",
            headers={"Content-Disposition":
                     f'attachment; filename="{task_id}_results.csv"'},
        )

    # -- cross-task comparison ----------------------------------------------

    @app.post("/api/v1/compare")
    def compare(task_a: str, task_b: str) -> CompareResponse:
        res_a = _completed_results(task_a)
        res_b = _completed_results(task_b)
        by_gene_b = {r.gene: r for r in res_b}
        common = sorted({r.gene for r in res_a} & set(by_gene_b))
        a_map = {r.gene: r for r in res_a}
        both_sig, concordant = 0, 0
        xs, ys = [], []
        for g in common:
            ra, rb = a_map[g], by_gene_b[g]
            xs.append(ra.effect_size)
            ys.append(rb.effect_size)
            if ra.significant and rb.significant:
                both_sig += 1
                if (ra.effect_size > 0) == (rb.effect_size > 0) or (
                    ra.effect_size == 0 and rb.effect_size == 0
                ):
                    concordant += 1
        return CompareResponse(
            task_a=task_a,
            task_b=task_b,
            n_common_genes=len(common),
            n_both_significant=both_sig,
            n_concordant_sign=concordant,
            effect_size_correlation=_pearson(xs, ys),
            common_genes=common,
        )

    return app


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    r = sxy / math.sqrt(sxx * syy)
    return max(-1.0, min(1.0, r))
