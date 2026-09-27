"""FastAPI application: HTTP layer for the analysis service."""
from __future__ import annotations

import csv
import io
import json
import math
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from typing import Any

from fastapi import FastAPI, File, Form, Query, Request, UploadFile
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, PlainTextResponse, Response

from .config import Settings
from .errors import ServiceError
from .models import (
    CreateTaskRequest,
    ExpressionMatrix,
    TaskStatus,
)
from .persistence import Store
from .results import compare, filter_rows, to_csv, to_json
from .runner import TaskRunner
from .validation import validate_request


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


class Service:
    def __init__(self, settings: Settings | None = None) -> None:
        self.settings = settings or Settings.from_env()
        self.store = Store(self.settings)
        self.store.load()
        self.runner = TaskRunner(self.settings, self.store)

    def create_task(self, request: CreateTaskRequest) -> dict:
        # Validate BEFORE any persistence; rejects never create a task.
        validated = validate_request(request)

        pending = sum(
            1
            for t in self.store.list_tasks()
            if t.status in (TaskStatus.PENDING, TaskStatus.RUNNING)
        )
        # queue_size bounds the total number of unfinished tasks
        # (running + pending). With queue_size = N and concurrency C, at
        # most N-C tasks wait in the queue; 0 disables the limit.
        if self.settings.queue_size and pending >= self.settings.queue_size:
            raise ServiceError(
                "task queue is full, retry later",
                code="QUEUE_FULL",
                http_status=429,
                details={"limit": self.settings.queue_size, "active": pending},
            )

        task_id = uuid.uuid4().hex
        from .models import TaskDetail

        task = TaskDetail(
            id=task_id,
            name=request.name,
            status=TaskStatus.PENDING,
            created_at=_utcnow(),
            control=validated.control_group,
            treatment=validated.treatment_group,
            design=request.design,
            n_genes=len(validated.gene_ids),
            n_samples=len(validated.sample_ids),
        )
        self.store.save_task(task)
        self.runner.submit(task_id, request, validated)
        return task.model_dump()


def create_app(settings: Settings | None = None) -> FastAPI:
    service = Service(settings)

    @asynccontextmanager
    async def lifespan(app: FastAPI) -> Any:
        service.runner.start_recovered()
        try:
            yield
        finally:
            service.runner.shutdown()

    app = FastAPI(
        title="Differential Expression Analysis Service",
        version="0.1.0",
        lifespan=lifespan,
    )
    app.state.service = service

    # --------------------------------------------------------------- errors
    @app.exception_handler(ServiceError)
    async def _service_error(request: Request, exc: ServiceError) -> JSONResponse:
        return JSONResponse(status_code=exc.http_status, content=exc.to_dict())

    @app.exception_handler(RequestValidationError)
    async def _request_validation(
        request: Request, exc: RequestValidationError
    ) -> JSONResponse:
        return JSONResponse(
            status_code=422,
            content={
                "error": {
                    "code": "INVALID_REQUEST",
                    "message": "request payload failed schema validation",
                    "details": {"errors": exc.errors()},
                }
            },
        )

    @app.get("/health")
    async def health() -> dict:
        return {"status": "ok"}

    # --------------------------------------------------------------- tasks
    @app.post("/tasks", status_code=201)
    async def create_task(request: CreateTaskRequest) -> dict:
        return service.create_task(request)

    @app.post("/tasks/from-file", status_code=201)
    async def create_task_from_file(
        matrix: UploadFile = File(...),
        metadata: str = Form(...),
    ) -> dict:
        try:
            meta = json.loads(metadata)
        except json.JSONDecodeError as exc:
            raise ServiceError(
                "metadata must be valid JSON",
                code="INVALID_METADATA",
                details={"position": exc.pos},
            )
        text = (await matrix.read()).decode("utf-8")
        expr = _parse_matrix_csv(text)
        payload = {**meta, "matrix": expr}
        request = CreateTaskRequest.model_validate(payload)
        return service.create_task(request)

    @app.get("/tasks")
    async def list_tasks(status: TaskStatus | None = None) -> dict:
        tasks = service.store.list_tasks()
        if status is not None:
            tasks = [t for t in tasks if t.status == status]
        return {"tasks": [t.model_dump() for t in tasks]}

    @app.get("/tasks/{task_id}")
    async def get_task(task_id: str) -> dict:
        return service.store.get_task(task_id).model_dump()

    @app.post("/tasks/{task_id}/cancel")
    async def cancel_task(task_id: str) -> dict:
        return service.runner.cancel(task_id).model_dump()

    # ------------------------------------------------------------- results
    @app.get("/tasks/{task_id}/result")
    async def get_result(
        task_id: str,
        gene: str | None = None,
        gene_contains: str | None = None,
        padj_max: float | None = None,
        pvalue_max: float | None = None,
        log2fc_min_abs: float | None = None,
        direction: str | None = None,
        top: int | None = Query(default=None),
    ) -> dict:
        result = service.store.load_result(task_id)
        rows = filter_rows(
            result,
            gene=gene,
            gene_contains=gene_contains,
            padj_max=padj_max,
            pvalue_max=pvalue_max,
            log2fc_min_abs=log2fc_min_abs,
            direction=direction,
            top=top,
        )
        payload = result.model_dump()
        payload["rows"] = [r.model_dump() for r in rows]
        payload["n_returned"] = len(rows)
        return payload

    @app.get("/tasks/{task_id}/export")
    async def export_result(
        task_id: str,
        format: str = "csv",
        gene_contains: str | None = None,
        padj_max: float | None = None,
        top: int | None = None,
    ) -> Response:
        result = service.store.load_result(task_id)
        rows = filter_rows(
            result, gene_contains=gene_contains, padj_max=padj_max, top=top
        )
        if format == "csv":
            return PlainTextResponse(
                to_csv(rows),
                media_type="text/csv",
                headers={
                    "Content-Disposition": f'attachment; filename="{task_id}.csv"'
                },
            )
        if format == "json":
            return PlainTextResponse(
                to_json(result, rows),
                media_type="application/json",
                headers={
                    "Content-Disposition": f'attachment; filename="{task_id}.json"'
                },
            )
        raise ServiceError(
            "format must be 'csv' or 'json'", code="INVALID_FORMAT"
        )

    @app.get("/results/search")
    async def search_results(
        gene: str | None = None,
        gene_contains: str | None = None,
        padj_max: float | None = None,
        direction: str | None = None,
        status: TaskStatus = TaskStatus.SUCCEEDED,
    ) -> dict:
        """Cross-task gene/condition search over every persisted result.

        Each completed task is snapshotted once per request, so concurrent
        task updates cannot interleave within a single response.
        """
        hits: list[dict] = []
        for task in service.store.list_tasks():
            if task.status != status:
                continue
            result = service.store.load_result(task.id)
            rows = filter_rows(
                result,
                gene=gene,
                gene_contains=gene_contains,
                padj_max=padj_max,
                direction=direction,
            )
            for row in rows:
                hits.append(
                    {
                        "task_id": task.id,
                        "task_name": task.name,
                        "control": task.control,
                        "treatment": task.treatment,
                        "row": row.model_dump(),
                    }
                )
        # deterministic order: task id then rank
        hits.sort(key=lambda h: (h["task_id"], h["row"]["rank"]))
        return {"n_hits": len(hits), "hits": hits}

    @app.post("/compare")
    async def compare_tasks(body: dict) -> list:
        task_ids = body.get("task_ids")
        if not isinstance(task_ids, list) or not task_ids:
            raise ServiceError(
                "task_ids must be a non-empty list", code="INVALID_COMPARE_REQUEST"
            )
        results = {tid: service.store.load_result(tid) for tid in task_ids}
        return compare(results, padj_max=body.get("padj_max"))

    return app


def _parse_matrix_csv(text: str) -> dict:
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)
    if len(rows) < 2:
        raise ServiceError(
            "matrix CSV needs a header row and at least one gene row",
            code="INVALID_MATRIX_FILE",
        )
    header, *body_rows = rows
    sample_ids = [h.strip() for h in header[1:]]
    gene_ids: list[str] = []
    values: list[list[float]] = []
    for line_no, row in enumerate(body_rows, start=2):
        if not row or all(c.strip() == "" for c in row):
            continue
        if len(row) != len(header):
            raise ServiceError(
                f"CSV line {line_no} has {len(row)} columns, expected {len(header)}",
                code="INVALID_MATRIX_FILE",
                details={"line": line_no},
            )
        gene_ids.append(row[0].strip())
        try:
            numeric = [float(v) for v in row[1:]]
        except ValueError:
            raise ServiceError(
                f"non-numeric value on CSV line {line_no}",
                code="INVALID_MATRIX_FILE",
                details={"line": line_no},
            )
        if any(not math.isfinite(v) for v in numeric):
            raise ServiceError(
                f"non-finite (NaN/inf) value on CSV line {line_no}",
                code="INVALID_MATRIX_FILE",
                details={"line": line_no},
            )
        values.append(numeric)
    if not gene_ids:
        raise ServiceError(
            "matrix CSV contains no gene rows", code="INVALID_MATRIX_FILE"
        )
    return ExpressionMatrix(
        gene_ids=gene_ids, sample_ids=sample_ids, values=values
    ).model_dump()
