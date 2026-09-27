"""Task manager: bounded worker pool, cancellation, timeout, recovery."""

from __future__ import annotations

import threading
import uuid
from concurrent.futures import ThreadPoolExecutor

from app.config import RecoveryPolicy, Settings
from app.errors import (
    QueueFullError,
    ServiceError,
    TaskConflictError,
    TaskNotFoundError,
)
from app.models import AnalysisRequest, AnalysisSummary, TaskSnapshot, TaskStatus
from app.stats.pipeline import Cancelled, run_analysis
from app.store import TaskStore
from app.validation import group_column_indices, validate_request


class TaskManager:
    """Owns the worker pool and the store. All state transitions go here."""

    def __init__(self, settings: Settings, store: TaskStore):
        self._settings = settings
        self._store = store
        self._pool: ThreadPoolExecutor | None = None
        self._submit_lock = threading.Lock()
        self._inflight = 0  # queued + running; bounds total accepted work
        self._tokens: dict[str, threading.Event] = {}
        self._tokens_lock = threading.Lock()
        self._timers: dict[str, threading.Timer] = {}
        self._progress: dict[str, tuple[int, int]] = {}

    # -- lifecycle ----------------------------------------------------------

    def start(self) -> None:
        """Open the pool and recover unfinished tasks per recovery policy."""
        self._pool = ThreadPoolExecutor(
            max_workers=self._settings.max_workers, thread_name_prefix="dea"
        )
        for task_id in self._store.unfinished_task_ids():
            if self._settings.recovery_policy is RecoveryPolicy.RESUME:
                self._requeue(task_id)
            else:
                self._store.mark_terminal(
                    task_id,
                    TaskStatus.INTERRUPTED,
                    error={
                        "code": "interrupted",
                        "message": "service stopped before the task finished",
                        "details": {},
                    },
                )

    def _requeue(self, task_id: str) -> None:
        with self._tokens_lock:
            self._tokens[task_id] = threading.Event()
        with self._submit_lock:
            self._inflight += 1
        assert self._pool is not None
        self._pool.submit(self._run_recovered, task_id)

    def _run_recovered(self, task_id: str) -> None:
        """Rows from a previous life are reset to PENDING, then re-run."""
        try:
            self._store.reset_for_recovery(task_id)
            self._run(task_id)
        finally:
            with self._submit_lock:
                self._inflight -= 1

    def shutdown(self, wait: bool = False) -> None:
        with self._tokens_lock:
            timers = list(self._timers.values())
        for t in timers:
            t.cancel()
        if self._pool is not None:
            self._pool.shutdown(wait=wait, cancel_futures=True)

    # -- submission ---------------------------------------------------------

    def submit(self, req: AnalysisRequest) -> TaskSnapshot:
        """Validate, persist, enqueue. Idempotent on client-supplied task_id."""
        validate_request(req, self._settings)
        a_idx, b_idx = group_column_indices(req)

        if req.task_id is not None:
            snap = self._store.get_snapshot(req.task_id)
            if snap is not None:
                stored_req, _, _ = self._store.load_request(req.task_id)
                if stored_req.model_dump() == req.model_dump():
                    return snap  # identical re-submission: no duplicate work
                raise TaskConflictError(
                    "task_id already used with a different request",
                    {"task_id": req.task_id},
                )
            task_id = req.task_id
        else:
            task_id = f"dea-{uuid.uuid4().hex[:16]}"

        with self._submit_lock:
            if self._inflight >= self._settings.max_queue_size + self._settings.max_workers:
                raise QueueFullError(
                    "too many pending tasks",
                    {
                        "max_queue_size": self._settings.max_queue_size,
                        "max_workers": self._settings.max_workers,
                    },
                )
            created = self._store.create_task(task_id, req, a_idx, b_idx)
            if not created:  # lost a race on a client-supplied id
                snap = self._store.get_snapshot(task_id)
                assert snap is not None
                return snap
            self._inflight += 1
            with self._tokens_lock:
                self._tokens[task_id] = threading.Event()
            assert self._pool is not None
            self._pool.submit(self._run, task_id)
        snap = self._store.get_snapshot(task_id)
        assert snap is not None
        return snap

    # -- execution ------------------------------------------------------------

    def _timeout_seconds(self, req: AnalysisRequest) -> float:
        t = req.config.timeout_seconds
        if t is None:
            t = self._settings.default_timeout_seconds
        return min(t, self._settings.max_timeout_seconds)

    def _run(self, task_id: str) -> None:
        timer: threading.Timer | None = None
        try:
            if not self._store.mark_running(task_id):
                return  # cancelled (or otherwise transitioned) while queued
            req, a_idx, b_idx = self._store.load_request(task_id)
            with self._tokens_lock:
                token = self._tokens.setdefault(task_id, threading.Event())

            timeout = self._timeout_seconds(req)
            timer = threading.Timer(timeout, self._on_timeout, args=(task_id,))
            timer.daemon = True
            with self._tokens_lock:
                self._timers[task_id] = timer
            timer.start()

            results = run_analysis(
                req.genes,
                req.matrix,
                a_idx,
                b_idx,
                req.config,
                token,
                progress=lambda done, total: self._progress.__setitem__(
                    task_id, (done, total)
                ),
                sleep_per_gene=self._settings.slow_task_sleep_per_gene,
            )
            summary = AnalysisSummary(
                n_genes_input=len(req.genes),
                n_genes_tested=len(results),
                n_significant=sum(1 for r in results if r.significant),
                group_a=req.groups[0].name,
                group_b=req.groups[1].name,
            )
            # False => task was cancelled/timed out mid-run: discard results.
            self._store.save_results(task_id, req.config, results, summary)
        except Cancelled:
            self._store.mark_terminal(task_id, TaskStatus.CANCELLED)
        except ServiceError as e:
            self._store.mark_terminal(task_id, TaskStatus.FAILED, error=e.to_body()["error"])
        except Exception as e:  # never leak a worker crash as a stuck task
            self._store.mark_terminal(
                task_id,
                TaskStatus.FAILED,
                error={"code": "internal_error", "message": str(e), "details": {}},
            )
        finally:
            if timer is not None:
                timer.cancel()
            with self._tokens_lock:
                self._timers.pop(task_id, None)
            self._progress.pop(task_id, None)
            with self._submit_lock:
                self._inflight -= 1

    def _on_timeout(self, task_id: str) -> None:
        with self._tokens_lock:
            token = self._tokens.get(task_id)
        if token is not None:
            token.set()  # ask the pipeline to abort cooperatively
        self._store.mark_terminal(
            task_id,
            TaskStatus.TIMED_OUT,
            error={
                "code": "timed_out",
                "message": "analysis exceeded its time budget",
                "details": {},
            },
            from_statuses=(TaskStatus.RUNNING,),
        )

    # -- queries --------------------------------------------------------------

    def get(self, task_id: str) -> TaskSnapshot:
        snap = self._store.get_snapshot(task_id)
        if snap is None:
            raise TaskNotFoundError("unknown task", {"task_id": task_id})
        return snap

    def list(self, status: TaskStatus | None = None) -> list[TaskSnapshot]:
        return self._store.list_snapshots(status)

    def get_progress(self, task_id: str) -> tuple[int, int] | None:
        return self._progress.get(task_id)

    def cancel(self, task_id: str) -> TaskSnapshot:
        """Cooperative cancel: pending tasks never run, running tasks stop."""
        snap = self.get(task_id)
        if snap.status.terminal:
            return snap  # idempotent: cancelling a finished task is a no-op
        with self._tokens_lock:
            token = self._tokens.get(task_id)
        if token is not None:
            token.set()
        self._store.mark_terminal(task_id, TaskStatus.CANCELLED)
        return self.get(task_id)
