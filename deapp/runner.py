"""Concurrent task orchestration.

* a bounded ``ThreadPoolExecutor`` enforces the max-concurrency resource
  limit, so one slow task never blocks more than its own worker;
* every task owns a :class:`threading.Event` used for cooperative cancel and
  timeout (a watchdog timer trips the same event);
* the state transitions pending -> running -> succeeded/failed/cancelled are
  persisted atomically; cancellation before completion leaves *no* result
  file and a task is marked succeeded only after the result is durably
  written, so finished results can never be overwritten with a failure.
"""
from __future__ import annotations

import threading
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

from .config import Settings
from .engine import TaskCancelled, run_analysis
from .errors import TaskStateError
from .models import CreateTaskRequest, TaskDetail, TaskStatus
from .persistence import Store
from .validation import ValidatedInput, validate_request


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


class TaskRunner:
    def __init__(self, settings: Settings, store: Store) -> None:
        self.settings = settings
        self.store = store
        self._pool = ThreadPoolExecutor(
            max_workers=settings.max_concurrency,
            thread_name_prefix="dea-task",
        )
        self._cancel_events: dict[str, threading.Event] = {}
        self._lock = threading.Lock()

    # ------------------------------------------------------------------
    def start_recovered(self) -> None:
        """Re-submit tasks that were pending/running before a restart."""
        for task in self.store.list_tasks():
            if task.status == TaskStatus.PENDING:
                self._submit(task.id)

    def submit(
        self,
        task_id: str,
        request: CreateTaskRequest,
        validated: ValidatedInput,
    ) -> None:
        with self._lock:
            event = threading.Event()
            self._cancel_events[task_id] = event
        self.store.save_request(task_id, request.model_dump(mode="json"))
        self._submit(task_id)

    def _submit(self, task_id: str) -> None:
        with self._lock:
            event = self._cancel_events.setdefault(task_id, threading.Event())
        self._pool.submit(self._run, task_id, event)

    # ------------------------------------------------------------------
    def cancel(self, task_id: str) -> TaskDetail:
        task = self.store.get_task(task_id)
        if task.status in (TaskStatus.SUCCEEDED, TaskStatus.FAILED, TaskStatus.CANCELLED):
            raise TaskStateError(
                f"task is already {task.status.value}",
                details={"status": task.status.value},
            )
        with self._lock:
            event = self._cancel_events.get(task_id)
        if event is not None:
            event.set()
        # A still-queued task may not reach the worker promptly; mark it
        # cancelled here and the worker treats the event as authoritative.
        if task.status == TaskStatus.PENDING:
            task.status = TaskStatus.CANCELLED
            task.finished_at = _utcnow()
            task.error = {
                "code": "TASK_CANCELLED",
                "message": "task cancelled before execution",
                "details": {"timeout": False},
            }
            self.store.save_task(task)
        return self.store.get_task(task_id)

    # ------------------------------------------------------------------
    def shutdown(self) -> None:
        self._pool.shutdown(wait=False, cancel_futures=True)

    # ------------------------------------------------------------------
    def _run(self, task_id: str, event: threading.Event) -> None:
        try:
            task = self.store.get_task(task_id)
            if task.status == TaskStatus.CANCELLED or event.is_set():
                return
            if task.status == TaskStatus.SUCCEEDED and self.store.has_result(task_id):
                return  # idempotent: finished result is never overwritten

            task.status = TaskStatus.RUNNING
            task.started_at = _utcnow()
            task.error = None
            self.store.save_task(task)

            request = CreateTaskRequest.model_validate(self.store.load_request(task_id))
            validated = validate_request(request)

            timeout = (
                request.config.timeout_seconds
                if request.config.timeout_seconds is not None
                else (self.settings.default_timeout or None)
            )
            timed_out = False
            timer: threading.Timer | None = None
            if timeout is not None:
                def _on_timeout() -> None:
                    nonlocal timed_out
                    timed_out = True
                    event.set()

                timer = threading.Timer(timeout, _on_timeout)
                timer.daemon = True
                timer.start()
            try:
                result = run_analysis(
                    task_id, validated, request.config, is_cancelled=event.is_set
                )
                if event.is_set():
                    raise TaskCancelled()
            finally:
                if timer is not None:
                    timer.cancel()

            # ---- commit point: result first, status after --------------
            self.store.save_result(result)
            fresh = self.store.get_task(task_id)
            fresh.status = TaskStatus.SUCCEEDED
            fresh.finished_at = _utcnow()
            fresh.started_at = task.started_at
            fresh.strategy = result.strategy
            fresh.n_genes = result.n_genes
            fresh.n_samples = len(validated.sample_ids)
            fresh.error = None
            self.store.save_task(fresh)
        except TaskCancelled:
            self._mark_terminal(task_id, TaskStatus.CANCELLED, "TASK_CANCELLED",
                                "task cancelled or timed out", timeout=timed_out)
        except Exception as exc:  # analysis failure (should not happen post-validation)
            self._mark_terminal(
                task_id,
                TaskStatus.FAILED,
                "ANALYSIS_FAILED",
                str(exc),
                traceback=traceback.format_exc(),
            )
        finally:
            with self._lock:
                self._cancel_events.pop(task_id, None)

    def _mark_terminal(
        self,
        task_id: str,
        status: TaskStatus,
        code: str,
        message: str,
        *,
        traceback: str | None = None,
        timeout: bool = False,
    ) -> None:
        task = self.store.get_task(task_id)
        # Never downgrade a succeeded task.
        if task.status == TaskStatus.SUCCEEDED and self.store.has_result(task_id):
            return
        task.status = status
        task.finished_at = _utcnow()
        details: dict = {"timeout": timeout} if timeout else {}
        if traceback is not None:
            details["traceback"] = traceback
        task.error = {"code": code, "message": message, "details": details}
        self.store.save_task(task)
