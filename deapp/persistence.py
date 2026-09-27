"""File-backed task/result persistence with crash-safe atomic writes.

On-disk layout::

    data_dir/
      tasks/<task_id>.json          # task record (atomic temp+rename)
      results/<task_id>.json        # full result (atomic temp+rename)
      requests/<task_id>.json       # original validated request snapshot

A task is only ever marked ``succeeded`` *after* its result file has been
fully written and fsynced, so:

* a crash while running leaves the task in running/pending and no result;
* a succeeded task always has a readable result;
* cancellation/timeout never produces a partial result file.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from .config import Settings
from .errors import ResultNotFoundError, TaskNotFoundError
from .models import ResultResponse, TaskDetail, TaskStatus


def _atomic_write(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-", suffix=".json")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _read_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


class Store:
    """Persistent store; keeps task records cached in memory."""

    def __init__(self, settings: Settings) -> None:
        self.settings = settings
        self.root = Path(settings.data_dir)
        self.tasks_dir = self.root / "tasks"
        self.results_dir = self.root / "results"
        self.requests_dir = self.root / "requests"
        for d in (self.tasks_dir, self.results_dir, self.requests_dir):
            d.mkdir(parents=True, exist_ok=True)
        self._tasks: dict[str, TaskDetail] = {}

    # ----- load -----------------------------------------------------------
    def load(self) -> None:
        """Read every task record and reconcile unfinished ones."""
        for path in sorted(self.tasks_dir.glob("*.json")):
            try:
                data = _read_json(path)
            except (json.JSONDecodeError, OSError):
                # A torn task record: nothing was ever finished, drop it.
                try:
                    path.unlink()
                except OSError:
                    pass
                continue
            task = TaskDetail.model_validate(data)
            self._tasks[task.id] = task

        for task in list(self._tasks.values()):
            result_path = self._result_path(task.id)
            if task.status == TaskStatus.SUCCEEDED:
                # A succeeded record without a result file means the
                # result-write crashed before the rename; it cannot be trusted.
                if not result_path.exists():
                    self._reconcile_unfinished(task)
            elif task.status in (TaskStatus.PENDING, TaskStatus.RUNNING):
                self._reconcile_unfinished(task)

    def _reconcile_unfinished(self, task: TaskDetail) -> None:
        if self.settings.recovery_policy == "cancel":
            task.status = TaskStatus.CANCELLED
            task.finished_at = task.finished_at
            task.error = {
                "code": "INTERRUPTED",
                "message": "task interrupted by service restart and marked cancelled",
                "details": {"policy": "cancel"},
            }
        else:
            # Re-queue from scratch; stale temp result files are ignored.
            task.status = TaskStatus.PENDING
            task.started_at = None
            task.error = None
        self.save_task(task)

    # ----- paths ----------------------------------------------------------
    def _task_path(self, task_id: str) -> Path:
        return self.tasks_dir / f"{task_id}.json"

    def _result_path(self, task_id: str) -> Path:
        return self.results_dir / f"{task_id}.json"

    def _request_path(self, task_id: str) -> Path:
        return self.requests_dir / f"{task_id}.json"

    # ----- task CRUD ------------------------------------------------------
    def save_task(self, task: TaskDetail) -> None:
        self._tasks[task.id] = task
        _atomic_write(self._task_path(task.id), task.model_dump_json(indent=2))

    def get_task(self, task_id: str) -> TaskDetail:
        task = self._tasks.get(task_id)
        if task is None:
            raise TaskNotFoundError(f"task {task_id!r} does not exist")
        return task

    def list_tasks(self) -> list[TaskDetail]:
        return sorted(self._tasks.values(), key=lambda t: (t.created_at, t.id))

    def exists(self, task_id: str) -> bool:
        return task_id in self._tasks

    def save_request(self, task_id: str, payload: dict) -> None:
        _atomic_write(
            self._request_path(task_id),
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        )

    def load_request(self, task_id: str) -> dict:
        path = self._request_path(task_id)
        if not path.exists():
            raise TaskNotFoundError(f"stored request for {task_id!r} is missing")
        return _read_json(path)

    # ----- results --------------------------------------------------------
    def save_result(self, result: ResultResponse) -> None:
        # Write the result FIRST, then flip the task to succeeded.
        _atomic_write(
            self._result_path(result.task_id),
            result.model_dump_json(indent=2),
        )

    def load_result(self, task_id: str) -> ResultResponse:
        path = self._result_path(task_id)
        if not path.exists():
            raise ResultNotFoundError(
                f"no completed result for task {task_id!r}",
                details={"status": self.get_task(task_id).status.value
                         if self.exists(task_id) else None},
            )
        return ResultResponse.model_validate(_read_json(path))

    def has_result(self, task_id: str) -> bool:
        return self._result_path(task_id).exists()
