"""SQLite persistence: task records and results survive restarts."""

from __future__ import annotations

import threading

from app.models import AnalysisConfig, AnalysisRequest, GeneResult, TaskSnapshot, TaskStatus


class TaskStore:
    """Thread-safe SQLite-backed store. One row per task; results as JSON."""

    def __init__(self, db_path: str):
        self._db_path = db_path
        self._lock = threading.Lock()
        raise NotImplementedError

    def close(self) -> None:
        raise NotImplementedError

    def create_task(
        self, task_id: str, req: AnalysisRequest, group_a_idx: list[int], group_b_idx: list[int]
    ) -> bool:
        """Insert a PENDING task. Return False if task_id already exists."""
        raise NotImplementedError

    def get_task_id_if_exists(self, task_id: str) -> bool:
        raise NotImplementedError

    def load_request(self, task_id: str) -> tuple[AnalysisRequest, list[int], list[int]]:
        """Return the stored request and resolved group column indices."""
        raise NotImplementedError

    def set_status(
        self,
        task_id: str,
        status: TaskStatus,
        *,
        started: bool = False,
        finished: bool = False,
        error: dict | None = None,
    ) -> None:
        raise NotImplementedError

    def save_results(
        self, task_id: str, config: AnalysisConfig, results: list[GeneResult]
    ) -> None:
        """Atomically store results and mark the task COMPLETED."""
        raise NotImplementedError

    def load_results(self, task_id: str) -> list[GeneResult]:
        raise NotImplementedError

    def get_snapshot(self, task_id: str) -> TaskSnapshot | None:
        raise NotImplementedError

    def list_snapshots(self, status: TaskStatus | None = None) -> list[TaskSnapshot]:
        raise NotImplementedError

    def unfinished_task_ids(self) -> list[str]:
        """Tasks in PENDING/RUNNING — used by crash recovery at startup."""
        raise NotImplementedError
