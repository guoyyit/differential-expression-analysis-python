"""Task manager: bounded worker pool, cancellation, timeout, recovery."""

from __future__ import annotations

import threading

from app.config import Settings
from app.models import AnalysisRequest, TaskSnapshot, TaskStatus
from app.store import TaskStore


class TaskManager:
    """Owns the worker pool and the store. All state transitions go here."""

    def __init__(self, settings: Settings, store: TaskStore):
        self._settings = settings
        self._store = store
        self._lock = threading.Lock()
        raise NotImplementedError

    def start(self) -> None:
        """Open the pool and recover unfinished tasks per recovery policy."""
        raise NotImplementedError

    def shutdown(self, wait: bool = False) -> None:
        raise NotImplementedError

    def submit(self, req: AnalysisRequest) -> TaskSnapshot:
        """Validate, persist, enqueue. Idempotent on client-supplied task_id."""
        raise NotImplementedError

    def get(self, task_id: str) -> TaskSnapshot:
        raise NotImplementedError

    def list(self, status: TaskStatus | None = None) -> list[TaskSnapshot]:
        raise NotImplementedError

    def cancel(self, task_id: str) -> TaskSnapshot:
        """Cooperative cancel: pending tasks never run, running tasks stop."""
        raise NotImplementedError
