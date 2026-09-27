"""SQLite persistence: task records and results survive restarts.

All state transitions are conditional UPDATEs guarded by the expected
current status, so a late worker can never overwrite a terminal state
(cancelled/timed-out) and a completed task can never be re-marked failed.
"""

from __future__ import annotations

import json
import sqlite3
import threading
from datetime import datetime, timezone

from app.models import (
    AnalysisConfig,
    AnalysisRequest,
    AnalysisSummary,
    GeneResult,
    TaskSnapshot,
    TaskStatus,
)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS tasks (
    task_id      TEXT PRIMARY KEY,
    status       TEXT NOT NULL,
    created_at   TEXT NOT NULL,
    started_at   TEXT,
    finished_at  TEXT,
    config_json  TEXT NOT NULL,
    request_json TEXT NOT NULL,
    group_a_idx  TEXT NOT NULL,
    group_b_idx  TEXT NOT NULL,
    error_json   TEXT,
    summary_json TEXT,
    results_json TEXT
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class TaskStore:
    """Thread-safe SQLite-backed store. One row per task; results as JSON."""

    def __init__(self, db_path: str):
        self._db_path = db_path
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(db_path, check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL")
        with self._lock, self._conn:
            self._conn.executescript(_SCHEMA)

    def close(self) -> None:
        with self._lock:
            self._conn.close()

    # -- creation ---------------------------------------------------------

    def create_task(
        self, task_id: str, req: AnalysisRequest, group_a_idx: list[int], group_b_idx: list[int]
    ) -> bool:
        """Insert a PENDING task. Return False if task_id already exists."""
        with self._lock, self._conn:
            cur = self._conn.execute(
                "INSERT OR IGNORE INTO tasks "
                "(task_id, status, created_at, config_json, request_json, "
                " group_a_idx, group_b_idx) VALUES (?,?,?,?,?,?,?)",
                (
                    task_id,
                    TaskStatus.PENDING.value,
                    _now(),
                    req.config.model_dump_json(),
                    req.model_dump_json(),
                    json.dumps(group_a_idx),
                    json.dumps(group_b_idx),
                ),
            )
            return cur.rowcount == 1

    def exists(self, task_id: str) -> bool:
        with self._lock:
            row = self._conn.execute(
                "SELECT 1 FROM tasks WHERE task_id = ?", (task_id,)
            ).fetchone()
            return row is not None

    def load_request(self, task_id: str) -> tuple[AnalysisRequest, list[int], list[int]]:
        """Return the stored request and resolved group column indices."""
        with self._lock:
            row = self._conn.execute(
                "SELECT request_json, group_a_idx, group_b_idx FROM tasks WHERE task_id = ?",
                (task_id,),
            ).fetchone()
        if row is None:
            raise KeyError(task_id)
        return (
            AnalysisRequest.model_validate_json(row["request_json"]),
            json.loads(row["group_a_idx"]),
            json.loads(row["group_b_idx"]),
        )

    # -- guarded transitions ----------------------------------------------

    def mark_running(self, task_id: str) -> bool:
        """PENDING -> RUNNING. False if the task left PENDING meanwhile."""
        with self._lock, self._conn:
            cur = self._conn.execute(
                "UPDATE tasks SET status = ?, started_at = ? "
                "WHERE task_id = ? AND status = ?",
                (TaskStatus.RUNNING.value, _now(), task_id, TaskStatus.PENDING.value),
            )
            return cur.rowcount == 1

    def save_results(
        self,
        task_id: str,
        config: AnalysisConfig,
        results: list[GeneResult],
        summary: AnalysisSummary,
    ) -> bool:
        """Atomically store results and mark COMPLETED — only from RUNNING.

        Returns False (results discarded) if the task was cancelled or
        timed out while the worker was computing.
        """
        payload = json.dumps([r.model_dump() for r in results])
        with self._lock, self._conn:
            cur = self._conn.execute(
                "UPDATE tasks SET status = ?, finished_at = ?, results_json = ?, "
                "summary_json = ? WHERE task_id = ? AND status = ?",
                (
                    TaskStatus.COMPLETED.value,
                    _now(),
                    payload,
                    summary.model_dump_json(),
                    task_id,
                    TaskStatus.RUNNING.value,
                ),
            )
            return cur.rowcount == 1

    def mark_terminal(
        self,
        task_id: str,
        status: TaskStatus,
        *,
        error: dict | None = None,
        from_statuses: tuple[TaskStatus, ...] = (TaskStatus.PENDING, TaskStatus.RUNNING),
    ) -> bool:
        """Move to a terminal state from an allowed set. Never touches
        COMPLETED tasks, so finished results can't be mis-marked failed."""
        if not status.terminal:
            raise ValueError(f"not a terminal status: {status}")
        placeholders = ",".join("?" for _ in from_statuses)
        with self._lock, self._conn:
            cur = self._conn.execute(
                f"UPDATE tasks SET status = ?, finished_at = ?, error_json = ? "
                f"WHERE task_id = ? AND status IN ({placeholders})",
                (
                    status.value,
                    _now(),
                    json.dumps(error) if error is not None else None,
                    task_id,
                    *(s.value for s in from_statuses),
                ),
            )
            return cur.rowcount == 1

    # -- reads --------------------------------------------------------------

    def load_results(self, task_id: str) -> list[GeneResult]:
        with self._lock:
            row = self._conn.execute(
                "SELECT results_json FROM tasks WHERE task_id = ?", (task_id,)
            ).fetchone()
        if row is None:
            raise KeyError(task_id)
        if row["results_json"] is None:
            return []
        return [GeneResult.model_validate(r) for r in json.loads(row["results_json"])]

    def get_snapshot(self, task_id: str) -> TaskSnapshot | None:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE task_id = ?", (task_id,)
            ).fetchone()
        return self._to_snapshot(row) if row else None

    def list_snapshots(self, status: TaskStatus | None = None) -> list[TaskSnapshot]:
        with self._lock:
            if status is None:
                rows = self._conn.execute(
                    "SELECT * FROM tasks ORDER BY created_at, task_id"
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT * FROM tasks WHERE status = ? ORDER BY created_at, task_id",
                    (status.value,),
                ).fetchall()
        return [self._to_snapshot(r) for r in rows]

    def unfinished_task_ids(self) -> list[str]:
        """Tasks in PENDING/RUNNING — used by crash recovery at startup."""
        with self._lock:
            rows = self._conn.execute(
                "SELECT task_id FROM tasks WHERE status IN (?, ?) ORDER BY created_at, task_id",
                (TaskStatus.PENDING.value, TaskStatus.RUNNING.value),
            ).fetchall()
        return [r["task_id"] for r in rows]

    @staticmethod
    def _to_snapshot(row: sqlite3.Row) -> TaskSnapshot:
        summary = row["summary_json"]
        return TaskSnapshot(
            task_id=row["task_id"],
            status=TaskStatus(row["status"]),
            created_at=row["created_at"],
            started_at=row["started_at"],
            finished_at=row["finished_at"],
            config=AnalysisConfig.model_validate_json(row["config_json"]),
            error=json.loads(row["error_json"]) if row["error_json"] else None,
            summary=AnalysisSummary.model_validate_json(summary) if summary else None,
        )
