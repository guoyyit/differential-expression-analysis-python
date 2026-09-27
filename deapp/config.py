"""Service configuration, overridable via environment variables.

Environment variables (prefix ``DEA_``):

``DEA_DATA_DIR``         directory for task/result persistence (default ``./data``)
``DEA_MAX_CONCURRENCY``  at most this many tasks run at once
``DEA_QUEUE_SIZE``       max number of queued (pending) tasks; 0 = unbounded
``DEA_DEFAULT_TIMEOUT``  per-task timeout in seconds; 0 disables the default
"""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class Settings:
    data_dir: str = "data"
    max_concurrency: int = 2
    queue_size: int = 0
    default_timeout: float = 30.0
    # Strategy applied to RUNNING/PENDING tasks found on disk at startup.
    # "recover" -> re-queue; "cancel" -> mark cancelled.
    recovery_policy: str = "recover"

    @staticmethod
    def from_env() -> "Settings":
        def _int(name: str, default: int) -> int:
            raw = os.environ.get(name)
            return int(raw) if raw not in (None, "") else default

        def _float(name: str, default: float) -> float:
            raw = os.environ.get(name)
            return float(raw) if raw not in (None, "") else default

        return Settings(
            data_dir=os.environ.get("DEA_DATA_DIR", "data"),
            max_concurrency=max(1, _int("DEA_MAX_CONCURRENCY", 2)),
            queue_size=max(0, _int("DEA_QUEUE_SIZE", 0)),
            default_timeout=_float("DEA_DEFAULT_TIMEOUT", 30.0),
            recovery_policy=os.environ.get("DEA_RECOVERY_POLICY", "recover"),
        )
