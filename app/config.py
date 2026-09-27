"""Server configuration for the differential expression analysis service."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class RecoveryPolicy(str, Enum):
    """What to do with tasks that were unfinished when the service stopped."""

    RESUME = "resume"  # re-queue unfinished tasks and run them again
    FAIL = "fail"  # mark unfinished tasks as interrupted failures


@dataclass(frozen=True)
class Settings:
    """Runtime settings. All limits are enforced by the task manager."""

    db_path: str = "dea_service.db"
    max_workers: int = 2  # bound on concurrently executing analyses
    max_queue_size: int = 16  # bound on pending tasks; overflow is rejected
    default_timeout_seconds: float = 300.0
    max_timeout_seconds: float = 3600.0
    recovery_policy: RecoveryPolicy = RecoveryPolicy.RESUME
    max_genes: int = 200_000
    max_samples: int = 10_000
    slow_task_sleep_per_gene: float = 0.0  # test hook: artificial per-gene delay
