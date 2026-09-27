"""Service entrypoint: ``uv run uvicorn main:app`` (or ``fastapi run``)."""

from __future__ import annotations

import os

from app.config import RecoveryPolicy, Settings
from app.service import create_app

settings = Settings(
    db_path=os.environ.get("DEA_DB_PATH", "dea_service.db"),
    max_workers=int(os.environ.get("DEA_MAX_WORKERS", "2")),
    max_queue_size=int(os.environ.get("DEA_MAX_QUEUE_SIZE", "16")),
    default_timeout_seconds=float(os.environ.get("DEA_DEFAULT_TIMEOUT_SECONDS", "300")),
    max_timeout_seconds=float(os.environ.get("DEA_MAX_TIMEOUT_SECONDS", "3600")),
    recovery_policy=RecoveryPolicy(os.environ.get("DEA_RECOVERY_POLICY", "resume")),
)

app = create_app(settings)
