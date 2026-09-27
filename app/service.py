"""FastAPI wiring: HTTP surface over the task manager."""

from __future__ import annotations

from fastapi import FastAPI

from app.config import Settings
from app.manager import TaskManager
from app.store import TaskStore


def create_app(settings: Settings | None = None) -> FastAPI:
    """Build the app with its own manager/store (test-friendly factory)."""
    raise NotImplementedError
