"""Structured errors. Every client-facing failure carries a stable code."""

from __future__ import annotations

from typing import Any


class ServiceError(Exception):
    """Base class. ``code`` is stable API surface; ``details`` is diagnostic."""

    code = "internal_error"
    http_status = 500

    def __init__(self, message: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.message = message
        self.details = details or {}

    def to_body(self) -> dict[str, Any]:
        return {
            "error": {
                "code": self.code,
                "message": self.message,
                "details": self.details,
            }
        }


class ValidationError(ServiceError):
    code = "validation_error"
    http_status = 422


class TaskNotFoundError(ServiceError):
    code = "task_not_found"
    http_status = 404


class TaskConflictError(ServiceError):
    code = "task_conflict"
    http_status = 409


class QueueFullError(ServiceError):
    code = "queue_full"
    http_status = 429


class ResultNotReadyError(ServiceError):
    code = "result_not_ready"
    http_status = 409
