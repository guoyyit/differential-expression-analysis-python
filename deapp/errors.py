"""Structured error types used across the service.

Every error that reaches the HTTP layer is rendered as::

    {
        "error": {
            "code": "INVALID_SAMPLE_COUNT",
            "message": "human readable text",
            "details": {"...": "..."}
        }
    }
"""
from __future__ import annotations

from typing import Any


class ServiceError(Exception):
    """Base class for errors carrying a stable machine code and details."""

    code: str = "SERVICE_ERROR"
    http_status: int = 400

    def __init__(
        self,
        message: str,
        *,
        details: dict[str, Any] | None = None,
        code: str | None = None,
        http_status: int | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.details: dict[str, Any] = details or {}
        if code is not None:
            self.code = code
        if http_status is not None:
            self.http_status = http_status

    def to_dict(self) -> dict[str, Any]:
        return {
            "error": {
                "code": self.code,
                "message": self.message,
                "details": self.details,
            }
        }


class ValidationError(ServiceError):
    """The request payload failed semantic validation."""

    code = "VALIDATION_ERROR"
    http_status = 422


class TaskNotFoundError(ServiceError):
    """Referenced task id does not exist."""

    code = "TASK_NOT_FOUND"
    http_status = 404


class TaskStateError(ServiceError):
    """Operation not allowed in the current task state."""

    code = "TASK_STATE_ERROR"
    http_status = 409


class ResultNotFoundError(ServiceError):
    """No result is available (yet) for the task."""

    code = "RESULT_NOT_FOUND"
    http_status = 404


class InvalidFilterError(ServiceError):
    """A query filter value is malformed."""

    code = "INVALID_FILTER"
    http_status = 422
