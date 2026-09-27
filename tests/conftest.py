from __future__ import annotations

import time

import pytest
from fastapi.testclient import TestClient

from app.config import Settings
from app.service import create_app


@pytest.fixture()
def settings(tmp_path):
    return Settings(db_path=str(tmp_path / "test.db"), max_workers=2, max_queue_size=4)


@pytest.fixture()
def client(settings):
    app = create_app(settings)
    with TestClient(app) as c:
        yield c


@pytest.fixture()
def payload():
    """Deterministic two-group dataset; g0 is strongly upregulated in trt."""
    return {
        "genes": ["g0", "g1", "g2", "g3", "g4"],
        "samples": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "matrix": [
            [10, 12, 11, 60, 62, 61],
            [20, 22, 21, 20, 21, 22],
            [100, 101, 99, 100, 100, 101],
            [5, 6, 5, 30, 32, 31],
            [50, 49, 51, 10, 11, 9],
        ],
        "groups": [
            {"name": "ctrl", "samples": ["s1", "s2", "s3"]},
            {"name": "trt", "samples": ["s4", "s5", "s6"]},
        ],
    }


def wait_terminal(client: TestClient, task_id: str, seconds: float = 30.0) -> dict:
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        snap = client.get(f"/api/v1/analyses/{task_id}").json()
        if snap["status"] in ("completed", "failed", "cancelled", "timed_out",
                              "interrupted"):
            return snap
        time.sleep(0.02)
    raise AssertionError(f"task {task_id} did not reach a terminal state")


def submit_and_wait(client: TestClient, payload: dict, **over) -> dict:
    body = {**payload, **over}
    task_id = client.post("/api/v1/analyses", json=body).json()["task_id"]
    return wait_terminal(client, task_id)


@pytest.fixture()
def make_client(tmp_path):
    """Factory for clients with custom settings (slow tasks, pool sizes...)."""
    created = []

    def _make(**over):
        over.setdefault("db_path", str(tmp_path / f"db{len(created)}.db"))
        settings = Settings(**over)
        tc = TestClient(create_app(settings))
        tc.__enter__()
        created.append(tc)
        return tc, settings

    yield _make

    for tc in created:
        tc.__exit__(None, None, None)
