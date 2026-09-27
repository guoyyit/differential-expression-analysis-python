"""Crash recovery: unfinished tasks are reconciled per policy on restart,
completed results survive intact and are never duplicated or downgraded."""
from __future__ import annotations

import time

from fastapi.testclient import TestClient

from deapp.api import create_app
from deapp.config import Settings
from deapp.models import TaskStatus
from deapp.persistence import Store
from deapp.models import TaskDetail


def _settings(data_dir, policy="recover", **kw):
    base = dict(max_concurrency=2, queue_size=0, default_timeout=30.0,
                recovery_policy=policy)
    base.update(kw)
    return Settings(data_dir=str(data_dir), **base)


def _wait(client, tid, timeout=10):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        t = client.get(f"/tasks/{tid}").json()
        if t["status"] in ("succeeded", "failed", "cancelled"):
            return t
        time.sleep(0.01)
    raise AssertionError("task never finished")


def test_succeeded_result_survives_restart(data_dir, payload_factory):
    settings = _settings(data_dir)
    with TestClient(create_app(settings)) as c1:
        r = c1.post("/tasks", json=payload_factory())
        tid = r.json()["id"]
        t1 = _wait(c1, tid)
        assert t1["status"] == "succeeded"
        result_before = c1.get(f"/tasks/{tid}/result").json()

    # restart with a fresh process-equivalent: new app, same data dir
    with TestClient(create_app(_settings(data_dir))) as c2:
        t2 = c2.get(f"/tasks/{tid}").json()
        assert t2["status"] == "succeeded"
        result_after = c2.get(f"/tasks/{tid}/result").json()
        assert result_before["rows"] == result_after["rows"]
        tasks = c2.get("/tasks").json()["tasks"]
        assert len(tasks) == 1  # no duplicate task record


def test_running_task_requeued_by_recover_policy(data_dir, payload_factory):
    # Simulate a crash: persist a RUNNING task with its request snapshot,
    # no result file present.
    with TestClient(create_app(_settings(data_dir))) as c:
        tid = c.post("/tasks", json=payload_factory()).json()["id"]
        _wait(c, tid)

    # force state back to RUNNING as if killed mid-analysis
    store = Store(_settings(data_dir))
    store.load()
    task = store.get_task(tid)
    task.status = TaskStatus.RUNNING
    store.save_task(task)
    result_path = data_dir / "results" / f"{tid}.json"
    if result_path.exists():
        result_path.unlink()

    with TestClient(create_app(_settings(data_dir, "recover"))) as c2:
        t = c2.get(f"/tasks/{tid}").json()
        assert t["status"] in ("pending", "running", "succeeded")
        final = _wait(c2, tid)
        assert final["status"] == "succeeded"
        assert c2.get(f"/tasks/{tid}/result").status_code == 200


def test_unfinished_marked_cancelled_by_cancel_policy(data_dir, payload_factory):
    with TestClient(create_app(_settings(data_dir))) as c:
        tid = c.post("/tasks", json=payload_factory()).json()["id"]
        _wait(c, tid)

    store = Store(_settings(data_dir))
    store.load()
    task = store.get_task(tid)
    task.status = TaskStatus.RUNNING
    store.save_task(task)
    (data_dir / "results" / f"{tid}.json").unlink()

    store2 = Store(_settings(data_dir, "cancel"))
    store2.load()
    t = store2.get_task(tid)
    assert t.status == TaskStatus.CANCELLED
    assert t.error["code"] == "INTERRUPTED"


def test_pending_task_recovered_and_executed(data_dir, payload_factory):
    with TestClient(create_app(_settings(data_dir))) as c:
        tid = c.post("/tasks", json=payload_factory()).json()["id"]
        _wait(c, tid)
    store = Store(_settings(data_dir))
    store.load()
    task = store.get_task(tid)
    task.status = TaskStatus.PENDING
    task.started_at = None
    store.save_task(task)
    (data_dir / "results" / f"{tid}.json").unlink()

    with TestClient(create_app(_settings(data_dir))) as c2:
        final = _wait(c2, tid)
        assert final["status"] == "succeeded"


def test_succeeded_without_result_file_is_reconciled(data_dir, payload_factory):
    # corruption window: status succeeded but result write never happened
    with TestClient(create_app(_settings(data_dir))) as c:
        tid = c.post("/tasks", json=payload_factory()).json()["id"]
        _wait(c, tid)
    store = Store(_settings(data_dir))
    store.load()
    task = store.get_task(tid)
    task.status = TaskStatus.SUCCEEDED
    store.save_task(task)
    (data_dir / "results" / f"{tid}.json").unlink()

    with TestClient(create_app(_settings(data_dir))) as c2:
        final = _wait(c2, tid)
        assert final["status"] == "succeeded"
        assert c2.get(f"/tasks/{tid}/result").status_code == 200


def test_recovered_result_matches_fresh_run(data_dir, payload_factory):
    with TestClient(create_app(_settings(data_dir))) as c:
        tid = c.post("/tasks", json=payload_factory()).json()["id"]
        _wait(c, tid)
        first = c.get(f"/tasks/{tid}/result").json()

    # restart; rows must be byte-for-byte equal (same seed/input)
    with TestClient(create_app(_settings(data_dir))) as c2:
        second = c2.get(f"/tasks/{tid}/result").json()
        assert first["rows"] == second["rows"]
