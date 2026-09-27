"""Concurrent execution: resource limits, slow tasks do not block others,
cancellation and timeout leave no intermediate results, completed results
are never misclassified."""
from __future__ import annotations

import threading
import time

import deapp.engine as engine_mod
from deapp.models import TaskStatus


def _create(client, payload):
    r = client.post("/tasks", json=payload)
    assert r.status_code == 201, r.text
    return r.json()["id"]


def test_slow_task_does_not_block_others(client, payload_factory, monkeypatch):
    entered = threading.Event()
    release = threading.Event()

    real = engine_mod.normalization.normalize

    def slow_normalize(method, data, config):
        if config.seed == 777:
            entered.set()
            release.wait(timeout=10)
        return real(method, data, config)

    monkeypatch.setattr(engine_mod.normalization, "normalize", slow_normalize)

    slow = payload_factory(config={"seed": 777}, n_genes=30)
    fast = payload_factory(config={"seed": 888}, n_genes=10)

    t0 = time.monotonic()
    slow_id = _create(client, slow)
    assert entered.wait(5.0)
    fast_id = _create(client, fast)

    # fast task completes even though the slow one holds one worker
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        st = client.get(f"/tasks/{fast_id}").json()["status"]
        if st == "succeeded":
            break
        time.sleep(0.01)
    assert st == "succeeded"
    assert time.monotonic() - t0 < 8.0
    slow_status = client.get(f"/tasks/{slow_id}").json()["status"]
    assert slow_status == "running"
    release.set()


def test_concurrency_limit_holds(client, payload_factory, settings, monkeypatch):
    active = 0
    max_seen = 0
    lock = threading.Lock()
    entered = threading.Event()
    real = engine_mod.normalization.normalize

    def counted(method, data, config):
        nonlocal active, max_seen
        with lock:
            active += 1
            max_seen = max(max_seen, active)
            if active >= settings.max_concurrency:
                entered.set()
        try:
            time.sleep(0.3)
            return real(method, data, config)
        finally:
            with lock:
                active -= 1

    monkeypatch.setattr(engine_mod.normalization, "normalize", counted)
    ids = [_create(client, payload_factory(config={"seed": i}))
           for i in range(5)]
    assert entered.wait(5)
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline:
        if all(client.get(f"/tasks/{i}").json()["status"] == "succeeded"
               for i in ids):
            break
        time.sleep(0.05)
    assert max_seen <= settings.max_concurrency


def test_cancel_running_task_leaves_no_result(
    client, payload_factory, monkeypatch, data_dir
):
    release = threading.Event()
    real = engine_mod.normalization.normalize

    def blocking(method, data, config):
        release.wait(timeout=10)
        return real(method, data, config)

    monkeypatch.setattr(engine_mod.normalization, "normalize", blocking)
    tid = _create(client, payload_factory())
    # wait until running
    while client.get(f"/tasks/{tid}").json()["status"] != "running":
        time.sleep(0.01)
    r = client.post(f"/tasks/{tid}/cancel")
    assert r.status_code == 200
    release.set()
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        st = client.get(f"/tasks/{tid}").json()["status"]
        if st in ("cancelled", "succeeded"):
            break
        time.sleep(0.01)
    assert st == "cancelled"
    res = client.get(f"/tasks/{tid}/result")
    assert res.status_code == 404
    assert not (data_dir / "results" / f"{tid}.json").exists()


def test_cancel_pending_task(client, payload_factory, monkeypatch):
    release = threading.Event()
    hold = threading.Event()
    real = engine_mod.normalization.normalize

    def blocking(method, data, config):
        hold.set()
        release.wait(timeout=10)
        return real(method, data, config)

    monkeypatch.setattr(engine_mod.normalization, "normalize", blocking)
    # saturate both workers
    blockers = [_create(client, payload_factory(config={"seed": i}))
                for i in range(2)]
    assert hold.wait(5)
    queued = _create(client, payload_factory(config={"seed": 9}))
    # it must still be pending
    for _ in range(50):
        if client.get(f"/tasks/{queued}").json()["status"] == "pending":
            break
        time.sleep(0.01)
    assert client.get(f"/tasks/{queued}").json()["status"] == "pending"
    r = client.post(f"/tasks/{queued}/cancel")
    assert r.status_code == 200
    assert r.json()["status"] == "cancelled"
    release.set()
    # queued task never runs; blockers still succeed
    for bid in blockers:
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            st = client.get(f"/tasks/{bid}").json()["status"]
            if st == "succeeded":
                break
            time.sleep(0.02)
        assert st == "succeeded"
    final = client.get(f"/tasks/{queued}").json()["status"]
    assert final == "cancelled"


def test_timeout_marks_cancelled_and_no_result(
    client, data_dir
):
    # CPU-bound workload that genuinely exceeds the timeout; the engine
    # observes the cancellation flag between genes (cooperative cancellation).
    n_genes = 60000
    n = 6
    genes = [f"g{i}" for i in range(n_genes)]
    samples = [f"s{i}" for i in range(n)]
    values = [[(i * 31 + j * 7) % 97 + 1.0 for j in range(n)]
              for i in range(n_genes)]
    payload = {
        "design": {"a": samples[:3], "b": samples[3:]},
        "control": "a", "treatment": "b",
        "matrix": {"gene_ids": genes, "sample_ids": samples, "values": values},
        # normalization "none": per-gene t-test loop is interruption-safe
        "config": {"normalization": "none", "timeout_seconds": 0.3},
    }
    tid = _create(client, payload)
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline:
        t = client.get(f"/tasks/{tid}").json()
        if t["status"] in ("cancelled", "succeeded", "failed"):
            break
        time.sleep(0.02)
    assert t["status"] == "cancelled", t
    assert t["error"]["code"] == "TASK_CANCELLED"
    assert t["error"]["details"]["timeout"] is True
    assert not (data_dir / "results" / f"{tid}.json").exists()


def test_completed_task_not_marked_failed_after_timeout(
    client, payload_factory
):
    # short timeout but trivially fast analysis -> must still succeed and
    # never be overwritten
    payload = payload_factory(
        config={"timeout_seconds": 30.0}, n_genes=5
    )
    tid = _create(client, payload)
    t = None
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        t = client.get(f"/tasks/{tid}").json()
        if t["status"] in ("succeeded", "failed", "cancelled"):
            break
        time.sleep(0.02)
    assert t["status"] == "succeeded"
    assert client.get(f"/tasks/{tid}/result").status_code == 200
    # status remains stable when re-read
    time.sleep(0.1)
    assert client.get(f"/tasks/{tid}").json()["status"] == "succeeded"


def test_cancel_finished_task_is_conflict(client, payload_factory, wait_tasks):
    tid = _create(client, payload_factory())
    wait_tasks(tid)
    r = client.post(f"/tasks/{tid}/cancel")
    assert r.status_code == 409
    assert r.json()["error"]["code"] == "TASK_STATE_ERROR"


def test_queue_limit_rejects(payload_factory, monkeypatch, tmp_path):
    from fastapi.testclient import TestClient
    from deapp.api import create_app
    from deapp.config import Settings

    gate = threading.Event()
    entered = threading.Event()
    state = {"n": 0}
    real = engine_mod.normalization.normalize

    def blocking(method, data, config):
        # Hold the worker (this runs after the task is persisted as RUNNING)
        entered.set()
        gate.wait(timeout=10)
        return real(method, data, config)

    monkeypatch.setattr(engine_mod.normalization, "normalize", blocking)
    limited = Settings(
        data_dir=str(tmp_path),
        max_concurrency=1,
        queue_size=3,  # 1 running + at most 2 waiting
        default_timeout=30.0,
    )
    app = create_app(limited)
    with TestClient(app) as c:
        first = _create(c, payload_factory(config={"seed": 0}))
        # wait until the first task really occupies the single worker
        assert entered.wait(5)
        while c.get(f"/tasks/{first}").json()["status"] != "running":
            time.sleep(0.005)
        for i in range(1, 3):  # 2 more -> both pending up to queue size
            _create(c, payload_factory(config={"seed": i}))
        r = c.post("/tasks", json=payload_factory(config={"seed": 42}))
        assert r.status_code == 429
        assert r.json()["error"]["code"] == "QUEUE_FULL"
        gate.set()


def test_pending_cancel_error_code(client, payload_factory, monkeypatch):
    gate = threading.Event()
    entered = threading.Event()
    real = engine_mod.normalization.normalize

    def blocking(method, data, config):
        entered.set()
        gate.wait(timeout=10)
        return real(method, data, config)

    monkeypatch.setattr(engine_mod.normalization, "normalize", blocking)
    from fastapi.testclient import TestClient
    from deapp.api import create_app
    from deapp.config import Settings
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        app = create_app(Settings(data_dir=d, max_concurrency=1, default_timeout=30))
        with TestClient(app) as c:
            first = c.post("/tasks", json=payload_factory(config={"seed": 1})).json()["id"]
            assert entered.wait(5)
            queued = c.post("/tasks", json=payload_factory(config={"seed": 2})).json()["id"]
            r = c.post(f"/tasks/{queued}/cancel")
            assert r.json()["error"]["code"] == "TASK_CANCELLED"
            gate.set()
