"""Concurrency: isolation, back-pressure, no blocking across tasks."""
from __future__ import annotations

import threading

from .conftest import submit_and_wait, wait_terminal


def _slow_payload(n=300):
    genes = [f"g{i}" for i in range(n)]
    matrix = [[float(10 + i + j) for j in range(6)] for i in range(n)]
    return {
        "genes": genes,
        "samples": [f"s{j}" for j in range(6)],
        "matrix": matrix,
        "groups": [{"name": "ctrl", "samples": ["s0", "s1", "s2"]},
                   {"name": "trt", "samples": ["s3", "s4", "s5"]}],
    }


def test_slow_task_does_not_block_others(make_client, payload):
    c, _ = make_client(max_workers=2, slow_task_sleep_per_gene=0.05)
    slow_id = c.post("/api/v1/analyses", json=_slow_payload(300)).json()["task_id"]
    fast_id = c.post("/api/v1/analyses", json=payload).json()["task_id"]
    fast = wait_terminal(c, fast_id, seconds=10)
    assert fast["status"] == "completed"
    slow_snap = c.get(f"/api/v1/analyses/{slow_id}").json()
    assert slow_snap["status"] in ("pending", "running")  # still chewing
    assert wait_terminal(c, slow_id)["status"] == "completed"


def test_queue_full_is_rejected(make_client):
    c, _ = make_client(max_workers=1, max_queue_size=1, slow_task_sleep_per_gene=0.05)
    # capacity = 1 worker + 1 queued
    c.post("/api/v1/analyses", json=_slow_payload(200))
    c.post("/api/v1/analyses", json=_slow_payload(200))
    r = c.post("/api/v1/analyses", json=_slow_payload(200))
    assert r.status_code == 429
    assert r.json()["error"]["code"] == "queue_full"


def test_cancel_running_task(make_client):
    c, _ = make_client(max_workers=1, slow_task_sleep_per_gene=0.05)
    tid = c.post("/api/v1/analyses", json=_slow_payload(500)).json()["task_id"]
    import time
    time.sleep(0.15)  # let it start
    r = c.post(f"/api/v1/analyses/{tid}/cancel")
    assert r.status_code == 200
    snap = wait_terminal(c, tid, seconds=10)
    assert snap["status"] == "cancelled"
    results = c.get(f"/api/v1/analyses/{tid}/results")
    assert results.status_code == 409  # no partial results ever exposed
    # the freed worker accepts new work
    assert wait_terminal(
        c, c.post("/api/v1/analyses", json=_slow_payload(10)).json()["task_id"]
    )["status"] == "completed"


def test_timeout_leaves_no_results(make_client):
    c, _ = make_client(max_workers=1, slow_task_sleep_per_gene=0.05)
    body = _slow_payload(500)
    body["config"] = {"timeout_seconds": 0.2}
    tid = c.post("/api/v1/analyses", json=body).json()["task_id"]
    snap = wait_terminal(c, tid, seconds=10)
    assert snap["status"] == "timed_out"
    results = c.get(f"/api/v1/analyses/{tid}/results")
    assert results.status_code == 409
    # after the slow worker finally exits, status stays timed_out
    import time
    time.sleep(1.0)
    assert c.get(f"/api/v1/analyses/{tid}").json()["status"] == "timed_out"


def test_completed_task_cannot_be_failed(make_client, payload):
    c, _ = make_client()
    snap = submit_and_wait(c, payload)
    assert snap["status"] == "completed"
    again = c.post(f"/api/v1/analyses/{snap['task_id']}/cancel").json()
    assert again["status"] == "completed"  # cancel is a no-op on terminal


def test_pagination_stable_while_other_tasks_update(make_client, payload):
    c, _ = make_client(max_workers=2, slow_task_sleep_per_gene=0.02)
    done = submit_and_wait(c, payload)
    c.post("/api/v1/analyses", json=_slow_payload(200))  # churn in background
    pages = []
    for offset in (0, 2, 4):
        p = c.get(f"/api/v1/analyses/{done['task_id']}/results",
                  params={"offset": offset, "limit": 2, "sort": "rank"}).json()
        pages.append([r["rank"] for r in p["results"]])
    flat = [r for page in pages for r in page]
    assert flat == [1, 2, 3, 4, 5]  # disjoint, ordered, complete coverage


def test_many_parallel_tasks_all_complete(make_client, payload):
    c, _ = make_client(max_workers=4, max_queue_size=8)
    ids = [c.post("/api/v1/analyses", json=payload).json()["task_id"]
           for _ in range(10)]
    statuses = [wait_terminal(c, t)["status"] for t in ids]
    assert statuses == ["completed"] * 10
