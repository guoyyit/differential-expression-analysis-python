"""Persistence across restarts: resume / fail policies, no loss or dupes."""
from __future__ import annotations

import time

from fastapi.testclient import TestClient

from app.config import RecoveryPolicy, Settings
from app.manager import TaskManager
from app.models import AnalysisRequest
from app.service import create_app
from app.store import TaskStore

from .conftest import wait_terminal


def _payload(n: int = 5) -> dict:
    genes = [f"g{i}" for i in range(n)]
    matrix = [[float(10 * (i + 1) + j) for j in range(6)] for i in range(n)]
    if n:
        matrix[0] = [10, 12, 11, 60, 62, 61]
    return {
        "genes": genes,
        "samples": [f"s{j}" for j in range(6)],
        "matrix": matrix,
        "groups": [{"name": "ctrl", "samples": ["s0", "s1", "s2"]},
                   {"name": "trt", "samples": ["s3", "s4", "s5"]}],
    }


def _req(body: dict) -> AnalysisRequest:
    return AnalysisRequest.model_validate(body)


def _spawn(db_path, **over):
    store = TaskStore(db_path)
    mgr = TaskManager(Settings(db_path=db_path, **over), store)
    mgr.start()
    return mgr, store


def _wait(mgr, tid, seconds=20.0):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if mgr.get(tid).status.terminal:
            return mgr.get(tid)
        time.sleep(0.05)
    raise AssertionError(f"task {tid} did not finish")


def test_running_task_resumed_after_restart(tmp_path):
    db = str(tmp_path / "r.db")
    m1, s1 = _spawn(db, max_workers=1, slow_task_sleep_per_gene=0.3)
    tid = m1.submit(_req(_payload())).task_id
    time.sleep(0.2)
    assert m1.get(tid).status.value == "running"
    m1.shutdown(wait=False)  # hard stop: worker abandoned mid-gene
    s1.close()

    m2, s2 = _spawn(db, max_workers=1)
    snap = _wait(m2, tid)
    assert snap.status.value == "completed"
    assert snap.summary.n_genes_tested == 5
    # exactly one row per gene: recovered run never duplicates output
    rows = s2.load_results(tid)
    assert sorted(r.gene for r in rows) == sorted(_payload()["genes"])
    assert len({r.gene for r in rows}) == len(rows)
    m2.shutdown(wait=True)
    s2.close()


def test_pending_task_resumed_after_restart(tmp_path):
    db = str(tmp_path / "p.db")
    m1, s1 = _spawn(db, max_workers=1, slow_task_sleep_per_gene=0.3)
    m1.submit(_req(_payload(n=200)))  # occupies the single worker
    time.sleep(0.1)
    pending = m1.submit(_req(_payload()))
    m1.shutdown(wait=False)
    s1.close()

    m2, s2 = _spawn(db, max_workers=2)
    assert _wait(m2, pending.task_id).status.value == "completed"
    m2.shutdown(wait=True)
    s2.close()


def test_fail_policy_marks_interrupted(tmp_path):
    db = str(tmp_path / "f.db")
    m1, s1 = _spawn(db, max_workers=1, slow_task_sleep_per_gene=0.3)
    tid = m1.submit(_req(_payload())).task_id
    time.sleep(0.2)
    m1.shutdown(wait=False)
    s1.close()

    m2, s2 = _spawn(db, max_workers=1, recovery_policy=RecoveryPolicy.FAIL)
    snap = m2.get(tid)
    assert snap.status.value == "interrupted"
    assert snap.error["code"] == "interrupted"
    assert s2.load_results(tid) == []  # interrupted tasks expose no results
    m2.shutdown(wait=True)
    s2.close()


def test_completed_results_survive_restart_unchanged(tmp_path):
    db = str(tmp_path / "c.db")
    m1, s1 = _spawn(db, max_workers=1)
    tid = m1.submit(_req(_payload())).task_id
    _wait(m1, tid)
    before = [r.model_dump() for r in s1.load_results(tid)]
    m1.shutdown(wait=True)
    s1.close()

    m2, s2 = _spawn(db, max_workers=1)  # recovery must not re-run finished task
    assert m2.get(tid).status.value == "completed"
    after = [r.model_dump() for r in s2.load_results(tid)]
    assert before == after
    m2.shutdown(wait=True)
    s2.close()


def test_recovery_through_http_layer(tmp_path):
    db = str(tmp_path / "http.db")
    settings = Settings(db_path=db, max_workers=1, slow_task_sleep_per_gene=0.3)
    with TestClient(create_app(settings)) as c1:
        tid = c1.post("/api/v1/analyses", json=_payload()).json()["task_id"]
        time.sleep(0.2)
        assert c1.get(f"/api/v1/analyses/{tid}").json()["status"] == "running"
    # process restarts with the same db file
    settings2 = Settings(db_path=db, max_workers=1)
    with TestClient(create_app(settings2)) as c2:
        snap = wait_terminal(c2, tid, seconds=20)
        assert snap["status"] == "completed"
        rows = c2.get(f"/api/v1/analyses/{tid}/results",
                      params={"limit": 100}).json()["results"]
        assert len(rows) == 5
