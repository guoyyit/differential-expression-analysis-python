"""Every advertised normalization/correction combination works end to end,
and results are reproducible across fresh service instances (simulating
separate processes, not just separate requests)."""
from __future__ import annotations

import time

from fastapi.testclient import TestClient

from deapp.api import create_app
from deapp.config import Settings


def _wait(c, tid, timeout=20):
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        t = c.get(f"/tasks/{tid}").json()
        if t["status"] in ("succeeded", "failed", "cancelled"):
            return t
        time.sleep(0.02)
    raise AssertionError("timeout")


def _payload(normalization, correction="benjamini-hochberg", test="welch-ttest"):
    genes = [f"g{i:03d}" for i in range(20)]
    samples = [f"s{i}" for i in range(8)]
    values = [
        [((i + 1) * 13 + (j + 3) * 7) % 50 + 1.0 for j in range(8)]
        for i in range(20)
    ]
    cfg = {"normalization": normalization, "correction": correction,
           "test": test, "seed": 7}
    if normalization == "tpm":
        cfg["gene_lengths"] = {g: 100.0 + i * 10 for i, g in enumerate(genes)}
    return {
        "design": {"ctl": samples[:4], "trt": samples[4:]},
        "control": "ctl", "treatment": "trt",
        "matrix": {"gene_ids": genes, "sample_ids": samples, "values": values},
        "config": cfg,
    }


def test_all_normalizations_run(tmp_path):
    for norm in ("none", "cpm", "tpm", "quantile", "tmm"):
        app = create_app(Settings(data_dir=str(tmp_path / norm)))
        with TestClient(app) as c:
            tid = c.post("/tasks", json=_payload(norm)).json()["id"]
            t = _wait(c, tid)
            assert t["status"] == "succeeded", (norm, t.get("error"))
            rows = c.get(f"/tasks/{tid}/result").json()["rows"]
            assert len(rows) == 20


def test_all_corrections_run(tmp_path):
    for corr in ("bonferroni", "benjamini-hochberg", "benjamini-yekutieli", "holm"):
        d = tmp_path / corr
        app = create_app(Settings(data_dir=str(d)))
        with TestClient(app) as c:
            tid = c.post("/tasks", json=_payload("cpm", corr)).json()["id"]
            t = _wait(c, tid)
            assert t["status"] == "succeeded", t.get("error")
            assert t["strategy"]["correction"] == corr


def test_mann_whitney_pipeline(tmp_path):
    app = create_app(Settings(data_dir=str(tmp_path / "mw")))
    with TestClient(app) as c:
        tid = c.post("/tasks", json=_payload("none", test="mann-whitney")).json()["id"]
        t = _wait(c, tid)
        assert t["status"] == "succeeded"
        rows = c.get(f"/tasks/{tid}/result").json()["rows"]
        assert all(0.0 <= r["p_value_adjusted"] <= 1.0 for r in rows)


def test_reproducibility_across_service_instances(tmp_path):
    payload = _payload("tmm", "benjamini-yekutieli")

    def run_once():
        app = create_app(Settings(data_dir=str(tmp_path / "a")))
        with TestClient(app) as c:
            tid = c.post("/tasks", json=payload).json()["id"]
            _wait(c, tid)
            return c.get(f"/tasks/{tid}/result").json()

    first = run_once()
    # a genuinely separate instance with its own data dir
    payload2 = _payload("tmm", "benjamini-yekutieli")
    app = create_app(Settings(data_dir=str(tmp_path / "b")))
    with TestClient(app) as c:
        tid = c.post("/tasks", json=payload2).json()["id"]
        _wait(c, tid)
        second = c.get(f"/tasks/{tid}/result").json()
    assert first["rows"] == second["rows"]


def test_interleaved_concurrent_tasks_remain_consistent(tmp_path):
    # submit several tasks then interleave result reads; each task's output
    # must remain self-consistent despite concurrent store activity
    payload = _payload("quantile")
    app = create_app(Settings(data_dir=str(tmp_path / "x"), max_concurrency=4))
    with TestClient(app) as c:
        ids = [c.post("/tasks", json=_payload("quantile", corr)).json()["id"]
               for corr in ("bonferroni", "benjamini-hochberg", "holm",
                            "benjamini-yekutieli")]
        snapshots = {}
        for tid in ids:
            _wait(c, tid)
            snapshots[tid] = c.get(f"/tasks/{tid}/result").json()["rows"]
        # repeat reads must be identical (no mutation after completion)
        for tid in ids:
            assert c.get(f"/tasks/{tid}/result").json()["rows"] == snapshots[tid]


def test_concurrent_load_results_are_individually_reproducible(tmp_path):
    # 8 different configs executed together; each result must equal the
    # result of the same payload run alone in a pristine data dir.
    configs = [
        ("none", "benjamini-hochberg", "welch-ttest"),
        ("cpm", "bonferroni", "student-ttest"),
        ("tpm", "holm", "welch-ttest"),
        ("quantile", "benjamini-yekutieli", "welch-ttest"),
        ("tmm", "benjamini-hochberg", "mann-whitney"),
        ("none", "bonferroni", "paired-ttest"),
        ("cpm", "benjamini-yekutieli", "welch-ttest"),
        ("quantile", "holm", "student-ttest"),
    ]
    payloads = [_payload(n, c, t) for n, c, t in configs]
    app = create_app(Settings(data_dir=str(tmp_path / "parallel"), max_concurrency=4))
    parallel_rows: dict[int, list] = {}
    with TestClient(app) as c:
        ids = [c.post("/tasks", json=p).json()["id"] for p in payloads]
        for i, tid in enumerate(ids):
            _wait(c, tid)
            parallel_rows[i] = c.get(f"/tasks/{tid}/result").json()["rows"]
    for i, p in enumerate(payloads):
        d = tmp_path / f"solo_{i}"
        with TestClient(create_app(Settings(data_dir=str(d)))) as c:
            tid = c.post("/tasks", json=p).json()["id"]
            _wait(c, tid)
            solo = c.get(f"/tasks/{tid}/result").json()["rows"]
        assert solo == parallel_rows[i], configs[i]
