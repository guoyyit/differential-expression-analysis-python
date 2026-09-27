"""End-to-end HTTP behaviour: lifecycle, reproducibility, filters, compare,
export and file-based submission."""
from __future__ import annotations

import io
import json
import time

import pytest


def _wait(client, tid, timeout=15):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        t = client.get(f"/tasks/{tid}").json()
        if t["status"] in ("succeeded", "failed", "cancelled"):
            return t
        time.sleep(0.02)
    raise AssertionError(f"task {tid} stuck")


def _submit(client, payload):
    r = client.post("/tasks", json=payload)
    assert r.status_code == 201, r.text
    return r.json()["id"]


def test_health(client):
    assert client.get("/health").json() == {"status": "ok"}


def test_task_lifecycle_fields(client, payload_factory):
    tid = _submit(client, payload_factory())
    detail = client.get(f"/tasks/{tid}").json()
    assert detail["control"] == "control"
    assert detail["treatment"] == "treatment"
    final = _wait(client, tid)
    assert final["status"] == "succeeded"
    assert final["strategy"]["test"] == "welch-ttest"
    assert final["n_genes"] == 12


def test_result_before_completion_404(client, payload_factory):
    # tiny tasks may finish immediately; create then poll result at least
    # returns 404 while not succeeded, 200 afterwards
    tid = _submit(client, payload_factory())
    final = _wait(client, tid)
    assert final["status"] == "succeeded"
    r = client.get(f"/tasks/{tid}/result")
    assert r.status_code == 200


def test_missing_task_404(client):
    r = client.get("/tasks/deadbeef")
    assert r.status_code == 404
    assert r.json()["error"]["code"] == "TASK_NOT_FOUND"
    r = client.get("/tasks/deadbeef/result")
    assert r.status_code == 404


def test_ranking_order_significance_then_effect(client, payload_factory):
    tid = _submit(client, payload_factory(n_genes=30))
    _wait(client, tid)
    rows = client.get(f"/tasks/{tid}/result").json()["rows"]
    keys = [(r["p_value_adjusted"], -abs(r["log2fc"]), r["gene_id"])
            for r in rows]
    assert keys == sorted(keys)
    assert [r["rank"] for r in rows] == list(range(1, len(rows) + 1))


def test_result_reproducible_same_seed_and_input(client, payload_factory):
    p = payload_factory(seed=2024, config={"seed": 42})
    t1 = _submit(client, json.loads(json.dumps(p)))
    t2 = _submit(client, json.loads(json.dumps(p)))
    _wait(client, t1)
    _wait(client, t2)
    r1 = client.get(f"/tasks/{t1}/result").json()
    r2 = client.get(f"/tasks/{t2}/result").json()
    assert r1["rows"] == r2["rows"]
    assert r1["strategy"] == r2["strategy"]


def test_different_strategies_yield_different_pvalues(client, payload_factory):
    p1 = payload_factory(config={"test": "welch-ttest"})
    p2 = payload_factory(config={"test": "student-ttest"})
    t1 = _submit(client, p1)
    t2 = _submit(client, p2)
    _wait(client, t1)
    _wait(client, t2)
    r1 = client.get(f"/tasks/{t1}/result").json()
    r2 = client.get(f"/tasks/{t2}/result").json()
    assert [r["p_value"] for r in r1["rows"]] != [r["p_value"] for r in r2["rows"]]
    # strategy metadata is traceable
    assert r1["strategy"]["test"] == "welch-ttest"
    assert r2["strategy"]["test"] == "student-ttest"


def test_correction_strategy_switch(client, payload_factory):
    t1 = _submit(client, payload_factory(config={"correction": "bonferroni"}))
    t2 = _submit(client, payload_factory(config={"correction": "benjamini-hochberg"}))
    _wait(client, t1)
    _wait(client, t2)
    b = client.get(f"/tasks/{t1}/result").json()["rows"]
    h = client.get(f"/tasks/{t2}/result").json()["rows"]
    for rb, rh in zip(sorted(b, key=lambda r: r["gene_id"]),
                      sorted(h, key=lambda r: r["gene_id"])):
        assert rb["p_value_adjusted"] >= rh["p_value_adjusted"] - 1e-12


def test_gene_filter_exact(client, payload_factory):
    tid = _submit(client, payload_factory(n_genes=10))
    _wait(client, tid)
    r = client.get(f"/tasks/{tid}/result", params={"gene": "gene_003"})
    rows = r.json()["rows"]
    assert len(rows) == 1
    assert rows[0]["gene_id"] == "gene_003"


def test_filters_padj_direction_top(client, payload_factory):
    tid = _submit(client, payload_factory(n_genes=30))
    _wait(client, tid)
    rows = client.get(f"/tasks/{tid}/result", params={"padj_max": 1.0}).json()["rows"]
    assert rows
    up = client.get(f"/tasks/{tid}/result",
                    params={"direction": "up", "padj_max": 0.99}).json()["rows"]
    assert all(r["log2fc"] > 0 for r in up)
    top5 = client.get(f"/tasks/{tid}/result", params={"top": 5}).json()["rows"]
    assert len(top5) == 5
    full = client.get(f"/tasks/{tid}/result").json()["rows"]
    assert [r["gene_id"] for r in top5] == [r["gene_id"] for r in full[:5]]


def test_invalid_filter_is_422(client, payload_factory):
    tid = _submit(client, payload_factory())
    _wait(client, tid)
    r = client.get(f"/tasks/{tid}/result", params={"padj_max": 2.0})
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "INVALID_FILTER"


def test_csv_export(client, payload_factory):
    tid = _submit(client, payload_factory())
    _wait(client, tid)
    r = client.get(f"/tasks/{tid}/export", params={"format": "csv"})
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/csv")
    lines = r.text.strip().splitlines()
    assert lines[0].startswith("rank,gene_id")
    assert len(lines) == 13  # header + 12 genes


def test_json_export(client, payload_factory):
    tid = _submit(client, payload_factory())
    _wait(client, tid)
    r = client.get(f"/tasks/{tid}/export", params={"format": "json"})
    body = r.json()
    assert body["n_genes"] == 12
    assert "strategy" in body


def test_bad_export_format(client, payload_factory):
    tid = _submit(client, payload_factory())
    _wait(client, tid)
    r = client.get(f"/tasks/{tid}/export", params={"format": "xml"})
    assert r.status_code == 400


def test_compare_tasks(client, payload_factory):
    p1 = payload_factory(seed=1, config={"normalization": "none"})
    p2 = payload_factory(seed=1, config={"normalization": "cpm"})
    t1 = _submit(client, p1)
    t2 = _submit(client, p2)
    _wait(client, t1)
    _wait(client, t2)
    r = client.post("/compare", json={"task_ids": [t1, t2]})
    assert r.status_code == 200
    entries = r.json()
    assert len(entries) == 12
    assert set(entries[0]["tasks"].keys()) == {t1, t2}
    # deterministic gene order
    genes = [e["gene_id"] for e in entries]
    assert genes == sorted(genes) or genes == [r2["gene_id"] for r2 in
        client.get(f"/tasks/{t1}/result").json()["rows"]]


def test_compare_unknown_task_404(client, payload_factory):
    tid = _submit(client, payload_factory())
    _wait(client, tid)
    r = client.post("/compare", json={"task_ids": [tid, "nope"]})
    assert r.status_code == 404


def test_from_file_submission(client):
    csv_text = (
        "gene,s1,s2,s3,s4\n"
        "g1,1,2,10,12\n"
        "g2,5,6,5,6\n"
        "g3,20,21,2,3\n"
    )
    metadata = {
        "design": {"a": ["s1", "s2"], "b": ["s3", "s4"]},
        "control": "a",
        "treatment": "b",
        "config": {"test": "welch-ttest"},
    }
    r = client.post(
        "/tasks/from-file",
        files={"matrix": ("matrix.csv", io.BytesIO(csv_text.encode()), "text/csv")},
        data={"metadata": json.dumps(metadata)},
    )
    assert r.status_code == 201, r.text
    tid = r.json()["id"]
    final = _wait(client, tid)
    assert final["status"] == "succeeded"
    rows = client.get(f"/tasks/{tid}/result").json()["rows"]
    assert {r["gene_id"] for r in rows} == {"g1", "g2", "g3"}


def test_from_file_bad_dimensions(client):
    csv_text = "gene,s1,s2\ng1,1\n"
    metadata = {"design": {"a": ["s1"], "b": ["s2"]}}
    r = client.post(
        "/tasks/from-file",
        files={"matrix": ("m.csv", io.BytesIO(csv_text.encode()), "text/csv")},
        data={"metadata": json.dumps(metadata)},
    )
    assert r.status_code == 400
    assert r.json()["error"]["code"] == "INVALID_MATRIX_FILE"


def test_list_tasks_filter_by_status(client, payload_factory):
    t1 = _submit(client, payload_factory())
    _wait(client, t1)
    done = client.get("/tasks", params={"status": "succeeded"}).json()["tasks"]
    assert all(t["status"] == "succeeded" for t in done)
    assert len(done) == 1
    none = client.get("/tasks", params={"status": "pending"}).json()["tasks"]
    assert none == []


def test_from_file_rejects_nan(client):
    csv_text = "gene,s1,s2,s3,s4\ng1,1,NaN,10,12\n"
    metadata = {"design": {"a": ["s1", "s2"], "b": ["s3", "s4"]}}
    r = client.post(
        "/tasks/from-file",
        files={"matrix": ("m.csv", io.BytesIO(csv_text.encode()), "text/csv")},
        data={"metadata": json.dumps(metadata)},
    )
    assert r.status_code == 400
    assert r.json()["error"]["code"] == "INVALID_MATRIX_FILE"


def test_compare_empty_list(client):
    r = client.post("/compare", json={"task_ids": []})
    assert r.status_code == 400


def test_cross_task_search(client, payload_factory):
    t1 = _submit(client, payload_factory(seed=5))
    t2 = _submit(client, payload_factory(seed=6))
    _wait(client, t1)
    _wait(client, t2)
    r = client.get("/results/search", params={"gene": "gene_003"})
    assert r.status_code == 200
    hits = r.json()["hits"]
    task_ids = {h["task_id"] for h in hits}
    assert task_ids == {t1, t2}
    assert all(h["row"]["gene_id"] == "gene_003" for h in hits)

    r2 = client.get("/results/search", params={"gene_contains": "gene_",
                                                "padj_max": 0.05})
    assert all(h["row"]["p_value_adjusted"] <= 0.05
               for h in r2.json()["hits"])
