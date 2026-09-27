"""Degenerate-but-valid data: constant groups, infinite separation and
zero fold change must round-trip through JSON/persistence without NaN."""
from __future__ import annotations

import math
import time


def _wait(c, tid):
    end = time.monotonic() + 10
    while time.monotonic() < end:
        t = c.get(f"/tasks/{tid}").json()
        if t["status"] in ("succeeded", "failed", "cancelled"):
            return t
        time.sleep(0.02)
    raise AssertionError("timeout")


PAYLOAD = {
    "design": {"a": ["s1", "s2", "s3"], "b": ["s4", "s5", "s6"]},
    "control": "a",
    "treatment": "b",
    "matrix": {
        "gene_ids": ["sep_up", "sep_down", "identical"],
        "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "values": [
            [1, 1, 1, 10, 10, 10],
            [10, 10, 10, 1, 1, 1],
            [5, 5, 5, 5, 5, 5],
        ],
    },
}


def test_constant_groups_round_trip(tmp_path):
    from fastapi.testclient import TestClient
    from deapp.api import create_app
    from deapp.config import Settings

    with TestClient(create_app(Settings(data_dir=str(tmp_path)))) as c:
        for test in ("welch-ttest", "student-ttest"):
            p = {**PAYLOAD, "config": {"test": test}}
            tid = c.post("/tasks", json=p).json()["id"]
            t = _wait(c, tid)
            assert t["status"] == "succeeded", t.get("error")
            r = c.get(f"/tasks/{tid}/result")
            assert r.status_code == 200
            raw = r.text
            assert "NaN" not in raw and "Infinity" not in raw
            rows = {row["gene_id"]: row for row in r.json()["rows"]}
            assert rows["sep_up"]["p_value"] == 0.0
            assert rows["sep_up"]["statistic"] >= 1e308 - 1e294
            assert rows["sep_up"]["log2fc"] > 0
            assert rows["sep_down"]["statistic"] <= -1e308 + 1e294
            assert rows["identical"]["p_value"] == 1.0
            assert rows["identical"]["statistic"] == 0.0
            assert rows["identical"]["log2fc"] == 0.0
            for row in rows.values():
                assert math.isfinite(row["statistic"])
                assert 0.0 <= row["p_value_adjusted"] <= 1.0


def test_paired_constant_difference(tmp_path):
    from fastapi.testclient import TestClient
    from deapp.api import create_app
    from deapp.config import Settings

    p = {**PAYLOAD, "config": {"test": "paired-ttest"}}
    with TestClient(create_app(Settings(data_dir=str(tmp_path)))) as c:
        tid = c.post("/tasks", json=p).json()["id"]
        t = _wait(c, tid)
        assert t["status"] == "succeeded", t.get("error")
        rows = {r["gene_id"]: r for r in c.get(f"/tasks/{tid}/result").json()["rows"]}
        assert rows["sep_up"]["p_value"] == 0.0
        assert math.isfinite(rows["sep_up"]["statistic"])
