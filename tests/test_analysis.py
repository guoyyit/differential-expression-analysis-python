"""Statistical semantics, strategy switching, ranking and reproducibility."""
from __future__ import annotations

import threading

import numpy as np
import pytest
from scipy import stats as ss

from app.models import AnalysisConfig, CorrectionMethod, NormalizationMethod, TestMethod
from app.stats.correction import adjust_pvalues
from app.stats.pipeline import run_analysis
from app.stats.tests import mann_whitney, welch_t

from .conftest import submit_and_wait


# -- unit-level statistics ---------------------------------------------------

def test_welch_matches_scipy():
    rng = np.random.default_rng(7)
    a = rng.normal(0, 1, 9)
    b = rng.normal(0.8, 1.3, 11)
    t, p = welch_t(a, b)
    ref = ss.ttest_ind(b, a, equal_var=False)
    assert t == pytest.approx(ref.statistic)
    assert p == pytest.approx(ref.pvalue)


def test_mann_whitney_matches_scipy():
    rng = np.random.default_rng(11)
    a = rng.normal(0, 1, 10)
    b = rng.normal(1.0, 1, 12)
    _, p = mann_whitney(a, b)
    ref = ss.mannwhitneyu(b, a, alternative="two-sided", method="asymptotic")
    assert p == pytest.approx(ref.pvalue)


def test_constant_groups_are_defined():
    assert welch_t(np.ones(4), np.ones(4)) == (0.0, 1.0)
    assert mann_whitney(np.ones(4), np.ones(4)) == (0.0, 1.0)


def test_bh_is_monotone_and_bounded():
    p = np.array([0.9, 0.001, 0.5, 0.01, 0.04])
    q = adjust_pvalues(p, CorrectionMethod.BH)
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= -1e-15)
    assert np.all((q >= 0) & (q <= 1))
    assert q[1] == pytest.approx(0.005)


# -- end-to-end through the API ------------------------------------------------

def test_result_fields_and_ranking(client, payload):
    snap = submit_and_wait(client, payload)
    assert snap["status"] == "completed"
    assert snap["summary"]["n_genes_tested"] == 5
    page = client.get(f"/api/v1/analyses/{snap['task_id']}/results").json()
    rows = page["results"]
    assert [r["rank"] for r in rows] == [1, 2, 3, 4, 5]
    for r in rows:
        assert set(r) >= {"statistic", "p_value", "adjusted_p_value",
                          "effect_size", "cohens_d", "significant"}
        assert 0.0 <= r["p_value"] <= 1.0
        assert 0.0 <= r["adjusted_p_value"] <= 1.0
    # sorted by significance then effect magnitude
    adj = [r["adjusted_p_value"] for r in rows]
    assert adj == sorted(adj)
    # g0 (strongly up in trt) must be significant with positive effect
    g0 = next(r for r in rows if r["gene"] == "g0")
    assert g0["significant"] and g0["effect_size"] > 0
    # g4 is down-regulated in trt
    g4 = next(r for r in rows if r["gene"] == "g4")
    assert g4["effect_size"] < 0


def test_strategy_switch_changes_output_and_is_traceable(client, payload):
    base = submit_and_wait(client, payload)
    tid_base = base["task_id"]
    alt = submit_and_wait(client, payload, config={
        "normalization": "median_ratio",
        "test": "mann_whitney",
        "correction": "bonferroni",
    })
    tid_alt = alt["task_id"]
    snap = client.get(f"/api/v1/analyses/{tid_alt}").json()
    assert snap["config"]["normalization"] == "median_ratio"
    assert snap["config"]["test"] == "mann_whitney"
    assert snap["config"]["correction"] == "bonferroni"
    r_base = client.get(f"/api/v1/analyses/{tid_base}/results").json()["results"]
    r_alt = client.get(f"/api/v1/analyses/{tid_alt}/results").json()["results"]
    # different strategy => different statistics somewhere
    assert any(
        a["statistic"] != b["statistic"] or a["adjusted_p_value"] != b["adjusted_p_value"]
        for a, b in zip(sorted(r_base, key=lambda r: r["gene"]),
                        sorted(r_alt, key=lambda r: r["gene"]))
    )


def test_all_strategy_combinations_complete(client, payload):
    for norm in NormalizationMethod:
        for test in TestMethod:
            for corr in CorrectionMethod:
                snap = submit_and_wait(client, payload, config={
                    "normalization": norm.value,
                    "test": test.value,
                    "correction": corr.value,
                })
                assert snap["status"] == "completed", (norm, test, corr, snap)


def test_reproducible_same_input_same_output(client, payload):
    a = submit_and_wait(client, payload)
    b = submit_and_wait(client, payload)
    ra = client.get(f"/api/v1/analyses/{a['task_id']}/results",
                    params={"limit": 1000}).json()["results"]
    rb = client.get(f"/api/v1/analyses/{b['task_id']}/results",
                    params={"limit": 1000}).json()["results"]
    assert ra == rb  # row-by-row identical, including floats and ranks


def test_reproducible_under_concurrency(client, payload):
    """Concurrent scheduling must not perturb any task's output."""
    snap = submit_and_wait(client, payload)
    solo = client.get(f"/api/v1/analyses/{snap['task_id']}/results",
                      params={"limit": 1000}).json()["results"]
    # run several copies at once, then compare each against the solo run
    ids = [client.post("/api/v1/analyses", json=payload).json()["task_id"]
           for _ in range(4)]
    from .conftest import wait_terminal
    for tid in ids:
        assert wait_terminal(client, tid)["status"] == "completed"
        rows = client.get(f"/api/v1/analyses/{tid}/results",
                          params={"limit": 1000}).json()["results"]
        assert rows == solo


def test_pipeline_level_determinism(payload):
    cfg = AnalysisConfig()
    r1 = run_analysis(payload["genes"], payload["matrix"], [0, 1, 2], [3, 4, 5],
                      cfg, threading.Event())
    r2 = run_analysis(payload["genes"], payload["matrix"], [0, 1, 2], [3, 4, 5],
                      cfg, threading.Event())
    assert [r.model_dump() for r in r1] == [r.model_dump() for r in r2]


def test_min_total_count_prefilter(client, payload):
    p = dict(payload)
    p["genes"] = payload["genes"] + ["zero"]
    p["matrix"] = payload["matrix"] + [[0, 0, 0, 0, 0, 0]]
    snap = submit_and_wait(client, p)
    assert snap["summary"]["n_genes_input"] == 6
    assert snap["summary"]["n_genes_tested"] == 5  # all-zero gene dropped
