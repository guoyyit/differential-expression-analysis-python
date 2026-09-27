"""HTTP API behavior: filters, sorts, pagination, export, compare, idempotency."""
from __future__ import annotations

from .conftest import submit_and_wait, wait_terminal


def _completed_id(client, payload):
    return submit_and_wait(client, payload)["task_id"]


def test_list_filter_by_status(client, payload):
    _completed_id(client, payload)
    all_rows = client.get("/api/v1/analyses").json()
    assert len(all_rows) == 1
    completed = client.get("/api/v1/analyses", params={"status": "completed"}).json()
    assert {t["task_id"] for t in completed} == {all_rows[0]["task_id"]}
    assert client.get("/api/v1/analyses", params={"status": "pending"}).json() == []


def test_unknown_task_404(client):
    assert client.get("/api/v1/analyses/deadbeef").status_code == 404
    assert client.post("/api/v1/analyses/deadbeef/cancel").status_code == 404


def test_results_not_ready(client, payload):
    tid = client.post("/api/v1/analyses", json=payload).json()["task_id"]
    # the task is likely not completed yet; endpoint must be consistent
    snap = client.get(f"/api/v1/analyses/{tid}").json()
    r = client.get(f"/api/v1/analyses/{tid}/results")
    if snap["status"] == "completed":
        assert r.status_code == 200
    else:
        assert r.status_code == 409
        assert r.json()["error"]["code"] == "result_not_ready"
    wait_terminal(client, tid)


def test_filter_by_gene_list(client, payload):
    tid = _completed_id(client, payload)
    r = client.get(f"/api/v1/analyses/{tid}/results",
                   params=[("gene", "g0"), ("gene", "g4")]).json()
    assert r["total"] == 2
    assert {row["gene"] for row in r["results"]} == {"g0", "g4"}


def test_significance_filter(client, payload):
    tid = _completed_id(client, payload)
    all_rows = client.get(f"/api/v1/analyses/{tid}/results",
                          params={"limit": 1000}).json()["results"]
    sig = client.get(f"/api/v1/analyses/{tid}/results",
                     params={"significant_only": True, "limit": 1000}).json()["results"]
    assert sig == [r for r in all_rows if r["significant"]]
    thr = client.get(f"/api/v1/analyses/{tid}/results",
                     params={"max_adjusted_p": 0.01, "limit": 1000}).json()["results"]
    assert all(r["adjusted_p_value"] <= 0.01 for r in thr)
    eff = client.get(f"/api/v1/analyses/{tid}/results",
                     params={"min_abs_effect_size": 1.0, "limit": 1000}).json()["results"]
    assert all(abs(r["effect_size"]) >= 1.0 for r in eff)


def test_sorting_and_pagination(client, payload):
    tid = _completed_id(client, payload)
    rows = client.get(f"/api/v1/analyses/{tid}/results",
                      params={"sort": "effect_size", "order": "desc",
                              "limit": 1000}).json()["results"]
    effects = [r["effect_size"] for r in rows]
    assert effects == sorted(effects, reverse=True)
    page1 = client.get(f"/api/v1/analyses/{tid}/results",
                       params={"offset": 0, "limit": 2}).json()
    page2 = client.get(f"/api/v1/analyses/{tid}/results",
                       params={"offset": 2, "limit": 2}).json()
    assert page1["total"] == 5 and page2["total"] == 5
    assert {r["gene"] for r in page1["results"]}.isdisjoint(
        r["gene"] for r in page2["results"]
    )


def test_export_csv_and_json(client, payload):
    tid = _completed_id(client, payload)
    csv_r = client.get(f"/api/v1/analyses/{tid}/results/export")
    assert csv_r.status_code == 200
    assert csv_r.headers["content-type"].startswith("text/csv")
    lines = csv_r.text.splitlines()
    assert lines[0].split(",")[:5] == ["rank", "gene", "statistic",
                                       "p_value", "adjusted_p_value"]
    assert len(lines) == 6  # header + 5 genes
    json_r = client.get(f"/api/v1/analyses/{tid}/results/export",
                        params={"format": "json"})
    assert json_r.status_code == 200
    assert len(json_r.json()["results"]) == 5
    filtered = client.get(f"/api/v1/analyses/{tid}/results/export",
                          params={"significant_only": True})
    assert len(filtered.text.splitlines()) <= 6


def test_compare_tasks(client, payload):
    a = _completed_id(client, payload)
    body = {**payload, "config": {"correction": "bonferroni"}}
    b = _completed_id(client, body)
    cmp = client.post("/api/v1/compare", params={"task_a": a, "task_b": b}).json()
    assert cmp["n_common_genes"] == 5
    assert cmp["task_a"] == a and cmp["task_b"] == b
    assert cmp["effect_size_correlation"] == 1.0  # same data => same effects
    # same task twice
    same = client.post("/api/v1/compare", params={"task_a": a, "task_b": a}).json()
    assert same["n_both_significant"] == same["n_common_genes"]
    # unfinished/unknown task in compare
    assert client.post("/api/v1/compare",
                       params={"task_a": a, "task_b": "nope"}).status_code == 404


def test_idempotent_client_task_id(client, payload):
    body = {**payload, "task_id": "fixed-id"}
    r1 = client.post("/api/v1/analyses", json=body)
    r2 = client.post("/api/v1/analyses", json=body)
    assert r1.json()["task_id"] == "fixed-id"
    assert r2.json()["task_id"] == "fixed-id"
    assert len(client.get("/api/v1/analyses").json()) == 1
    other = {**body, "genes": ["x", "y", "p", "q", "r"]}
    other["matrix"] = [[float(v) for v in row] for row in body["matrix"]]
    r3 = client.post("/api/v1/analyses", json=other)
    assert r3.status_code == 409
    assert r3.json()["error"]["code"] == "task_conflict"
    assert wait_terminal(client, "fixed-id")["status"] == "completed"


def test_generated_ids_are_distinct(client, payload):
    r1 = client.post("/api/v1/analyses", json=payload).json()
    r2 = client.post("/api/v1/analyses", json=payload).json()
    assert r1["task_id"] != r2["task_id"]


def test_status_enum_validation(client):
    r = client.get("/api/v1/analyses", params={"status": "bogus"})
    assert r.status_code == 422
