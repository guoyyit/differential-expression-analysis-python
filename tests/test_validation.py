"""Illegal input must be rejected structurally without creating a task."""
from __future__ import annotations


def _base():
    return {
        "design": {"a": ["s1", "s2"], "b": ["s3", "s4"]},
        "control": "a",
        "treatment": "b",
        "matrix": {
            "gene_ids": ["g1", "g2"],
            "sample_ids": ["s1", "s2", "s3", "s4"],
            "values": [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
            ],
        },
        "config": {},
    }


def test_valid_request_accepted(client, wait_tasks):
    r = client.post("/tasks", json=_base())
    assert r.status_code == 201
    wait_tasks(r.json()["id"])


def test_too_few_samples_per_group(client):
    p = _base()
    p["design"] = {"a": ["s1"], "b": ["s3", "s4"]}
    r = client.post("/tasks", json=p)
    assert r.status_code == 422
    body = r.json()
    assert body["error"]["code"] == "VALIDATION_ERROR"
    assert body["error"]["details"]["reason"] == "not_enough_samples"


def test_group_count_must_be_two(client):
    p = _base()
    p["design"] = {"a": ["s1", "s2"], "b": ["s3", "s4"], "c": ["s5", "s6"]}
    r = client.post("/tasks", json=p)
    assert r.status_code == 422
    assert r.json()["error"]["details"]["reason"] == "invalid_group_count"


def test_empty_group_rejected(client):
    p = _base()
    p["design"] = {"a": ["s1", "s2"], "b": []}
    r = client.post("/tasks", json=p)
    assert r.status_code == 422


def test_sample_in_multiple_groups(client):
    p = _base()
    p["design"] = {"a": ["s1", "s2"], "b": ["s1", "s4"]}
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "sample_in_multiple_groups"


def test_design_sample_not_in_matrix(client):
    p = _base()
    p["design"] = {"a": ["s1", "s2"], "b": ["sX", "s4"]}
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "unknown_sample"


def test_duplicate_samples_in_matrix(client):
    p = _base()
    p["matrix"]["sample_ids"] = ["s1", "s1", "s3", "s4"]
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "duplicate_sample"


def test_matrix_row_dimension_mismatch(client):
    p = _base()
    p["matrix"]["values"] = [[1, 2, 3], [5, 6, 7, 8]]
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "dimension_mismatch"


def test_matrix_row_count_mismatch(client):
    p = _base()
    p["matrix"]["values"] = [[1, 2, 3, 4]]
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "dimension_mismatch"


def _post_raw(client, p):
    import json

    return client.post(
        "/tasks",
        content=json.dumps(p),
        headers={"content-type": "application/json"},
    )


def test_missing_value_rejected(client):
    p = _base()
    p["matrix"]["values"][0][2] = float("nan")
    r = _post_raw(client, p)
    assert r.status_code == 422
    assert r.json()["error"]["details"]["reason"] == "missing_value"


def test_infinite_value_rejected(client):
    p = _base()
    p["matrix"]["values"][0][2] = float("inf")
    r = _post_raw(client, p)
    assert r.status_code == 422
    assert r.json()["error"]["details"]["reason"] == "non_finite_value"


def test_duplicate_gene_rejected(client):
    p = _base()
    p["matrix"]["gene_ids"] = ["g1", "g1"]
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "duplicate_gene"


def test_paired_requires_equal_groups(client):
    p = _base()
    p["config"] = {"test": "paired-ttest"}
    p["design"] = {"a": ["s1", "s2"], "b": ["s3"]}
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "unequal_paired_groups"


def test_paired_equal_groups_ok(client, wait_tasks):
    p = _base()
    p["config"] = {"test": "paired-ttest"}
    p["design"] = {"a": ["s1", "s2"], "b": ["s3", "s4"]}
    r = client.post("/tasks", json=p)
    assert r.status_code == 201
    wait_tasks(r.json()["id"])


def test_tpm_requires_gene_lengths(client):
    p = _base()
    p["config"] = {"normalization": "tpm"}
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "missing_gene_lengths"


def test_tpm_incomplete_gene_lengths(client):
    p = _base()
    p["config"] = {"normalization": "tpm", "gene_lengths": {"g1": 100}}
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "incomplete_gene_lengths"


def test_bad_enum_value_is_structured(client):
    p = _base()
    p["config"] = {"test": "nonsense"}
    r = client.post("/tasks", json=p)
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "INVALID_REQUEST"


def test_rejected_input_creates_no_task(client):
    p = _base()
    p["design"] = {"a": ["s1"], "b": ["s3", "s4"]}
    r = client.post("/tasks", json=p)
    assert r.status_code == 422
    listed = client.get("/tasks").json()["tasks"]
    assert listed == []


def test_unknown_comparison_group(client):
    p = _base()
    p["control"] = "zzz"
    r = client.post("/tasks", json=p)
    assert r.json()["error"]["details"]["reason"] == "unknown_group"
