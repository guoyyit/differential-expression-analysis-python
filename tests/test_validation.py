"""Invalid submissions: structured errors, and nothing is persisted."""
from __future__ import annotations

import copy
import json

import pytest


def _post(client, payload, **over):
    body = copy.deepcopy(payload)
    body.update(over)
    return client.post("/api/v1/analyses", json=body)


def _post_raw(client, payload, **over):
    """Post through a permissive encoder so NaN/Infinity reach the server."""
    body = copy.deepcopy(payload)
    body.update(over)
    return client.post(
        "/api/v1/analyses",
        content=json.dumps(body),  # allow_nan: emits NaN/Infinity literally
        headers={"content-type": "application/json"},
    )


def test_matrix_row_count_mismatch(client, payload):
    r = _post(client, payload, matrix=payload["matrix"][:-1])
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "validation_error"


def test_matrix_col_count_mismatch(client, payload):
    m = [row[:-1] for row in payload["matrix"]]
    r = _post(client, payload, matrix=m)
    assert r.status_code == 422
    assert r.json()["error"]["details"]["field"] == "matrix"


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_matrix_non_finite_values(client, payload, bad):
    m = copy.deepcopy(payload["matrix"])
    m[1][2] = bad
    r = _post_raw(client, payload, matrix=m)
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "validation_error"
    assert r.json()["error"]["details"]["field"] == "matrix"


def test_matrix_negative_values(client, payload):
    m = copy.deepcopy(payload["matrix"])
    m[1][2] = -1.0
    r = _post(client, payload, matrix=m)
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "validation_error"


def test_duplicate_genes_rejected(client, payload):
    r = _post(client, payload, genes=["g0", "g0", "g2", "g3", "g4"])
    assert r.status_code == 422


def test_duplicate_samples_rejected(client, payload):
    r = _post(client, payload,
              samples=["s1", "s1", "s3", "s4", "s5", "s6"])
    assert r.status_code == 422


def test_unknown_group_sample(client, payload):
    groups = [{"name": "ctrl", "samples": ["s1", "sX"]},
              {"name": "trt", "samples": ["s4"]}]
    r = _post(client, payload, groups=groups)
    assert r.status_code == 422
    assert r.json()["error"]["details"]["field"] == "groups"


def test_sample_in_two_groups(client, payload):
    groups = [{"name": "ctrl", "samples": ["s1", "s2"]},
              {"name": "trt", "samples": ["s2", "s3"]}]
    r = _post(client, payload, groups=groups)
    assert r.status_code == 422


def test_wrong_group_count(client, payload):
    r = _post(client, payload, groups=[{"name": "only", "samples": ["s1"]}])
    assert r.status_code == 422  # schema: exactly two groups required


def test_timeout_above_limit(client, payload, settings):
    cfg = {"timeout_seconds": settings.max_timeout_seconds + 1}
    r = _post(client, payload, config=cfg)
    assert r.status_code == 422


def test_invalid_submission_leaves_no_trace(client, payload):
    before = client.get("/api/v1/analyses").json()
    m = copy.deepcopy(payload["matrix"])
    m[0][0] = float("nan")
    assert _post_raw(client, payload, matrix=m).status_code == 422
    after = client.get("/api/v1/analyses").json()
    assert before == after  # no partial task, no partial results


def test_malformed_body_is_structured(client):
    r = client.post("/api/v1/analyses", json={"genes": []})
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "validation_error"
