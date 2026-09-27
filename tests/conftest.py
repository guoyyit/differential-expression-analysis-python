"""Shared fixtures: isolated data dir, deterministic clients, payloads."""
from __future__ import annotations

import os
import random
import tempfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from deapp.api import create_app, Service
from deapp.config import Settings


@pytest.fixture
def data_dir():
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def settings(data_dir):
    return Settings(
        data_dir=str(data_dir),
        max_concurrency=2,
        queue_size=0,
        default_timeout=30.0,
        recovery_policy="recover",
    )


@pytest.fixture
def app(settings):
    return create_app(settings)


@pytest.fixture
def client(app):
    with TestClient(app) as c:
        yield c


@pytest.fixture
def service(settings):
    s = Service(settings)
    yield s
    s.runner.shutdown()


@pytest.fixture
def payload_factory():
    def make(n_control=4, n_treatment=4, n_genes=12, seed=123, effect=5.0,
             config=None, design_type="dict"):
        rng = random.Random(seed)
        c = [f"c{i}" for i in range(n_control)]
        t = [f"t{i}" for i in range(n_treatment)]
        samples = c + t
        values = []
        genes = [f"gene_{i:03d}" for i in range(n_genes)]
        for g in range(n_genes):
            base = 20 + g
            row = [max(0.1, base + rng.gauss(0, 1.0)) for _ in c]
            shift = effect if g % 3 == 0 else (0.0 if g % 3 == 1 else -effect)
            row += [max(0.1, base + shift + rng.gauss(0, 1.0)) for _ in t]
            values.append(row)
        design = {"control": c, "treatment": t}
        payload = {
            "name": f"task-seed-{seed}",
            "design": design,
            "control": "control",
            "treatment": "treatment",
            "matrix": {
                "gene_ids": genes,
                "sample_ids": samples,
                "values": values,
            },
            "config": config or {},
        }
        return payload
    return make


import time


@pytest.fixture
def wait_tasks(client):
    def _wait(task_id=None, timeout=10.0):
        deadline = time.monotonic() + timeout
        if task_id is not None:
            while time.monotonic() < deadline:
                t = client.get(f"/tasks/{task_id}").json()
                if t["status"] in ("succeeded", "failed", "cancelled"):
                    return t
                time.sleep(0.01)
            raise AssertionError(f"task {task_id} did not finish")
        while time.monotonic() < deadline:
            pending = [
                t for t in client.get("/tasks").json()["tasks"]
                if t["status"] in ("pending", "running")
            ]
            if not pending:
                return
            time.sleep(0.01)
        raise AssertionError("tasks did not settle")
    return _wait
