"""Minimal client walk-through. Run the server first:

    uv run uvicorn main:app --host 127.0.0.1 --port 8000

then:  uv run python examples/example_client.py
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import httpx

BASE = os.environ.get("DEA_BASE_URL", "http://127.0.0.1:8000")


def main() -> None:
    here = Path(__file__).parent
    metadata = {
        "name": "example",
        "design": {"control": ["s1", "s2", "s3"], "treatment": ["s4", "s5", "s6"]},
        "control": "control",
        "treatment": "treatment",
        "config": {
            "normalization": "cpm",
            "test": "welch-ttest",
            "correction": "benjamini-hochberg",
            "seed": 42,
        },
    }
    with open(here / "example_matrix.csv", "rb") as fh, httpx.Client() as http:
        r = http.post(
            f"{BASE}/tasks/from-file",
            files={"matrix": ("example_matrix.csv", fh, "text/csv")},
            data={"metadata": json.dumps(metadata)},
            timeout=30,
        )
        r.raise_for_status()
        task_id = r.json()["id"]
        print("created task:", task_id)

        while True:
            t = http.get(f"{BASE}/tasks/{task_id}").json()
            if t["status"] in ("succeeded", "failed", "cancelled"):
                break
            time.sleep(0.2)
        print("status:", t["status"], "| strategy:", t["strategy"])

        rows = http.get(f"{BASE}/tasks/{task_id}/result").json()["rows"]
        for row in rows:
            print(
                f"#{row['rank']} {row['gene_id']:8s} "
                f"log2fc={row['log2fc']:+.3f} "
                f"p={row['p_value']:.3g} padj={row['p_value_adjusted']:.3g}"
            )

        csv_export = http.get(
            f"{BASE}/tasks/{task_id}/export", params={"format": "csv"}
        )
        print("\nCSV export (first 3 lines):")
        print("\n".join(csv_export.text.splitlines()[:3]))


if __name__ == "__main__":
    main()
