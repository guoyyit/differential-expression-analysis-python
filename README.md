# Local Differential Expression Analysis Service

A self-contained HTTP service for transcriptome differential-expression
analysis. Submit an expression matrix plus a two-group design; the service
normalizes the data, runs a per-gene statistical test, applies multiple
testing correction, and returns results ranked by significance and effect
size. Everything runs against the local loopback and the local filesystem —
no external services are contacted.

All computations are implemented in pure Python with deterministic summation
(`math.fsum`), stable tie-breaking and recorded strategy metadata: **the same
input, seed and configuration always produce byte-identical results**,
independent of thread scheduling.

## 1. Running

```bash
uv sync
uv run uvicorn main:app --host 127.0.0.1 --port 8000
# or
uv run python -m deapp
```

OpenAPI / Swagger UI: <http://127.0.0.1:8000/docs>

### Configuration (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `DEA_DATA_DIR` | `data` | Directory for persisted tasks, requests and results |
| `DEA_MAX_CONCURRENCY` | `2` | At most this many tasks execute simultaneously |
| `DEA_QUEUE_SIZE` | `0` | Max number of **unfinished** tasks (running + queued); `0` = unlimited. A full queue answers `429 QUEUE_FULL` |
| `DEA_DEFAULT_TIMEOUT` | `30` | Default per-task timeout in seconds; `0` disables it |
| `DEA_RECOVERY_POLICY` | `recover` | `recover` re-queues interrupted tasks on startup; `cancel` marks them `cancelled` with error `INTERRUPTED` |

## 2. API overview

| Method & path | Purpose |
|---|---|
| `POST /tasks` | Create an analysis task from a JSON payload |
| `POST /tasks/from-file` | Create a task uploading a CSV matrix (`multipart/form-data`) |
| `GET /tasks` | List tasks (optional `?status=succeeded`) |
| `GET /tasks/{id}` | Task status, strategy metadata and error details |
| `POST /tasks/{id}/cancel` | Cooperatively cancel a queued/running task |
| `GET /tasks/{id}/result` | Completed result with filter/query parameters |
| `GET /tasks/{id}/export?format=csv\|json` | Download filtered results |
| `GET /results/search` | Cross-task search by gene / significance / direction |
| `POST /compare` | Join results of several tasks gene by gene |
| `GET /health` | Liveness probe |

### Creating a task

```json
{
  "name": "experiment-1",
  "design": {"control": ["s1", "s2", "s3"], "treatment": ["s4", "s5", "s6"]},
  "control": "control",
  "treatment": "treatment",
  "matrix": {
    "gene_ids": ["GENE_A", "GENE_B"],
    "sample_ids": ["s1", "s2", "s3", "s4", "s5", "s6"],
    "values": [[10, 11, 9, 80, 85, 75], [5, 6, 5, 5, 6, 5]]
  },
  "config": {
    "normalization": "cpm",
    "test": "welch-ttest",
    "correction": "benjamini-hochberg",
    "seed": 42,
    "timeout_seconds": 60
  }
}
```

`control`/`treatment` may be omitted when exactly two groups exist; they then
default to the alphabetically sorted group names (first = control). File
submission uses `metadata` (JSON string, same fields minus `matrix`) plus a
CSV whose first column is `gene` and remaining columns are sample ids.

### Task states

`pending → running → succeeded | failed | cancelled`

Validation errors are returned **before** any task or partial result is
created. Errors use a stable envelope:

```json
{"error": {"code": "VALIDATION_ERROR",
           "message": "...",
           "details": {"reason": "not_enough_samples", "group": "control"}}}
```

Error codes include `INVALID_REQUEST` (schema), `VALIDATION_ERROR`
(semantic, with a machine `reason`), `TASK_NOT_FOUND`, `TASK_STATE_ERROR`,
`RESULT_NOT_FOUND`, `INVALID_FILTER`, `QUEUE_FULL`, `INVALID_MATRIX_FILE`.

## 3. Statistical semantics

### Normalization (`config.normalization`)

| Value | Method |
|---|---|
| `none` | Values used as submitted |
| `cpm` | Counts per million: each sample scaled to sum 1,000,000 |
| `tpm` | Transcripts per million: divide by `gene_lengths` (required, positive), then scale to 1,000,000 |
| `quantile` | Quantile normalization: each column is replaced by the per-rank mean across samples; tied values receive averaged rank means |
| `tmm` | edgeR-style trimmed mean of M values: upper-quartile reference sample, 30% M-trim / 5% A-trim, asymptotic variance weights, factors rescaled to geometric mean 1, then CPM-style scaling |

### Tests (`config.test`), treatment vs control

| Value | Statistic reported | Notes |
|---|---|---|
| `welch-ttest` | t | Welch–Satterthwaite degrees of freedom; default |
| `student-ttest` | t | Pooled variance |
| `paired-ttest` | t | Requires equal group sizes; tests mean of paired differences |
| `mann-whitney` | signed z | Wilcoxon rank-sum, tie-corrected normal approximation with continuity correction |

p-values are two-sided. The t-distribution is evaluated through the regularized
incomplete beta function (Lentz continued fraction, 1e-14 tolerance) and
matches SciPy to ~1e-10 precision. Zero-variance cases (both groups constant)
return a signed infinite statistic in memory, serialized as the finite JSON
sentinel **±1e308** (paired with `p = 0`) when the groups differ, or
`statistic = 0, p = 1` when they coincide; no `NaN`/`Infinity` tokens ever
appear in an API response or persisted file.

### Effect size

`log2fc = log2((mean_treatment + eps)/(mean_control + eps))` computed on the
**normalized** group means. `eps = 1e-300` only keeps the value finite for an
all-zero group; for ordinary positive expression it is the plain mean ratio.

### Multiple testing correction (`config.correction`)

* `bonferroni` — `min(1, m·p)`
* `holm` — Holm step-down FWER
* `benjamini-hochberg` — BH FDR (default)
* `benjamini-yekutieli` — BY FDR (dependence-robust, multiplies by the m-th
  harmonic number)

All procedures enforce monotonicity and break ties by original gene order;
they were validated against `statsmodels.stats.multitest`.

### Ranking

Rows are ranked by `(adjusted p-value ascending, |log2fc| descending,
|statistic| descending, gene_id ascending)`, giving a total, reproducible
order even with ties.

The chosen methods, seed, pseudocount and group sizes are echoed back in
`strategy` of every task/result for full traceability.

## 4. Filtering, comparison and export

`GET /tasks/{id}/result` accepts: `gene` (exact), `gene_contains`,
`padj_max`, `pvalue_max`, `log2fc_min_abs`, `direction=up|down`, `top`.
Filtering never reorders rows — rank order is preserved, and results are
immutable once completed, so reads stay consistent while other tasks update
the store.

`POST /compare` body `{"task_ids": ["...", "..."], "padj_max": 0.05}` joins
per-gene statistics across tasks in deterministic order.

Export supports `csv` (header `rank,gene_id,mean_control,mean_treatment,
log2fc,statistic,p_value,p_value_adjusted`) and `json`.

## 5. Concurrency, cancellation and timeouts

* A bounded thread pool (`DEA_MAX_CONCURRENCY`) runs tasks; a slow task only
  occupies its own worker.
* Cancellation is cooperative: a task checks its cancel flag between genes.
  Cancellation (manual or by timeout watchdog) marks the task `cancelled` and
  **never writes a result file**.
* The commit point is *write result file (atomic temp + fsync + rename), then
  flip status to `succeeded`*. A succeeded task always has a readable result
  and is never downgraded by a late timer or exception.

## 6. Persistence and crash recovery

```
data/
  tasks/<id>.json      # task record (atomic writes)
  results/<id>.json    # full result
  requests/<id>.json   # original request snapshot
```

On startup the store scans task records and reconciles unfinished work:

* `pending`/`running` tasks are re-queued (`recover`) or marked `cancelled`
  with `INTERRUPTED` (`cancel`);
* a `succeeded` record whose result file is missing (crash inside the commit
  window) is re-queued;
* completed results are read back unchanged and never duplicated.

## 7. Reproducibility guarantees

* Pure-Python arithmetic with `math.fsum`; deterministic tie-breaking in
  ranking, quantile normalization and M/A trimming;
* the seed is part of the configuration (used by any future randomized
  method) and is recorded in the strategy metadata;
* timestamps are confined to task bookkeeping and never enter the numerical
  result;
* the same payload submitted twice (or replayed after a restart) yields the
  same rows field-for-field.

## 8. Tests

```bash
uv run pytest
```

The suite (67+ tests) covers: illegal inputs (sample counts, group
categories, dimensions, missing/non-finite values), concurrency and resource
limits, cancellation and timeout leaving no intermediate results,
restart/recovery policies, strategy switching cross-checked against
SciPy/statsmodels, and end-to-end bit-level reproducibility.
