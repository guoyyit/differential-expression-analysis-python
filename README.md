# Differential Expression Analysis Service

本地转录组差异表达分析服务。客户端通过 HTTP 提交表达量矩阵与分组设计，
服务完成归一化、两组差异检验与多重检验校正，返回按显著性与效应量排序的
结果表。全部计算在本进程内完成，仅依赖本地 SQLite 持久化，不访问外部服务。

## 快速开始

```bash
uv sync
uv run uvicorn main:app --host 127.0.0.1 --port 8000
```

提交一个分析（也可把 JSON 存入文件后用 `curl --data @file.json`）：

```bash
curl -s -X POST http://127.0.0.1:8000/api/v1/analyses \
  -H 'content-type: application/json' \
  --data @example_request.json
```

`example_request.json`：

```json
{
  "genes": ["g0", "g1"],
  "samples": ["s1", "s2", "s3", "s4"],
  "matrix": [[10, 12, 60, 62], [20, 21, 20, 22]],
  "groups": [
    {"name": "ctrl", "samples": ["s1", "s2"]},
    {"name": "trt",  "samples": ["s3", "s4"]}
  ],
  "config": {"normalization": "log_cpm", "test": "welch", "correction": "bh"},
  "task_id": "my-analysis-1"
}
```

## HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/api/v1/health` | 健康检查 |
| POST | `/api/v1/analyses` | 提交分析，返回 `202` 与任务快照 |
| GET | `/api/v1/analyses?status=` | 列出任务（可按状态过滤） |
| GET | `/api/v1/analyses/{id}` | 任务状态、配置回显、错误信息、汇总 |
| POST | `/api/v1/analyses/{id}/cancel` | 取消（终态任务上为幂等空操作） |
| GET | `/api/v1/analyses/{id}/results` | 过滤/排序/分页查询结果 |
| GET | `/api/v1/analyses/{id}/results/export?format=csv\|json` | 导出结果 |
| POST | `/api/v1/compare?task_a=&task_b=` | 跨任务结果对比 |

结果查询参数：`gene`（可重复，精确匹配）、`significant_only`、
`max_adjusted_p`、`min_abs_effect_size`、`sort`（`rank` /
`adjusted_p_value` / `p_value` / `effect_size` / `abs_effect_size` /
`statistic` / `gene`）、`order`、`offset`、`limit`。

错误统一为结构化 JSON：`{"error": {"code", "message", "details"}}`，
状态码：`422` 输入非法、`404` 任务不存在、`409` 冲突/结果未就绪、
`429` 队列已满。

## 统计语义

- **归一化**（`config.normalization`）
  - `none`：原始计数；
  - `cpm`：每样本按文库大小缩放为 counts-per-million；
  - `log_cpm`（默认）：`log2(CPM + 1)`；
  - `median_ratio`：DESeq 风格 median-of-ratios 尺寸因子校正后取
    `log2(x + 1)`；无全阳性基因时回退为文库大小因子（行为确定）。
- **预过滤**：`config.min_total_count`（默认 1.0）丢弃总原始计数低于阈值
  的基因，被丢弃基因不参与检验与校正。
- **差异检验**（`config.test`），零假设为两组分布相同：
  - `welch`（默认）：Welch 双样本 t 检验（不要求方差齐性），统计量为 t，
    正号表示 B 组高于 A 组；常数组等退化情形返回定义好的值（均值相等
    时 `t=0, p=1`），绝不产生 NaN；
  - `mann_whitney`：秩和检验，正态近似 + 并列校正 + 连续性校正，
    统计量为 z。
- **多重检验校正**（`config.correction`）：`bh`（Benjamini-Hochberg FDR，
  默认）、`bonferroni`、`none`。
- **效应量**：`effect_size = log2(mean_B + 1) - log2(mean_A + 1)`
  （基于归一化后均值，伪计数保证有限）；`cohens_d` 为合并标准差标准化
  的均值差， pooled SD 为 0 或样本量不足时为 `null`。
- **显著性**：`adjusted_p_value <= config.alpha`（默认 0.05）。
- **排序**：`rank` 按 `(校正后 p, 原始 p, -|效应量|, 基因名)` 全序排列，
  任何并列都有确定次序。

## 任务语义

- 状态机：`pending → running → completed | failed | cancelled | timed_out
  | interrupted`。终态不可逆：已完成任务永远不会被改写为失败，取消/超时
  后迟到的计算结果会被丢弃，绝不留下中间结果。
- **幂等**：客户端提供 `task_id` 时，相同请求体重提交返回原任务；不同
  请求体复用同一 `task_id` 返回 `409`。
- **并发与资源**：固定大小工作池（`max_workers`）+ 有界等待队列
  （`max_queue_size`），超限返回 `429`；慢任务只占用自己的 worker。
- **超时**：`config.timeout_seconds`（上限 `max_timeout_seconds`），
  到期任务标记 `timed_out` 并协作式中断计算。
- **持久化与恢复**：任务与结果存于 SQLite（WAL）。重启时
  `recovery_policy=resume` 重跑未完成任务（结果不丢失不重复），
  `fail` 则标记为 `interrupted`。
- **可复现性**：管线不使用随机数；每个基因的统计量独立计算，校正与
  排序为确定性全序。相同输入（无论并发调度如何）产生逐条一致的结果。

## 配置项

服务级配置通过环境变量（见 `main.py`）或 `Settings` 注入：

| 环境变量 | 默认 | 说明 |
| --- | --- | --- |
| `DEA_DB_PATH` | `dea_service.db` | SQLite 文件路径 |
| `DEA_MAX_WORKERS` | `2` | 并行执行的分析数上限 |
| `DEA_MAX_QUEUE_SIZE` | `16` | 排队任务上限（超出拒绝） |
| `DEA_DEFAULT_TIMEOUT_SECONDS` | `300` | 任务默认超时 |
| `DEA_MAX_TIMEOUT_SECONDS` | `3600` | 任务超时上限 |
| `DEA_RECOVERY_POLICY` | `resume` | `resume` / `fail` |

分析级配置见 `AnalysisConfig`（`normalization`、`test`、`correction`、
`alpha`、`min_total_count`、`timeout_seconds`、`seed`）。每次任务的配置
随快照回显回显，策略选择完全可追溯。

## 测试

```bash
uv run pytest
```

覆盖：非法输入（矩阵维度/缺失值/分组设计）、统计正确性（对照 scipy）、
策略切换、并发与队列背压、取消与超时、崩溃恢复（resume/fail）、
幂等提交、过滤排序分页导出、跨任务对比与结果可复现性。
