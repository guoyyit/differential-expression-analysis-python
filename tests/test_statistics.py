"""Strategy correctness, switching and determinism at the computation layer."""
from __future__ import annotations

import math

import pytest

from deapp import multipletest, normalization, statistics
from deapp.engine import run_analysis
from deapp.models import AnalysisConfig, CorrectionMethod, NormalizationMethod, TestMethod
from deapp.statistics import gene_stat
from deapp.validation import ValidatedInput


def _validated(c=(1.0, 2, 3, 4), t=(10.0, 12, 11, 13)):
    return ValidatedInput(
        gene_ids=("g1",),
        sample_ids=tuple(f"c{i}" for i in c) + tuple(f"t{i}" for i in t),
        control_group="c",
        treatment_group="t",
        control_samples=tuple(f"c{i}" for i in c),
        treatment_samples=tuple(f"t{i}" for i in t),
        values=((tuple(c), tuple(t)),),
    )


def test_welch_significant():
    gs = gene_stat(TestMethod.WELCH, "g1", (1, 2, 3, 4), (10, 12, 11, 13))
    assert gs.log2fc > 2.0
    assert gs.p_value < 0.001
    assert gs.statistic > 0


def test_negative_effect_direction():
    gs = gene_stat(TestMethod.WELCH, "g1", (10, 12, 11, 13), (1, 2, 3, 4))
    assert gs.log2fc < -2.0
    assert gs.statistic < 0


def test_identical_groups_pvalue_one():
    gs = gene_stat(TestMethod.WELCH, "g1", (5, 6, 7), (5, 6, 7))
    assert gs.p_value == 1.0
    assert gs.statistic == 0.0


def test_student_and_welch_differ_under_unequal_variance():
    c = (1.0, 2, 1.5, 2.5)
    t = (10.0, 50, 8, 55)
    w = gene_stat(TestMethod.WELCH, "g", c, t)
    s = gene_stat(TestMethod.STUDENT, "g", c, t)
    assert w.p_value != s.p_value


def test_paired_matches_difference_design():
    gs = gene_stat(TestMethod.PAIRED, "g", (1, 2, 3), (4, 5, 6))
    assert gs.p_value < 0.05


def test_mann_whitney_handles_ties():
    gs = gene_stat(TestMethod.MANN_WHITNEY, "g", (1, 1, 2, 2), (3, 3, 4, 4))
    assert gs.p_value < 0.05


def test_p_values_bounded():
    for method in TestMethod:
        if method == TestMethod.PAIRED:
            c, t = (1.0, 2, 3), (1.0, 2, 3)
        else:
            c, t = (1.0, 1, 2), (1.0, 1, 2)
        gs = gene_stat(method, "g", c, t)
        assert 0.0 <= gs.p_value <= 1.0


def test_bh_monotone_by_rank():
    pvals = [0.001, 0.5, 0.02, 0.3, 0.0]
    adj = multipletest.benjamini_hochberg(pvals)
    # smallest p gets smallest adjusted value
    assert adj[4] <= adj[0] <= adj[2] <= adj[3] <= adj[1]
    assert all(0.0 <= a <= 1.0 for a in adj)


def test_correction_ordering_of_methods():
    pvals = [0.01, 0.02, 0.03, 0.7]
    bonf = multipletest.bonferroni(pvals)
    bh = multipletest.benjamini_hochberg(pvals)
    for a, b in zip(bonf, bh):
        assert a >= b - 1e-12


def test_normalizations_preserve_shape():
    data = _validated()
    for method in NormalizationMethod:
        cfg = AnalysisConfig(normalization=method,
                             gene_lengths={"g1": 100.0} if method == NormalizationMethod.TPM else None)
        out = normalization.normalize(method, data, cfg)
        assert len(out) == 1
        assert len(out[0][0]) == 4 and len(out[0][1]) == 4


def test_quantile_makes_distributions_equal():
    data = _validated(c=(1.0, 2, 3, 4), t=(100.0, 200, 300, 400))
    out = normalization.normalize(NormalizationMethod.QUANTILE, data, AnalysisConfig())
    all_cols = sorted(out[0][0]) + sorted(out[0][1])
    # after quantile normalization corresponding ranks are identical
    assert sorted(out[0][0]) == sorted(out[0][1])


def test_cpm_scales_to_library_size():
    data = ValidatedInput(
        gene_ids=("g1", "g2"),
        sample_ids=("c1", "c2", "t1", "t2"),
        control_group="c", treatment_group="t",
        control_samples=("c1", "c2"), treatment_samples=("t1", "t2"),
        values=(((1.0, 3.0), (2.0, 4.0)),   # g1: c=1,3  t=2,4
                ((9.0, 7.0), (8.0, 6.0))),  # g2
    )
    out = normalization.normalize(NormalizationMethod.CPM, data, AnalysisConfig())
    # out[gene] -> (control values, treatment values); each *sample* sums to 1e6
    sample_sums = [
        out[0][0][0] + out[1][0][0],  # c1
        out[0][0][1] + out[1][0][1],  # c2
        out[0][1][0] + out[1][1][0],  # t1
        out[0][1][1] + out[1][1][1],  # t2
    ]
    for total in sample_sums:
        assert math.isclose(total, 1_000_000.0, rel_tol=1e-12)


def test_engine_result_is_sorted_and_ranked():
    data = _validated()
    res = run_analysis("t", data, AnalysisConfig())
    assert [r.rank for r in res.rows] == [1]
    assert res.strategy.test == TestMethod.WELCH.value


def test_engine_strategy_traceability():
    data = _validated()
    cfg = AnalysisConfig(
        normalization=NormalizationMethod.CPM,
        test=TestMethod.STUDENT,
        correction=CorrectionMethod.BONFERRONI,
        seed=99,
    )
    res = run_analysis("t", data, cfg)
    assert res.strategy.normalization == "cpm"
    assert res.strategy.test == "student-ttest"
    assert res.strategy.correction == "bonferroni"
    assert res.strategy.seed == 99


def test_engine_bit_identical_across_calls():
    data = _validated(c=tuple(float(i) for i in range(2, 8)),
                      t=tuple(float(i) * 3 + 1 for i in range(2, 8)))
    cfg = AnalysisConfig(normalization=NormalizationMethod.TMM,
                         test=TestMethod.WELCH, correction=CorrectionMethod.BH,
                         seed=7)
    r1 = run_analysis("a", data, cfg)
    r2 = run_analysis("b", data, cfg)
    a = [r.model_dump() for r in r1.rows]
    b = [r.model_dump() for r in r2.rows]
    assert a == b
