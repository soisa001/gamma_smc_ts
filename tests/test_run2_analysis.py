"""Statistic, null and p-value tests for the run2 study.

These run on synthetic tables so the arithmetic is checked independently of any
simulation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gamma_smc_aou import run2_analysis as analysis
from gamma_smc_aou.run2_config import CHB_ARM, SIGNIFICANCE_LEVEL


def _statistics(selected: list[float], neutral: list[float], threshold: float = 4500.0):
    """Build a minimal observed-statistics frame for one threshold."""
    rows = []
    for mode, values in (("selected", selected), ("neutral", neutral)):
        for index, value in enumerate(values):
            for genotype_class in ("overall", "hom_carrier"):
                rows.append(
                    {
                        "arm_id": CHB_ARM.arm_id,
                        "mode": mode,
                        "unit_id": f"{mode}-{index}",
                        "replicate_index": index,
                        "genotype_class": genotype_class,
                        "n_pairs": 100,
                        "threshold_years": threshold,
                        "threshold_generations": threshold / 25.0,
                        "p_tmrca_lt_threshold": value
                        if genotype_class == "overall"
                        else value * 0.5,
                    }
                )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# p-value arithmetic
# ---------------------------------------------------------------------------


def test_upper_tail_pvalue_counts_ties_as_exceedances():
    null_values = np.array([0.0, 0.1, 0.2, 0.3])
    assert analysis.upper_tail_pvalue(0.2, null_values) == pytest.approx(3 / 5)
    assert analysis.upper_tail_pvalue(0.31, null_values) == pytest.approx(1 / 5)


def test_upper_tail_pvalue_never_reaches_zero():
    null_values = np.zeros(100)
    assert analysis.upper_tail_pvalue(1.0, null_values) == pytest.approx(1 / 101)


def test_upper_tail_pvalue_is_nan_without_a_null():
    assert np.isnan(analysis.upper_tail_pvalue(0.5, np.array([])))
    assert np.isnan(analysis.upper_tail_pvalue(np.nan, np.array([0.1])))


# ---------------------------------------------------------------------------
# Table construction
# ---------------------------------------------------------------------------


def test_pvalues_use_the_overall_class_only():
    statistics = _statistics(selected=[0.9] * 5, neutral=[0.0] * 100)
    table = analysis.pvalues(statistics, CHB_ARM)
    assert len(table) == 5
    assert set(table["n_null"]) == {100}
    assert table["p_value"].unique() == pytest.approx([1 / 101])
    assert table["significant"].all()


def test_pvalues_do_not_fire_when_the_signal_is_absent():
    rng = np.random.default_rng(0)
    neutral = list(rng.uniform(0, 0.2, size=100))
    selected = list(rng.uniform(0, 0.2, size=50))
    table = analysis.pvalues(_statistics(selected, neutral), CHB_ARM)
    assert 0.0 <= table["significant"].mean() <= 0.25


def test_power_summary_reports_the_detection_rate():
    statistics = _statistics(selected=[0.9, 0.9, 0.0, 0.0], neutral=[0.0] * 100)
    table = analysis.pvalues(statistics, CHB_ARM)
    power = analysis.power_summary(table, CHB_ARM)
    assert len(power) == 1
    record = power.iloc[0]
    assert record["n_selected"] == 4
    assert record["n_significant"] == 2
    assert record["power"] == pytest.approx(0.5)
    assert record["p_value_floor"] == pytest.approx(1 / 101)
    assert record["significance_level"] == SIGNIFICANCE_LEVEL


def test_null_distribution_reports_quantiles():
    neutral = list(np.linspace(0.0, 1.0, 101))
    table = analysis.null_distribution(_statistics([0.5], neutral), CHB_ARM)
    record = table.iloc[0]
    assert record["n_null"] == 101
    assert record["min"] == pytest.approx(0.0)
    assert record["max"] == pytest.approx(1.0)
    assert record["q50"] == pytest.approx(0.5)
    assert record["q95"] == pytest.approx(0.95)


def test_null_calibration_is_near_the_nominal_level():
    rng = np.random.default_rng(7)
    neutral = list(rng.normal(0.1, 0.02, size=200))
    table = analysis.null_calibration(_statistics([0.5], neutral), CHB_ARM)
    record = table.iloc[0]
    assert record["n_null"] == 200
    # Leave-one-out p-values of an exchangeable null are uniform, so the false
    # positive rate must sit near alpha rather than at 0 or 1.
    assert 0.0 <= record["false_positive_rate"] <= 0.12
    assert record["median_p_value"] == pytest.approx(0.5, abs=0.1)


def test_analyse_arm_requires_replicates(tmp_path):
    with pytest.raises(RuntimeError, match="no completed replicates"):
        analysis.analyse_arm(tmp_path, CHB_ARM)
