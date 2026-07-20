from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import kstest


STAT_COLUMN = "mean_p_tmrca_lt_threshold"


def monte_carlo_pvalue(observed: float, null: Sequence[float]) -> float:
    """One-sided upper-tail Monte Carlo p-value with the +1 correction."""
    values = np.asarray(null, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or not np.isfinite(observed):
        return np.nan
    return float((1 + np.count_nonzero(values >= observed)) / (values.size + 1))


def randomized_rank_pvalue(observed: float, null: Sequence[float], uniform: float) -> float:
    """Randomized upper-tail rank used only to diagnose discreteness/ties."""
    values = np.asarray(null, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or not np.isfinite(observed):
        return np.nan
    greater = np.count_nonzero(values > observed)
    tied_rank_slots = np.count_nonzero(values == observed) + 1  # include observed
    return float((greater + uniform * tied_rank_slots) / (values.size + 1))


def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    if not vals.size:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    adjusted = ranked * vals.size / np.arange(1, vals.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[ok] = restored
    return out


def _relative_position(frame: pd.DataFrame, sequence_length: float) -> np.ndarray:
    position = frame["position_0based"].to_numpy(dtype=float)
    return np.clip(position / sequence_length, 0.0, 1.0)


def calibrate_sites(
    observed: pd.DataFrame,
    simulations: Sequence[pd.DataFrame],
    *,
    observed_length: float,
    simulation_lengths: Sequence[float] | None = None,
    match: str = "nearest",
    stat_column: str = STAT_COLUMN,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Match one null statistic per replicate to every observed segregating site.

    ``nearest`` matches relative position within the region. ``region_max`` uses
    the maximum statistic in each replicate and therefore controls the per-region
    scan more conservatively. Replicate sites are never pooled as independent.
    """
    if match not in {"nearest", "region_max"}:
        raise ValueError("match must be 'nearest' or 'region_max'")
    if simulation_lengths is None:
        simulation_lengths = [observed_length] * len(simulations)
    if len(simulation_lengths) != len(simulations):
        raise ValueError("simulation_lengths must match simulations")

    obs_rel = _relative_position(observed, observed_length)
    null = np.full((len(observed), len(simulations)), np.nan)
    for replicate, (sim, length) in enumerate(zip(simulations, simulation_lengths)):
        if sim.empty:
            continue
        sim_stats = sim[stat_column].to_numpy(dtype=float)
        if match == "region_max":
            null[:, replicate] = np.nanmax(sim_stats)
            continue
        sim_rel = _relative_position(sim, length)
        order = np.argsort(sim_rel)
        sim_rel, sim_stats = sim_rel[order], sim_stats[order]
        insertion = np.searchsorted(sim_rel, obs_rel)
        right = np.clip(insertion, 0, len(sim_rel) - 1)
        left = np.clip(insertion - 1, 0, len(sim_rel) - 1)
        choose_right = np.abs(sim_rel[right] - obs_rel) < np.abs(sim_rel[left] - obs_rel)
        nearest = np.where(choose_right, right, left)
        null[:, replicate] = sim_stats[nearest]

    result = observed.copy()
    result["mc_p_upper"] = [
        monte_carlo_pvalue(obs, row)
        for obs, row in zip(result[stat_column].to_numpy(dtype=float), null)
    ]
    result["bh_q"] = bh_fdr(result["mc_p_upper"])
    result["null_mean"] = np.nanmean(null, axis=1)
    result["null_sd"] = np.nanstd(null, axis=1, ddof=1)
    result["n_null"] = np.sum(np.isfinite(null), axis=1)
    return result, null


def calibration_metrics(pvalues: Sequence[float]) -> dict[str, float]:
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if not p.size:
        return {"n": 0, "mean_p": np.nan, "ks_uniform_p": np.nan, "fraction_p_lt_0_05": np.nan}
    return {
        "n": int(p.size),
        "mean_p": float(np.mean(p)),
        "ks_uniform_p": float(kstest(p, "uniform").pvalue),
        "fraction_p_lt_0_05": float(np.mean(p < 0.05)),
    }
