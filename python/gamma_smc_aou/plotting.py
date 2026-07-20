from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .calibration import calibration_metrics


def plot_scan(results: pd.DataFrame, null: np.ndarray, output_dir: str | Path) -> dict[str, float]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    x = results["position_1based"]
    axes[0, 0].plot(x, results["mean_p_tmrca_lt_threshold"], lw=0.7)
    axes[0, 0].plot(x, results["null_mean"], lw=0.7, color="grey", alpha=0.8)
    axes[0, 0].set(xlabel="Segregating-site position", ylabel="P(T < threshold)", title="Observed and null mean")
    axes[0, 1].scatter(x, -np.log10(np.maximum(results["mc_p_upper"], np.finfo(float).tiny)), s=7)
    axes[0, 1].axhline(-np.log10(0.05), color="firebrick", ls="--", lw=0.8)
    axes[0, 1].set(xlabel="Segregating-site position", ylabel="-log10 Monte Carlo p", title="Selection scan")
    lead = int(np.nanargmin(results["mc_p_upper"])) if len(results) else None
    if lead is not None:
        axes[1, 0].hist(null[lead, np.isfinite(null[lead])], bins=30, alpha=0.8)
        axes[1, 0].axvline(results.iloc[lead]["mean_p_tmrca_lt_threshold"], color="firebrick")
    axes[1, 0].set(xlabel="Null P(T < threshold)", ylabel="Replicates", title="Lead-site null")
    p = np.sort(results["mc_p_upper"].dropna().to_numpy())
    expected = (np.arange(len(p)) + 0.5) / max(len(p), 1)
    axes[1, 1].scatter(expected, p, s=7)
    axes[1, 1].plot([0, 1], [0, 1], color="grey", lw=0.8)
    axes[1, 1].set(xlabel="Expected uniform quantile", ylabel="Observed p", title="Calibration QQ")
    fig.savefig(output_dir / "scan_verification.png", dpi=160)
    plt.close(fig)
    metrics = calibration_metrics(p)
    with (output_dir / "calibration_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    return metrics
