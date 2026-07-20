from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .calibration import calibration_metrics


def plot_truth_tmrca_relationship(
    summaries: list[Path],
    *,
    sequence_length: float,
    effective_size: float,
    threshold_generations: float,
    output_dir: str | Path,
    relative_position: float = 0.5,
) -> dict[str, float]:
    """Plot segregating-site truth against an unascertained random-locus reference."""
    rows = []
    for replicate, path in enumerate(summaries):
        frame = pd.read_csv(path, sep="\t")
        if frame.empty:
            continue
        relative = frame["position_0based"].to_numpy(dtype=float) / sequence_length
        row = frame.iloc[int(np.argmin(np.abs(relative - relative_position)))]
        rows.append({
            "replicate": replicate,
            "position_0based": int(row["position_0based"]),
            "mean_tmrca_generations": float(row["mean_tmrca_generations"]),
            "fraction_tmrca_lt_threshold": float(row["mean_p_tmrca_lt_threshold"]),
        })
    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError("no nonempty truth summaries")
    expected_mean = 2.0 * effective_size
    expected_recent = 1.0 - np.exp(-threshold_generations / expected_mean)
    rho, rho_p = spearmanr(
        table["mean_tmrca_generations"], table["fraction_tmrca_lt_threshold"]
    )
    metrics = {
        "n_replicates": int(len(table)),
        "random_locus_expected_mean_tmrca_generations": float(expected_mean),
        "simulated_mean_tmrca_generations": float(table["mean_tmrca_generations"].mean()),
        "random_locus_expected_fraction_recent": float(expected_recent),
        "simulated_mean_fraction_recent": float(table["fraction_tmrca_lt_threshold"].mean()),
        "spearman_rho": float(rho),
        "spearman_p": float(rho_p),
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "truth_tmrca_relationship.tsv", sep="\t", index=False)
    with (output_dir / "truth_tmrca_relationship.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].scatter(
        table["mean_tmrca_generations"], table["fraction_tmrca_lt_threshold"],
        s=22, alpha=0.65, label="simulated replicates",
    )
    x_min = max(1.0, float(table["mean_tmrca_generations"].min()) * 0.9)
    x_max = float(table["mean_tmrca_generations"].max()) * 1.1
    x_curve = np.linspace(x_min, x_max, 300)
    axes[0].plot(
        x_curve, 1.0 - np.exp(-threshold_generations / x_curve),
        color="grey", lw=1.2, label="random-locus exponential reference",
    )
    axes[0].scatter([expected_mean], [expected_recent], marker="*", s=130, color="firebrick", label="random-locus 2Ne reference")
    axes[0].set(
        xlabel="Simulated mean within-pair TMRCA (generations)",
        ylabel=f"True fraction TMRCA < {threshold_generations:g} generations",
        title="Truth versus simulated TMRCA",
    )
    axes[0].legend(fontsize=8)
    axes[1].hist(table["mean_tmrca_generations"], bins=20, alpha=0.8)
    axes[1].axvline(expected_mean, color="firebrick", lw=1.2, label="random-locus 2Ne")
    axes[1].axvline(table["mean_tmrca_generations"].mean(), color="black", ls="--", lw=1.0, label="simulation mean")
    axes[1].set(
        xlabel="Mean within-pair TMRCA (generations)", ylabel="Replicates",
        title="2,000-diploid midpoint truth",
    )
    axes[1].legend(fontsize=8)
    fig.savefig(output_dir / "truth_vs_simulated_tmrca.png", dpi=180)
    plt.close(fig)
    return metrics


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
