from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STAT = "mean_p_tmrca_lt_threshold"


def align_nearest(truth: pd.DataFrame, decoded: pd.DataFrame) -> pd.DataFrame:
    if truth.empty or decoded.empty:
        return pd.DataFrame()
    tx = truth["position_0based"].to_numpy(dtype=float)
    dx = decoded["position_0based"].to_numpy(dtype=float)
    order = np.argsort(tx)
    tx = tx[order]
    insertion = np.searchsorted(tx, dx)
    right = np.clip(insertion, 0, len(tx) - 1)
    left = np.clip(insertion - 1, 0, len(tx) - 1)
    nearest = np.where(np.abs(tx[right] - dx) < np.abs(tx[left] - dx), right, left)
    truth_rows = truth.iloc[order[nearest]].reset_index(drop=True)
    return pd.DataFrame({
        "position_0based": decoded["position_0based"].to_numpy(),
        "truth": truth_rows[STAT].to_numpy(dtype=float),
        "decoded": decoded[STAT].to_numpy(dtype=float),
        "truth_mean_tmrca_generations": truth_rows["mean_tmrca_generations"].to_numpy(dtype=float),
        "decoded_mean_tmrca_generations": decoded["mean_tmrca_generations"].to_numpy(dtype=float),
    })


def error_metrics(aligned: pd.DataFrame) -> dict[str, float]:
    ok = np.isfinite(aligned["truth"]) & np.isfinite(aligned["decoded"])
    truth = aligned.loc[ok, "truth"].to_numpy()
    predicted = aligned.loc[ok, "decoded"].to_numpy()
    error = predicted - truth
    return {
        "n": int(len(error)),
        "mae_recent_probability": float(np.mean(np.abs(error))) if len(error) else np.nan,
        "rmse_recent_probability": float(np.sqrt(np.mean(error ** 2))) if len(error) else np.nan,
        "bias_recent_probability": float(np.mean(error)) if len(error) else np.nan,
        "pearson_recent_probability": float(np.corrcoef(truth, predicted)[0, 1]) if len(error) > 1 else np.nan,
    }


def evaluate_pairs(pairs: list[tuple[Path, Path]], output_dir: str | Path) -> dict[str, float]:
    frames = []
    for replicate, (truth_path, decoded_path) in enumerate(pairs):
        aligned = align_nearest(
            pd.read_csv(truth_path, sep="\t"), pd.read_csv(decoded_path, sep="\t")
        )
        aligned["replicate"] = replicate
        frames.append(aligned)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["truth", "decoded"])
    metrics = error_metrics(combined)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_dir / "decoder_errors.tsv", sep="\t", index=False)
    with (output_dir / "decoder_error_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    if not combined.empty:
        axes[0].scatter(combined["truth"], combined["decoded"], s=5, alpha=0.25)
        axes[0].plot([0, 1], [0, 1], color="grey")
        residual = combined["decoded"] - combined["truth"]
        axes[1].hist(residual, bins=50)
        first = combined[combined["replicate"] == combined["replicate"].min()]
        axes[2].plot(first["position_0based"], first["truth"], label="true", lw=0.8)
        axes[2].plot(first["position_0based"], first["decoded"], label="Gamma-SMC", lw=0.8)
        axes[2].legend()
    axes[0].set(xlabel="True fraction recent", ylabel="Decoded posterior fraction", title="Decoder accuracy")
    axes[1].set(xlabel="Decoded - true", ylabel="Sites", title="Residual loss")
    axes[2].set(xlabel="Segregating-site position", ylabel="Fraction recent", title="First-coalescence branch")
    fig.savefig(output_dir / "decoder_verification.png", dpi=160)
    plt.close(fig)
    return metrics
