#!/usr/bin/env python3
"""Audit and plot a completed high-frequency rejection screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


THRESHOLDS = (0.0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate(
    attempts: pd.DataFrame,
    contract: dict,
    *,
    expected_attempts: int,
) -> None:
    required = {
        "attempt",
        "seed",
        "accepted",
        "focal_allele_outcome",
        "population_allele_frequency",
    }
    missing = required.difference(attempts.columns)
    if missing:
        raise ValueError(f"attempt table is missing columns: {sorted(missing)}")
    attempts.sort_values("attempt", inplace=True)
    observed = attempts["attempt"].to_numpy(dtype=int)
    if not np.array_equal(observed, np.arange(expected_attempts, dtype=int)):
        raise ValueError("attempt IDs are not contiguous over the expected range")
    expected_seeds = int(contract["seed"]) + observed * 10
    if not np.array_equal(
        attempts["seed"].to_numpy(dtype=int),
        expected_seeds,
    ):
        raise ValueError("attempt seeds do not match the recorded contract")
    af = attempts["population_allele_frequency"].to_numpy(dtype=float)
    if np.any((af < 0) | (af > 1)):
        raise ValueError("allele frequencies must be in [0, 1]")
    expected_accepted = (
        af >= float(contract["minimum_population_af"])
    ) & attempts["focal_allele_outcome"].eq("segregating").to_numpy()
    observed_accepted = (
        attempts["accepted"].astype(str).str.lower().eq("true").to_numpy()
    )
    if not np.array_equal(observed_accepted, expected_accepted):
        raise ValueError("saved acceptance decisions do not match the contract")


def _plot(
    attempts: pd.DataFrame,
    contract: dict,
    tail: pd.DataFrame,
    output: Path,
) -> None:
    af = attempts["population_allele_frequency"].to_numpy(dtype=float)
    positive = af[af > 0]
    target = float(contract["minimum_population_af"])
    figure, axes = plt.subplots(1, 2, figsize=(15, 6.2))

    lower = min(float(positive.min()), 1e-4)
    bins = np.geomspace(lower, max(target * 1.05, positive.max() * 1.05), 34)
    axes[0].hist(
        positive,
        bins=bins,
        color="#3A6EA5",
        edgecolor="white",
        linewidth=0.5,
    )
    axes[0].axvline(
        positive.max(),
        color="#D17C0B",
        linewidth=2,
        label=f"maximum = {positive.max():.4%}",
    )
    axes[0].axvline(
        target,
        color="#A23E48",
        linestyle="--",
        linewidth=2.2,
        label=f"target = {target:.0%}",
    )
    axes[0].set(
        xscale="log",
        xlabel="Present-day selected-allele frequency",
        ylabel="Exact SLiM trajectories",
        title=f"Non-lost trajectories ({len(positive):,}/{len(af):,})",
    )
    axes[0].legend(frameon=False)

    plotting_probability = (
        tail["n_exceeding"].to_numpy(dtype=float) + 1
    ) / (len(attempts) + 1)
    nonzero = tail["n_exceeding"].to_numpy(dtype=int) > 0
    threshold_percent = 100 * tail["threshold"].to_numpy(dtype=float)
    axes[1].plot(
        threshold_percent[nonzero],
        plotting_probability[nonzero],
        color="#3A6EA5",
        marker="o",
        linewidth=2,
        label="observed exceedance",
    )
    axes[1].scatter(
        threshold_percent[~nonzero],
        plotting_probability[~nonzero],
        facecolors="none",
        edgecolors="#A23E48",
        marker="v",
        s=90,
        linewidth=2,
        label="zero observed; 1/(n+1) plotting position",
    )
    axes[1].axvline(
        100 * target,
        color="#A23E48",
        linestyle="--",
        linewidth=2.2,
    )
    for row, probability in zip(
        tail.itertuples(index=False),
        plotting_probability,
        strict=True,
    ):
        if row.n_exceeding == 0:
            axes[1].annotate(
                "0",
                (100 * row.threshold, probability),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                color="#A23E48",
            )
    axes[1].set(
        yscale="log",
        xlabel="AF threshold (%)",
        ylabel="Fraction of all trajectories",
        title="Present-day AF exceedance tail",
    )
    axes[1].legend(frameon=False)
    figure.suptitle(
        "s=0.01 single-copy sweep, 180 generations: AF>10% feasibility",
        fontsize=17,
        y=1.01,
    )
    figure.tight_layout()
    figure.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--attempts", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-attempts", type=int, default=20_000)
    args = parser.parse_args()
    if args.expected_attempts < 1:
        parser.error("--expected-attempts must be positive")

    attempts_path = args.attempts.resolve()
    contract_path = args.contract.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    attempts = pd.read_csv(attempts_path, sep="\t")
    with contract_path.open(encoding="utf-8") as handle:
        contract = json.load(handle)
    _validate(attempts, contract, expected_attempts=args.expected_attempts)

    af = attempts["population_allele_frequency"].to_numpy(dtype=float)
    tail = pd.DataFrame(
        {
            "threshold": THRESHOLDS,
            "n_exceeding": [int(np.count_nonzero(af > x)) for x in THRESHOLDS],
        }
    )
    tail["fraction_exceeding"] = tail["n_exceeding"] / len(attempts)
    tail.to_csv(output_dir / "af_exceedance_counts.tsv", sep="\t", index=False)
    attempts.nlargest(20, "population_allele_frequency").to_csv(
        output_dir / "top_20_trajectories.tsv",
        sep="\t",
        index=False,
    )
    _plot(
        attempts,
        contract,
        tail,
        output_dir / "s0p01_af10_feasibility.png",
    )

    positive = af[af > 0]
    target = float(contract["minimum_population_af"])
    n_at_or_above = int(np.count_nonzero(af >= target))
    n_strictly_above = int(np.count_nonzero(af > target))
    summary = {
        "result": (
            "qualifying_trajectory_found"
            if n_at_or_above
            else "no_qualifying_trajectory"
        ),
        "exact_slim_attempts": int(len(attempts)),
        "selection_coefficient": float(contract["selection_coefficient"]),
        "variant_age_generations": int(contract["variant_age_generations"]),
        "minimum_population_af": target,
        "n_at_or_above_target": n_at_or_above,
        "n_strictly_above_target": n_strictly_above,
        "n_lost": int(np.count_nonzero(af == 0)),
        "n_segregating": int(np.count_nonzero(af > 0)),
        "survival_fraction": float(np.mean(af > 0)),
        "maximum_population_af": float(af.max()),
        "maximum_attempt": int(attempts.iloc[int(np.argmax(af))]["attempt"]),
        "maximum_seed": int(attempts.iloc[int(np.argmax(af))]["seed"]),
        "positive_af_median": float(np.median(positive)),
        "positive_af_mean": float(np.mean(positive)),
        "positive_af_q95": float(np.quantile(positive, 0.95)),
        "zero_event_one_sided_95_upper_probability": float(
            1 - 0.05 ** (1 / len(attempts))
        ),
        "decoder_run": False,
        "decoder_reason": (
            "No selected genealogy met the requested present-day AF "
            "conditioning."
        ),
        "input_sha256": {
            attempts_path.name: _sha256(attempts_path),
            contract_path.name: _sha256(contract_path),
        },
        "output_files": [
            "af_exceedance_counts.tsv",
            "top_20_trajectories.tsv",
            "s0p01_af10_feasibility.png",
        ],
    }
    with (output_dir / "screen_summary.json").open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
