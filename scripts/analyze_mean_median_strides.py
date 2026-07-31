#!/usr/bin/env python3
"""Compare posterior mean and median hard calls at 1 kb and 10 kb."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CALLED_COLUMN = "frac_recent_4500"
SOFT_COLUMNS = [
    "position_0based",
    "position_1based",
    "n_pairs",
    "mean_p_tmrca_lt_threshold",
    "mean_tmrca_generations",
    "mean_p_lt_4500",
]


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if np.std(left) == 0 or np.std(right) == 0:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def _assert_soft_identity(first: pd.DataFrame, second: pd.DataFrame) -> None:
    if len(first) != len(second):
        raise RuntimeError("paired call-rule profiles have different lengths")
    for column in SOFT_COLUMNS:
        if not np.allclose(
            first[column].to_numpy(),
            second[column].to_numpy(),
            rtol=0,
            atol=1e-12,
            equal_nan=True,
        ):
            raise RuntimeError(
                f"mean/median call choice changed soft column {column!r}"
            )


def _read_selected(result_dir: Path, call: str) -> pd.DataFrame:
    return pd.read_csv(
        result_dir / f"selected_decoded_{call}_call_profile.tsv",
        sep="\t",
    )


def _read_neutral(result_dir: Path, call: str) -> pd.DataFrame:
    return pd.read_csv(
        result_dir / f"neutral_decoded_{call}_call_profiles.tsv.gz",
        sep="\t",
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)


def _center(frame: pd.DataFrame, position: int = 5_000_000) -> pd.Series:
    rows = frame.loc[frame["position_0based"].eq(position)]
    if len(rows) != 1:
        raise RuntimeError(f"expected one selected center row; observed {len(rows)}")
    return rows.iloc[0]


def _shared_stride_metrics(
    stride_1kb: pd.DataFrame,
    stride_10kb: pd.DataFrame,
    *,
    keys: list[str],
    value: str,
) -> dict:
    downsampled = stride_1kb.loc[
        stride_1kb["position_0based"].mod(10_000).eq(0),
        [*keys, value],
    ]
    merged = downsampled.merge(
        stride_10kb[[*keys, value]],
        on=keys,
        how="inner",
        suffixes=("_1kb", "_10kb"),
        validate="one_to_one",
    )
    if len(merged) != len(stride_10kb):
        raise RuntimeError("1 kb and 10 kb profiles do not share the expected grid")
    left = merged[f"{value}_1kb"].to_numpy(dtype=float)
    right = merged[f"{value}_10kb"].to_numpy(dtype=float)
    difference = left - right
    return {
        "shared_rows": int(len(merged)),
        "pearson_correlation": _correlation(left, right),
        "mean_difference_1kb_minus_10kb": float(np.mean(difference)),
        "mean_absolute_difference": float(np.mean(np.abs(difference))),
        "maximum_absolute_difference": float(np.max(np.abs(difference))),
        "identical_rows": int(np.count_nonzero(difference == 0)),
    }


def _comparison_plot(
    summary: pd.DataFrame,
    selected: dict[tuple[str, str], pd.DataFrame],
    neutral: dict[tuple[str, str], pd.DataFrame],
    truth_1kb: pd.DataFrame,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    colors = {"mean": "#4976b8", "median": "#d45d48"}
    hard = summary.loc[
        summary["statistic"].isin(
            ["posterior_mean_call", "posterior_median_call"]
        )
    ].copy()
    hard["call"] = hard["recent_call"]
    hard["label"] = (
        hard["stride_label"] + "\n" + hard["call"].str.capitalize()
    )

    axis = axes[0, 0]
    x = np.arange(len(hard))
    lower = hard["center_null_ci95_lower"].to_numpy(dtype=float)
    upper = hard["center_null_ci95_upper"].to_numpy(dtype=float)
    axis.vlines(x, lower, upper, color="#888888", linewidth=7, alpha=0.45)
    axis.scatter(
        x,
        hard["center_null_mean"],
        color="black",
        s=55,
        zorder=3,
        label="neutral mean and 95% range",
    )
    for index, row in enumerate(hard.itertuples()):
        axis.scatter(
            index,
            row.center_selected,
            color=colors[row.call],
            s=105,
            zorder=4,
        )
        axis.annotate(
            f"p={row.center_monte_carlo_p_upper:.4f}",
            (index, row.center_selected),
            xytext=(0, -22),
            textcoords="offset points",
            ha="center",
        )
    axis.axhline(
        0.098,
        color="#188977",
        linestyle="--",
        linewidth=1.8,
        label="tree truth at selected center",
    )
    axis.set_xticks(x, hard["label"])
    axis.margins(x=0.08)
    axis.set_ylabel("Fraction of pairs called recent")
    axis.set_title("A. Both hard-call rules recover the selected center")
    axis.legend(frameon=False, loc="lower right")

    axis = axes[0, 1]
    truth_local = truth_1kb.loc[
        truth_1kb["position_0based"].between(4_000_000, 6_000_000)
    ]
    axis.plot(
        (truth_local["position_0based"] - 5_000_000) / 1e3,
        truth_local["truth_fraction_recent"],
        color="#188977",
        linewidth=2.0,
        label="tree truth",
    )
    for stride_label, linestyle, linewidth in [
        ("1 kb", "-", 1.7),
        ("10 kb", "--", 2.0),
    ]:
        for call in ["mean", "median"]:
            frame = selected[(stride_label, call)]
            local = frame.loc[
                frame["position_0based"].between(4_000_000, 6_000_000)
            ]
            axis.plot(
                (local["position_0based"] - 5_000_000) / 1e3,
                local[CALLED_COLUMN],
                color=colors[call],
                linestyle=linestyle,
                linewidth=linewidth,
                label=f"{call.capitalize()} call, {stride_label}",
            )
    axis.axvline(0, color="#999999", linewidth=1)
    axis.set_xlabel("Distance from selected site (kb)")
    axis.set_ylabel("Fraction of pairs called recent")
    axis.set_title("B. Median is closer at the peak; mean is conservative")
    axis.legend(frameon=False, ncol=2, fontsize=9)

    axis = axes[1, 0]
    limits = [0.0, 0.105]
    axis.plot(limits, limits, color="#777777", linestyle=":", linewidth=1.5)
    for call in ["mean", "median"]:
        one = selected[("1 kb", call)]
        ten = selected[("10 kb", call)]
        downsampled = one.loc[
            one["position_0based"].mod(10_000).eq(0),
            ["position_0based", CALLED_COLUMN],
        ]
        paired = downsampled.merge(
            ten[["position_0based", CALLED_COLUMN]],
            on="position_0based",
            suffixes=("_1kb", "_10kb"),
            validate="one_to_one",
        )
        left = paired[f"{CALLED_COLUMN}_1kb"].to_numpy()
        right = paired[f"{CALLED_COLUMN}_10kb"].to_numpy()
        axis.scatter(
            left,
            right,
            color=colors[call],
            s=13,
            alpha=0.45,
            label=(
                f"{call.capitalize()}: "
                f"r={_correlation(left, right):.6f}"
            ),
        )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("1 kb result at shared positions")
    axis.set_ylabel("Independent 10 kb result")
    axis.set_title("C. Hard-call profiles are stable across stride")
    axis.legend(frameon=False)

    axis = axes[1, 1]
    markers = {"1 kb": "o", "10 kb": "s"}
    for stride_label in ["1 kb", "10 kb"]:
        mean_center = neutral[(stride_label, "mean")].loc[
            lambda frame: frame["position_0based"].eq(5_000_000),
            ["replicate", CALLED_COLUMN],
        ]
        median_center = neutral[(stride_label, "median")].loc[
            lambda frame: frame["position_0based"].eq(5_000_000),
            ["replicate", CALLED_COLUMN],
        ]
        paired = mean_center.merge(
            median_center,
            on="replicate",
            suffixes=("_mean", "_median"),
            validate="one_to_one",
        )
        axis.scatter(
            paired[f"{CALLED_COLUMN}_mean"],
            paired[f"{CALLED_COLUMN}_median"],
            marker=markers[stride_label],
            s=38,
            alpha=0.5,
            label=f"neutral, {stride_label}",
        )
    axis.scatter(
        0.05,
        0.0685,
        marker="*",
        s=230,
        color="#7a3db8",
        edgecolor="black",
        linewidth=0.6,
        label="selected center",
        zorder=5,
    )
    limits = [0.002, 0.072]
    axis.plot(limits, limits, color="#777777", linestyle=":", linewidth=1.5)
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_xlabel("Fraction recent: posterior-mean rule")
    axis.set_ylabel("Fraction recent: posterior-median rule")
    axis.set_title("D. Median calls more pairs recent in selected and null data")
    axis.legend(frameon=False)

    fig.suptitle(
        "Gamma-SMC posterior mean versus median at "
        r"$\mu=1.29\times10^{-8}$",
        fontsize=21,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    source = args.source.resolve()
    analysis = root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    result_dirs = {
        "1 kb": root / "stride1kb",
        "10 kb": root / "stride10kb",
    }
    source_metrics = _load_json(source / "metrics.json")
    mutation_rate = float(source_metrics["mutation_rate"])
    if not np.isclose(mutation_rate, 1.29e-8, rtol=0, atol=1e-20):
        raise RuntimeError(f"expected mutation rate 1.29e-8; observed {mutation_rate}")

    selected: dict[tuple[str, str], pd.DataFrame] = {}
    neutral: dict[tuple[str, str], pd.DataFrame] = {}
    summary_parts = []
    soft_identity = {}
    run_contract = {}
    for stride_label, expected_stride in [("1 kb", 1_000), ("10 kb", 10_000)]:
        result_dir = result_dirs[stride_label]
        metrics = _load_json(result_dir / "decoded_study_metrics.json")
        required = {
            "stride_bp": expected_stride,
            "neutral_replicates": 100,
            "neutral_workers": 12,
            "recent_call": "mean",
            "comparison_recent_call": "median",
            "cache_size_bp": 1_000,
        }
        for key, expected in required.items():
            if metrics[key] != expected:
                raise RuntimeError(
                    f"{stride_label} {key} is {metrics[key]!r}; "
                    f"expected {expected!r}"
                )
        run_contract[stride_label] = {
            **required,
            "observed_output_positions": int(metrics["observed_output_positions"]),
            "elapsed_seconds": float(metrics["elapsed_seconds"]),
        }
        part = pd.read_csv(
            result_dir / "posterior_summary_rule_center_comparison.tsv",
            sep="\t",
        )
        part.insert(0, "stride_label", stride_label)
        part.insert(1, "stride_bp", expected_stride)
        summary_parts.append(part)

        for call in ["mean", "median"]:
            selected[(stride_label, call)] = _read_selected(result_dir, call)
            neutral[(stride_label, call)] = _read_neutral(result_dir, call)
        _assert_soft_identity(
            selected[(stride_label, "mean")],
            selected[(stride_label, "median")],
        )
        _assert_soft_identity(
            neutral[(stride_label, "mean")],
            neutral[(stride_label, "median")],
        )
        soft_identity[stride_label] = {
            "selected_rows": int(len(selected[(stride_label, "mean")])),
            "neutral_rows": int(len(neutral[(stride_label, "mean")])),
            "all_soft_columns_identical_at_tolerance_1e-12": True,
        }

    summary = pd.concat(summary_parts, ignore_index=True)
    summary.to_csv(
        analysis / "mean_median_stride_center_summary.tsv",
        sep="\t",
        index=False,
    )

    truth_1kb = pd.read_csv(
        result_dirs["1 kb"] / "selected_decoded_vs_truth.tsv.gz",
        sep="\t",
        usecols=["position_0based", "truth_fraction_recent"],
    )
    comparisons = {}
    for call in ["mean", "median"]:
        comparisons[call] = {
            "selected_stride": _shared_stride_metrics(
                selected[("1 kb", call)],
                selected[("10 kb", call)],
                keys=["position_0based"],
                value=CALLED_COLUMN,
            ),
            "neutral_stride": _shared_stride_metrics(
                neutral[("1 kb", call)],
                neutral[("10 kb", call)],
                keys=["replicate", "position_0based"],
                value=CALLED_COLUMN,
            ),
        }

    rule_effect = {}
    paired_center_rows = []
    for stride_label in ["1 kb", "10 kb"]:
        selected_mean = _center(selected[(stride_label, "mean")])
        selected_median = _center(selected[(stride_label, "median")])
        mean_null_center = neutral[(stride_label, "mean")].loc[
            lambda frame: frame["position_0based"].eq(5_000_000),
            ["replicate", CALLED_COLUMN],
        ]
        median_null_center = neutral[(stride_label, "median")].loc[
            lambda frame: frame["position_0based"].eq(5_000_000),
            ["replicate", CALLED_COLUMN],
        ]
        paired = mean_null_center.merge(
            median_null_center,
            on="replicate",
            suffixes=("_mean", "_median"),
            validate="one_to_one",
        )
        paired.insert(0, "stride_label", stride_label)
        paired_center_rows.append(paired)
        mean_values = paired[f"{CALLED_COLUMN}_mean"].to_numpy(dtype=float)
        median_values = paired[f"{CALLED_COLUMN}_median"].to_numpy(dtype=float)
        rule_effect[stride_label] = {
            "selected_center_mean_call": float(selected_mean[CALLED_COLUMN]),
            "selected_center_median_call": float(selected_median[CALLED_COLUMN]),
            "selected_center_median_minus_mean": float(
                selected_median[CALLED_COLUMN] - selected_mean[CALLED_COLUMN]
            ),
            "selected_center_truth": 0.098,
            "neutral_center_mean_call_average": float(np.mean(mean_values)),
            "neutral_center_median_call_average": float(np.mean(median_values)),
            "neutral_center_median_minus_mean_average": float(
                np.mean(median_values - mean_values)
            ),
            "neutral_center_rule_correlation": _correlation(
                mean_values,
                median_values,
            ),
        }
    pd.concat(paired_center_rows, ignore_index=True).to_csv(
        analysis / "paired_null_center_calls.tsv",
        sep="\t",
        index=False,
    )

    metrics = {
        "mutation_rate": mutation_rate,
        "interpretation": (
            "recent-call mean versus median changes only hard calls; "
            "soft posterior probabilities are invariant"
        ),
        "run_contract": run_contract,
        "soft_posterior_identity": soft_identity,
        "stride_comparison": comparisons,
        "call_rule_effect": rule_effect,
    }
    (analysis / "mean_median_stride_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n",
        encoding="utf-8",
    )
    _comparison_plot(
        summary,
        selected,
        neutral,
        truth_1kb,
        analysis / "mu1p29e8_mean_vs_median_1kb_10kb.png",
    )


if __name__ == "__main__":
    main()
