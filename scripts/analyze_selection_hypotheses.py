#!/usr/bin/env python3
"""Analyze matched stride, mutation-rate, and carrier-class sweep tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit
from scipy.stats import beta

from gamma_smc_aou.decoder import run_within_decoder


SOFT_COLUMNS = [
    "mean_p_tmrca_lt_threshold",
    "mean_p_lt_4500",
    "mean_tmrca_generations",
]


def _at_stride(frame: pd.DataFrame, stride: int) -> pd.DataFrame:
    positions = frame["position_0based"].to_numpy(dtype=float)
    keep = np.isclose(np.mod(positions, stride), 0.0, atol=1e-8, rtol=0)
    return frame.loc[keep].reset_index(drop=True)


def _center_row(frame: pd.DataFrame, center: float) -> pd.Series:
    index = int(
        np.argmin(
            np.abs(frame["position_0based"].to_numpy(dtype=float) - center)
        )
    )
    return frame.iloc[index]


def _metric(metrics: dict, current: str, legacy: str) -> float:
    value = metrics.get(current, metrics.get(legacy))
    if value is None:
        raise KeyError(f"metrics contain neither {current!r} nor {legacy!r}")
    return float(value)


def _write_pair_files(carriers: pd.DataFrame, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for pair_class in ["hom_ref", "heterozygous", "hom_alt"]:
        subset = carriers.loc[carriers["focal_genotype_class"] == pair_class]
        path = output_dir / f"{pair_class}.pairs.tsv"
        lines = [
            "# focal-genotype-stratified within-individual haplotype pairs",
            *[
                f"{int(row.gamma_smc_haplotype_0)}\t"
                f"{int(row.gamma_smc_haplotype_1)}"
                for row in subset.itertuples()
            ],
        ]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        paths[pair_class] = path
    return paths


def _compare_low_stride_outputs(low_1kb: Path, low_10kb: Path) -> dict:
    selected_1 = _at_stride(
        pd.read_csv(
            low_1kb / "selected_decoded_recent_probability_profile.tsv",
            sep="\t",
        ),
        10_000,
    )
    selected_10 = pd.read_csv(
        low_10kb / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
    )
    neutral_1 = _at_stride(
        pd.read_csv(
            low_1kb / "neutral_decoded_recent_probability_profiles.tsv.gz",
            sep="\t",
        ),
        10_000,
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    neutral_10 = pd.read_csv(
        low_10kb / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    if len(selected_1) != len(selected_10) or len(neutral_1) != len(neutral_10):
        raise RuntimeError("1 kb and 10 kb profiles do not share the expected grid")
    result = {}
    for label, left, right in [
        ("selected", selected_1, selected_10),
        ("neutral", neutral_1, neutral_10),
    ]:
        for column in SOFT_COLUMNS:
            difference = np.abs(
                left[column].to_numpy(dtype=float)
                - right[column].to_numpy(dtype=float)
            )
            result[f"{label}_{column}_max_abs_difference"] = float(
                difference.max(initial=0.0)
            )
    result["selected_shared_positions"] = int(len(selected_1))
    result["neutral_shared_rows"] = int(len(neutral_1))
    return result


def _decode_high_selected_at_10kb(
    high_source: Path,
    executable: Path,
    output_dir: Path,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "selected_decoded_recent_probability_profile.tsv"
    run_path = path.with_suffix(path.suffix + ".run.json")
    if path.exists() and run_path.exists():
        return json.loads(run_path.read_text(encoding="utf-8"))
    run = run_within_decoder(
        executable,
        high_source / "selected_s0p05_af30.vcf.gz",
        path,
        scaled_mutation_rate=0.0005,
        recombination_to_mutation_ratio=0.8,
        mutation_rate=1.25e-8,
        threshold_years=4_500,
        generation_time=25,
        input_format="vcf",
        output_at_stride=10_000,
        output_at_hets=False,
        only_within=True,
        recent_call="mean",
        cache_size=1_000,
        threads=1,
    )
    (output_dir / "selected_only_metrics.json").write_text(
        json.dumps({
            "mutation_rate": 1.25e-8,
            "stride_bp": 10_000,
            "neutral_replicates": 0,
            "interpretation": (
                "Independent selected decode used only for the stride "
                "sensitivity comparison; no 10 kb high-mutation null p-value."
            ),
            **run,
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    return run


def _decode_carrier_classes(
    sources: dict[str, Path],
    pair_paths: dict[str, Path],
    executable: Path,
    output_dir: Path,
) -> pd.DataFrame:
    rows = []
    configurations = {
        "low_mu": {
            "mutation_rate": 1.29e-9,
            "theta": 5.16e-5,
            "rho_over_theta": 1e-8 / 1.29e-9,
        },
        "high_mu": {
            "mutation_rate": 1.25e-8,
            "theta": 0.0005,
            "rho_over_theta": 0.8,
        },
    }
    for configuration, source in sources.items():
        rates = configurations[configuration]
        for pair_class, pairs_path in pair_paths.items():
            summary_path = output_dir / f"{configuration}_{pair_class}_stride10kb.tsv"
            run_path = summary_path.with_suffix(summary_path.suffix + ".run.json")
            if summary_path.exists() and run_path.exists():
                run = json.loads(run_path.read_text(encoding="utf-8"))
            else:
                run = run_within_decoder(
                    executable,
                    source / "selected_s0p05_af30.vcf.gz",
                    summary_path,
                    scaled_mutation_rate=rates["theta"],
                    recombination_to_mutation_ratio=rates["rho_over_theta"],
                    mutation_rate=rates["mutation_rate"],
                    threshold_years=4_500,
                    generation_time=25,
                    input_format="vcf",
                    output_at_stride=10_000,
                    output_at_hets=False,
                    only_within=False,
                    pairs_file=pairs_path,
                    recent_call="mean",
                    cache_size=1_000,
                    threads=1,
                )
            profile = pd.read_csv(summary_path, sep="\t")
            center = _center_row(profile, 5_000_000)
            rows.append({
                "configuration": configuration,
                "mutation_rate": rates["mutation_rate"],
                "pair_class": pair_class,
                "n_pairs": int(center["n_pairs"]),
                "center_mean_p_recent": float(
                    center["mean_p_tmrca_lt_threshold"]
                ),
                "center_mean_tmrca_generations_decoded": float(
                    center["mean_tmrca_generations"]
                ),
                "center_called_fraction_recent_posterior_mean_rule": float(
                    center["frac_recent_4500"]
                ),
                "decode_seconds": float(run["decode_seconds"]),
            })
    return pd.DataFrame(rows)


def _truth_carrier_summary(
    source: Path,
    carriers: pd.DataFrame,
    neutral_truth: Path,
) -> tuple[pd.DataFrame, dict]:
    ts = tskit.load(source / "selected_s0p05_af30.trees")
    tree = ts.at(5_000_000)
    frame = carriers.copy()
    frame["tmrca_generations"] = [
        tree.tmrca(int(row.sample_node_0), int(row.sample_node_1))
        for row in frame.itertuples()
    ]
    frame["tmrca_lt_180_generations"] = frame["tmrca_generations"] < 180
    rows = []
    for pair_class, subset in [
        ("all_within_individual", frame),
        *list(frame.groupby("focal_genotype_class", sort=False)),
    ]:
        rows.append({
            "pair_class": pair_class,
            "n_pairs": int(len(subset)),
            "center_fraction_recent_truth": float(
                subset["tmrca_lt_180_generations"].mean()
            ),
            "center_mean_tmrca_generations_truth": float(
                subset["tmrca_generations"].mean()
            ),
        })
    neutral = pd.read_csv(neutral_truth / "neutral_statistics.tsv", sep="\t")
    summary = {
        "neutral_replicates": int(len(neutral)),
        "neutral_mean_fraction_recent": float(
            neutral["center_fraction_recent"].mean()
        ),
        "neutral_fraction_recent_ci95": [
            float(neutral["center_fraction_recent"].quantile(0.025)),
            float(neutral["center_fraction_recent"].quantile(0.975)),
        ],
        "neutral_mean_tmrca_generations": float(
            neutral["center_mean_tmrca_generations"].mean()
        ),
        "neutral_mean_tmrca_ci95": [
            float(neutral["center_mean_tmrca_generations"].quantile(0.025)),
            float(neutral["center_mean_tmrca_generations"].quantile(0.975)),
        ],
    }
    return pd.DataFrame(rows), summary


def _factor_summary(
    cells: list[tuple[str, float, int, Path, int, bool]],
) -> pd.DataFrame:
    rows = []
    for label, mutation_rate, stride, path, sites, derived in cells:
        metrics = json.loads((path / "decoded_study_metrics.json").read_text())
        comparison = metrics.get("posterior_summary_rule_comparison") or {}
        soft = comparison.get("soft_probability")
        if soft is None:
            selected = _metric(
                metrics,
                "center_observed_mean_p_recent",
                "center_observed_statistic",
            )
            neutral_mean = _metric(
                metrics,
                "center_null_mean_p_recent",
                "center_null_mean_statistic",
            )
            neutral_lower = float(metrics["center_null_ci95_lower"])
            neutral_upper = float(metrics["center_null_ci95_upper"])
            exceedances = int(metrics["center_neutral_exceedances"])
            pvalue = float(metrics["center_monte_carlo_p_upper"])
        else:
            selected = float(soft["center_selected"])
            neutral_mean = float(soft["center_null_mean"])
            neutral_lower = float(soft["center_null_ci95_lower"])
            neutral_upper = float(soft["center_null_ci95_upper"])
            exceedances = int(soft["center_neutral_exceedances"])
            pvalue = float(soft["center_monte_carlo_p_upper"])
        rows.append({
            "configuration": label,
            "mutation_rate": mutation_rate,
            "retained_sites_selected": sites,
            "stride_bp": stride,
            "center_truth_fraction_recent": 0.098,
            "center_selected_mean_p_recent": selected,
            "center_null_mean_p_recent": neutral_mean,
            "center_null_ci95_lower": neutral_lower,
            "center_null_ci95_upper": neutral_upper,
            "center_neutral_exceedances": exceedances,
            "center_monte_carlo_p_upper": pvalue,
            "derived_by_downsampling_1kb": derived,
        })
    return pd.DataFrame(rows)


def _plot_summary(
    factor: pd.DataFrame,
    truth: pd.DataFrame,
    truth_null: dict,
    carrier: pd.DataFrame,
    low_1kb: Path,
    low_10kb: Path,
    high_1kb: Path,
    high_10kb: Path,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))
    colors = {"low_mu": "#4976b8", "high_mu": "#d45d48"}

    axis = axes[0, 0]
    x = np.arange(len(factor))
    y = factor["center_selected_mean_p_recent"].to_numpy()
    lower = factor["center_null_ci95_lower"].to_numpy()
    upper = factor["center_null_ci95_upper"].to_numpy()
    axis.vlines(x, lower, upper, color="#777777", linewidth=6, alpha=0.5)
    axis.scatter(x, factor["center_null_mean_p_recent"], color="black", s=55, label="neutral mean")
    low = factor["mutation_rate"].to_numpy() < 1e-8
    axis.scatter(
        x[low], y[low], color=colors["low_mu"], s=95,
        label="selected, low μ", zorder=3,
    )
    axis.scatter(
        x[~low], y[~low], color=colors["high_mu"], s=95,
        label="selected, high μ", zorder=3,
    )
    axis.annotate(
        "selected only",
        (x[-1], y[-1]),
        xytext=(0, -13),
        textcoords="offset points",
        ha="center",
        va="top",
        fontsize=10,
    )
    axis.set_xticks(x, [f"{row.stride_bp // 1000} kb\n$\\mu$={row.mutation_rate:.2g}" for row in factor.itertuples()])
    axis.set_ylabel("Mean posterior P(TMRCA < 180 generations)")
    axis.set_title("A. Center statistic: mutation information, not stride")
    axis.legend(frameon=False)

    axis = axes[0, 1]
    profiles = [
        ("low μ, 1 kb", low_1kb, colors["low_mu"], "-"),
        ("low μ, 10 kb", low_10kb, colors["low_mu"], "--"),
        ("high μ, 1 kb", high_1kb, colors["high_mu"], "-"),
        ("high μ, 10 kb", high_10kb, colors["high_mu"], "--"),
    ]
    for label, path, color, linestyle in profiles:
        frame = pd.read_csv(
            path / "selected_decoded_recent_probability_profile.tsv",
            sep="\t",
        )
        local = frame.loc[
            frame["position_0based"].between(4_500_000, 5_500_000)
        ]
        axis.plot(
            (local["position_0based"] - 5_000_000) / 1_000,
            local["mean_p_tmrca_lt_threshold"],
            color=color,
            linestyle=linestyle,
            linewidth=1.8,
            label=label,
        )
    axis.axvline(0, color="black", linewidth=1, alpha=0.5)
    axis.set_xlabel("Distance from selected site (kb)")
    axis.set_ylabel("Mean posterior P(recent)")
    axis.set_title("B. Local decoded profiles")
    axis.legend(frameon=False, ncol=2)

    axis = axes[1, 0]
    order = ["hom_alt", "heterozygous", "hom_ref", "all_within_individual"]
    truth_indexed = truth.set_index("pair_class").loc[order]
    tx = np.arange(len(order))
    axis.bar(
        tx,
        truth_indexed["center_mean_tmrca_generations_truth"],
        color=["#5aa469", "#d9a441", "#8b6bb1", "#777777"],
    )
    neutral_mean = truth_null["neutral_mean_tmrca_generations"]
    neutral_ci = truth_null["neutral_mean_tmrca_ci95"]
    axis.axhline(neutral_mean, color="black", linestyle="--", label="neutral mean")
    axis.axhspan(neutral_ci[0], neutral_ci[1], color="black", alpha=0.08, label="neutral 95% range")
    axis.set_yscale("log")
    axis.set_xticks(tx, ["hom-alt", "heterozygous", "hom-ref", "all pairs"])
    axis.set_ylabel("True mean TMRCA (generations, log scale)")
    axis.set_title("C. Carrier contrast is real—but hom-ref is not neutral")
    axis.legend(frameon=False)

    axis = axes[1, 1]
    class_order = ["hom_alt", "heterozygous", "hom_ref"]
    width = 0.34
    cx = np.arange(len(class_order))
    for offset, configuration in [(-width / 2, "low_mu"), (width / 2, "high_mu")]:
        subset = (
            carrier.loc[carrier["configuration"] == configuration]
            .set_index("pair_class")
            .loc[class_order]
        )
        axis.bar(
            cx + offset,
            subset["center_mean_p_recent"],
            width=width,
            color=colors[configuration],
            label=configuration.replace("_", " "),
        )
    axis.axhline(
        factor.loc[factor["mutation_rate"] < 1e-8, "center_null_mean_p_recent"].iloc[0],
        color="black",
        linestyle=":",
        label="low-μ neutral mean",
    )
    axis.set_yscale("log")
    axis.set_xticks(cx, ["hom-alt", "heterozygous", "hom-ref"])
    axis.set_ylabel("Decoded mean posterior P(recent), log scale")
    axis.set_title("D. Low mutation rate erases the hom-alt signal")
    axis.legend(frameon=False)

    fig.suptitle(
        "Why the selected sweep was recovered before: matched hypothesis tests",
        fontsize=22,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--existing-low-10kb", type=Path, required=True)
    parser.add_argument("--neutral-truth", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    low_source = root / "source_low_mu"
    high_source = root / "source_high_mu"
    low_1kb = root / "low_mu_stride1kb"
    high_1kb = root / "high_mu_stride1kb"
    high_10kb = root / "high_mu_stride10kb_selected_only"
    low_10kb = args.existing_low_10kb.resolve()
    analysis_dir = root / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    stride_validation = _compare_low_stride_outputs(low_1kb, low_10kb)
    _decode_high_selected_at_10kb(high_source, args.executable, high_10kb)
    high_independent = pd.read_csv(
        high_10kb / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
    )
    high_1kb_at_10kb = _at_stride(
        pd.read_csv(
            high_1kb / "selected_decoded_recent_probability_profile.tsv",
            sep="\t",
        ),
        10_000,
    )
    stride_validation["high_mu_selected_mean_p_recent_max_abs_difference"] = float(
        np.max(
            np.abs(
                high_independent["mean_p_tmrca_lt_threshold"].to_numpy()
                - high_1kb_at_10kb["mean_p_tmrca_lt_threshold"].to_numpy()
            )
        )
    )

    carriers = pd.read_csv(
        low_source / "selected_focal_carrier_pairs.tsv",
        sep="\t",
    )
    pair_paths = _write_pair_files(carriers, analysis_dir / "pairs")
    carrier = _decode_carrier_classes(
        {"low_mu": low_source, "high_mu": high_source},
        pair_paths,
        args.executable,
        analysis_dir,
    )
    carrier.to_csv(
        analysis_dir / "carrier_class_decoded_center.tsv",
        sep="\t",
        index=False,
    )

    truth, truth_null = _truth_carrier_summary(
        low_source,
        carriers,
        args.neutral_truth.resolve(),
    )
    truth.to_csv(
        analysis_dir / "carrier_class_truth_center.tsv",
        sep="\t",
        index=False,
    )
    (analysis_dir / "truth_neutral_center_summary.json").write_text(
        json.dumps(truth_null, indent=2) + "\n",
        encoding="utf-8",
    )

    low_sites = tskit.load(low_source / "selected_s0p05_af30.trees").num_sites
    high_sites = tskit.load(high_source / "selected_s0p05_af30.trees").num_sites
    factor = _factor_summary([
        ("low_mu_stride1kb", 1.29e-9, 1_000, low_1kb, low_sites, False),
        ("low_mu_stride10kb", 1.29e-9, 10_000, low_10kb, low_sites, False),
        ("high_mu_stride1kb", 1.25e-8, 1_000, high_1kb, high_sites, False),
    ])
    high_10_center = _center_row(high_independent, 5_000_000)
    factor = pd.concat([
        factor,
        pd.DataFrame([{
            "configuration": "high_mu_stride10kb_selected_only",
            "mutation_rate": 1.25e-8,
            "retained_sites_selected": high_sites,
            "stride_bp": 10_000,
            "center_truth_fraction_recent": 0.098,
            "center_selected_mean_p_recent": float(
                high_10_center["mean_p_tmrca_lt_threshold"]
            ),
            "center_null_mean_p_recent": np.nan,
            "center_null_ci95_lower": np.nan,
            "center_null_ci95_upper": np.nan,
            "center_neutral_exceedances": np.nan,
            "center_monte_carlo_p_upper": np.nan,
            "derived_by_downsampling_1kb": False,
        }]),
    ], ignore_index=True)
    factor.to_csv(
        analysis_dir / "stride_mutation_factorial_center.tsv",
        sep="\t",
        index=False,
    )

    low_exceedances = int(
        factor.loc[
            factor["configuration"] == "low_mu_stride10kb",
            "center_neutral_exceedances",
        ].iloc[0]
    )
    n_null = 100
    exceedance_ci = [
        float(beta.ppf(0.025, low_exceedances, n_null - low_exceedances + 1))
        if low_exceedances
        else 0.0,
        float(beta.ppf(0.975, low_exceedances + 1, n_null - low_exceedances))
        if low_exceedances < n_null
        else 1.0,
    ]
    summary = {
        "primary_statistic": "mean across pairs of posterior P(TMRCA < 180 generations)",
        "posterior_mean_call_is_not_primary_statistic": True,
        "stride_validation": stride_validation,
        "low_mu_null_exceedance_probability_clopper_pearson_95": exceedance_ci,
        "high_mu_10kb_calibration": (
            "selected-only independent decode; the 100-null p-value is from "
            "the complete high-mutation 1 kb run"
        ),
        "truth_null_summary": truth_null,
    }
    (analysis_dir / "hypothesis_test_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    _plot_summary(
        factor,
        truth,
        truth_null,
        carrier,
        low_1kb,
        low_10kb,
        high_1kb,
        high_10kb,
        analysis_dir / "selection_recovery_hypothesis_tests.png",
    )


if __name__ == "__main__":
    main()
