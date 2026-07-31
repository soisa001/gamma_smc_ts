#!/usr/bin/env python3
"""Analyze a matched two-rate, two-stride selected-sweep experiment."""

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


def _compare_stride_outputs(stride_1kb: Path, stride_10kb: Path) -> dict:
    selected_1 = _at_stride(
        pd.read_csv(
            stride_1kb / "selected_decoded_recent_probability_profile.tsv",
            sep="\t",
        ),
        10_000,
    )
    selected_10 = pd.read_csv(
        stride_10kb / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
    )
    neutral_1 = _at_stride(
        pd.read_csv(
            stride_1kb / "neutral_decoded_recent_probability_profiles.tsv.gz",
            sep="\t",
        ),
        10_000,
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    neutral_10 = pd.read_csv(
        stride_10kb / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    if len(selected_1) != len(selected_10) or len(neutral_1) != len(neutral_10):
        raise RuntimeError("1 kb and 10 kb profiles do not share the expected grid")
    if not np.array_equal(
        selected_1["position_0based"].to_numpy(),
        selected_10["position_0based"].to_numpy(),
    ):
        raise RuntimeError("selected profiles have different 10 kb positions")
    for column in ["replicate", "position_0based"]:
        if not np.array_equal(
            neutral_1[column].to_numpy(),
            neutral_10[column].to_numpy(),
        ):
            raise RuntimeError(
                f"neutral profiles differ in aligned {column!r} values"
            )
    result = {}
    for label, left, right in [
        ("selected", selected_1, selected_10),
        ("neutral", neutral_1, neutral_10),
    ]:
        for column in SOFT_COLUMNS:
            left_values = left[column].to_numpy(dtype=float)
            right_values = right[column].to_numpy(dtype=float)
            difference = np.abs(
                left_values
                - right_values
            )
            result[f"{label}_{column}_max_abs_difference"] = float(
                difference.max(initial=0.0)
            )
            result[f"{label}_{column}_pearson_correlation"] = float(
                np.corrcoef(left_values, right_values)[0, 1]
            )
    result["selected_shared_positions"] = int(len(selected_1))
    result["neutral_shared_rows"] = int(len(neutral_1))
    return result


def _compare_rate_outputs(rate_a: Path, rate_b: Path) -> dict:
    selected_a = pd.read_csv(
        rate_a / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
    )
    selected_b = pd.read_csv(
        rate_b / "selected_decoded_recent_probability_profile.tsv",
        sep="\t",
    )
    neutral_a = pd.read_csv(
        rate_a / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    neutral_b = pd.read_csv(
        rate_b / "neutral_decoded_recent_probability_profiles.tsv.gz",
        sep="\t",
    ).sort_values(["replicate", "position_0based"]).reset_index(drop=True)
    if len(selected_a) != len(selected_b) or len(neutral_a) != len(neutral_b):
        raise RuntimeError("mutation-rate profiles have different dimensions")
    if not np.array_equal(
        selected_a["position_0based"].to_numpy(),
        selected_b["position_0based"].to_numpy(),
    ):
        raise RuntimeError("selected mutation-rate profiles use different grids")
    for column in ["replicate", "position_0based"]:
        if not np.array_equal(
            neutral_a[column].to_numpy(),
            neutral_b[column].to_numpy(),
        ):
            raise RuntimeError(
                f"neutral mutation-rate profiles differ in {column!r}"
            )

    result = {}
    for label, left, right in [
        ("selected", selected_a, selected_b),
        ("neutral", neutral_a, neutral_b),
    ]:
        for column in SOFT_COLUMNS:
            left_values = left[column].to_numpy(dtype=float)
            right_values = right[column].to_numpy(dtype=float)
            difference = left_values - right_values
            result[f"{label}_{column}_mean_difference_a_minus_b"] = float(
                np.mean(difference)
            )
            result[f"{label}_{column}_max_abs_difference"] = float(
                np.max(np.abs(difference), initial=0.0)
            )
            result[f"{label}_{column}_pearson_correlation"] = float(
                np.corrcoef(left_values, right_values)[0, 1]
            )
    result["selected_shared_positions"] = int(len(selected_a))
    result["neutral_shared_rows"] = int(len(neutral_a))
    return result


def _decode_carrier_classes(
    sources: dict[str, Path],
    pair_paths: dict[str, Path],
    executable: Path,
    output_dir: Path,
) -> pd.DataFrame:
    rows = []
    for configuration, source in sources.items():
        metrics = json.loads((source / "metrics.json").read_text())
        mutation_rate = float(metrics["mutation_rate"])
        ancestral_size = float(
            metrics["demography"]["ancestral_population_size"]
        )
        recombination_rate = float(metrics["recombination_rate"])
        theta = 4 * ancestral_size * mutation_rate
        rho_over_theta = recombination_rate / mutation_rate
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
                    scaled_mutation_rate=theta,
                    recombination_to_mutation_ratio=rho_over_theta,
                    mutation_rate=mutation_rate,
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
                "mutation_rate": mutation_rate,
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
    rate_a_1kb: Path,
    rate_a_10kb: Path,
    rate_b_1kb: Path,
    rate_b_10kb: Path,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))
    colors = {"rate_a": "#4976b8", "rate_b": "#d45d48"}
    rate_a = float(
        factor.loc[
            factor["configuration"].str.startswith("rate_a"),
            "mutation_rate",
        ].iloc[0]
    )
    rate_b = float(
        factor.loc[
            factor["configuration"].str.startswith("rate_b"),
            "mutation_rate",
        ].iloc[0]
    )

    axis = axes[0, 0]
    x = np.arange(len(factor))
    y = factor["center_selected_mean_p_recent"].to_numpy()
    lower = factor["center_null_ci95_lower"].to_numpy()
    upper = factor["center_null_ci95_upper"].to_numpy()
    axis.vlines(x, lower, upper, color="#777777", linewidth=6, alpha=0.5)
    axis.scatter(x, factor["center_null_mean_p_recent"], color="black", s=55, label="neutral mean")
    is_rate_a = factor["configuration"].str.startswith("rate_a").to_numpy()
    axis.scatter(
        x[is_rate_a], y[is_rate_a], color=colors["rate_a"], s=95,
        label=f"selected, μ={rate_a:.3g}", zorder=3,
    )
    axis.scatter(
        x[~is_rate_a], y[~is_rate_a], color=colors["rate_b"], s=95,
        label=f"selected, μ={rate_b:.3g}", zorder=3,
    )
    for index, row in enumerate(factor.itertuples()):
        offset = 10 if row.center_selected_mean_p_recent < 0.04 else -20
        axis.annotate(
            f"p={row.center_monte_carlo_p_upper:.3g}",
            (index, row.center_selected_mean_p_recent),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va="bottom" if offset > 0 else "top",
            fontsize=10,
        )
    axis.set_xticks(x, [f"{row.stride_bp // 1000} kb\n$\\mu$={row.mutation_rate:.2g}" for row in factor.itertuples()])
    axis.set_ylabel("Mean posterior P(TMRCA < 180 generations)")
    axis.set_title("A. Both corrected mutation rates recover the sweep")
    axis.legend(frameon=False)

    axis = axes[0, 1]
    profiles = [
        (f"μ={rate_a:.3g}, 1 kb", rate_a_1kb, colors["rate_a"], "-"),
        (f"μ={rate_a:.3g}, 10 kb", rate_a_10kb, colors["rate_a"], "--"),
        (f"μ={rate_b:.3g}, 1 kb", rate_b_1kb, colors["rate_b"], "-"),
        (f"μ={rate_b:.3g}, 10 kb", rate_b_10kb, colors["rate_b"], "--"),
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
    for offset, configuration in [
        (-width / 2, "rate_a"),
        (width / 2, "rate_b"),
    ]:
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
            label=(
                f"μ={rate_a:.3g}"
                if configuration == "rate_a"
                else f"μ={rate_b:.3g}"
            ),
        )
    axis.axhline(
        factor.iloc[0]["center_null_mean_p_recent"],
        color="black",
        linestyle=":",
        label=f"μ={rate_a:.3g} neutral mean",
    )
    axis.set_yscale("log")
    axis.set_xticks(cx, ["hom-alt", "heterozygous", "hom-ref"])
    axis.set_ylabel("Decoded mean posterior P(recent), log scale")
    axis.set_title("D. Both rates recover the hom-alt carrier signal")
    axis.legend(frameon=False)

    fig.suptitle(
        "Corrected mutation-rate comparison on matched selected and null simulations",
        fontsize=22,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--rate-a-source", type=Path, required=True)
    parser.add_argument("--rate-a-1kb", type=Path, required=True)
    parser.add_argument("--rate-a-10kb", type=Path, required=True)
    parser.add_argument("--rate-b-source", type=Path, required=True)
    parser.add_argument("--rate-b-1kb", type=Path, required=True)
    parser.add_argument("--rate-b-10kb", type=Path, required=True)
    parser.add_argument("--expected-rate-a", type=float, default=1.29e-8)
    parser.add_argument("--expected-rate-b", type=float, default=1.25e-8)
    parser.add_argument("--neutral-truth", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    rate_a_source = args.rate_a_source.resolve()
    rate_a_1kb = args.rate_a_1kb.resolve()
    rate_a_10kb = args.rate_a_10kb.resolve()
    rate_b_source = args.rate_b_source.resolve()
    rate_b_1kb = args.rate_b_1kb.resolve()
    rate_b_10kb = args.rate_b_10kb.resolve()
    analysis_dir = root / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    def source_rate(path: Path) -> float:
        metrics = json.loads((path / "metrics.json").read_text())
        return float(metrics["mutation_rate"])

    rate_a = source_rate(rate_a_source)
    rate_b = source_rate(rate_b_source)
    if not np.isclose(rate_a, args.expected_rate_a, rtol=0, atol=1e-20):
        raise RuntimeError(
            f"rate A is {rate_a:g}; expected {args.expected_rate_a:g}"
        )
    if not np.isclose(rate_b, args.expected_rate_b, rtol=0, atol=1e-20):
        raise RuntimeError(
            f"rate B is {rate_b:g}; expected {args.expected_rate_b:g}"
        )

    stride_validation = {}
    for prefix, stride_1kb, stride_10kb in [
        ("rate_a", rate_a_1kb, rate_a_10kb),
        ("rate_b", rate_b_1kb, rate_b_10kb),
    ]:
        stride_validation.update({
            f"{prefix}_{key}": value
            for key, value in _compare_stride_outputs(
                stride_1kb,
                stride_10kb,
            ).items()
        })
    fixed_stride_rate_validation = {
        "stride_1kb": _compare_rate_outputs(rate_a_1kb, rate_b_1kb),
        "stride_10kb": _compare_rate_outputs(rate_a_10kb, rate_b_10kb),
    }

    carriers = pd.read_csv(
        rate_a_source / "selected_focal_carrier_pairs.tsv",
        sep="\t",
    )
    pair_paths = _write_pair_files(carriers, analysis_dir / "pairs")
    carrier = _decode_carrier_classes(
        {"rate_a": rate_a_source, "rate_b": rate_b_source},
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
        rate_a_source,
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

    rate_a_sites = tskit.load(
        rate_a_source / "selected_s0p05_af30.trees"
    ).num_sites
    rate_b_sites = tskit.load(
        rate_b_source / "selected_s0p05_af30.trees"
    ).num_sites
    factor = _factor_summary([
        ("rate_a_stride1kb", rate_a, 1_000, rate_a_1kb, rate_a_sites, False),
        ("rate_a_stride10kb", rate_a, 10_000, rate_a_10kb, rate_a_sites, False),
        ("rate_b_stride1kb", rate_b, 1_000, rate_b_1kb, rate_b_sites, False),
        ("rate_b_stride10kb", rate_b, 10_000, rate_b_10kb, rate_b_sites, False),
    ])
    factor.to_csv(
        analysis_dir / "stride_mutation_factorial_center.tsv",
        sep="\t",
        index=False,
    )
    factor_indexed = factor.set_index("configuration")
    center_rate_comparison = {}
    for stride_label in ["1kb", "10kb"]:
        row_a = factor_indexed.loc[f"rate_a_stride{stride_label}"]
        row_b = factor_indexed.loc[f"rate_b_stride{stride_label}"]
        selected_difference = (
            row_a["center_selected_mean_p_recent"]
            - row_b["center_selected_mean_p_recent"]
        )
        center_rate_comparison[f"stride_{stride_label}"] = {
            "selected_a_minus_b": float(selected_difference),
            "selected_relative_difference_a_vs_b": float(
                row_a["center_selected_mean_p_recent"]
                / row_b["center_selected_mean_p_recent"]
                - 1
            ),
            "null_mean_a_minus_b": float(
                row_a["center_null_mean_p_recent"]
                - row_b["center_null_mean_p_recent"]
            ),
            "rate_a_p_upper": float(row_a["center_monte_carlo_p_upper"]),
            "rate_b_p_upper": float(row_b["center_monte_carlo_p_upper"]),
        }

    rate_a_exceedances = int(
        factor.loc[
            factor["configuration"] == "rate_a_stride10kb",
            "center_neutral_exceedances",
        ].iloc[0]
    )
    n_null = 100
    exceedance_ci = [
        float(
            beta.ppf(
                0.025,
                rate_a_exceedances,
                n_null - rate_a_exceedances + 1,
            )
        )
        if rate_a_exceedances
        else 0.0,
        float(
            beta.ppf(
                0.975,
                rate_a_exceedances + 1,
                n_null - rate_a_exceedances,
            )
        )
        if rate_a_exceedances < n_null
        else 1.0,
    ]
    summary = {
        "primary_statistic": "mean across pairs of posterior P(TMRCA < 180 generations)",
        "posterior_mean_call_is_not_primary_statistic": True,
        "mutation_rates": {
            "rate_a": rate_a,
            "rate_b": rate_b,
            "relative_difference_rate_a_vs_rate_b": rate_a / rate_b - 1,
        },
        "stride_validation": stride_validation,
        "fixed_stride_rate_validation": fixed_stride_rate_validation,
        "center_rate_comparison": center_rate_comparison,
        "rate_a_10kb_null_exceedance_probability_clopper_pearson_95": exceedance_ci,
        "all_four_cells_are_independent_complete_100_null_calibrations": True,
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
        rate_a_1kb,
        rate_a_10kb,
        rate_b_1kb,
        rate_b_10kb,
        analysis_dir / "corrected_mutation_rate_comparison.png",
    )


if __name__ == "__main__":
    main()
