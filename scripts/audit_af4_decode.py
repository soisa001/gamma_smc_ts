#!/usr/bin/env python3
"""Independently audit the exact s=0.01, AF approximately 4% decode."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from gamma_smc_aou.calibration import bh_fdr


CALLED_COLUMN = "frac_recent_4500"
DEFINITIONS = (
    (
        "soft_probability",
        "mean_p_tmrca_lt_threshold",
        "selected_decoded_mean_call_profile.tsv",
        "neutral_decoded_mean_call_profiles.tsv.gz",
        "comparison_soft_probability_spatial_calibration.tsv",
    ),
    (
        "posterior_mean_call",
        CALLED_COLUMN,
        "selected_decoded_mean_call_profile.tsv",
        "neutral_decoded_mean_call_profiles.tsv.gz",
        "comparison_posterior_mean_call_spatial_calibration.tsv",
    ),
    (
        "posterior_median_call",
        CALLED_COLUMN,
        "selected_decoded_median_call_profile.tsv",
        "neutral_decoded_median_call_profiles.tsv.gz",
        "comparison_posterior_median_call_spatial_calibration.tsv",
    ),
)


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _center(frame: pd.DataFrame, position: float = 5_000_000) -> pd.Series:
    return frame.iloc[
        int(
            np.argmin(
                np.abs(frame["position_0based"].to_numpy(dtype=float) - position)
            )
        )
    ]


def _density_counts(
    scan: pd.DataFrame,
    neutral: pd.DataFrame,
    statistic_column: str,
    *,
    center: float = 5_000_000,
    local_half_width: float = 500_000,
) -> dict:
    positions = scan["position_0based"].to_numpy(dtype=float)
    pivot = (
        neutral.loc[neutral["position_0based"].isin(positions)]
        .pivot(
            index="replicate",
            columns="position_0based",
            values=statistic_column,
        )
        .reindex(columns=positions)
    )
    if pivot.isna().any().any():
        raise RuntimeError("neutral profiles do not contain a complete scan grid")
    n_null = len(pivot)
    ranks = pivot.rank(axis=0, method="min").to_numpy(dtype=float)
    leave_one_out_p = (1 + n_null - ranks) / n_null
    local = np.abs(positions - center) <= local_half_width
    observed_p = scan["p_upper"].to_numpy(dtype=float)
    rows = {}
    for label, mask in (("global", np.ones(len(positions), dtype=bool)), ("local", local)):
        observed_count = int(np.count_nonzero(observed_p[mask] < 0.05))
        null_counts = np.count_nonzero(leave_one_out_p[:, mask] < 0.05, axis=1)
        rows[label] = {
            "n_windows": int(np.count_nonzero(mask)),
            "selected_positive_windows": observed_count,
            "null_mean_positive_windows": float(null_counts.mean()),
            "null_median_positive_windows": float(np.median(null_counts)),
            "null_exceedances": int(np.count_nonzero(null_counts >= observed_count)),
            "monte_carlo_p_upper": float(
                (1 + np.count_nonzero(null_counts >= observed_count))
                / (1 + n_null)
            ),
        }
    rows["local"]["half_width_bp"] = int(local_half_width)
    return rows


def _assert_contract(source_metrics: dict, decode_metrics: dict) -> None:
    expected_source = {
        "selection_coefficient": 0.01,
        "mutation_rate": 1.29e-8,
        "minimum_population_allele_frequency": 0.04,
        "sample_diploids": 2_000,
        "sequence_length": 10_000_000,
        "variant_age_generations": 180,
    }
    for key, expected in expected_source.items():
        observed = source_metrics[key]
        if isinstance(expected, float):
            if not np.isclose(float(observed), expected, rtol=0, atol=1e-20):
                raise RuntimeError(f"source {key}={observed!r}; expected {expected!r}")
        elif observed != expected:
            raise RuntimeError(f"source {key}={observed!r}; expected {expected!r}")
    rejection = source_metrics["selected_rejection"]
    if int(rejection["accepted_attempt_zero_based"]) != 4_141:
        raise RuntimeError("source is not deterministic selected attempt 4141")
    if not np.isclose(
        float(rejection["population_allele_frequency"]),
        0.040225,
        rtol=0,
        atol=1e-12,
    ):
        raise RuntimeError("source population AF is not 0.040225")

    expected_decode = {
        "stride_bp": 10_000,
        "neutral_replicates": 100,
        "neutral_workers": 12,
        "recent_call": "mean",
        "comparison_recent_call": "median",
        "calibration_statistic": "called-fraction",
        "cache_size_bp": 1_000,
        "decoder_threads_per_replicate": 1,
    }
    for key, expected in expected_decode.items():
        if decode_metrics[key] != expected:
            raise RuntimeError(
                f"decode {key}={decode_metrics[key]!r}; expected {expected!r}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    source = root / "source"
    decoded = root / "stride10kb"
    output = root / "analysis"
    output.mkdir(parents=True, exist_ok=True)

    source_metrics = _read_json(source / "metrics.json")
    decode_metrics = _read_json(decoded / "decoded_study_metrics.json")
    _assert_contract(source_metrics, decode_metrics)

    selected_mean = pd.read_csv(
        decoded / "selected_decoded_mean_call_profile.tsv", sep="\t"
    )
    selected_median = pd.read_csv(
        decoded / "selected_decoded_median_call_profile.tsv", sep="\t"
    )
    neutral_mean = pd.read_csv(
        decoded / "neutral_decoded_mean_call_profiles.tsv.gz", sep="\t"
    )
    neutral_median = pd.read_csv(
        decoded / "neutral_decoded_median_call_profiles.tsv.gz", sep="\t"
    )
    soft_columns = [
        "position_0based",
        "n_pairs",
        "mean_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
    ]
    if not np.allclose(
        selected_mean[soft_columns].to_numpy(dtype=float),
        selected_median[soft_columns].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise RuntimeError("selected soft posterior outputs differ by call rule")
    neutral_sort = ["replicate", "position_0based"]
    neutral_mean = neutral_mean.sort_values(neutral_sort).reset_index(drop=True)
    neutral_median = neutral_median.sort_values(neutral_sort).reset_index(drop=True)
    if not np.allclose(
        neutral_mean[[*neutral_sort, *soft_columns[1:]]].to_numpy(dtype=float),
        neutral_median[[*neutral_sort, *soft_columns[1:]]].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise RuntimeError("neutral soft posterior outputs differ by call rule")

    frames = {
        "selected_decoded_mean_call_profile.tsv": selected_mean,
        "selected_decoded_median_call_profile.tsv": selected_median,
        "neutral_decoded_mean_call_profiles.tsv.gz": neutral_mean,
        "neutral_decoded_median_call_profiles.tsv.gz": neutral_median,
    }
    rows = []
    density_results = {}
    for slug, statistic, selected_name, neutral_name, scan_name in DEFINITIONS:
        selected = frames[selected_name]
        neutral = frames[neutral_name]
        scan = pd.read_csv(decoded / scan_name, sep="\t")
        selected_center = _center(selected)
        neutral_center = neutral.loc[
            neutral["position_0based"].eq(selected_center["position_0based"]),
            statistic,
        ].to_numpy(dtype=float)
        exceedances = int(
            np.count_nonzero(neutral_center >= float(selected_center[statistic]))
        )
        p_center = float((1 + exceedances) / (1 + len(neutral_center)))
        scan_center = _center(scan)
        if not np.isclose(p_center, float(scan_center["p_upper"]), rtol=0, atol=1e-12):
            raise RuntimeError(f"{slug} center p-value does not reproduce")
        q_recomputed = bh_fdr(scan["p_upper"].to_numpy(dtype=float))
        if not np.allclose(
            q_recomputed,
            scan["q_bh"].to_numpy(dtype=float),
            rtol=0,
            atol=1e-12,
        ):
            raise RuntimeError(f"{slug} BH q-values do not reproduce")
        density = _density_counts(scan, neutral, statistic)
        density_results[slug] = density
        rows.append(
            {
                "statistic": slug,
                "center_selected": float(selected_center[statistic]),
                "center_null_mean": float(neutral_center.mean()),
                "center_null_ci95_lower": float(np.quantile(neutral_center, 0.025)),
                "center_null_ci95_upper": float(np.quantile(neutral_center, 0.975)),
                "center_neutral_exceedances": exceedances,
                "center_monte_carlo_p_upper": p_center,
                "n_pointwise_p_lt_0_05_windows": int(
                    np.count_nonzero(scan["p_upper"] < 0.05)
                ),
                "n_bh_q_lt_0_05_windows": int(
                    np.count_nonzero(q_recomputed < 0.05)
                ),
                "global_selected_positive_windows": density["global"][
                    "selected_positive_windows"
                ],
                "global_density_monte_carlo_p_upper": density["global"][
                    "monte_carlo_p_upper"
                ],
                "local_selected_positive_windows": density["local"][
                    "selected_positive_windows"
                ],
                "local_density_monte_carlo_p_upper": density["local"][
                    "monte_carlo_p_upper"
                ],
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(output / "significance_summary.tsv", sep="\t", index=False)
    truth = source_metrics["selected_rejection"]
    audit = {
        "study": "exact s=0.01 AF=0.040225 selected attempt 4141",
        "contract_validated": True,
        "soft_outputs_identical_between_call_rules": True,
        "selected_rows": int(len(selected_mean)),
        "neutral_rows": int(len(neutral_mean)),
        "neutral_replicates": int(neutral_mean["replicate"].nunique()),
        "truth_center_fraction_recent": float(truth["truth_center_fraction_recent"]),
        "truth_center_monte_carlo_p_upper": float(
            truth["truth_center_monte_carlo_p_upper"]
        ),
        "density": density_results,
        "hashes": {
            str(path.relative_to(root)): _sha256(path)
            for path in (
                source / "selected.trees",
                source / "selected.vcf.gz",
                decoded / "selected_decoded_mean_call_profile.tsv",
                decoded / "selected_decoded_median_call_profile.tsv",
                decoded / "neutral_decoded_mean_call_profiles.tsv.gz",
                decoded / "neutral_decoded_median_call_profiles.tsv.gz",
            )
        },
    }
    (output / "audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
