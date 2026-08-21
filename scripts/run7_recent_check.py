#!/usr/bin/env python
"""P(TMRCA < x) over a fine, recent x grid for the 5 kya sweep arms.

run7's standard grid starts at 1,000 years and jumps to 4,500, which is far too
coarse for a sweep that began 200 generations ago. This scores the recent arms
against the same shared neutral null over a grid that actually brackets the
sweep, and prints the introgressed and 55 kya de novo arms alongside so the
regimes can be compared on identical cutoffs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run6_scan import build_matched_null  # noqa: E402
from gamma_smc_aou.run7_analysis import (  # noqa: E402
    ALPHA,
    MIN_CLASS_PAIRS,
    load_arm,
    load_mode,
)
from gamma_smc_aou.run7_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
)

FINE_THRESHOLDS_YEARS = (
    500.0,
    1_000.0,
    2_000.0,
    3_000.0,
    5_000.0,
    7_500.0,
    10_000.0,
    15_000.0,
    20_000.0,
    30_000.0,
    50_000.0,
)

ARMS = (
    ("recent 5 kya, s=0.05", "recent5k_s0p05"),
    ("recent 5 kya, s=0.1", "recent5k_s0p1"),
    ("de novo 55 kya, s=0.05", "denovo_s0p05"),
    ("introgressed 55 kya, s=0.003", "eas_introgressed_s0p003"),
)


def tail_p(observed: float, null: np.ndarray) -> float:
    if null.size == 0:
        return float("nan")
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", default=str(REPO / "sim_results_run7"))
    parser.add_argument("--genotype-class", default="unconditional",
                        choices=("unconditional", "hom_carrier"))
    args = parser.parse_args(argv)

    neutral = load_mode(args.study_root, None)
    rng = np.random.default_rng(20261302)
    null_cache: dict[tuple[float, int], np.ndarray] = {}

    rows = []
    for label, arm in ARMS:
        replicates = load_arm(args.study_root, arm)
        if not replicates:
            print(f"  (no replicates for {arm}, skipping)")
            continue
        for replicate in replicates:
            focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
            if args.genotype_class == "unconditional" or replicate.genotypes is None:
                mask = np.ones(replicate.tmrca.shape[1], dtype=bool)
            else:
                mask = replicate.genotypes == 2
            size = int(mask.sum())
            if size < MIN_CLASS_PAIRS:
                continue
            for years in FINE_THRESHOLDS_YEARS:
                generations = float(years) / GENERATION_TIME_YEARS
                value = float(np.mean(replicate.tmrca[focal, mask] < generations))
                key = (float(years), size)
                if key not in null_cache:
                    null_cache[key] = build_matched_null(
                        neutral, generations, size, rng=rng
                    )
                null = null_cache[key]
                rows.append(
                    {
                        "arm": label,
                        "arm_id": arm,
                        "replicate_index": replicate.replicate_index,
                        "sample_af": replicate.sample_af,
                        "threshold_years": float(years),
                        "statistic": value,
                        "null_mean": float(null.mean()) if null.size else np.nan,
                        "percentile": float(np.mean(null < value))
                        + 0.5 * float(np.mean(null == value)),
                        "p_value": tail_p(value, null),
                    }
                )

    table = pd.DataFrame(rows)
    out = Path(args.study_root) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(
        out / f"recent_check_{args.genotype_class}.tsv", sep="\t", index=False
    )

    summary = (
        table.groupby(["arm", "threshold_years"])
        .agg(
            n=("p_value", "size"),
            mean_af=("sample_af", "mean"),
            observed=("statistic", "mean"),
            null=("null_mean", "mean"),
            auc=("percentile", "mean"),
            power=("p_value", lambda s: float(np.mean(s < ALPHA))),
        )
        .reset_index()
    )
    summary.to_csv(
        out / f"recent_check_{args.genotype_class}_summary.tsv", sep="\t", index=False
    )

    pd.set_option("display.width", 250)
    print(f"=== genotype class: {args.genotype_class} ===\n")
    for metric in ("observed", "auc", "power"):
        print(f"--- {metric} ---")
        print(
            summary.pivot_table(
                index="threshold_years", columns="arm", values=metric
            ).to_string(float_format=lambda x: f"{x:.3f}")
        )
        print()
    print("--- mean final AF ---")
    print(
        summary.groupby("arm")["mean_af"].first().to_string(
            float_format=lambda x: f"{x:.3f}"
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
