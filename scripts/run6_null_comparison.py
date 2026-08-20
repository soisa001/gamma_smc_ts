#!/usr/bin/env python
"""Score one observation against four candidate nulls.

The hom-carrier statistic is only as good as what it is compared against, so
this puts the four plausible nulls side by side on the *same* observation --
run5's CHB selected replicates, hom-carrier ``P(TMRCA < x)`` at the focal base:

1. ``chb_unconditional``   -- CHB's own demography, no allele, unstratified.
2. ``eas_unconditional``   -- the EAS demography, no allele, unstratified. This
   is the "borrow another population's fitted model" proposal, and it carries a
   demographic mismatch that has nothing to do with selection.
3. ``snp_hom_matched``     -- hom for a neutral SNP of matched frequency. Correct
   in spirit but the wrong ascertainment: every copy of a SNP descends from one
   mutation, so this null is biased *young*.
4. ``introgressed_hom_matched`` -- hom for a neutrally introgressed tract at
   matched introgression frequency. This is the ascertainment the real analysis
   actually performs, and it is the null the claim has to clear.

Each replicate is scored against the null bin matching its own final frequency,
so a replicate that swept to 80% is compared with neutral loci that drifted to
80% rather than with the genome-wide average.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run6_scan import build_matched_null, load_replicates  # noqa: E402

THRESHOLDS_YEARS = (10_000.0, 20_000.0, 30_000.0, 50_000.0)
FREQUENCY_EDGES = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.01]
MIN_CLASS_PAIRS = 10
ALPHA = 0.05


def observed_hom_carrier(arm_dir: Path, generation_time: float) -> pd.DataFrame:
    """Per replicate: final frequency and hom-carrier P(TMRCA<x) at the focal base."""
    rows = []
    for replicate in load_replicates(arm_dir, "selected"):
        if replicate.genotypes is None:
            continue
        mask = replicate.genotypes == 2
        n_hom = int(mask.sum())
        if n_hom < MIN_CLASS_PAIRS:
            continue
        focal = int(np.argmin(np.abs(replicate.positions - 5_000_000)))
        frequency = float(replicate.genotypes.sum()) / (2 * replicate.genotypes.size)
        for years in THRESHOLDS_YEARS:
            generations = years / generation_time
            rows.append(
                {
                    "replicate_index": replicate.replicate_index,
                    "sample_frequency": frequency,
                    "n_hom_pairs": n_hom,
                    "threshold_years": years,
                    "observed": float(
                        np.mean(replicate.tmrca[focal, mask] < generations)
                    ),
                }
            )
    return pd.DataFrame(rows)


def unconditional_null(arm_dir: Path, generation_time: float, n_pairs: int, rng):
    out = {}
    neutral = load_replicates(arm_dir, "neutral")
    for years in THRESHOLDS_YEARS:
        out[years] = build_matched_null(
            neutral, years / generation_time, n_pairs, rng=rng
        )
    return out


def tail_p(observed: float, null: np.ndarray) -> float:
    """Upper-tail p-value with the conventional +1 correction."""
    if null.size == 0:
        return float("nan")
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(REPO / "sim_results_run6" / "null_comparison"))
    args = parser.parse_args(argv)

    rng = np.random.default_rng(20260825)
    chb = Path("sim_results_run5/chb_from_introgression")
    eas = Path("sim_results_run6/eas_standing")

    observed = observed_hom_carrier(chb, 29.0)
    observed["frequency_bin"] = pd.cut(
        observed["sample_frequency"], FREQUENCY_EDGES, right=False
    )

    # Null 1 and 2: unstratified neutral, drawn at each replicate's own pair count.
    sizes = sorted(observed["n_hom_pairs"].unique())
    chb_null = {n: unconditional_null(chb, 29.0, n, rng) for n in sizes}
    eas_null = {n: unconditional_null(eas, 25.0, n, rng) for n in sizes}

    # Null 3: hom for a matched-frequency neutral SNP (EAS demography).
    snp = pd.read_csv(
        "sim_results_run6/conditioned_null/conditioned_null_raw.tsv", sep="\t"
    )
    snp["frequency_bin"] = pd.cut(
        snp["derived_frequency"], FREQUENCY_EDGES, right=False
    )

    # Null 4: hom for a neutrally introgressed tract at matched frequency (CHB).
    intro = pd.read_csv(
        "sim_results_run6/introgressed_null/introgressed_null_raw.tsv", sep="\t"
    )
    intro["frequency_bin"] = pd.cut(
        intro["introgression_frequency"], FREQUENCY_EDGES, right=False
    )

    records = []
    for _, row in observed.iterrows():
        years = row["threshold_years"]
        n = int(row["n_hom_pairs"])
        value = row["observed"]
        candidates = {
            "chb_unconditional": chb_null[n][years],
            "eas_unconditional": eas_null[n][years],
            "snp_hom_matched": snp.loc[
                (snp["threshold_years"] == years)
                & (snp["frequency_bin"] == row["frequency_bin"]),
                "p_hom_derived",
            ].to_numpy(),
            "introgressed_hom_matched": intro.loc[
                (intro["threshold_years"] == years)
                & (intro["frequency_bin"] == row["frequency_bin"]),
                "p_hom_introgressed",
            ].to_numpy(),
        }
        for name, null in candidates.items():
            records.append(
                {
                    "replicate_index": row["replicate_index"],
                    "threshold_years": years,
                    "sample_frequency": row["sample_frequency"],
                    "null": name,
                    "null_size": int(null.size),
                    "null_mean": float(np.mean(null)) if null.size else np.nan,
                    "observed": value,
                    "p_value": tail_p(value, np.asarray(null, dtype=float)),
                }
            )
    table = pd.DataFrame(records)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    observed.to_csv(out / "observed_hom_carrier.tsv", sep="\t", index=False)
    table.to_csv(out / "null_comparison_raw.tsv", sep="\t", index=False)

    summary = (
        table.groupby(["threshold_years", "null"])
        .agg(
            n_replicates=("p_value", "size"),
            median_null_size=("null_size", "median"),
            mean_null=("null_mean", "mean"),
            mean_observed=("observed", "mean"),
            power=("p_value", lambda s: float(np.mean(s < ALPHA))),
            median_p=("p_value", "median"),
        )
        .reset_index()
    )
    summary.to_csv(out / "null_comparison_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(json.dumps({"observed_replicates": int(observed["replicate_index"].nunique())}, indent=2))
    print()
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
