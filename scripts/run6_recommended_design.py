#!/usr/bin/env python
"""Score the design that the study can actually run.

The intended analysis compares *empirical* windows against a *simulated* neutral
null drawn from the EAS history inferred on those same samples. That null model
contains no archaic admixture, so it cannot produce a homozygous-introgressed
class at all -- which means the hom-carrier statistic has an ascertainment the
null cannot match, and any comparison against an unconditional null credits the
ascertainment as though it were selection.

Two escapes exist, and this script scores both on the same replicates:

``conditioned``
    Keep the hom-carrier statistic and pay for it with a matched-frequency
    conditioned null. Correct only to the extent the null's conditioning matches
    the real one; ``run6_conditioned_null`` and ``run6_introgressed_null``
    measure how large that residual gap is.

``unconditional``
    Drop the conditioning from the *observation* as well: score each window by
    the within-individual ``P(TMRCA < x)`` over all individuals, and compare it
    with the same quantity in neutral simulations. No ascertainment, so no
    conditioning is owed. This works because a swept locus drags the whole
    window's distribution, not merely the carriers'.
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
ALPHA = 0.05
FOCAL = 5_000_000


def tail_p(observed, null):
    if null.size == 0:
        return np.nan
    return float((np.sum(null >= observed) + 1) / (null.size + 1))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(REPO / "sim_results_run6" / "recommended"))
    args = parser.parse_args(argv)

    rng = np.random.default_rng(20260826)
    arms = {
        "chb_from_introgression": (Path("sim_results_run5/chb_from_introgression"), 29.0),
        "eas_standing": (Path("sim_results_run6/eas_standing"), 25.0),
    }
    eas_neutral = load_replicates(arms["eas_standing"][0], "neutral")

    records = []
    for arm_id, (arm_dir, generation_time) in arms.items():
        selected = load_replicates(arm_dir, "selected")
        own_neutral = load_replicates(arm_dir, "neutral")

        for statistic in ("unconditional", "hom_carrier"):
            for replicate in selected:
                if replicate.genotypes is None:
                    continue
                if statistic == "unconditional":
                    mask = np.ones(replicate.genotypes.size, dtype=bool)
                else:
                    mask = replicate.genotypes == 2
                size = int(mask.sum())
                if size < 10:
                    continue
                focal = int(np.argmin(np.abs(replicate.positions - FOCAL)))
                frequency = float(replicate.genotypes.sum()) / (2 * replicate.genotypes.size)

                for years in THRESHOLDS_YEARS:
                    value = float(
                        np.mean(replicate.tmrca[focal, mask] < years / generation_time)
                    )
                    for null_name, neutral, null_gt in (
                        ("own_demography", own_neutral, generation_time),
                        ("eas_simulated", eas_neutral, 25.0),
                    ):
                        null = build_matched_null(
                            neutral, years / null_gt, size, rng=rng
                        )
                        records.append(
                            {
                                "arm": arm_id,
                                "statistic": statistic,
                                "null": null_name,
                                "threshold_years": years,
                                "sample_frequency": frequency,
                                "observed": value,
                                "null_mean": float(null.mean()) if null.size else np.nan,
                                "p_value": tail_p(value, null),
                            }
                        )

    table = pd.DataFrame(records)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "recommended_design_raw.tsv", sep="\t", index=False)

    summary = (
        table.groupby(["arm", "statistic", "null", "threshold_years"])
        .agg(
            n=("p_value", "size"),
            null_mean=("null_mean", "mean"),
            observed_mean=("observed", "mean"),
            power=("p_value", lambda s: float(np.mean(s < ALPHA))),
        )
        .reset_index()
    )
    summary.to_csv(out / "recommended_design_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(json.dumps({"rows": int(len(table))}, indent=2))
    print()
    for arm in summary["arm"].unique():
        print(f"=== {arm} ===")
        block = summary[summary["arm"] == arm]
        print(
            block.pivot_table(
                index="threshold_years",
                columns=["statistic", "null"],
                values="power",
            ).to_string(float_format=lambda x: f"{x:.3f}")
        )
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
