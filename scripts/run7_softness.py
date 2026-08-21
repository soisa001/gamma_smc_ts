#!/usr/bin/env python
"""Why stronger selection gives a *weaker* recent-TMRCA signal.

The suspicion is that the effect is about sweep softness, not sweep strength. An
allele under weak selection lingers at low frequency for many generations, and
while the carrier pool is small its lineages coalesce with each other, so the
sample of carriers ends up descending from few founding archaic haplotypes.
Strong selection pulls the allele out of the low-frequency regime quickly, so
little of that coalescence happens and many founders survive into the present --
a softer sweep, whose carriers meet only back in the archaic branch.

Two measurements test that directly, both on hom-carrier pairs at the focal base:

* the fraction coalescing **after the pulse**, which is the fraction that traces
  to a shared founder rather than to two different archaic lineages;
* the median TMRCA, which says where the mass actually sits.

If softness is the mechanism, the post-pulse fraction should fall as s rises.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run7_analysis import ARM_SPECS, load_arm, load_mode  # noqa: E402
from gamma_smc_aou.run7_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    PULSE_GENERATIONS,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", default=str(REPO / "sim_results_run7"))
    args = parser.parse_args(argv)

    rows = []
    specs = list(ARM_SPECS) + [("neutral", 0.0, None)]
    for origin, coefficient, arm in specs:
        replicates = (
            load_mode(args.study_root, None)
            if arm is None
            else load_arm(args.study_root, arm)
        )
        for replicate in replicates:
            focal = int(np.argmin(np.abs(replicate.positions - FOCAL_POSITION_BP)))
            if replicate.genotypes is None or origin == "neutral":
                mask = np.ones(replicate.tmrca.shape[1], dtype=bool)
            else:
                mask = replicate.genotypes == 2
            if mask.sum() < 10:
                continue
            tmrca = replicate.tmrca[focal, mask].astype(float)
            rows.append(
                {
                    "origin": origin,
                    "selection_coefficient": coefficient,
                    "replicate_index": replicate.replicate_index,
                    "sample_af": replicate.sample_af,
                    "n_pairs": int(mask.sum()),
                    "post_pulse_fraction": float(np.mean(tmrca < PULSE_GENERATIONS)),
                    "median_tmrca_gen": float(np.median(tmrca)),
                }
            )

    table = pd.DataFrame(rows)
    out = Path(args.study_root) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "sweep_softness.tsv", sep="\t", index=False)

    summary = (
        table.groupby(["origin", "selection_coefficient"])
        .agg(
            n=("post_pulse_fraction", "size"),
            mean_af=("sample_af", "mean"),
            mean_pairs=("n_pairs", "mean"),
            post_pulse_fraction=("post_pulse_fraction", "mean"),
            median_tmrca_gen=("median_tmrca_gen", "median"),
        )
        .reset_index()
    )
    summary["median_tmrca_ky"] = (
        summary["median_tmrca_gen"] * GENERATION_TIME_YEARS / 1000
    )
    summary.to_csv(out / "sweep_softness_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 220)
    print(
        f"pulse = {PULSE_GENERATIONS:,.0f} generations "
        f"({PULSE_GENERATIONS * GENERATION_TIME_YEARS / 1000:.0f} kya)\n"
    )
    print(
        summary[
            [
                "origin",
                "selection_coefficient",
                "n",
                "mean_af",
                "mean_pairs",
                "post_pulse_fraction",
                "median_tmrca_ky",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
