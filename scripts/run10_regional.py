#!/usr/bin/env python
"""Does averaging over a region rescue detection at realistic frequencies?

The variance diagnostic showed the null's width at 30-50 kya is genealogical
rather than sampling noise, so a larger panel cannot help. Regional averaging
attacks the *other* term: trees decorrelate along the sequence on a scale of
roughly ``1 / r``, so a window spanning many recombination-independent
genealogies averages over independent draws of exactly the variance that
dominates.

Two levers are measured on the same replicates:

* **window width** -- the statistic averaged over 10 kb (a single stride) up to
  4 Mb, centred on the focal site.
* **carrier conditioning** -- overall against carrier-carrier pairs.

Everything is scored against the neutral arm treated identically, so a wider
window is compared with a wider-window null and the comparison stays honest.

The frequency band matters more than the arm here, so replicates are pooled and
binned on final allele frequency, with the 0.25-0.45 band called out: that is
what real introgressed regions actually reach.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from gamma_smc_aou.run8_analysis import load_arm  # noqa: E402
from gamma_smc_aou.run9_config import (  # noqa: E402
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    SCENARIOS,
)

THRESHOLDS = (10_000.0, 30_000.0, 50_000.0)
WINDOWS_BP = (10_000, 100_000, 500_000, 1_000_000, 2_000_000, 4_000_000)
CLASSES = ("overall", "carrier_carrier")
AF_EDGES = (0.0, 0.15, 0.25, 0.45, 0.60, 1.01)
MIN_PAIRS = 50


def windowed(profile: np.ndarray, positions: np.ndarray, width: int) -> float:
    """Mean of the per-position statistic over a window centred on the focal site."""
    half = width / 2.0
    mask = np.abs(positions - FOCAL_POSITION_BP) <= half
    if not mask.any():
        mask = np.argmin(np.abs(positions - FOCAL_POSITION_BP))
        return float(profile[mask])
    return float(profile[mask].mean())


def collect(study_root: str | Path) -> pd.DataFrame:
    rows = []
    arms = [("neutral", None)] + [(s.scenario_id, s.scenario_id) for s in SCENARIOS]
    for arm_label, scenario_id in arms:
        for replicate in load_arm(study_root, scenario_id):
            positions = replicate.positions
            for years in THRESHOLDS:
                generations = years / GENERATION_TIME_YEARS
                below = replicate.tmrca < generations  # (positions, pairs)
                for name in CLASSES:
                    mask = replicate.class_mask(name)
                    if int(mask.sum()) < MIN_PAIRS:
                        continue
                    profile = below[:, mask].mean(axis=1)
                    for width in WINDOWS_BP:
                        rows.append(
                            {
                                "arm": arm_label,
                                "replicate_index": replicate.replicate_index,
                                "sample_af": replicate.sample_af,
                                "threshold_years": years,
                                "genotype_class": name,
                                "window_bp": width,
                                "statistic": windowed(profile, positions, width),
                            }
                        )
    return pd.DataFrame(rows)


def score(table: pd.DataFrame) -> pd.DataFrame:
    """AUC and power against the neutral arm scored the same way."""
    neutral = table[table["arm"] == "neutral"]
    selected = table[table["arm"] != "neutral"].copy()
    selected["af_bin"] = pd.cut(selected["sample_af"], AF_EDGES, right=False)

    rows = []
    for (years, name, width, af_bin), block in selected.groupby(
        ["threshold_years", "genotype_class", "window_bp", "af_bin"], observed=True
    ):
        # The null always uses the *overall* class: the neutral arm has no
        # carriers, so a carrier-matched null does not exist here. Carrier rows
        # are therefore optimistic and bound headroom.
        null = neutral[
            (neutral["threshold_years"] == years)
            & (neutral["window_bp"] == width)
            & (neutral["genotype_class"] == "overall")
        ]["statistic"].to_numpy(dtype=float)
        if null.size == 0 or block.empty:
            continue
        values = block["statistic"].to_numpy(dtype=float)
        exceed = np.array([int(np.sum(null >= v)) for v in values])
        percentile = np.array(
            [float(np.mean(null < v)) + 0.5 * float(np.mean(null == v)) for v in values]
        )
        rows.append(
            {
                "threshold_years": years,
                "genotype_class": name,
                "window_bp": width,
                "af_bin": str(af_bin),
                "n": int(values.size),
                "mean_statistic": float(values.mean()),
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)),
                "auc": float(percentile.mean()),
                "power": float(np.mean(exceed <= 1)),
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", default=str(REPO / "sim_results_run9"))
    args = parser.parse_args(argv)

    table = collect(args.study_root)
    out = Path(args.study_root) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "regional_raw.tsv", sep="\t", index=False)
    summary = score(table)
    summary.to_csv(out / "regional_summary.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    band = "[0.25, 0.45)"
    print("#" * 90)
    print(f"REALISTIC BAND  AF in {band}  — what real introgressed regions reach")
    print("#" * 90)
    for name in CLASSES:
        part = summary[
            (summary["af_bin"] == band) & (summary["genotype_class"] == name)
        ]
        if part.empty:
            continue
        print(f"\n=== {name} : power ===")
        print(
            part.pivot_table(index="window_bp", columns="threshold_years", values="power")
            .to_string(float_format=lambda x: f"{x:.3f}")
        )
        print(f"\n=== {name} : AUC ===")
        print(
            part.pivot_table(index="window_bp", columns="threshold_years", values="auc")
            .to_string(float_format=lambda x: f"{x:.3f}")
        )
    print()
    print("#" * 90)
    print("NULL SD BY WINDOW (does averaging shrink the genealogical variance?)")
    print("#" * 90)
    part = summary[summary["genotype_class"] == "overall"].drop_duplicates(
        ["threshold_years", "window_bp"]
    )
    print(
        part.pivot_table(index="window_bp", columns="threshold_years", values="null_sd")
        .to_string(float_format=lambda x: f"{x:.4f}")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
