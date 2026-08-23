#!/usr/bin/env python
"""Selected introgressed regions against the frequency-matched neutral null.

This is the comparison the whole series has been building toward. Earlier
carrier-carrier results were scored against a neutral arm with no carriers in it,
which cannot separate "selected" from "introgressed". Here the null is neutral
loci whose archaic ancestry drifted to the *same* frequency, scored on the *same*
statistic, so anything that survives is attributable to selection.

Only the pulse arms are compared. The de novo arms define their carrier class by
a modern allele rather than an archaic tract, so an archaic-ancestry-matched null
is not the right comparator for them -- theirs would be hom-derived pairs at a
matched-frequency neutral SNP, which is a different construction.

Matching is on final frequency because that is what an empirical analysis can
observe: for a pulse arm the selected allele rides archaic haplotypes, so its
sample frequency and the local archaic ancestry frequency are the same quantity.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

AF_EDGES = (0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.71)
THRESHOLDS = (10_000.0, 30_000.0, 50_000.0)
PULSE_ARMS = (
    "pulse_s0p003",
    "pulse_s0p004",
    "pulse_s0p005",
    "pulse_delayed_s0p005",
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run9", default=str(REPO / "sim_results_run9"))
    parser.add_argument("--run10", default=str(REPO / "sim_results_run10"))
    args = parser.parse_args(argv)

    null = pd.read_csv(Path(args.run10) / "tables" / "matched_null_raw.tsv", sep="\t")
    null["af_bin"] = pd.cut(null["archaic_frequency"], AF_EDGES, right=False)

    # run9's per-replicate carrier-carrier statistic at the single 10 kb stride.
    selected = pd.read_csv(Path(args.run9) / "tables" / "regional_raw.tsv", sep="\t")
    selected = selected[
        (selected["window_bp"] == 10_000)
        & (selected["genotype_class"] == "carrier_carrier")
        & (selected["arm"].isin(PULSE_ARMS))
    ].copy()
    selected["af_bin"] = pd.cut(selected["sample_af"], AF_EDGES, right=False)

    rows = []
    for (af_bin, years), block in selected.groupby(
        ["af_bin", "threshold_years"], observed=True
    ):
        reference = null[
            (null["af_bin"] == af_bin) & (null["threshold_years"] == years)
        ]["carrier_carrier"].to_numpy(dtype=float)
        if reference.size < 30 or block.empty:
            continue
        values = block["statistic"].to_numpy(dtype=float)
        exceed = np.array([int(np.sum(reference >= v)) for v in values])
        p = (exceed + 1.0) / (reference.size + 1.0)
        percentile = np.array(
            [
                float(np.mean(reference < v)) + 0.5 * float(np.mean(reference == v))
                for v in values
            ]
        )
        rows.append(
            {
                "af_bin": str(af_bin),
                "threshold_years": years,
                "n_selected": int(values.size),
                "n_null_loci": int(reference.size),
                "selected_mean": float(values.mean()),
                "null_mean": float(reference.mean()),
                "null_q95": float(np.quantile(reference, 0.95)),
                "auc": float(percentile.mean()),
                "power_p01": float(np.mean(p <= 0.01)),
                "power_p05": float(np.mean(p <= 0.05)),
                "median_p": float(np.median(p)),
            }
        )

    table = pd.DataFrame(rows)
    out = Path(args.run10) / "tables"
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "headtohead.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print("#" * 96)
    print("SELECTED PULSE ARMS vs FREQUENCY-MATCHED NEUTRAL INTROGRESSION")
    print("  carrier-carrier statistic, both sides; null = neutral loci at the same")
    print("  archaic frequency, so surviving signal is attributable to selection")
    print("#" * 96)
    print()
    print(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print()
    print("=" * 96)
    print("HEADLINE: AUC by frequency band and cutoff")
    print("=" * 96)
    print(
        table.pivot_table(index="af_bin", columns="threshold_years", values="auc")
        .to_string(float_format=lambda x: f"{x:.3f}")
    )
    print()
    print("=" * 96)
    print("HEADLINE: power at p <= 0.01")
    print("=" * 96)
    print(
        table.pivot_table(index="af_bin", columns="threshold_years", values="power_p01")
        .to_string(float_format=lambda x: f"{x:.3f}")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
