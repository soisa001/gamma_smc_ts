#!/usr/bin/env python
"""Run run9's analysis and write the tables."""
from __future__ import annotations
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np
import pandas as pd

from gamma_smc_aou.run9_analysis import (
    af_stratified_table, exceedance_table, focal_table, quantile_table, scan_table,
)
from gamma_smc_aou.run9_config import SCENARIOS

STUDY = REPO / "sim_results_run9"
TABLES = STUDY / "tables"
ORDER = [s.scenario_id for s in SCENARIOS]


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    focal = focal_table(STUDY); focal.to_csv(TABLES / "focal_table.tsv", sep="\t", index=False)
    exceed = exceedance_table(focal); exceed.to_csv(TABLES / "exceedance.tsv", sep="\t", index=False)
    quant = quantile_table(focal); quant.to_csv(TABLES / "quantiles.tsv", sep="\t", index=False)
    strat = af_stratified_table(focal); strat.to_csv(TABLES / "af_stratified.tsv", sep="\t", index=False)
    scan = scan_table(STUDY); scan.to_csv(TABLES / "scan.tsv", sep="\t", index=False)

    pd.set_option("display.width", 260)
    o = exceed[exceed["genotype_class"] == "overall"]

    for metric, label in (("power_rule", "POWER (<=1 of 300 neutral exceeds)"), ("auc", "AUC")):
        print("#" * 100)
        print(f"UNCONDITIONAL {label}")
        print("#" * 100)
        for m in ("truth", "decode"):
            b = o[o["method"] == m]
            if b.empty: continue
            print(f"\n-- {m} --")
            print(b.pivot_table(index="threshold_years", columns="arm", values=metric)
                   .reindex(columns=ORDER).to_string(float_format=lambda x: f"{x:.3f}"))
        print()

    print("#" * 100)
    print("QUANTILE CURVE: p-value by percentile of the statistic (truth, overall)")
    print("#" * 100)
    q = quant[(quant["method"] == "truth") & (quant["genotype_class"] == "overall")]
    for arm in ORDER:
        part = q[q["arm"] == arm]
        if part.empty: continue
        print(f"\n-- {arm} --")
        print(part.pivot_table(index="threshold_years", columns="quantile",
                               values="p_corrected").to_string(float_format=lambda x: f"{x:.4f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
