#!/usr/bin/env python
"""Run run8's three checks and write the tables."""
from __future__ import annotations
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))

import numpy as np
import pandas as pd

from gamma_smc_aou.run8_analysis import exceedance_table, focal_table, scan_table
from gamma_smc_aou.run8_config import SCENARIOS

STUDY = REPO / "sim_results_run8"
TABLES = STUDY / "tables"


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    focal = focal_table(STUDY)
    focal.to_csv(TABLES / "focal_table.tsv", sep="\t", index=False)
    exceed = exceedance_table(focal)
    exceed.to_csv(TABLES / "exceedance.tsv", sep="\t", index=False)
    scan = scan_table(STUDY)
    scan.to_csv(TABLES / "scan.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    order = [s.scenario_id for s in SCENARIOS]

    print("#" * 78)
    print("CHECK 1 + 3: UNCONDITIONAL, scored against the unconditioned null")
    print("#" * 78)
    for method in ("truth", "decode"):
        block = exceed[(exceed["method"] == method)
                       & (exceed["genotype_class"] == "overall")]
        if block.empty:
            continue
        print(f"\n--- {method} : AUC ---")
        print(block.pivot_table(index="threshold_years", columns="arm",
                                values="auc").reindex(columns=order)
              .to_string(float_format=lambda x: f"{x:.3f}"))
        print(f"\n--- {method} : power (at most 1 of 100 neutral exceeds) ---")
        print(block.pivot_table(index="threshold_years", columns="arm",
                                values="power_rule").reindex(columns=order)
              .to_string(float_format=lambda x: f"{x:.3f}"))

    print()
    print("#" * 78)
    print("CHECK 2: CARRIER-PAIR STRATIFICATION (truth only), AUC")
    print("#" * 78)
    block = exceed[exceed["method"] == "truth"]
    for arm in order:
        part = block[block["arm"] == arm]
        if part.empty:
            continue
        print(f"\n--- {arm} ---")
        print(part.pivot_table(index="threshold_years", columns="genotype_class",
                               values="auc").to_string(float_format=lambda x: f"{x:.3f}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
