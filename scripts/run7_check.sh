#!/bin/bash
cd /home/mew/AllOfUs_Phase2/gamma_smc_ts || exit 1
.venv/bin/python - <<'PY'
import pandas as pd, numpy as np
f = pd.read_csv("sim_results_run7/tables/focal_table.tsv", sep="\t")
n = f[f["genotype_class"] == "neutral"]
print("=== null calibration: neutral arm scored against its own null ===")
print("  (a calibrated test gives ~uniform p, so FPR at 0.05 should be ~0.05)")
for y in sorted(n["threshold_years"].unique()):
    p = n.loc[n["threshold_years"] == y, "p_value"]
    print("  %7.0f y : n=%3d  FPR@0.05 = %.3f   median p = %.3f   mean p = %.3f"
          % (y, p.size, (p < 0.05).mean(), p.median(), p.mean()))

print()
print("=== final AF and hom-carrier counts by arm ===")
sel = f[(f["genotype_class"] == "unconditional") & (f["threshold_years"] == 30000.0)]
for (origin, s), g in sel.groupby(["origin", "selection_coefficient"]):
    print("  %-14s s=%-6g n=%3d  AF %.3f  hom-carrier diploids %.1f"
          % (origin, s, len(g), g["sample_af"].mean(), g["hom_carrier_diploids"].mean()))

print()
print("=== mean statistic vs null, 30 ky, hom carriers ===")
h = f[(f["genotype_class"] == "hom_carrier") & (f["threshold_years"] == 30000.0)]
for (origin, s), g in h.groupby(["origin", "selection_coefficient"]):
    print("  %-14s s=%-6g  observed %.3f  null %.3f  ratio %.2f"
          % (origin, s, g["statistic"].mean(), g["null_mean"].mean(),
             g["statistic"].mean() / max(g["null_mean"].mean(), 1e-9)))
PY
