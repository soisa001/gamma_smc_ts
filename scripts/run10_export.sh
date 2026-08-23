#!/bin/bash
SRC=/home/mew/AllOfUs_Phase2/gamma_smc_ts
DST=/mnt/c/Users/Lenovo/AllOfUs_Phase2/gamma_smc_ts
mkdir -p "$DST/sim_results_run10/figures" "$DST/sim_results_run10/tables" "$DST/sim_results_run9/tables"
cp -f "$SRC"/sim_results_run10/figures/*.png "$DST/sim_results_run10/figures/" 2>/dev/null
for f in matched_null_summary headtohead; do
  cp -f "$SRC/sim_results_run10/tables/$f.tsv" "$DST/sim_results_run10/tables/" 2>/dev/null
done
# run9 diagnostics produced by run10 scripts
for f in variance_decomposition regional_summary; do
  cp -f "$SRC/sim_results_run9/tables/$f.tsv" "$DST/sim_results_run9/tables/" 2>/dev/null
done
du -sh "$DST/sim_results_run10"; find "$DST/sim_results_run10" -type f | wc -l
