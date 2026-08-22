#!/bin/bash
SRC=/home/mew/AllOfUs_Phase2/gamma_smc_ts/sim_results_run9
DST=/mnt/c/Users/Lenovo/AllOfUs_Phase2/gamma_smc_ts/sim_results_run9
mkdir -p "$DST/figures" "$DST/tables" "$DST/config" "$DST/logs"
cp -f "$SRC"/figures/*.png "$DST/figures/" 2>/dev/null
cp -f "$SRC"/config/*.json "$DST/config/" 2>/dev/null
cp -f "$SRC"/logs/*.tsv "$DST/logs/" 2>/dev/null
# scan.tsv is 35 MB of per-position profile; the rest are small enough to track.
for f in af_stratified exceedance quantiles; do
  cp -f "$SRC/tables/$f.tsv" "$DST/tables/" 2>/dev/null
done
du -sh "$DST"; find "$DST" -type f | wc -l
