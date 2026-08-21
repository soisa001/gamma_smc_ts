#!/bin/bash
SRC=/home/mew/AllOfUs_Phase2/gamma_smc_ts/sim_results_run8
DST=/mnt/c/Users/Lenovo/AllOfUs_Phase2/gamma_smc_ts/sim_results_run8
mkdir -p "$DST/figures" "$DST/tables" "$DST/config" "$DST/logs"
cp -f "$SRC"/figures/*.png "$DST/figures/" 2>/dev/null
cp -f "$SRC"/tables/*.tsv  "$DST/tables/"  2>/dev/null
cp -f "$SRC"/config/*.json "$DST/config/"  2>/dev/null
cp -f "$SRC"/logs/*.tsv    "$DST/logs/"    2>/dev/null
du -sh "$DST"; find "$DST" -type f | wc -l
