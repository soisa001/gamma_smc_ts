#!/usr/bin/env bash
set -euo pipefail
COALESCENT=${COSI2_COALESCENT:-../../external/cosi2/coalescent}
MAXATTEMPTS=${COSI_MAXATTEMPTS:-10000}
test -x "$COALESCENT"
RUN_ROOT=../cosi2_framework_runs
mkdir -p "$RUN_ROOT/smoke"
# Syntax smoke only: the narrow Han boundary gate may reject the one draw.
# This script does not implement the required outer rejection loop.
mkdir -p "$RUN_ROOT/smoke/eas_s0p01"
COSI_MAXATTEMPTS="$MAXATTEMPTS" COSI_SAVE_TRAJ="$RUN_ROOT/smoke/eas_s0p01/trajectory.tsv" "$COALESCENT" -p configs/eas_s0p01.par -n 1 -r 20240524 -u 1 -m -M -e > "$RUN_ROOT/smoke/eas_s0p01/simulation.ms" 2> "$RUN_ROOT/smoke/eas_s0p01/simulation.stderr.log"
uv run python ../../scripts/run_cosi2_framework.py --repo-root ../.. validate-trajectory --cell-id eas_s0p01 --trajectory "$RUN_ROOT/smoke/eas_s0p01/trajectory.tsv"
mkdir -p "$RUN_ROOT/smoke/han_introgressed_s0p01"
COSI_MAXATTEMPTS="$MAXATTEMPTS" COSI_SAVE_TRAJ="$RUN_ROOT/smoke/han_introgressed_s0p01/trajectory.tsv" "$COALESCENT" -p configs/han_introgressed_s0p01.par -n 1 -r 20240624 -u 1 -m -M -e > "$RUN_ROOT/smoke/han_introgressed_s0p01/simulation.ms" 2> "$RUN_ROOT/smoke/han_introgressed_s0p01/simulation.stderr.log"
uv run python ../../scripts/run_cosi2_framework.py --repo-root ../.. validate-trajectory --cell-id han_introgressed_s0p01 --trajectory "$RUN_ROOT/smoke/han_introgressed_s0p01/trajectory.tsv"
