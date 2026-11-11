#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PHASE="${1:-all}"
case "$PHASE" in all|merge|compare) ;; *) echo 'Phase must be all, merge, or compare' >&2; exit 2 ;; esac
bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
    if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
  done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
CURRENT="${SV_CURRENT_ROOT:-$HOME/chr1_popped_simple_merge}"
SAME="${SV_SAME_START_ROOT:-$HOME/chr1_popped_same_start}"
REPORT="${SV_SAME_START_REPORT:-$SAME/compare_current}"
test -f "$CURRENT/prepare/done.json" || { echo "Missing completed preparation in $CURRENT" >&2; exit 1; }
if [[ "$PHASE" == all || "$PHASE" == merge ]]; then
  "$UV" run --no-project --python 3.12 scripts/run_popped_simple_sv.py \
    --root "$SAME" --reuse-prepared "$CURRENT" --phase merge --same-start \
    --similarities 0.7 0.9 --threads 4 --memory-gb 26 --bcftools "$BCFTOOLS_BIN"
fi
if [[ "$PHASE" == all || "$PHASE" == compare ]]; then
  "$UV" run --no-project --python 3.12 scripts/compare_popped_merge_modes.py \
    --current-root "$CURRENT" --same-start-root "$SAME" --root "$REPORT" \
    --threads 4 --bcftools "$BCFTOOLS_BIN"
fi
