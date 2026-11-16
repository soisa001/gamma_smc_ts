#!/usr/bin/env bash
set -euo pipefail
umask 077
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
BATCH="${POPPED_FULL_ROOT:-$HOME/popped_full_same_original_pos_p70}"
OUT="${TLR_COMPARE_ROOT:-$BATCH/tlr_before_after}"
export TMPDIR="$OUT/runtime/tmp" UV_CACHE_DIR="$OUT/runtime/uv-cache"
mkdir -p "$TMPDIR" "$UV_CACHE_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Reuse installed binaries read-only; no changes to the running merge checkout.
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
    for base in "$HOME/popped_full_merge_code/gamma_smc_ts" "$HOME/gamma_smc_ts_reports" "$HOME/gamma_smc_ts"; do
        for candidate in "$base/.native/bin/bcftools" "$base/.tools/singleton-native/bin/bcftools"; do
            if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break 2; fi
        done
    done
fi
bash scripts/bootstrap_singleton_tools.sh
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
    for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
        if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
    done
fi
"$REPO/.tools/uv-bin/uv" run --no-project --python 3.12 scripts/compare_tlr_popped_merged.py \
    --batch-root "$BATCH" --out "$OUT" --bcftools "${BCFTOOLS_BIN:-bcftools}" \
    --plot-dir "${TLR_PLOT_DIR:-$HOME/ashg_revised/23_shortlist_followup/custom_loci/TLR_popped/internal_only}" "$@"
