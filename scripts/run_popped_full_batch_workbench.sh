#!/usr/bin/env bash
set -euo pipefail
umask 077
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
ROOT="${POPPED_FULL_ROOT:-$HOME/popped_full_same_original_pos_p70}"
export TMPDIR="$ROOT/runtime/tmp"
export UV_CACHE_DIR="$ROOT/runtime/uv-cache"
mkdir -p "$TMPDIR" "$UV_CACHE_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
    for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
        if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
    done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
ARGS=(--root "$ROOT" --jobs "${POPPED_FULL_JOBS:-4}" --memory-gb "${POPPED_FULL_MEMORY_GB:-26}" --bcftools "$BCFTOOLS_BIN")
if [[ -n "${POPPED_INPUT_TEMPLATE:-}" ]]; then
    ARGS+=(--input-template "$POPPED_INPUT_TEMPLATE")
fi
if [[ -n "${POPPED_BILLING_PROJECT:-}" ]]; then
    ARGS+=(--billing-project "$POPPED_BILLING_PROJECT")
fi
echo 'Source: gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.{chrom}.popped.bcf'
echo "Output: $ROOT; jobs=${POPPED_FULL_JOBS:-4}; total memory budget=${POPPED_FULL_MEMORY_GB:-26} GB"
"$UV" run --no-project --python 3.12 scripts/run_popped_full_batch.py "${ARGS[@]}" "$@"
