#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PHASE="${1:-plan}"
case "$PHASE" in plan|prepare|merge|all) ;; *) echo "Unknown phase: $PHASE" >&2; exit 2 ;; esac
bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
    if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
  done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ "$PHASE" == all || "$PHASE" == merge ]]; then
  SV_PILOT_INTEGRATION=1 "$UV" run --no-project --python 3.12 \
    --with 'truvari==5.5.0' --with 'google-crc32c==1.8.0' \
    python -m unittest discover -s tests -p test_sv_phased_batch.py -v
fi
ARGS=(--phase "$PHASE" --root "${SV_PHASED_ROOT:-$HOME/sv_merge_chr1_21_phased}" \
  --jobs "${SV_CHROM_JOBS:-4}" --memory-gb "${SV_MEMORY_GB:-26}" --bcftools "$BCFTOOLS_BIN")
if [[ -n "${SV_LOCAL_INPUT_TEMPLATE:-}" ]]; then ARGS+=(--input-template "$SV_LOCAL_INPUT_TEMPLATE"); fi
if [[ "${SV_KEEP_STAGED_BCF:-0}" == 1 ]]; then ARGS+=(--keep-staged-bcf); fi
"$UV" run --no-project --python 3.12 --with 'truvari==5.5.0' scripts/run_sv_phased_batch.py "${ARGS[@]}"
