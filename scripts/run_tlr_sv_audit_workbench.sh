#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
    if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
  done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT="${TLR_AUDIT_ROOT:-$HOME/tlr_popped_sv_audit}"
INPUT="${TLR_POPPED_INPUT:-}"
if [[ -z "$INPUT" ]]; then
  for candidate in \
    "$ROOT/inputs/aou_lr_phase2_v1.chr4.popped.bcf" \
    "$HOME/chr4_popped_representation/inputs/aou_lr_phase2_v1.chr4.popped.bcf" \
    "$HOME/panel_popped_vcf/aou_lr_phase2_v1.chr4.popped.bcf"; do
    if [[ -s "$candidate" && -s "$candidate.csi" ]]; then INPUT="$candidate"; break; fi
  done
fi
ARGS=(--root "$ROOT" --region "${TLR_REGION:-chr4:38660000-38960000}" --threads 4 --memory-gb 26 --bcftools "$BCFTOOLS_BIN")
if [[ -n "$INPUT" ]]; then
  echo "Using local popped panel: $INPUT"
  ARGS+=(--input-bcf "$INPUT")
else
  echo 'No indexed local chr4 popped BCF found; staging full chr4 BCF/index once. This download can be large.'
  echo 'Source: gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.chr4.popped.bcf'
fi
"$UV" run --no-project --python 3.12 scripts/audit_tlr_sv.py "${ARGS[@]}" "$@"
