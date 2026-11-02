#!/usr/bin/env bash
set -euo pipefail
# Run ONLY in an authorized data VM. Outputs remain on local VM disk.
# Public dependency bootstrap and synthetic tests precede any cloud access.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PHASE="${1:-plan}"
case "$PHASE" in plan|stage|prepare|merge|compare|all) ;; *) echo "Unknown phase: $PHASE" >&2; exit 2 ;; esac
ROOT="${SV_MERGE_ROOT:-$HOME/sv_merge_chr1_p70}"
SOURCE_URI="${SV_PANEL_URI:-gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf}"
THREADS="${SV_THREADS:-4}"
MEMORY_GB="${SV_MEMORY_GB:-26}"
[[ "$ROOT" = /* ]] || { echo 'SV_MERGE_ROOT must be absolute' >&2; exit 2; }

bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
    if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
  done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ "$PHASE" == all || "$PHASE" == merge ]]; then
  # These tests use generated data only. A failure prevents the data run.
  SV_PILOT_INTEGRATION=1 "$UV" run --no-project --python 3.12 \
    --with 'truvari==5.5.0' --with 'matplotlib==3.10.8' --with 'google-crc32c==1.8.0' \
    python -m unittest discover -s tests -p test_chr1_sv_merge.py -v
fi

ARGS=(--root "$ROOT" --phase "$PHASE" --source-uri "$SOURCE_URI" \
  --chrom chr1 --threshold 0.7 --refdist 500 --min-size 50 --max-size 100000 \
  --threads "$THREADS" --memory-gb "$MEMORY_GB" --bcftools "$BCFTOOLS_BIN")
if [[ "$PHASE" == all || "$PHASE" == stage ]]; then
  : "${GOOGLE_PROJECT:?Set GOOGLE_PROJECT to the requester-pays billing project}"
  # The default full-panel URI was explicitly confirmed by the user.
  ARGS+=(--confirm-source-uri "$SOURCE_URI" --billing-project "$GOOGLE_PROJECT")
fi
"$UV" run --no-project --python 3.12 --with 'truvari==5.5.0' \
  scripts/run_chr1_sv_merge.py "${ARGS[@]}"
