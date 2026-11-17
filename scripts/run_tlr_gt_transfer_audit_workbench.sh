#!/usr/bin/env bash
set -euo pipefail
umask 077
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
BATCH="${POPPED_FULL_ROOT:-$HOME/popped_full_same_original_pos_p70}"
OUT="${TLR_GT_AUDIT_ROOT:-$BATCH/tlr_gt_transfer_audit_38754377}"
export TMPDIR="$OUT/runtime/tmp" UV_CACHE_DIR="$OUT/runtime/uv-cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$UV_CACHE_DIR" "$REPO/.tools/uv-bin"
UV="$REPO/.tools/uv-bin/uv"
if [[ ! -x "$UV" ]]; then
    curl --proto '=https' --tlsv1.2 -fsSL \
        https://releases.astral.sh/github/uv/releases/download/0.11.16/uv-installer.sh \
        -o "$REPO/.tools/tlr-audit-uv-installer.sh"
    UV_INSTALL_DIR="$REPO/.tools/uv-bin" UV_NO_MODIFY_PATH=1 \
        sh "$REPO/.tools/tlr-audit-uv-installer.sh"
fi
"$UV" run --no-project --python 3.12 scripts/audit_tlr_gt_transfer.py \
    --batch-root "$BATCH" --out "$OUT" "$@"
