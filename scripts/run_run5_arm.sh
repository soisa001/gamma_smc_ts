#!/usr/bin/env bash
set -euo pipefail

# Space-constrained runner for one run5 arm.
#
# Simulates and decodes a single arm, then deletes that arm's tree sequences.
# Peak disk use stays near one arm's worth (~2.6 GB) instead of the whole
# study's (~8 GB), which matters when the host disk is nearly full. The raw
# per-pair TMRCA matrices are kept -- those are the artifact worth preserving,
# since every threshold and stratification is computed from them. Re-decoding
# would require re-simulating, which the recorded seeds make reproducible.
#
# Usage: bash scripts/run_run5_arm.sh <arm_id> [workers]

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARM="${1:?usage: run_run5_arm.sh <arm_id> [workers]}"
WORKERS="${2:-8}"

UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
[[ -z "$UV" ]] && UV="$(command -v uv || true)"
[[ -n "$UV" ]] || { echo "uv is unavailable" >&2; exit 2; }

export MPLBACKEND=Agg
SLIM="${SLIM_BIN:-$REPO/.native-stdpopsim/bin/slim}"
DECODER="${GAMMA_SMC_BIN:-$REPO/bin/gamma_smc}"

run_phase() {
    "$UV" run --project "$REPO" --frozen --all-extras --no-sync \
        python "$REPO/scripts/run_run5.py" "$@"
}

echo "=== ${ARM}: simulate ==="
run_phase simulate --slim-bin "$SLIM" --arm "$ARM" --workers "$WORKERS" 2>&1 | tail -6

echo "=== ${ARM}: decode ==="
run_phase decode --decoder-bin "$DECODER" --arm "$ARM" --workers "$WORKERS" 2>&1 | tail -6

echo "=== ${ARM}: dropping tree sequences ==="
find "$REPO/sim_results_run5/$ARM" -name simulation.trees -delete
df -h "$REPO" | tail -1
