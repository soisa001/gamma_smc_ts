#!/usr/bin/env bash
set -euo pipefail

# Phase runner for the run4 study.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKERS="${RUN4_WORKERS:-4}"
THREADS="${RUN4_THREADS:-1}"

usage() {
    cat <<'EOF'
Usage: bash sim_results_run4/run.sh PHASE [options]

PHASE: validate | smoke | simulate | decode | analyze | plot | all

Environment overrides:
  RUN4_WORKERS   parallel replicates (default: 4)
  RUN4_THREADS   decoder threads per replicate (default: 1)
  SLIM_BIN       SLiM 4.2.2 (default: .native-stdpopsim/bin/slim)
  GAMMA_SMC_BIN  Gamma-SMC (default: bin/gamma_smc)
EOF
}

PHASE="${1:-}"
case "$PHASE" in
    validate|smoke|simulate|decode|analyze|plot|all) shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
esac

UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
[[ -z "$UV" ]] && UV="$(command -v uv || true)"
[[ -n "$UV" ]] || { echo "uv is unavailable; run scripts/bootstrap_run2.sh first." >&2; exit 2; }

export MPLBACKEND=Agg
SLIM="${SLIM_BIN:-$REPO/.native-stdpopsim/bin/slim}"
DECODER="${GAMMA_SMC_BIN:-$REPO/bin/gamma_smc}"

uv_args=(run --project "$REPO" --frozen --all-extras --no-sync)
run_phase() {
    local phase="$1"; shift
    "$UV" "${uv_args[@]}" python "$REPO/scripts/run_run4.py" "$phase" "$@"
}

case "$PHASE" in
    smoke)
        run_phase simulate --slim-bin "$SLIM" --workers "$WORKERS" --index 0 --index 1 --index 2
        run_phase decode --decoder-bin "$DECODER" --workers "$WORKERS" --threads "$THREADS" --index 0 --index 1 --index 2
        run_phase analyze
        run_phase plot
        ;;
    all)
        run_phase validate --slim-bin "$SLIM" --decoder-bin "$DECODER"
        run_phase simulate --slim-bin "$SLIM" --workers "$WORKERS" "$@"
        run_phase decode --decoder-bin "$DECODER" --workers "$WORKERS" --threads "$THREADS" "$@"
        run_phase analyze
        run_phase plot
        ;;
    simulate) run_phase simulate --slim-bin "$SLIM" --workers "$WORKERS" "$@" ;;
    decode)   run_phase decode --decoder-bin "$DECODER" --workers "$WORKERS" --threads "$THREADS" "$@" ;;
    *)        run_phase "$PHASE" "$@" ;;
esac
