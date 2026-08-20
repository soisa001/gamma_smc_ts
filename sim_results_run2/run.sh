#!/usr/bin/env bash
set -euo pipefail

# Phase runner for the run2 study, mirroring introgression_EAS_sim/run.sh.
# Everything goes through the locked uv environment.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKERS="${RUN2_WORKERS:-4}"

usage() {
    cat <<'EOF'
Usage: bash sim_results_run2/run.sh PHASE [runner options]

PHASE is one of:
  validate   generate and assert the SLiM schedule (no SLiM binary required)
  smoke      three replicates per arm and mode, then analyze and plot
  simulate   all replicates
  analyze    result tables from whatever replicates exist
  plot       figures from the result tables
  all        validate, simulate, analyze, plot

Environment overrides:
  RUN2_WORKERS   parallel replicates (default: 4)
  SLIM_BIN       SLiM 4.2.2 executable (default: .native-stdpopsim/bin/slim)
  RUN2_UV_NO_SYNC=1  pass --no-sync to uv run

Additional options are forwarded to scripts/run_run2.py, for example:
  bash sim_results_run2/run.sh simulate --arm chb_ancient_eurasia
EOF
}

PHASE="${1:-}"
case "$PHASE" in
    validate|smoke|simulate|analyze|plot|all) shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
esac

if [[ ! "$WORKERS" =~ ^[1-9][0-9]*$ ]]; then
    echo "RUN2_WORKERS must be a positive integer, got: $WORKERS" >&2
    exit 2
fi

UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
[[ -z "$UV" && -x "${HOME}/.local/bin/uv" ]] && UV="${HOME}/.local/bin/uv"
[[ -z "$UV" ]] && UV="$(command -v uv || true)"
[[ -n "$UV" ]] || {
    echo "uv is unavailable; run bash scripts/bootstrap_run2.sh first." >&2
    exit 2
}

export MPLBACKEND=Agg

needs_slim=0
case "$PHASE" in
    smoke|simulate|all) needs_slim=1 ;;
esac

runner_args=()
if [[ "$needs_slim" -eq 1 ]]; then
    SLIM="${SLIM_BIN:-$REPO/.native-stdpopsim/bin/slim}"
    [[ -x "$SLIM" ]] || {
        echo "SLiM 4.2.2 is unavailable at $SLIM; run scripts/bootstrap_run2.sh first." >&2
        exit 2
    }
    slim_version="$("$SLIM" -v 2>&1)"
    if ! grep -Eq '(^|[^0-9])4\.2\.2([^0-9]|$)' <<<"$slim_version"; then
        echo "run2 requires SLiM 4.2.2, but $SLIM reported:" >&2
        echo "$slim_version" >&2
        exit 2
    fi
    runner_args+=(--slim-bin "$SLIM" --workers "$WORKERS")
elif [[ -n "${SLIM_BIN:-}" ]]; then
    runner_args+=(--slim-bin "$SLIM_BIN")
fi

uv_args=(run --project "$REPO" --frozen --all-extras)
if [[ "${RUN2_UV_NO_SYNC:-0}" == "1" ]]; then
    uv_args+=(--no-sync)
fi

run_phase() {
    local phase="$1"
    shift
    printf 'Running:'
    printf ' %q' "$UV" "${uv_args[@]}" python "$REPO/scripts/run_run2.py" "$phase" "$@"
    printf '\n'
    "$UV" "${uv_args[@]}" python "$REPO/scripts/run_run2.py" "$phase" "$@"
}

case "$PHASE" in
    smoke)
        run_phase simulate "${runner_args[@]}" --index 0 --index 1 --index 2 "$@"
        run_phase analyze
        run_phase plot
        ;;
    all)
        run_phase validate "${runner_args[@]}"
        run_phase simulate "${runner_args[@]}" "$@"
        run_phase analyze
        run_phase plot
        ;;
    *)
        run_phase "$PHASE" "${runner_args[@]}" "$@"
        ;;
esac
