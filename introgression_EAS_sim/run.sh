#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCENARIO="introgression"
WORKERS="${EAS_SWEEP_WORKERS:-4}"
THREADS="${EAS_SWEEP_THREADS:-4}"

usage() {
    cat <<'EOF'
Usage: bash introgression_EAS_sim/run.sh PHASE [runner options]

PHASE is one of: plan, simulate, decode, plot, all

Environment overrides:
  EAS_SWEEP_WORKERS   concurrent SLiM specifications (default: 4)
  EAS_SWEEP_THREADS   Gamma-SMC threads (default: 4)
  SLIM_BIN            SLiM 4.2.2 executable
  GAMMA_SMC_BIN       Gamma-SMC executable
  EAS_SWEEP_UV_NO_SYNC=1  pass --no-sync to uv run

Additional options are passed to scripts/run_eas_sweep_study.py. For example:
  bash introgression_EAS_sim/run.sh simulate --spec s0p010_af50
EOF
}

PHASE="${1:-}"
case "$PHASE" in
    plan|simulate|decode|plot|all) shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
esac

for value_name in WORKERS THREADS; do
    value="${!value_name}"
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        echo "$value_name must be a positive integer, got: $value" >&2
        exit 2
    fi
done

UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
[[ -z "$UV" && -x "${HOME}/.local/bin/uv" ]] && UV="${HOME}/.local/bin/uv"
[[ -z "$UV" ]] && UV="$(command -v uv || true)"
[[ -n "$UV" ]] || {
    echo "uv is unavailable; run bash scripts/bootstrap_eas_sweep_study.sh first." >&2
    exit 2
}

export MPLBACKEND=Agg
if [[ -d "$REPO/.native/lib" ]]; then
    export LD_LIBRARY_PATH="$REPO/.native/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

runner_args=(
    "$PHASE"
    --scenario "$SCENARIO"
    --repo-root "$REPO"
    --workers "$WORKERS"
    --threads "$THREADS"
)

if [[ "$PHASE" == "simulate" || "$PHASE" == "all" ]]; then
    SLIM="${SLIM_BIN:-$REPO/.native-stdpopsim/bin/slim}"
    [[ -x "$SLIM" ]] || {
        echo "SLiM 4.2.2 is unavailable at $SLIM; run the study bootstrap first." >&2
        exit 2
    }
    slim_version="$("$SLIM" -v 2>&1)"
    if ! grep -Eq '(^|[^0-9])4\.2\.2([^0-9]|$)' <<<"$slim_version"; then
        echo "The EAS study requires SLiM 4.2.2, but $SLIM reported:" >&2
        echo "$slim_version" >&2
        exit 2
    fi
    runner_args+=(--slim-bin "$SLIM")
fi

if [[ "$PHASE" == "decode" || "$PHASE" == "all" ]]; then
    DECODER="${GAMMA_SMC_BIN:-$REPO/bin/gamma_smc}"
    [[ -x "$DECODER" ]] || {
        echo "Gamma-SMC is unavailable at $DECODER; run the full study bootstrap first." >&2
        exit 2
    }
    runner_args+=(--decoder-bin "$DECODER")
elif [[ "$PHASE" == "plot" && -n "${GAMMA_SMC_BIN:-}" ]]; then
    DECODER="$GAMMA_SMC_BIN"
    [[ -x "$DECODER" ]] || {
        echo "Explicit GAMMA_SMC_BIN is unavailable at $DECODER." >&2
        exit 2
    }
    runner_args+=(--decoder-bin "$DECODER")
fi

runner_args+=("$@")
uv_args=(run --project "$REPO" --frozen --all-extras)
if [[ "${EAS_SWEEP_UV_NO_SYNC:-0}" == "1" ]]; then
    uv_args+=(--no-sync)
fi

printf 'Running:'
printf ' %q' "$UV" "${uv_args[@]}" python "$REPO/scripts/run_eas_sweep_study.py" "${runner_args[@]}"
printf '\n'
exec "$UV" "${uv_args[@]}" python \
    "$REPO/scripts/run_eas_sweep_study.py" "${runner_args[@]}"
