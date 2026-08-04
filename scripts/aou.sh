#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLBACKEND=Agg
UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
[[ -z "$UV" && -x "$HOME/.local/bin/uv" ]] && UV="$HOME/.local/bin/uv"
[[ -z "$UV" ]] && UV="$(command -v uv || true)"
if [[ -z "$UV" ]]; then
    echo "uv is not installed; run scripts/bootstrap_uv.sh first." >&2
    exit 2
fi

[[ -z "${SLIM_BIN:-}" && -x "$REPO/.native/bin/slim" ]] && \
    export SLIM_BIN="$REPO/.native/bin/slim"
[[ -z "${GAMMA_SMC_BIN:-}" && -x "$REPO/bin/gamma_smc" ]] && \
    export GAMMA_SMC_BIN="$REPO/bin/gamma_smc"
if [[ -d "$REPO/.native/lib" ]]; then
    export LD_LIBRARY_PATH="$REPO/.native/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

if [[ "${1:-}" == "--sync-only" ]]; then
    exec "$UV" sync --project "$REPO" --frozen --all-extras
fi

run_args=(run --project "$REPO" --frozen --all-extras)
if [[ "${AOU_UV_NO_SYNC:-0}" == "1" ]]; then
    run_args+=(--no-sync)
fi
exec "$UV" "${run_args[@]}" gamma-smc-aou "$@"
