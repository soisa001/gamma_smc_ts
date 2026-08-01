#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UV=""
[[ -x "$REPO/.tools/uv-bin/uv" ]] && UV="$REPO/.tools/uv-bin/uv"
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

exec "$UV" run --project "$REPO" --frozen --all-extras gamma-smc-aou "$@"
