#!/usr/bin/env bash
# Run the three phases in one restartable Linux job.
set -euo pipefail
cd "$(dirname "$0")/.."
raw="${EAS_CAUSAL_RAW:-$HOME/sim/eas_causal_50k_20261123}"
workers="${EAS_WORKERS:-4}"
mkdir -p "$raw"
exec 9>"$raw/workflow.lock"
flock -n 9 || { echo "A workflow is already active in $raw"; exit 1; }
if [[ -n "${EAS_PYTHON:-}" ]]; then
    run=(uv --no-config run --no-project --python "$EAS_PYTHON" python)
else
    run=(uv run --frozen --extra selection python)
fi
"${run[@]}" scripts/run_eas_causal_50k.py --raw "$raw" --workers "$workers"
"${run[@]}" scripts/report_eas_causal_50k.py prepare --raw "$raw"
"${run[@]}" scripts/report_eas_causal_50k.py report
