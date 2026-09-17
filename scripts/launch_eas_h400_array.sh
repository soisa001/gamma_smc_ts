#!/usr/bin/env bash
# Generate or audit saved inputs only. Decoding is a separate future task.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-archive}"
case "${phase}" in
  archive)
    bash "${BASH_SOURCE[0]}" simulate-only
    exec bash "${BASH_SOURCE[0]}" audit-simulations
    ;;
  simulate-only|simulation-status|audit-simulations) ;;
  *) echo "Allowed phases: archive, simulate-only, simulation-status, audit-simulations" >&2; exit 2 ;;
esac
export SIM_CONFIG_FILE="${repo_dir}/configs/eas_q02_50k_h400_array.json"
export SIM_OUTPUT_DIR="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_q02_h400}"
exec bash "${repo_dir}/scripts/launch_fresh_eas.sh" "${phase}"
