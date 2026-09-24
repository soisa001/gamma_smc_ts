#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export SIM_CONFIG_FILE="${SIM_CONFIG_FILE:-${repo_dir}/configs/segregating_introgression.json}"
export SIM_OUTPUT_DIR="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/segregating_introgression_h400}"
export SLIM_BIN="${SLIM_BIN:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native/bin/slim}"
exec bash "${repo_dir}/scripts/launch_origin_onset.sh" "$@"
