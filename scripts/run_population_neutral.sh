#!/usr/bin/env bash
# Explicit, restartable phases: plan, simulate, decode, report, or run.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
root="${NEUTRAL_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/all_population_neutral_20261001}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1 UV_NO_CONFIG=1 UV_LINK_MODE=copy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
phase="${1:-plan}"
if [ "$#" -gt 0 ]; then shift; fi
cd "$repo_dir"
exec "${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project \
  --python "${ORIGIN_PYTHON:-${repo_dir}/.venv-origin/bin/python}" \
  python -m gamma_smc_aou.population_neutral \
  --config "${NEUTRAL_CONFIG_FILE:-${repo_dir}/configs/all_population_neutral.json}" \
  --out "$root" --slim "${SLIM_BIN:-${native_dir}/bin/slim}" --phase "$phase" "$@"
