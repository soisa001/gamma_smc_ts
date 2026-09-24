#!/usr/bin/env bash
# Explicit phases; run is restartable and executes simulation -> audit -> decode -> analyze.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-plan}"
if [ "$#" -gt 0 ]; then shift; fi
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
python_bin="${ORIGIN_PYTHON:-${repo_dir}/.venv-origin/bin/python}"
out="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/origin_onset_pilot}"
slim_bin="${SLIM_BIN:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native-stdpopsim/bin/slim}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1 UV_NO_CONFIG=1 UV_LINK_MODE=copy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$repo_dir"
exec "$uv_bin" --no-config run --no-project --python "$python_bin" python -m gamma_smc_aou.origin_onset \
  --config "${SIM_CONFIG_FILE:-${repo_dir}/configs/origin_onset_pilot.json}" \
  --out "$out" --slim "$slim_bin" --decoder "${repo_dir}/bin/gamma_smc" \
  --legacy-root "${LEGACY_SIM_ROOT:-/mnt/d/phase2simselection/sim/eas_q02_h400}" \
  --phase "$phase" "$@"
