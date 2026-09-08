#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
sim_dir="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_q02}"
config_file="${SIM_CONFIG_FILE:-${repo_dir}/configs/eas_q02_50k.json}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
slim_bin="${SLIM_BIN:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native-stdpopsim/bin/slim}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export UV_NO_CONFIG=1 UV_LINK_MODE=copy
export PYTHONUNBUFFERED=1
mkdir -p "${sim_dir}/logs" "${sim_dir}/tmp"
export TMPDIR="${sim_dir}/tmp"
cd "${repo_dir}"
exec "${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python -m gamma_smc_aou.fresh_power --phase "${phase}" \
  --config "${config_file}" --out "${sim_dir}" \
  --slim "${slim_bin}" --decoder "${repo_dir}/bin/gamma_smc"
