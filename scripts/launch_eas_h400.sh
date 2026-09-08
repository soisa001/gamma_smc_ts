#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
export SIM_CONFIG_FILE="${repo_dir}/configs/eas_q02_50k_h400.json"
export SIM_OUTPUT_DIR="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_q02_h400}"
baseline_dir="${SIM_BASELINE_DIR:-/mnt/d/phase2simselection/sim/eas_q02}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export UV_NO_CONFIG=1 PYTHONUNBUFFERED=1
mkdir -p "${SIM_OUTPUT_DIR}/logs"
cd "${repo_dir}"
case "${phase}" in
  smoke|status)
    exec bash scripts/launch_fresh_eas.sh "${phase}"
    ;;
  run)
    bash scripts/launch_fresh_eas.sh run
    ;;
  analyse)
    ;;
  *)
    printf 'Unknown phase: %s\n' "${phase}" >&2
    exit 2
    ;;
esac
"${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python -m gamma_smc_aou.paired_peak_grid \
  --study "${SIM_OUTPUT_DIR}" --workers 20
"${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python -m gamma_smc_aou.panel_size_comparison \
  --baseline "${baseline_dir}" --larger "${SIM_OUTPUT_DIR}" \
  --out "${SIM_OUTPUT_DIR}/analysis/panel_size_comparison" --workers 20
