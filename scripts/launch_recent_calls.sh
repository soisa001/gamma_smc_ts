#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
baseline_dir="${SIM_BASELINE_DIR:-/mnt/d/phase2simselection/sim/eas_q02_h400}"
later_dir="${SIM_LATER_DIR:-/mnt/d/phase2simselection/sim/eas_q02_onset10k_h400}"
output_dir="${SIM_COMPARISON_DIR:-/mnt/d/phase2simselection/sim/eas_recent_calls}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export UV_NO_CONFIG=1 PYTHONUNBUFFERED=1
cd "${repo_dir}"
mkdir -p "${output_dir}/logs"
case "${phase}" in
  simulate|run)
    SIM_CONFIG_FILE="${repo_dir}/configs/eas_q02_10k_h400.json" SIM_OUTPUT_DIR="${later_dir}" \
      bash scripts/launch_fresh_eas.sh simulate
    if [[ "${phase}" == simulate ]]; then exit 0; fi
    ;;
  smoke|decode|analyse) ;;
  *) printf 'Unknown phase: %s\n' "${phase}" >&2; exit 2 ;;
esac
exec "${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python -m gamma_smc_aou.recent_call_comparison \
  --baseline "${baseline_dir}" --later "${later_dir}" --out "${output_dir}" \
  --decoder "${repo_dir}/bin/gamma_smc" --workers 20 --phase "${phase}"
