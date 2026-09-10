#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
shift || true
out="${JOINT_SCAN_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_joint_scan}"
native="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export UV_NO_CONFIG=1 PYTHONUNBUFFERED=1
cd "${repo_dir}"
mkdir -p "${out}/logs"
common=(--baseline "${SIM_BASELINE_DIR:-/mnt/d/phase2simselection/sim/eas_q02_h400}"
        --later "${SIM_LATER_DIR:-/mnt/d/phase2simselection/sim/eas_q02_onset10k_h400}"
        --out "${out}" --decoder "${repo_dir}/bin/gamma_smc" --workers 20)
uv_run=("${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" python)
if [[ "${phase}" == labels || "${phase}" == run ]]; then
  "${uv_run[@]}" -m gamma_smc_aou.joint_scan_labels "${common[@]}" \
    --slim "${SLIM_BIN:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native-stdpopsim/bin/slim}" "$@" \
    > "${out}/logs/labels.log" 2>&1
fi
if [[ "${phase}" == profiles || "${phase}" == run ]]; then
  if [[ ! -f bin/joint_pair_tmrca.so || cpp/joint_pair_tmrca.cpp -nt bin/joint_pair_tmrca.so ]]; then
    c++ -std=c++17 -O3 -shared -fPIC cpp/joint_pair_tmrca.cpp -o bin/joint_pair_tmrca.so.tmp
    mv bin/joint_pair_tmrca.so.tmp bin/joint_pair_tmrca.so
  fi
  "${uv_run[@]}" -m gamma_smc_aou.joint_scan_profiles "${common[@]}" "$@" \
    > "${out}/logs/profiles.log" 2>&1
fi
if [[ "${phase}" == analyse || "${phase}" == run ]]; then
  "${uv_run[@]}" -m gamma_smc_aou.joint_scan_evaluation --out "${out}" "$@" \
    > "${out}/logs/evaluation.log" 2>&1
fi
