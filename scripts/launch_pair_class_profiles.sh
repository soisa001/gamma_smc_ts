#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
shift || true
output_dir="${PAIR_CLASS_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_pair_classes}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export UV_NO_CONFIG=1 PYTHONUNBUFFERED=1
cd "${repo_dir}"
mkdir -p "${output_dir}/logs"
"${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python -m gamma_smc_aou.pair_class_profiles \
  --baseline "${SIM_BASELINE_DIR:-/mnt/d/phase2simselection/sim/eas_q02_h400}" \
  --later "${SIM_LATER_DIR:-/mnt/d/phase2simselection/sim/eas_q02_onset10k_h400}" \
  --comparison "${SIM_COMPARISON_DIR:-/mnt/d/phase2simselection/sim/eas_recent_calls}" \
  --out "${output_dir}" --decoder "${repo_dir}/bin/gamma_smc" \
  --workers 20 --phase "${phase}" "$@"
if [[ "${phase}" == analyse || "${phase}" == run ]]; then
  "${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
    python scripts/audit_pair_class_profiles.py --out "${output_dir}"
fi
