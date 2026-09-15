#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-run}"
input_dir="${JOINT_SCAN_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_joint_scan}"
out="${ALLELE_ABLATION_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/eas_allele_class_ablation}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${repo_dir}"
uv_run=("${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" python)
mkdir -p "${out}/logs"
if [[ "${phase}" == score || "${phase}" == run ]]; then
  "${uv_run[@]}" -m gamma_smc_aou.allele_class_ablation --input "${input_dir}" --out "${out}" --workers 20 > "${out}/logs/analysis.log" 2>&1
fi
if [[ "${phase}" == audit || "${phase}" == run ]]; then
  "${uv_run[@]}" scripts/audit_allele_class_ablation.py --out "${out}" --export "${repo_dir}/docs/results/allele_class_ablation" > "${out}/logs/audit.log" 2>&1
fi
