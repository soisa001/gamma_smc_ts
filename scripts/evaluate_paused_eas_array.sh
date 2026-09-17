#!/usr/bin/env bash
# Evaluation only: never resumes simulations or starts the decoder.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
sim_root="${SIM_ROOT:-/mnt/d/phase2simselection/sim}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "${repo_dir}"
exec "${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
  python scripts/evaluate_paused_array.py \
  --study "${sim_root}/eas_q02_h400" \
  --reference "${sim_root}/eas_allele_class_ablation" \
  --focal-reference "${sim_root}/eas_allele_focal_comparison" \
  --profiles "${sim_root}/eas_joint_scan/profile_manifest.json" \
  --out "${sim_root}/eas_h400_pause_evaluation_20260917" \
  --workers 20
