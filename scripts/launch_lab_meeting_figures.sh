#!/usr/bin/env bash
# Summarize saved inputs and build figures; never launch simulations or decoding.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
phase="${1:-build}"
figure_dir="/mnt/d/phase2simselection/sim/eas_lab_meeting_20260921"
runtime="${repo_dir}/.venv/bin/python"
export UV_LINK_MODE=copy MPLBACKEND=Agg PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
case "${phase}" in
  build|refresh|package) ;;
  *) echo 'Allowed phases: build, refresh, package' >&2; exit 2 ;;
esac
mkdir -p "${figure_dir}/logs"
exec > >(tee -a "${figure_dir}/logs/${phase}.log") 2>&1
printf 'UTC: %s\nPhase: %s\nRepo: %s\nOutput: %s\n' "$(date -u +%FT%TZ)" "${phase}" "${repo_dir}" "${figure_dir}"
cd "${repo_dir}"
git rev-parse HEAD
if [ "${phase}" = package ]; then
  exec "${uv_bin}" --no-config run --no-project --with pypdf==6.10.0 --with pillow \
    python scripts/package_lab_meeting_figures.py --root "${figure_dir}"
fi
if [ ! -x "${runtime}" ]; then
  "${uv_bin}" --no-config venv "${repo_dir}/.venv" --python 3.11
  "${uv_bin}" --no-config pip install --python "${runtime}" -e '.[selection]'
fi
run_python() {
  "${uv_bin}" --no-config run --no-project --python "${runtime}" python "$@"
}
if [ "${phase}" = refresh ]; then
  # The archive must already be complete. Hash-verified summaries are reused.
  run_python scripts/prepare_lab_meeting.py
  run_python scripts/evaluate_positional_array.py \
    --out "${figure_dir}/analysis/positional" --focal-dir "${figure_dir}/analysis"
  run_python scripts/summarize_positional_evaluation.py --directory "${figure_dir}/analysis/positional"
  bash scripts/launch_ihs_benchmark.sh run
fi
run_python scripts/lab_meeting_figures.py
