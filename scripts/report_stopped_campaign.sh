#!/usr/bin/env bash
# Analyze a frozen completed subset. This script never launches simulations.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
phase="${1:-report}"
root="${SIM_OUTPUT_DIR:-/mnt/d/phase2simselection/sim/segregating_introgression_h400}"
out="${REPORT_OUTPUT_DIR:-${root}/report_20260928}"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
python_bin="${ORIGIN_PYTHON:-${repo_dir}/.venv-origin/bin/python}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1 UV_NO_CONFIG=1 UV_LINK_MODE=copy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
case "$phase" in
    all) phases=(truth decoded report) ;;
    truth|decoded|report) phases=("$phase") ;;
    *) printf 'Expected all, truth, decoded, or report\n' >&2; exit 2 ;;
esac
cd "$repo_dir"
for step in "${phases[@]}"; do
    "$uv_bin" --no-config run --no-project --python "$python_bin" \
      python -m gamma_smc_aou.stopped_campaign_report \
      --root "$root" --out "$out" --phase "$step" --workers "${REPORT_WORKERS:-4}"
done
