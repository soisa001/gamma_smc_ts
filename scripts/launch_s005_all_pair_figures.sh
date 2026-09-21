#!/usr/bin/env bash
# Render or package saved results only; no simulation or decoder phase.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
phase="${1:-build}"
out="/mnt/d/phase2simselection/sim/eas_lab_meeting_s005_all_pairs_20260921"
export UV_LINK_MODE=copy MPLBACKEND=Agg PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
case "${phase}" in build|package) ;; *) echo 'Allowed phases: build, package' >&2; exit 2 ;; esac
mkdir -p "${out}/logs"
exec > >(tee -a "${out}/logs/${phase}.log") 2>&1
printf 'UTC: %s\nPhase: %s\nRepo: %s\nOutput: %s\n' "$(date -u +%FT%TZ)" "${phase}" "${repo_dir}" "${out}"
cd "${repo_dir}"
git rev-parse HEAD
if [ "${phase}" = build ]; then
  exec "${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" \
    python scripts/s005_all_pair_figures.py
fi
exec "${uv_bin}" --no-config run --no-project --with pypdf==6.10.0 --with pillow \
  python scripts/package_lab_meeting_figures.py --root "${out}"
