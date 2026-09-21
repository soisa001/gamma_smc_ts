#!/usr/bin/env bash
# Saved-tree decoding, analysis, rendering, and packaging are separate phases.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
phase="${1:-build}"
root="/mnt/d/phase2simselection/sim"
out="${root}/eas_lab_meeting_s002_20260921"
all_out="${root}/eas_lab_meeting_s002_all_pairs_20260921"
profiles="${root}/eas_s002_saved_decoding/profiles"
export UV_LINK_MODE=copy MPLBACKEND=Agg PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
case "${phase}" in decode|analyze|build|package) ;; *) echo 'Allowed phases: decode, analyze, build, package' >&2; exit 2 ;; esac
mkdir -p "${out}/logs"
exec > >(tee -a "${out}/logs/${phase}.log") 2>&1
printf 'UTC: %s\nPhase: %s\nRepo: %s\nOutput: %s\n' "$(date -u +%FT%TZ)" "${phase}" "${repo_dir}" "${out}"
cd "${repo_dir}"
git rev-parse HEAD
run_python() { "${uv_bin}" --no-config run --no-project --python "${repo_dir}/.venv/bin/python" python "$@"; }
case "${phase}" in
  decode) run_python scripts/decode_saved_s002.py --workers 20 ;;
  analyze)
    run_python scripts/evaluate_positional_array.py --out "${out}/analysis/positional" \
      --focal-dir "${root}/eas_lab_meeting_20260921/analysis" --selected-s .002 --selected-profiles "${profiles}" --only-selected-s
    run_python scripts/summarize_positional_evaluation.py --directory "${out}/analysis/positional"
    run_python scripts/evaluate_positional_distance.py --out "${out}/analysis/distance" \
      --selected-s .002 --selected-profiles "${profiles}" --focal-dir "${out}/analysis/positional" --all-pairs-only
    ;;
  build) run_python scripts/s002_lab_meeting_figures.py ;;
  package)
    for folder in "${out}" "${all_out}"; do
      "${uv_bin}" --no-config run --no-project --with pypdf==6.10.0 --with pillow \
        python scripts/package_lab_meeting_figures.py --root "${folder}"
    done
    ;;
esac
