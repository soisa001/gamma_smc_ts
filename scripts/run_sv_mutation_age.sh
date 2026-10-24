#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1 UV_NO_CONFIG=1 UV_LINK_MODE=copy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
native_dir="${GAMMA_NATIVE_DIR:-$repo_dir/.native}"
if [[ ! -d "$native_dir/lib" && -d /home/mew/AllOfUs_Phase2/gamma_smc_ts/.native/lib ]]; then
  native_dir=/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native
fi
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
cd "$repo_dir"
uv_bin="${UV_BIN:-$(command -v uv || true)}"
if [[ -z "$uv_bin" ]]; then uv_bin="$repo_dir/.tools/uv-bin/uv"; fi
exec "$uv_bin" --no-config run --no-project --python "${ORIGIN_PYTHON:-$repo_dir/.venv-origin/bin/python}" \
  python -m gamma_smc_aou.sv_mutation_age --out "${SIM_OUTPUT_DIR:-$repo_dir/../sim/sv_mutation_age_20261024}" "$@"
