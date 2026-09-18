#!/usr/bin/env bash
# Only computes iHS from saved tree sequences; no SLiM or Gamma-SMC invocation.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
runtime="${repo_dir}/.venv-ihs"
export UV_LINK_MODE=copy PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ ! -x "${runtime}/bin/python" ]; then
  "${uv_bin}" --no-config venv "${runtime}" --python 3.11
fi
"${uv_bin}" --no-config pip sync --python "${runtime}/bin/python" "${repo_dir}/configs/requirements-ihs.txt"
cd "${repo_dir}"
exec "${uv_bin}" --no-config run --no-project --python "${runtime}/bin/python" \
  python scripts/benchmark_ihs.py --phase "${1:-run}"
