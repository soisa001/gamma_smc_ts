#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
CASE_ENV="${CASE_ENV:-/home/mew/.cache/case_001_env}"
if [[ ! -x "$CASE_ENV/bin/python" ]]; then uv venv "$CASE_ENV"; fi
uv pip install --python "$CASE_ENV/bin/python" -r requirements.txt
case "${1:-predict}" in
  explore)
    "$CASE_ENV/bin/python" explore.py > exploration.txt
    "$CASE_ENV/bin/python" diagnostics.py > diagnostics.txt
    "$CASE_ENV/bin/python" motif.py > motif.txt
    ;;
  validate)
    "$CASE_ENV/bin/python" model.py > model_output.txt
    "$CASE_ENV/bin/python" advanced.py > advanced_output.txt
    ;;
  predict)
    "$CASE_ENV/bin/python" finalize.py > finalization_output.txt
    ;;
  *) echo 'Use: bash run.sh {explore|validate|predict}' >&2; exit 2;;
esac
