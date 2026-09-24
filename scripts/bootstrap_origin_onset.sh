#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
uv_bin="${UV_BIN:-${HOME}/.local/bin/uv}"
export UV_NO_CONFIG=1 UV_LINK_MODE=copy
cd "$repo_dir"
if [ ! -x .venv-origin/bin/python ]; then
  "$uv_bin" --no-config venv .venv-origin --python 3.12
fi
"$uv_bin" --no-config pip install --python .venv-origin/bin/python \
  --constraint configs/origin_onset_constraints.txt --editable '.[selection,test]'
