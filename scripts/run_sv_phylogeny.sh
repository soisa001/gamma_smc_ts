#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"
export MPLBACKEND=Agg UV_NO_CONFIG=1 PYTHONUNBUFFERED=1
export IQTREE_BIN="${IQTREE_BIN:-$(find "$repo_dir/.tools/phylo/iqtree" -type f -name iqtree2 -print -quit)}"
export RELATE_BIN="${RELATE_BIN:-$(find "$repo_dir/.tools/phylo/relate-linux" -type f -name Relate -print -quit)}"
exec uv run --no-project --python "$repo_dir/.tools/phylo/venv/bin/python" \
  python scripts/sv_phylogeny_benchmark.py "$@"
