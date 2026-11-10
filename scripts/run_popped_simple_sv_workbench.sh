#!/usr/bin/env bash
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
bash scripts/bootstrap_singleton_tools.sh
UV="$REPO/.tools/uv-bin/uv"
if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  for candidate in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
    if [[ -x "$candidate" ]]; then export BCFTOOLS_BIN="$candidate"; break; fi
  done
fi
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=(--root "${SV_SIMPLE_ROOT:-$HOME/chr1_popped_simple_merge}" --threads 4 --memory-gb 26 --bcftools "$BCFTOOLS_BIN")
# Prefer the input named in the completed comparison, with exact popped basename
# and an index. Never substitute bubble/imputed inputs or search arbitrary files.
INPUT="${SV_POPPED_INPUT:-}"
if [[ -z "$INPUT" ]]; then
  INPUT="$("$UV" run --no-project --python 3.12 python - <<'PY'
import json
from pathlib import Path
home = Path.home()
name = 'aou_lr_phase2_v1.chr1.popped.bcf'
candidates = []
metadata = home / 'chr1_sv_representation_comparison/comparison_metadata.json'
if metadata.exists():
    candidates.append(Path(json.loads(metadata.read_text())['inputs']['popped']))
candidates.extend(home / root / 'inputs' / name for root in (
    'chr1_confirmed_popped_representation', 'chr1_popped_representation',
    'chr1_popped_simple_merge'))
for path in candidates:
    if path.name == name and path.is_file() and Path(str(path) + '.csi').is_file():
        print(path)
        break
PY
)"
fi
if [[ -n "$INPUT" ]]; then
  echo "Using staged popped panel: $INPUT"
  ARGS+=(--input-bcf "$INPUT")
else
  echo 'No indexed popped input found in known locations; staging the confirmed popped chr1 URI.'
fi
"$UV" run --no-project --python 3.12 scripts/run_popped_simple_sv.py "${ARGS[@]}" "$@"
