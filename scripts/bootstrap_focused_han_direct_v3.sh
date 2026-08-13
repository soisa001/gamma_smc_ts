#!/usr/bin/env bash
set -euo pipefail

# Exact, isolated runtime bootstrap for focused Han direct-v3 production.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE="$REPO/scripts/bootstrap_eas_sweep_study.sh"
MAMBA="$REPO/.tools/micromamba"
PREFIX="$REPO/.native-stdpopsim"
SLIM="$PREFIX/bin/slim"
PYTHON_VERSION="3.11.15"
SLIM_VERSION="4.2.2"
SLIM_BUILD="h5888daf_1"
SLIM_PACKAGE_SHA256="d9ba79f22fe916bb5f76eb66284e6ec45dcb0e8694319a9f0112e07f9145aec3"
SLIM_BINARY_SHA256="67409ff190808c2967c949cae6697bf235010dbe832389cb4976d9089b5f5e1d"

for argument in "$@"; do
    if [[ "$argument" == "--python" || "$argument" == --python=* ]]; then
        echo "Han direct-v3 pins Python $PYTHON_VERSION; do not pass --python." >&2
        exit 2
    fi
done

bash "$BASE" "$@" --python "$PYTHON_VERSION"

"$MAMBA" install -y -p "$PREFIX" \
    -c conda-forge -c bioconda --strict-channel-priority \
    "slim=$SLIM_VERSION=$SLIM_BUILD"

metadata="$PREFIX/conda-meta/slim-$SLIM_VERSION-$SLIM_BUILD.json"
[[ -f "$metadata" && -x "$SLIM" ]] || {
    echo "Pinned Han direct-v3 SLiM package is absent." >&2
    exit 2
}
actual_package_sha="$($REPO/.venv/bin/python - "$metadata" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1], encoding="utf-8"))["sha256"])
PY
)"
actual_binary_sha="$(sha256sum "$SLIM" | awk '{print $1}')"
actual_python="$($REPO/.venv/bin/python -c 'import platform; print(platform.python_version())')"
actual_python_packages="$($REPO/.venv/bin/python - <<'PY'
import json
import msprime
import numpy
import pandas
import pyslim
import stdpopsim
import tskit

print(json.dumps({
    "msprime": msprime.__version__,
    "numpy": numpy.__version__,
    "pandas": pandas.__version__,
    "pyslim": pyslim.__version__,
    "stdpopsim": stdpopsim.__version__,
    "tskit": tskit.__version__,
}, sort_keys=True))
PY
)"
expected_python_packages='{"msprime": "1.4.2", "numpy": "2.4.6", "pandas": "3.0.3", "pyslim": "1.1.1", "stdpopsim": "0.3.0", "tskit": "1.0.3"}'

[[ "$actual_package_sha" == "$SLIM_PACKAGE_SHA256" ]] || {
    echo "Pinned Han direct-v3 SLiM package SHA-256 differs." >&2
    exit 2
}
[[ "$actual_binary_sha" == "$SLIM_BINARY_SHA256" ]] || {
    echo "Pinned Han direct-v3 SLiM binary SHA-256 differs." >&2
    exit 2
}
[[ "$actual_python" == "$PYTHON_VERSION" ]] || {
    echo "Pinned Han direct-v3 Python version differs." >&2
    exit 2
}
[[ "$actual_python_packages" == "$expected_python_packages" ]] || {
    echo "Pinned Han direct-v3 Python package versions differ:" >&2
    echo "$actual_python_packages" >&2
    exit 2
}

PYTHONPATH="$REPO/python" "$REPO/.venv/bin/python" \
    "$REPO/scripts/run_focused_han_direct_v3.py" --help >/dev/null
printf 'Han direct-v3 runtime verified: Python %s, SLiM %s=%s\n' \
    "$PYTHON_VERSION" "$SLIM_VERSION" "$SLIM_BUILD"
