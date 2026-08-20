#!/usr/bin/env bash
set -euo pipefail

# Lean bootstrap for the run2 study.
#
# run2 is tree-truth only, so unlike scripts/bootstrap_eas_sweep_study.sh this
# does not build Gamma-SMC and does not install the SLiM 5.2 prefix. It installs
# exactly two things: the locked Python environment, and stdpopsim-compatible
# SLiM 4.2.2 in the same .native-stdpopsim prefix the other EAS studies use, so
# the two bootstraps stay compatible if both are run.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="$REPO/.tools"
STUDY_NATIVE="$REPO/.native-stdpopsim"
STUDY_SLIM="$STUDY_NATIVE/bin/slim"
UV_VERSION="0.11.16"
MICROMAMBA_VERSION="2.6.2-1"
STUDY_SLIM_VERSION="4.2.2"
PYTHON_VERSION="3.11"
RUN_TESTS=1
RUN_VALIDATE=1

usage() {
    cat <<'EOF'
Usage: scripts/bootstrap_run2.sh [options]

Installs the locked Python environment and SLiM 4.2.2, then runs the run2
validation phase to prove the study is executable before any simulation starts.

Options:
  --skip-tests     Do not run pytest.
  --skip-validate  Do not run the run2 validate phase.
  --python VER     Python version for uv (default: 3.11).
  -h, --help       Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-tests) RUN_TESTS=0; shift ;;
        --skip-validate) RUN_VALIDATE=0; shift ;;
        --python) PYTHON_VERSION="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$(uname -s)" in
    Linux|Darwin) ;;
    *) echo "run2 must be bootstrapped on Linux or macOS (try WSL)." >&2; exit 2 ;;
esac

mkdir -p "$TOOLS"

uv_numeric_version() {
    "$1" --version 2>/dev/null | awk '{print $2}'
}

UV=""
SYSTEM_UV="$(command -v uv || true)"
LOCAL_UV="$TOOLS/uv-bin/uv"
if [[ -n "$SYSTEM_UV" && "$(uv_numeric_version "$SYSTEM_UV")" == "$UV_VERSION" ]]; then
    UV="$SYSTEM_UV"
elif [[ -x "$LOCAL_UV" && "$(uv_numeric_version "$LOCAL_UV")" == "$UV_VERSION" ]]; then
    UV="$LOCAL_UV"
else
    if [[ -n "$SYSTEM_UV" ]]; then
        echo "Ignoring incompatible $("$SYSTEM_UV" --version 2>/dev/null || echo uv) at $SYSTEM_UV."
    fi
    echo "Installing repository-pinned uv $UV_VERSION..."
    curl --proto '=https' --tlsv1.2 -LsSf \
        "https://releases.astral.sh/github/uv/releases/download/$UV_VERSION/uv-installer.sh" \
        -o "$TOOLS/uv-installer.sh"
    UV_INSTALL_DIR="$TOOLS/uv-bin" UV_NO_MODIFY_PATH=1 sh "$TOOLS/uv-installer.sh"
    UV="$LOCAL_UV"
fi
[[ "$(uv_numeric_version "$UV")" == "$UV_VERSION" ]] || {
    echo "Failed to select uv $UV_VERSION: $UV" >&2
    exit 2
}

MACHINE="$(uname -m)"
case "$(uname -s):$MACHINE" in
    Linux:x86_64) MAMBA_ASSET="micromamba-linux-64" ;;
    Darwin:x86_64) MAMBA_ASSET="micromamba-osx-64" ;;
    Darwin:arm64) MAMBA_ASSET="micromamba-osx-arm64" ;;
    *) echo "Unsupported platform: $(uname -s) $MACHINE" >&2; exit 2 ;;
esac

MAMBA="$TOOLS/micromamba"
if [[ ! -x "$MAMBA" ]]; then
    echo "Installing micromamba $MICROMAMBA_VERSION..."
    curl --proto '=https' --tlsv1.2 -Lsf \
        "https://github.com/mamba-org/micromamba-releases/releases/download/$MICROMAMBA_VERSION/$MAMBA_ASSET" \
        -o "$MAMBA"
    chmod +x "$MAMBA"
fi

echo "Installing SLiM $STUDY_SLIM_VERSION into $STUDY_NATIVE..."
if [[ -d "$STUDY_NATIVE/conda-meta" ]]; then
    "$MAMBA" install -y -p "$STUDY_NATIVE" \
        -c conda-forge -c bioconda --strict-channel-priority "slim=$STUDY_SLIM_VERSION"
else
    "$MAMBA" create -y -p "$STUDY_NATIVE" \
        -c conda-forge -c bioconda --strict-channel-priority "slim=$STUDY_SLIM_VERSION"
fi

[[ -x "$STUDY_SLIM" ]] || {
    echo "SLiM was not installed at $STUDY_SLIM" >&2
    exit 2
}
slim_version="$("$STUDY_SLIM" -v 2>&1)"
if ! grep -Eq "(^|[^0-9])${STUDY_SLIM_VERSION//./\\.}([^0-9]|$)" <<<"$slim_version"; then
    echo "SLiM version mismatch; run2 requires $STUDY_SLIM_VERSION:" >&2
    echo "$slim_version" >&2
    exit 2
fi
printf 'SLiM: %s\n' "$(head -n 1 <<<"$slim_version")"

echo "Installing the locked Python environment..."
"$UV" python install "$PYTHON_VERSION"
"$UV" sync --project "$REPO" --frozen --all-extras --python "$PYTHON_VERSION"

export MPLBACKEND=Agg

if [[ "$RUN_TESTS" -eq 1 ]]; then
    echo "Running the run2 test suite..."
    "$UV" run --project "$REPO" --frozen --all-extras pytest -q \
        "$REPO/tests/test_run2_models.py" \
        "$REPO/tests/test_run2_analysis.py" \
        "$REPO/tests/test_run2_simulate.py"
fi

if [[ "$RUN_VALIDATE" -eq 1 ]]; then
    echo "Running the run2 validation phase..."
    "$UV" run --project "$REPO" --frozen --all-extras python \
        "$REPO/scripts/run_run2.py" validate --slim-bin "$STUDY_SLIM" >/dev/null
    echo "Validation report: $REPO/sim_results_run2/validation_report.json"
fi

cat <<EOF

run2 is ready.
  SLiM $STUDY_SLIM_VERSION: $STUDY_SLIM
  Locked Python:  $REPO/.venv

Smoke test three replicates per arm and mode:
  bash sim_results_run2/run.sh smoke

Then run the study:
  RUN2_WORKERS=12 bash sim_results_run2/run.sh all
EOF
