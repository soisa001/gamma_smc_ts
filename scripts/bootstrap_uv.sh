#!/usr/bin/env bash
set -euo pipefail

# Rootless, reproducible bootstrap for Linux. Python is locked by uv; SLiM and
# the C++ toolchain/libraries live in a repository-local micromamba prefix.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="$REPO/.tools"
NATIVE="$REPO/.native"
UV_VERSION="0.11.16"
MICROMAMBA_VERSION="2.6.2-1"
PYTHON_VERSION="3.11"
FULL=1
RUN_TESTS=1
INSTALL_NATIVE=1

usage() {
    cat <<'EOF'
Usage: scripts/bootstrap_uv.sh [options]

Options:
  --simulation-only  Install Python + SLiM, but do not build Gamma-SMC.
  --skip-native      Use SLIM_BIN/native tools already supplied by the caller.
  --skip-tests       Verify imports/binaries but do not run pytest.
  --python VERSION   Python version for uv (default: 3.11).
  -h, --help         Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --simulation-only) FULL=0; shift ;;
        --skip-native) INSTALL_NATIVE=0; shift ;;
        --skip-tests) RUN_TESTS=0; shift ;;
        --python) PYTHON_VERSION="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

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
    UV_INSTALL_DIR="$TOOLS/uv-bin" UV_NO_MODIFY_PATH=1 \
        sh "$TOOLS/uv-installer.sh"
    UV="$LOCAL_UV"
fi
[[ "$(uv_numeric_version "$UV")" == "$UV_VERSION" ]] || {
    echo "Failed to select uv $UV_VERSION: $UV" >&2
    exit 2
}

SLIM="${SLIM_BIN:-}"
if [[ "$INSTALL_NATIVE" -eq 1 ]]; then
    MACHINE="$(uname -m)"
    case "$(uname -s):$MACHINE" in
        Linux:x86_64) MAMBA_ASSET="micromamba-linux-64" ;;
        Darwin:x86_64) MAMBA_ASSET="micromamba-osx-64" ;;
        Darwin:arm64) MAMBA_ASSET="micromamba-osx-arm64" ;;
        *) echo "Unsupported native platform: $(uname -s) $MACHINE" >&2; exit 2 ;;
    esac
    MAMBA="$TOOLS/micromamba"
    if [[ ! -x "$MAMBA" ]]; then
        echo "Installing micromamba $MICROMAMBA_VERSION..."
        curl --proto '=https' --tlsv1.2 -Lsf \
            "https://github.com/mamba-org/micromamba-releases/releases/download/$MICROMAMBA_VERSION/$MAMBA_ASSET" \
            -o "$MAMBA"
        chmod +x "$MAMBA"
    fi

    native_packages=("slim=5.2")
    if [[ "$FULL" -eq 1 ]]; then
        if [[ "$(uname -s):$MACHINE" != "Linux:x86_64" ]]; then
            echo "Gamma-SMC decoding requires Linux x86_64/AVX2; use --simulation-only here." >&2
            exit 2
        fi
        native_packages+=("gxx_linux-64" "make" "boost-cpp" "htslib" "bcftools" "zstd")
    fi
    if [[ -d "$NATIVE/conda-meta" ]]; then
        "$MAMBA" install -y -p "$NATIVE" -c conda-forge -c bioconda \
            --strict-channel-priority \
            "${native_packages[@]}"
    else
        "$MAMBA" create -y -p "$NATIVE" -c conda-forge -c bioconda \
            --strict-channel-priority \
            "${native_packages[@]}"
    fi
    SLIM="$NATIVE/bin/slim"
fi

echo "Installing locked Python environment..."
"$UV" python install "$PYTHON_VERSION"
"$UV" sync --project "$REPO" --frozen --all-extras --python "$PYTHON_VERSION"

DECODER="${GAMMA_SMC_BIN:-}"
if [[ "$FULL" -eq 1 ]]; then
    if ! grep -qw avx2 /proc/cpuinfo; then
        echo "Gamma-SMC requires an AVX2-capable CPU." >&2
        exit 2
    fi
    if [[ "$INSTALL_NATIVE" -eq 1 ]]; then
        CXX="$NATIVE/bin/x86_64-conda-linux-gnu-c++"
        MAKE="$NATIVE/bin/make"
        if [[ ! -x "$CXX" ]]; then
            echo "Conda C++ compiler was not installed at $CXX" >&2
            exit 2
        fi
        export CPATH="$NATIVE/include${CPATH:+:$CPATH}"
        export LIBRARY_PATH="$NATIVE/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
        export LD_LIBRARY_PATH="$NATIVE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    else
        CXX="${CXX:-$(command -v g++)}"
        MAKE="${MAKE:-$(command -v make)}"
    fi
    echo "Building Gamma-SMC decoder..."
    "$MAKE" -C "$REPO" clean
    "$MAKE" -C "$REPO" CXX="$CXX" MARCH=x86-64
    DECODER="$REPO/bin/gamma_smc"
fi

export SLIM_BIN="$SLIM"
[[ -n "$DECODER" ]] && export GAMMA_SMC_BIN="$DECODER"
verify_args=(--require-slim)
[[ "$FULL" -eq 1 ]] && verify_args+=(--require-decoder)
"$UV" run --project "$REPO" --frozen --all-extras python \
    "$REPO/scripts/verify_install.py" "${verify_args[@]}"

if [[ "$RUN_TESTS" -eq 1 ]]; then
    echo "Running test suite..."
    "$UV" run --project "$REPO" --frozen --all-extras pytest -q
fi

echo
echo "Ready. Run commands through: scripts/aou.sh <gamma-smc-aou arguments>"
