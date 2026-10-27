#!/usr/bin/env bash
set -euo pipefail
# Public tooling only; no cloud data access, no credentials, no repository push.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$REPO/.tools"
if ! command -v uv >/dev/null 2>&1 && [[ ! -x "$REPO/.tools/uv-bin/uv" ]]; then
  curl --proto '=https' --tlsv1.2 -fsSL \
    https://releases.astral.sh/github/uv/releases/download/0.11.16/uv-installer.sh \
    -o "$REPO/.tools/singleton-uv-installer.sh"
  UV_INSTALL_DIR="$REPO/.tools/uv-bin" UV_NO_MODIFY_PATH=1 \
    sh "$REPO/.tools/singleton-uv-installer.sh"
fi
if [[ -n "${BCFTOOLS_BIN:-}" ]]; then
  "$BCFTOOLS_BIN" --version
  exit 0
fi
for executable in "$REPO/.native/bin/bcftools" "$REPO/.tools/singleton-native/bin/bcftools"; do
  if [[ -x "$executable" ]] && "$executable" --version >/dev/null 2>&1; then
    echo "bcftools available: $executable"
    exit 0
  fi
done
if command -v bcftools >/dev/null 2>&1; then
  bcftools --version
  exit 0
fi
[[ "$(uname -s):$(uname -m)" == Linux:x86_64 ]] || {
  echo 'Automatic bcftools installation supports Linux x86_64; otherwise supply BCFTOOLS_BIN.' >&2
  exit 2
}
MAMBA="$REPO/.tools/micromamba"
if [[ ! -x "$MAMBA" ]]; then
  curl --proto '=https' --tlsv1.2 -fsSL \
    https://github.com/mamba-org/micromamba-releases/releases/download/2.6.2-1/micromamba-linux-64 \
    -o "$MAMBA.partial"
  chmod +x "$MAMBA.partial"
  mv "$MAMBA.partial" "$MAMBA"
fi
PREFIX="$REPO/.tools/singleton-native"
ACTION=create
[[ ! -d "$PREFIX/conda-meta" ]] || ACTION=install
"$MAMBA" "$ACTION" -y -p "$PREFIX" -c conda-forge -c bioconda \
  --strict-channel-priority 'bcftools=1.22' 'htslib=1.22'
"$PREFIX/bin/bcftools" --version
