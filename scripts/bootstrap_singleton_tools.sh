#!/usr/bin/env bash
set -euo pipefail
# Public tooling only; no cloud data access, no credentials, no repository push.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$REPO/.tools"
LOCAL_UV="$REPO/.tools/uv-bin/uv"
# The repository requires >=0.11.16,<0.12. A uv on PATH may be older or
# newer than that range. Pin a private copy and have callers use this path.
uv_ready() {
  [[ -x "$LOCAL_UV" ]] || return 1
  local version
  version="$("$LOCAL_UV" --version 2>/dev/null)" || return 1
  [[ "$version" =~ ^uv\ 0\.11\.16($|[[:space:]]) ]]
}
if ! uv_ready; then
  curl --proto '=https' --tlsv1.2 -fsSL \
    https://releases.astral.sh/github/uv/releases/download/0.11.16/uv-installer.sh \
    -o "$REPO/.tools/singleton-uv-installer.sh"
  UV_INSTALL_DIR="$REPO/.tools/uv-bin" UV_NO_MODIFY_PATH=1 \
    sh "$REPO/.tools/singleton-uv-installer.sh"
fi
uv_ready || { echo "Expected uv 0.11.16 at $LOCAL_UV after installation" >&2; exit 2; }
echo "uv available: $LOCAL_UV ($("$LOCAL_UV" --version))"
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
