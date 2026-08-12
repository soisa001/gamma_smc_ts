#!/usr/bin/env bash
set -euo pipefail

# Bootstrap the normal locked environment and Gamma-SMC build first. Keep the
# repository's SLiM 5.2 prefix intact, then install stdpopsim-compatible SLiM
# 4.2.2 into a separate prefix used only by the EAS sweep studies.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_BOOTSTRAP="$REPO/scripts/bootstrap_uv.sh"
MAMBA="$REPO/.tools/micromamba"
BASE_NATIVE="$REPO/.native"
STUDY_NATIVE="$REPO/.native-stdpopsim"
BASE_SLIM="$BASE_NATIVE/bin/slim"
STUDY_SLIM="$STUDY_NATIVE/bin/slim"
BASE_SLIM_VERSION="5.2"
STUDY_SLIM_VERSION="4.2.2"

usage() {
    cat <<'EOF'
Usage: scripts/bootstrap_eas_sweep_study.sh [bootstrap_uv.sh options]

Runs the normal repository bootstrap first, then installs SLiM 4.2.2 into the
separate .native-stdpopsim prefix. The normal .native SLiM 5.2 prefix is kept.

Useful forwarded options:
  --skip-tests       Verify the installation without running pytest.
  --simulation-only  Skip the Gamma-SMC build (decode will remain unavailable).
  --python VERSION   Select the locked Python environment version.

--skip-native is intentionally unsupported because both SLiM prefixes must be
created and verified by this bootstrap.
EOF
}

for argument in "$@"; do
    case "$argument" in
        -h|--help)
            usage
            exit 0
            ;;
        --skip-native)
            echo "--skip-native is incompatible with this bootstrap: the separate" >&2
            echo "SLiM 5.2 and 4.2.2 prefixes must both be verified." >&2
            exit 2
            ;;
    esac
done

bash "$BASE_BOOTSTRAP" "$@"

[[ -x "$MAMBA" ]] || {
    echo "Pinned micromamba was not installed at $MAMBA" >&2
    exit 2
}
[[ -x "$BASE_SLIM" ]] || {
    echo "The normal bootstrap did not install SLiM $BASE_SLIM_VERSION at $BASE_SLIM" >&2
    exit 2
}

version_output() {
    "$1" -v 2>&1
}

require_version() {
    local executable="$1"
    local expected="$2"
    local label="$3"
    local output
    output="$(version_output "$executable")"
    if ! grep -Eq "(^|[^0-9])${expected//./\\.}([^0-9]|$)" <<<"$output"; then
        echo "$label version mismatch; expected $expected:" >&2
        echo "$output" >&2
        exit 2
    fi
    printf '%s: %s\n' "$label" "$(head -n 1 <<<"$output")"
}

require_version "$BASE_SLIM" "$BASE_SLIM_VERSION" "normal SLiM prefix"

study_packages=("slim=$STUDY_SLIM_VERSION")
if [[ -d "$STUDY_NATIVE/conda-meta" ]]; then
    "$MAMBA" install -y -p "$STUDY_NATIVE" \
        -c conda-forge -c bioconda --strict-channel-priority \
        "${study_packages[@]}"
else
    "$MAMBA" create -y -p "$STUDY_NATIVE" \
        -c conda-forge -c bioconda --strict-channel-priority \
        "${study_packages[@]}"
fi

[[ -x "$STUDY_SLIM" ]] || {
    echo "Study SLiM executable was not installed at $STUDY_SLIM" >&2
    exit 2
}
require_version "$STUDY_SLIM" "$STUDY_SLIM_VERSION" "EAS-study SLiM prefix"
require_version "$BASE_SLIM" "$BASE_SLIM_VERSION" "preserved normal SLiM prefix"

cat <<EOF

EAS sweep study environment is ready.
  Gamma-SMC and general simulation prefix: $BASE_NATIVE
  EAS stdpopsim SLiM prefix:              $STUDY_NATIVE

Run a phase with, for example:
  bash no_introgression_EAS_sim/run.sh plan
EOF
