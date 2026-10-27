"""Exercise the real Bash bootstrap with stub tools; no network or data access."""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

REPO = Path(__file__).resolve().parents[1]
BASH = str(Path("C:/Program Files/Git/bin/bash.exe")) if os.name == "nt" else shutil.which("bash")


def shell_path(path):
    value = path.resolve().as_posix()
    return f"/{value[0].lower()}{value[2:]}" if os.name == "nt" else value


def executable(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode())
    path.chmod(0o755)


@pytest.mark.parametrize("private_version", [None, "0.11.3", "0.12.0", "0.11.16"])
def test_pinned_uv_replaces_incompatible_local_and_ignores_old_path(tmp_path, private_version):
    executable(tmp_path / "scripts/bootstrap_singleton_tools.sh", (REPO / "scripts/bootstrap_singleton_tools.sh").read_text())
    executable(tmp_path / "bin/uv", "#!/bin/sh\necho 'uv 0.11.3'\n")
    executable(tmp_path / "bin/bcftools", "#!/bin/sh\necho 'bcftools 1.22'\n")
    private = tmp_path / ".tools/uv-bin/uv"
    if private_version:
        executable(private, f"#!/bin/sh\necho 'uv {private_version}'\n")
    executable(tmp_path / "bin/curl", """#!/bin/sh
set -eu
echo download >> "$MOCK_DOWNLOAD_LOG"
while [ "$#" -gt 0 ]; do
  if [ "$1" = -o ]; then shift; out="$1"; fi
  shift
done
cat > "$out" <<'INSTALL'
#!/bin/sh
set -eu
test "$UV_NO_MODIFY_PATH" = 1
mkdir -p "$UV_INSTALL_DIR"
printf '#!/bin/sh\\necho "uv 0.11.16"\\n' > "$UV_INSTALL_DIR/uv"
chmod +x "$UV_INSTALL_DIR/uv"
INSTALL
""")
    command = ['bash', shell_path(tmp_path / "scripts/bootstrap_singleton_tools.sh")]
    script = 'export PATH="$1/bin:$PATH" BCFTOOLS_BIN="$1/bin/bcftools" MOCK_DOWNLOAD_LOG="$1/downloads"; shift; "$@"'
    for _ in range(2):
        result = subprocess.run([BASH, "-c", script, "--", shell_path(tmp_path), *command], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert private.exists(), "Old uv on PATH incorrectly prevented installing compatible private uv"
        version = subprocess.check_output([BASH, shell_path(private), "--version"], text=True).strip()
        assert version == "uv 0.11.16"
    downloads = tmp_path / "downloads"
    assert (downloads.read_text().splitlines() if downloads.exists() else []) == ([] if private_version == "0.11.16" else ["download"])
