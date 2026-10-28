"""Download pinned tools into a local directory, preserving SHA256 receipts.

Run with Python 3.12 after installing uv. Java 17+ must already be installed.
Linux uses upstream Relate 1.2.4. Windows builds the same source algorithms
with small filesystem/resource-reporting portability changes (see build script).
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tools", type=Path, default=ROOT / ".tools/phylo")
    args = parser.parse_args()
    base = args.tools.resolve()
    base.mkdir(parents=True, exist_ok=True)
    manifest = {}

    def download(url, name):
        path = base / name
        record = path.with_suffix(path.suffix + ".download.json")
        valid = False
        if path.exists() and record.exists():
            old = json.loads(record.read_text())
            valid = (
                old["url"] == url
                and old["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
            )
        if not valid:
            tmp = path.with_suffix(path.suffix + ".tmp")
            urllib.request.urlretrieve(url, tmp)
            tmp.replace(path)
        info = {
            "url": url,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
        record.write_text(json.dumps(info, indent=2) + "\n")
        manifest[name] = info
        return path

    def unzip(path, dest):
        with zipfile.ZipFile(path) as z:
            assert z.testzip() is None
            for member in z.infolist():
                assert (dest / member.filename).resolve().is_relative_to(dest.resolve())
            z.extractall(dest)

    def untar(path, dest):
        with tarfile.open(path) as t:
            t.extractall(dest, filter="data")

    for component in ("base", "app"):
        p = download(
            f"https://github.com/CompEvol/beast2/releases/download/v2.7.8/BEAST.{component}.package.v2.7.8.zip",
            f"beast-{component}.zip",
        )
        unzip(p, base / f"beast-{component}")
    download(
        "https://github.com/CompEvol/beast2/releases/download/v2.7.7/launcher.jar",
        "launcher.jar",
    )
    uv = shutil.which("uv")
    if not uv or not shutil.which("java"):
        raise RuntimeError("Install uv and Java 17+ first")
    venv = base / "venv"
    python = venv / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not python.exists():
        subprocess.run([uv, "venv", "--python", "3.12", str(venv)], check=True)
    subprocess.run(
        [
            uv,
            "pip",
            "install",
            "--python",
            str(python),
            "-r",
            str(ROOT / "configs/sv_phylogeny_requirements.txt"),
        ],
        check=True,
    )
    if os.name == "nt":
        p = download(
            "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Windows.zip",
            "iqtree.zip",
        )
        unzip(p, base / "iqtree")
        p = download(
            "https://github.com/skeeto/w64devkit/releases/download/v2.10.1/w64devkit-x64-2.10.1.7z.exe",
            "w64devkit.exe",
        )
        subprocess.run([str(p), "-y", "-o" + str(base)], check=True)
        untar(
            download("https://zlib.net/fossils/zlib-1.3.1.tar.gz", "zlib.tar.gz"), base
        )
        repo = base / "relate-sparse"
        commit = "b54ede259cbb0be095bc9c9a8bd18cdaf7e88b74"
        if not (repo / ".git").exists():
            subprocess.run(
                [
                    "git",
                    "clone",
                    "--depth",
                    "1",
                    "--filter=blob:none",
                    "--sparse",
                    "https://github.com/MyersGroup/relate.git",
                    str(repo),
                ],
                check=True,
            )
        subprocess.run(
            ["git", "-C", str(repo), "sparse-checkout", "set", "include"], check=True
        )
        subprocess.run(
            ["git", "-C", str(repo), "fetch", "--depth", "1", "origin", commit],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(repo), "checkout", "--detach", commit], check=True
        )
        if base != ROOT / ".tools/phylo":
            raise ValueError(
                "Windows build currently requires the default tools directory"
            )
        subprocess.run(
            [str(python), str(ROOT / "scripts/build_relate_windows.py")], check=True
        )
    else:
        p = download(
            "https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz",
            "iqtree-linux.tar.gz",
        )
        (base / "iqtree").mkdir(exist_ok=True)
        untar(p, base / "iqtree")
        p = download(
            "https://myersgroup.github.io/relate/download/relate_v1.2.4_x86_64_static.tgz",
            "relate-linux.tgz",
        )
        (base / "relate-linux").mkdir(exist_ok=True)
        untar(p, base / "relate-linux")
    (base / "downloads.json").write_text(json.dumps(manifest, indent=2) + "\n")
    subprocess.run(
        [uv, "pip", "freeze", "--python", str(python)],
        stdout=(base / "environment.txt").open("w"),
        check=True,
    )


if __name__ == "__main__":
    main()
