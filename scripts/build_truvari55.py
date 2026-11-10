# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = []
# ///
"""Build Truvari 5.5 and native wheels with uv; no sudo or controlled data access.

Run with uv --no-config run --script scripts/build_truvari55.py --out PATH.
Then PATH/venv/bin/python scripts/compare_stringent_recollapse.py ...
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
from urllib.request import urlopen
import zipfile


# Published source distributions, verified against their PyPI SHA256 metadata.
SOURCES = {
    "truvari": ("https://files.pythonhosted.org/packages/29/c8/5ab1ba49a3a6eba0cb568a22f09da8d90a963582b7b1f7187c73109efc21/truvari-5.5.0.tar.gz",
                "58aa2dc02de37597ba57ae18a2b32c662157e0e03d65f584a5b75e525a402d50"),
    "pywfa": ("https://files.pythonhosted.org/packages/93/40/ec1c77237515eb618ba7407d4ac954b9dd127808a303cd85b0928f8dc12a/pywfa-0.5.1.tar.gz",
              "e972bf53f9e6d8957e9105ecc22cf704ac4bfad4d882d79c82f11fc260381483"),
    "pytabix": ("https://files.pythonhosted.org/packages/84/6a/520ecf75c2ada77492cb4ed21fb22aed178e791df434ca083b59fffadddd/pytabix-0.1.tar.gz",
               "0774f1687ebd41811fb07a0e50951b6be72d7cc7e22ed2b18972eaf7482eb7d1"),
    "zlib": ("https://zlib.net/fossils/zlib-1.3.1.tar.gz",
             "9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23"),
}


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def save(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def fetch(name, root):
    url, expected = SOURCES[name]
    path = root / url.rsplit("/", 1)[-1]
    if not path.exists() or digest(path) != expected:
        tmp = path.with_suffix(".download")
        with urlopen(url, timeout=120) as source, tmp.open("wb") as target:
            shutil.copyfileobj(source, target)
        if digest(tmp) != expected:
            raise RuntimeError(f"Download checksum mismatch: {name}")
        tmp.replace(path)
    return path


def execute(command, log, **kwargs):
    with log.open("a") as handle:
        handle.write("COMMAND " + json.dumps([str(x) for x in command]) + "\n")
        handle.flush()
        result = subprocess.run([str(x) for x in command], stdout=handle, stderr=handle, **kwargs)
    if result.returncode:
        raise RuntimeError(f"Command failed; see {log}\n" + log.read_text(errors="replace")[-6000:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12) or sys.platform != "linux":
        parser.error("Run this builder on Linux with uv-managed Python 3.12")
    if args.workers < 1:
        parser.error("workers must be positive")
    uv = shutil.which("uv")
    for command in ("uv", "cc", "make", "ar"):
        if not shutil.which(command):
            parser.error(f"Required build tool not found: {command}")
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    sources = root / "sources"
    sources.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    build_identity = dict(python=sys.version, machine=platform.machine(), host=platform.node(),
                          compiler=subprocess.check_output(["cc", "--version"], text=True).splitlines()[0])

    print("Phase 1: build/check a local static zlib (provides zlib.h; no sudo).", flush=True)
    prefix = root / "zlib"
    receipt = root / "zlib.json"
    zfiles = ("include/zlib.h", "include/zconf.h", "lib/libz.a")
    recipe = dict(source=SOURCES["zlib"], build=build_identity, flags="-O2 -fPIC --static")
    cached = json.loads(receipt.read_text()) if receipt.exists() else {}
    if cached.get("recipe") != json.loads(json.dumps(recipe)) or not all(
        (prefix / f).is_file() and digest(prefix / f) == cached.get("hashes", {}).get(f) for f in zfiles
    ):
        archive = fetch("zlib", sources)
        with tarfile.open(archive) as tar:
            tar.extractall(root / "native_source", filter="data")
        native = root / "native_source/zlib-1.3.1"
        zenv = dict(environment, CFLAGS="-O2 -fPIC")
        execute(["./configure", "--static", f"--prefix={prefix}"], root / "zlib.log", cwd=native, env=zenv)
        execute(["make", f"-j{args.workers}"], root / "zlib.log", cwd=native, env=zenv)
        execute(["make", "install"], root / "zlib.log", cwd=native, env=zenv)
        save(receipt, dict(recipe=recipe, hashes={f: digest(prefix / f) for f in zfiles}))
    for key, subdir in (("CPATH", "include"), ("LIBRARY_PATH", "lib")):
        environment[key] = str(prefix / subdir) + (os.pathsep + environment[key] if environment.get(key) else "")

    print("Phase 2: uv build --wheel for Truvari 5.5.0, pywfa, and pytabix.", flush=True)
    def build(name):
        folder = root / "wheels" / name
        folder.mkdir(parents=True, exist_ok=True)
        receipt = folder / "build.json"
        recipe = dict(source=SOURCES[name], build=build_identity, zlib=digest(prefix / "lib/libz.a"))
        cached = json.loads(receipt.read_text()) if receipt.exists() else {}
        wheel = folder / cached.get("wheel", "missing")
        if cached.get("recipe") == json.loads(json.dumps(recipe)) and wheel.is_file() and digest(wheel) == cached["sha256"]:
            print(f"Reusing verified {wheel.name}", flush=True)
            return wheel
        archive = fetch(name, sources)
        execute([uv, "--no-config", "build", "--wheel", "--python", sys.executable,
                 "--out-dir", folder, archive], folder / "build.log", env=environment)
        wheels = list(folder.glob("*.whl"))
        if len(wheels) != 1:
            raise RuntimeError(f"Expected exactly one built wheel in {folder}")
        wheel = wheels[0]
        with zipfile.ZipFile(wheel) as contents:
            if contents.testzip() is not None:
                raise RuntimeError(f"Corrupt wheel: {wheel}")
            if name == "pywfa" and not any("align.cpython-312" in f and f.endswith(".so") for f in contents.namelist()):
                raise RuntimeError("pywfa wheel lacks its Python 3.12 native extension")
        save(receipt, dict(recipe=recipe, wheel=wheel.name, sha256=digest(wheel)))
        return wheel

    with ThreadPoolExecutor(max_workers=min(args.workers, 3)) as pool:
        wheels = list(pool.map(build, ("truvari", "pywfa", "pytabix")))
    print("Phase 3: install built wheels with uv and check imports/CLI.", flush=True)
    venv = root / "venv"
    python = venv / "bin/python"
    if not python.exists():
        execute([uv, "--no-config", "venv", "--python", sys.executable, venv], root / "install.log")
    execute([uv, "--no-config", "pip", "install", "--python", python, *wheels], root / "install.log", env=environment)
    execute([python, "-c", "import truvari, tabix; from pywfa.align import WavefrontAligner; "
             "assert truvari.__version__ == '5.5.0'; print(truvari.VariantRecord); print(truvari.__version__)"],
            root / "install.log", env=environment)
    execute([python, "-m", "truvari", "collapse", "--help"], root / "collapse_help.txt", env=environment)
    installed = subprocess.check_output([uv, "--no-config", "pip", "freeze", "--python", python], text=True)
    (root / "installed.txt").write_text(installed)
    save(root / "build_manifest.json", dict(utc=datetime.now(timezone.utc).isoformat(), sources=SOURCES,
         build=build_identity, workers=args.workers, script_sha256=digest(__file__),
         uv=subprocess.check_output([uv, "--version"], text=True).strip(),
         wheels={str(w): digest(w) for w in wheels}, installed=installed))
    print(f"Ready: Truvari 5.5.0 with {python}", flush=True)


if __name__ == "__main__":
    main()
