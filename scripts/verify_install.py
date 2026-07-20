#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def executable(candidates: list[str | Path]) -> str | None:
    for candidate in candidates:
        if not candidate:
            continue
        path = shutil.which(str(candidate))
        if path:
            return str(Path(path).resolve())
        candidate_path = Path(candidate)
        if candidate_path.is_file():
            return str(candidate_path.resolve())
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify the AoU simulation/decoder installation")
    parser.add_argument("--require-slim", action="store_true")
    parser.add_argument("--require-decoder", action="store_true")
    args = parser.parse_args()

    packages = [
        "gamma-smc-aou",
        "matplotlib",
        "msprime",
        "numpy",
        "pandas",
        "pyslim",
        "scipy",
        "tskit",
        "tszip",
        "zstandard",
    ]
    report: dict[str, object] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {
            name: importlib.metadata.version(name) for name in packages
        },
    }

    slim = executable([
        os.environ.get("SLIM_BIN", ""),
        "slim",
        REPO / ".native" / "bin" / "slim",
        REPO / ".native" / "Library" / "bin" / "slim.exe",
    ])
    report["slim"] = None
    if slim:
        completed = subprocess.run(
            [slim, "-v"], check=True, capture_output=True, text=True
        )
        slim_version = (completed.stdout or completed.stderr).splitlines()[0]
        report["slim"] = {"path": slim, "version": slim_version}
        if "version 5.2" not in slim_version:
            raise RuntimeError(f"expected SLiM 5.2, found: {slim_version}")
    elif args.require_slim:
        raise RuntimeError("SLiM 5.2 not found; set SLIM_BIN or rerun the bootstrap")

    decoder = executable([
        os.environ.get("GAMMA_SMC_BIN", ""),
        REPO / "bin" / "gamma_smc",
        "gamma_smc",
    ])
    report["gamma_smc_decoder"] = decoder
    if args.require_decoder and not decoder:
        raise RuntimeError("compiled Gamma-SMC decoder not found")

    if platform.system() == "Linux":
        cpuinfo = Path("/proc/cpuinfo")
        report["cpu_has_avx2"] = bool(
            cpuinfo.exists() and "avx2" in cpuinfo.read_text(errors="ignore")
        )
        if args.require_decoder and not report["cpu_has_avx2"]:
            raise RuntimeError("Gamma-SMC decoding requires AVX2")

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
