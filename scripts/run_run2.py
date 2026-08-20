#!/usr/bin/env python
"""Command-line entry point for the run2 selection-detection study.

    python scripts/run_run2.py validate                 # no SLiM binary needed
    python scripts/run_run2.py simulate --slim-bin ...  # 400 replicates
    python scripts/run_run2.py analyze
    python scripts/run_run2.py plot

``validate`` generates the SLiM scripts stdpopsim would run and asserts the run2
invariants against them -- zero rejection sampling, the focal draw one tick after
the pulse, one fitness callback per carrying population, and a reserved focal
base -- without launching anything.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "python"))

from gamma_smc_aou.run2_config import ARMS, MODES, N_REPLICATES  # noqa: E402
from gamma_smc_aou.run2_workflow import PHASES, run_phase  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=PHASES)
    parser.add_argument(
        "--repo-root", default=str(REPO_ROOT), help="repository root (default: inferred)"
    )
    parser.add_argument(
        "--slim-bin",
        default=None,
        help="SLiM 4.2.2 executable; required for simulate and all",
    )
    parser.add_argument(
        "--arm",
        action="append",
        choices=[arm.arm_id for arm in ARMS],
        help="restrict to one arm (repeatable; default: both)",
    )
    parser.add_argument(
        "--mode",
        action="append",
        choices=list(MODES),
        help="restrict to one mode (repeatable; default: both)",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=N_REPLICATES,
        help=f"replicates per arm per mode (default: {N_REPLICATES})",
    )
    parser.add_argument(
        "--index",
        action="append",
        type=int,
        help="run only these replicate indices (repeatable; for smoke tests)",
    )
    parser.add_argument("--workers", type=int, default=1, help="parallel replicates")
    parser.add_argument(
        "--save-trees",
        action="store_true",
        help="keep simulation.trees per replicate (large; off by default)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = run_phase(
        args.phase,
        args.repo_root,
        slim_path=args.slim_bin,
        arm_ids=args.arm,
        modes=tuple(args.mode) if args.mode else MODES,
        replicates=args.replicates,
        workers=args.workers,
        save_trees=args.save_trees,
        indices=args.index,
    )
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
