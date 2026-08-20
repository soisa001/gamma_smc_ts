#!/usr/bin/env python
"""Command-line entry point for the run5 selection-detection study.

    python scripts/run_run5.py validate
    python scripts/run_run5.py simulate --slim-bin ... --workers 20
    python scripts/run_run5.py decode  --decoder-bin ... --workers 20
    python scripts/run_run5.py analyze
    python scripts/run_run5.py plot
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "python"))

from gamma_smc_aou.run5_config import ARMS, MODES, N_REPLICATES  # noqa: E402
from gamma_smc_aou.run5_workflow import PHASES, run_phase  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=PHASES)
    parser.add_argument("--repo-root", default=str(REPO_ROOT))
    parser.add_argument("--slim-bin", default=None, help="SLiM 4.2.2 executable")
    parser.add_argument("--decoder-bin", default=None, help="Gamma-SMC executable")
    parser.add_argument(
        "--arm", action="append", choices=[a.arm_id for a in ARMS],
        help="restrict to one arm (repeatable)",
    )
    parser.add_argument(
        "--mode", action="append", choices=list(MODES),
        help="restrict to one mode (repeatable)",
    )
    parser.add_argument("--replicates", type=int, default=N_REPLICATES)
    parser.add_argument("--index", action="append", type=int, help="specific replicates")
    parser.add_argument("--workers", type=int, default=1, help="parallel replicates")
    parser.add_argument("--threads", type=int, default=1, help="decoder threads each")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result = run_phase(
        args.phase,
        args.repo_root,
        slim_path=args.slim_bin,
        decoder_path=args.decoder_bin,
        arm_ids=args.arm,
        modes=tuple(args.mode) if args.mode else MODES,
        replicates=args.replicates,
        workers=args.workers,
        threads=args.threads,
        indices=args.index,
    )
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
