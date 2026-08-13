#!/usr/bin/env python3
"""Write the deterministic focused-campaign RUN_RESULTS.md report."""

from __future__ import annotations

import argparse
from pathlib import Path

from gamma_smc_aou.focused_selection_report import generate_run_report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate the focused selection campaign results report"
    )
    parser.add_argument(
        "--campaign-dir",
        type=Path,
        default=Path("focused_selection_EAS_sim"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    path = generate_run_report(args.campaign_dir, output_path=args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
