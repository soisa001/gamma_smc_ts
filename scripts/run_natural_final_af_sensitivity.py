#!/usr/bin/env python3
"""Run the additive natural final-AF selection sensitivity study."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from gamma_smc_aou.natural_final_af_sensitivity import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
