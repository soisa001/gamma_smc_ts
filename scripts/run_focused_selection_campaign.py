#!/usr/bin/env python3
"""CLI wrapper for the focused selected-versus-neutral EAS campaign."""

import sys

from gamma_smc_aou.focused_selection_campaign import main

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
