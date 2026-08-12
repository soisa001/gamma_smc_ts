#!/usr/bin/env python3
"""CLI wrapper for the immutable, provenance-logged replicate campaign."""

import sys

from gamma_smc_aou.eas_introgression_replicates import main

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
