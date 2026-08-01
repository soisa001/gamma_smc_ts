"""AoU within-individual Gamma-SMC workflow."""

import os

# Workbench notebooks export matplotlib's optional inline backend into child
# processes. This package only writes figures to files, so select the portable
# headless backend before any eagerly imported module imports matplotlib.
os.environ["MPLBACKEND"] = "Agg"

from .calibration import calibrate_sites, monte_carlo_pvalue
from .simulation import SimulationConfig, simulate_replicates

__all__ = [
    "SimulationConfig",
    "calibrate_sites",
    "monte_carlo_pvalue",
    "simulate_replicates",
]
