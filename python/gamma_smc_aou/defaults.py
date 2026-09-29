"""Shared scientific and scan defaults for every Python entry point."""

DEFAULT_SCALED_MUTATION_RATE = 0.00075
DEFAULT_RECOMBINATION_TO_MUTATION_RATIO = 0.8
DEFAULT_MUTATION_RATE = 1.25e-8
DEFAULT_OUTPUT_STRIDE = 10_000
DEFAULT_CACHE_SIZE = 1_000
DEFAULT_GENERATION_TIME = 25

# With DEFAULT_SCALED_MUTATION_RATE and DEFAULT_MUTATION_RATE above,
# 2Ne = theta / (2 * mu) = 30,000 generations.
DEFAULT_DIPLOID_COALESCENT_SCALE = (
    DEFAULT_SCALED_MUTATION_RATE / (2 * DEFAULT_MUTATION_RATE)
)

# The scan reports one block of statistics per threshold. Two defaults:
# 10,000 years is the EPAS1-like recent-adaptation scale, and 50,000 years
# places the signal at the onset of archaic introgression. At 25 years per
# generation and 2Ne = 30,000 these have neutral P(T < t) of 1.32% and 6.45%,
# so their candidate screens are not interchangeable.
DEFAULT_THRESHOLD_YEARS = (10_000.0, 50_000.0)

# One decode evaluates every threshold from the same posteriors, so the extra
# cost per threshold is one table lookup per pair per position. The cap keeps
# the bit matrix (one bit per pair/position/threshold) and the summary width
# bounded.
MAX_THRESHOLDS = 5
