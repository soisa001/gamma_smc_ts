"""Coalescent predictions for run8, derived and then checked against the sims.

Everything here is worked out from standard theory *before* looking at the
simulated numbers, so that the comparison is a test rather than a description.
Four predictions, in the order they constrain the study.

1. The neutral baseline
-----------------------
For a single population with piecewise-constant size ``N(tau)`` the pairwise
coalescent rate at time ``tau`` back is ``1 / (2 N(tau))``, so

    P(T < t) = 1 - exp( - integral_0^t dtau / (2 N(tau)) )

This is exact for the PHLASH trajectory, which *is* a step function, so the
integral is a finite sum with no approximation. It fixes the null the whole study
is scored against, and disagreement would mean the pipeline is wrong rather than
the biology being interesting.

2. Coalescence during a sweep -- why strong selection gives a *weaker* signal
----------------------------------------------------------------------------
Take an additive allele, so the odds grow as ``exp(alpha t)`` with
``alpha = h s = s / 2``. Writing the logistic trajectory in terms of frequency,

    dt = dp / (alpha p (1 - p))

Two carrier lineages coalesce at rate ``1 / (2 N p)`` while the carrier class sits
at frequency ``p`` -- the carriers are their own small subpopulation. The expected
number of pairwise coalescent events across the sweep is therefore

    I = integral dt / (2 N p) = 1/(2 N alpha) * integral_{p0}^{pf} dp / (p^2 (1-p))

and since ``integral dp / (p^2 (1-p)) = -1/p + log(p/(1-p))``, for small ``p0``
the lower limit dominates:

    I  ~=  1 / (2 N alpha p0)  =  1 / (N s p0)

so the probability that a carrier pair coalesces during the sweep is

    P_sweep  ~=  1 - exp( -1 / (N s p0) )

**This decreases in s.** A strongly selected allele leaves the low-frequency
regime too quickly for its carriers to find common ancestors, so more founding
haplotypes survive -- a softer sweep with deeper carrier TMRCA. That is the
mechanism behind run7's non-monotonic power, stated as a formula.

Note ``N s`` is invariant under SLiM's rescaling: ``(N/Q)(sQ) = N s``. So this
prediction is unaffected by ``Q = 5``, which makes it a fair check.

The same formula separates the two origins. Introgression starts at
``p0 = m = 0.025``; a de novo allele starts at ``p0 = 1/(2N)``, which is smaller
by three orders of magnitude, so ``I`` is enormous and essentially every carrier
pair coalesces during the sweep. De novo sweeps are hard; pulse-fed sweeps are
soft by construction.

3. Sweep duration and final frequency
-------------------------------------
Integrating the same logistic gives the time to travel from ``p0`` to ``pf``:

    t = (1 / alpha) * log[ (pf / (1 - pf)) * ((1 - p0) / p0) ]

Conditioning on survival matters: a lineage that escapes drift does so by
reaching roughly ``1 / alpha`` copies, so the trajectory behaves as though it
started at ``p_eff ~= 1 / (2 N alpha)`` rather than at a single copy. Both the
raw and establishment-corrected predictions are reported, because the gap between
them is itself informative.

4. Hard floors on P(TMRCA < x)
------------------------------
Some of these statistics are exactly zero, not merely small, and that is worth
stating as a prediction because it is easy to check.

* A **heterozygote** at an introgressed locus carries one archaic and one modern
  haplotype. They cannot meet more recently than the human-Neanderthal split, so
  ``P(TMRCA < x) = 0`` for every ``x < 700 kya``.
* A **heterozygote** at a de novo locus carries one derived and one ancestral
  haplotype, which cannot meet more recently than the mutation, so
  ``P(TMRCA < x) = 0`` for every ``x`` below the onset.
* A **homozygous carrier** of an introgressed allele has both haplotypes tracing
  through the pulse, so any pair that does not share a post-pulse ancestor is
  pushed back into the archaic branch. Below the pulse the statistic sees only
  sweep-driven coalescence -- the admixture floor.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .eas_sweep_models import build_eas_demography_models, load_phlash_eas_npz
from .run7_config import (
    ADMIXTURE_PROPORTION,
    DOMINANCE_COEFFICIENT,
    GENERATION_TIME_YEARS,
    PHLASH_EAS_RELATIVE_PATH,
    PHLASH_EAS_SHA256,
    PULSE_GENERATIONS,
    SPLIT_GENERATIONS,
)

__all__ = [
    "NeTrajectory",
    "load_trajectory",
    "neutral_p_tmrca_below",
    "sweep_coalescence_probability",
    "sweep_duration_generations",
    "establishment_frequency",
    "predicted_final_frequency",
]


@dataclass(frozen=True)
class NeTrajectory:
    """A piecewise-constant Ne, as PHLASH reports it."""

    time_generations: np.ndarray  # left edge of each epoch, ascending from 0
    ne: np.ndarray

    def size_at(self, generations: float | np.ndarray) -> np.ndarray:
        index = np.searchsorted(self.time_generations, generations, side="right") - 1
        index = np.clip(index, 0, self.ne.size - 1)
        return self.ne[index]

    def harmonic_mean_to(self, generations: float) -> float:
        """Harmonic mean of Ne over [0, t], which is what sets coalescence."""
        integral = self.coalescent_intensity(generations)
        return float(generations / (2.0 * integral)) if integral > 0 else float("inf")

    def coalescent_intensity(self, generations: float) -> float:
        """``integral_0^t dtau / (2 N(tau))``, computed exactly on the steps."""
        if generations <= 0:
            return 0.0
        edges = self.time_generations
        total = 0.0
        for i in range(edges.size):
            lo = edges[i]
            hi = edges[i + 1] if i + 1 < edges.size else np.inf
            if lo >= generations:
                break
            span = min(hi, generations) - lo
            if span > 0:
                total += span / (2.0 * self.ne[i])
        return float(total)


def load_trajectory(repo_root: str | Path, quantile: str = "median") -> NeTrajectory:
    artifact = load_phlash_eas_npz(
        Path(repo_root) / PHLASH_EAS_RELATIVE_PATH, expected_sha256=PHLASH_EAS_SHA256
    )
    model = build_eas_demography_models(artifact)[quantile]
    return NeTrajectory(
        time_generations=np.asarray(model.time_generations, dtype=float),
        ne=np.asarray(model.ne, dtype=float),
    )


def neutral_p_tmrca_below(
    trajectory: NeTrajectory,
    thresholds_years: Sequence[float],
    *,
    generation_time: float = GENERATION_TIME_YEARS,
) -> np.ndarray:
    """Exact ``P(T < x)`` for a pair under the piecewise-constant history."""
    out = []
    for years in thresholds_years:
        generations = float(years) / generation_time
        out.append(1.0 - np.exp(-trajectory.coalescent_intensity(generations)))
    return np.asarray(out, dtype=float)


def sweep_coalescence_probability(
    effective_size: float,
    selection_coefficient: float,
    initial_frequency: float,
    *,
    final_frequency: float = 0.99,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> float:
    """``1 - exp(-I)`` with ``I`` the pairwise coalescent intensity over the sweep.

    The integral is evaluated in full rather than truncated to ``1/(N s p0)``, so
    the approximation quoted in the module docstring can be checked against it.
    """
    alpha = dominance * selection_coefficient
    if alpha <= 0 or not 0 < initial_frequency < final_frequency < 1:
        raise ValueError("need 0 < p0 < pf < 1 and a positive selection coefficient")

    def antiderivative(p: float) -> float:
        return -1.0 / p + np.log(p / (1.0 - p))

    intensity = (
        antiderivative(final_frequency) - antiderivative(initial_frequency)
    ) / (2.0 * effective_size * alpha)
    return float(1.0 - np.exp(-max(intensity, 0.0)))


def sweep_duration_generations(
    selection_coefficient: float,
    initial_frequency: float,
    final_frequency: float,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> float:
    alpha = dominance * selection_coefficient
    odds_final = final_frequency / (1.0 - final_frequency)
    odds_initial = initial_frequency / (1.0 - initial_frequency)
    return float(np.log(odds_final / odds_initial) / alpha)


def establishment_frequency(
    effective_size: float,
    selection_coefficient: float,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> float:
    """``p_eff = 1 / (2 N alpha)``: the frequency a surviving lineage escapes at."""
    alpha = dominance * selection_coefficient
    return float(1.0 / (2.0 * effective_size * alpha))


def predicted_final_frequency(
    selection_coefficient: float,
    initial_frequency: float,
    elapsed_generations: float,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> float:
    """Logistic forward projection of the odds."""
    alpha = dominance * selection_coefficient
    odds = initial_frequency / (1.0 - initial_frequency)
    odds *= np.exp(alpha * elapsed_generations)
    return float(odds / (1.0 + odds))


#: Hard floors: below these, the statistic is exactly zero for that class.
FLOORS_GENERATIONS = {
    "introgressed_het": SPLIT_GENERATIONS,
    "introgressed_hom_carrier": PULSE_GENERATIONS,
    "admixture_proportion": ADMIXTURE_PROPORTION,
}


def logistic_frequency_forward(
    selection_coefficient: float,
    initial_frequency: float,
    generations_since_onset: np.ndarray,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> np.ndarray:
    """Deterministic additive trajectory, odds growing as ``exp(alpha tau)``."""
    alpha = dominance * selection_coefficient
    odds = initial_frequency / (1.0 - initial_frequency)
    odds = odds * np.exp(alpha * np.asarray(generations_since_onset, dtype=float))
    return odds / (1.0 + odds)


def sweep_coalescence_intensity(
    trajectory: NeTrajectory,
    selection_coefficient: float,
    initial_frequency: float,
    onset_generations: float,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
    steps: int = 20_000,
) -> float:
    """Pairwise coalescent intensity for carriers, over the real ``N(t)``.

    The closed form ``1/(N s p0)`` assumes a constant population. These histories
    grow by an order of magnitude across the sweep -- PHLASH puts EAS at 3,869
    around the pulse and 37,637 today -- and since the carrier coalescent rate is
    ``1 / (2 N(t) p(t))``, holding ``N`` at its onset value overstates coalescence
    badly. Integrating over the actual trajectory is what makes the prediction
    comparable to the simulated post-pulse fractions.

    Time runs backwards from the present to the onset; ``p`` is the carrier
    frequency then, so early in the sweep (near the onset) the carrier pool is
    small and contributes most of the intensity.
    """
    if onset_generations <= 0:
        return 0.0
    back = np.linspace(0.0, float(onset_generations), int(steps))
    forward = float(onset_generations) - back
    frequency = logistic_frequency_forward(
        selection_coefficient, initial_frequency, forward, dominance=dominance
    )
    frequency = np.clip(frequency, 1.0 / (2.0 * trajectory.ne.max()), 1.0)
    size = trajectory.size_at(back)
    integrand = 1.0 / (2.0 * size * frequency)
    return float(np.trapezoid(integrand, back))


def sweep_coalescence_probability_trajectory(
    trajectory: NeTrajectory,
    selection_coefficient: float,
    initial_frequency: float,
    onset_generations: float,
    *,
    dominance: float = DOMINANCE_COEFFICIENT,
) -> float:
    """``1 - exp(-I)`` using the real ``N(t)``; comparable to the measured fraction."""
    intensity = sweep_coalescence_intensity(
        trajectory,
        selection_coefficient,
        initial_frequency,
        onset_generations,
        dominance=dominance,
    )
    return float(1.0 - np.exp(-intensity))
