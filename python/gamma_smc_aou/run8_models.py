"""Model and extended events for run8's three scenarios.

The demography is run7's graft, imported unchanged: the inferred PHLASH EAS
history carrying the whole human lineage, a Neanderthal branch of Ne 3,600
joining it at 28,000 generations (700 kya), and a 2.5% pulse at 2,200
generations (55 kya). Nothing about the model varies between scenarios -- only
where the selected allele comes from, when selection starts, and how strong it
is.

Two origins are supported:

``de_novo``
    stdpopsim's default placement, a single new copy in EAS at the onset. Used
    for the recent LCT-like sweep and for the post-introgression hard sweep.

``introgressed_pulse``
    The allele is fixed in the archaic branch just after the split, so the pulse
    delivers it on genuine archaic haplotypes, and it is required to land in
    EAS at 2.0-3.0%.

Both are conditioned on surviving to the present.
"""

from __future__ import annotations

from typing import Any

import stdpopsim

from .run7_config import (
    ARCHAIC_POPULATION,
    GENERATION_TIME_YEARS,
    FOCAL_POSITION_BP,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    FOCAL_SITE_ID,
    SLIM_SCALING_FACTOR,
    SPLIT_GENERATIONS,
    TARGET_FREQUENCY_BAND,
)
from .run7_models import (  # noqa: F401  (re-exported so run8 has one import site)
    build_eas_with_archaic,
    scoped_slim_patch,
    snap_generations,
)
from .run8_config import MUTATION_RATE, Scenario

__all__ = [
    "build_contig",
    "build_events",
    "build_stdpopsim_model",
    "scoped_slim_patch",
    "snap_generations",
]


def build_stdpopsim_model(repo_root, **kwargs) -> stdpopsim.DemographicModel:
    """The run7 graft, declared at run8's mutation rate.

    run7's builder stamps the model with 1.25e-8. Left alone that disagrees with
    the contig's 1.29e-8 and stdpopsim warns on every simulate() call; more to the
    point, the recorded provenance would then contradict what was actually run.
    The demography itself is untouched.
    """
    kwargs.setdefault("add_census", False)
    grafted = build_eas_with_archaic(repo_root, **kwargs)
    return stdpopsim.DemographicModel(
        id="PhlashEASArchaicGraftRun8",
        description="PHLASH EAS median history with a grafted Neanderthal branch",
        long_description=(
            "The inferred EAS history carrying the whole human lineage, with a "
            "Neanderthal branch of constant size joining it at the split and a "
            "single admixture pulse into EAS. Identical to run7's graft; only "
            "the declared mutation rate differs."
        ),
        generation_time=GENERATION_TIME_YEARS,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
        model=grafted.demography,
    )


def build_contig() -> stdpopsim.Contig:
    """run8's contig, at its own mutation rate.

    run7 simulated at 1.25e-8; run8 uses 1.29e-8 so that the simulation matches
    the rate the empirical analysis assumes. Gamma-SMC never sees this value --
    it is handed theta directly -- so the two are independent settings.
    """
    contig = stdpopsim.get_species("HomSap").get_contig(
        length=SEQUENCE_LENGTH_BP,
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
    )
    contig.add_single_site(
        id=FOCAL_SITE_ID,
        coordinate=FOCAL_POSITION_BP,
        description="run8 focal selected site",
    )
    return contig


def _survival(check) -> Any:
    return stdpopsim.ConditionOnAlleleFrequency(
        start_time=check,
        end_time=0.0,
        single_site_id=FOCAL_SITE_ID,
        population=EAS_POPULATION,
        op=">",
        allele_frequency=0.0,
    )


def build_events(scenario: Scenario, mode: str) -> tuple[Any, ...]:
    """Extended events for one scenario. The neutral mode carries none."""
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")

    onset = snap_generations(scenario.onset_generations)
    fitness = stdpopsim.ChangeMutationFitness(
        start_time=onset,
        end_time=0.0,
        single_site_id=FOCAL_SITE_ID,
        population=EAS_POPULATION,
        selection_coeff=float(scenario.selection_coefficient),
        dominance_coeff=DOMINANCE_COEFFICIENT,
    )

    if scenario.origin == "de_novo":
        # A GenerationAfter check is fine here: the mutation exists from the
        # onset tick, so the window does not need to straddle a migration.
        check = stdpopsim.GenerationAfter(onset)
        return (
            stdpopsim.DrawMutation(
                time=onset, single_site_id=FOCAL_SITE_ID, population=EAS_POPULATION
            ),
            fitness,
            _survival(check),
        )

    if scenario.origin != "introgressed_pulse":
        raise ValueError(f"unknown origin {scenario.origin!r}")

    # The frequency check must sit one tick *after* the pulse, as a plain time:
    # the allele only enters EAS when the migration fires, and a GenerationAfter
    # on both ends of the window collapses it so the check never runs.
    check = onset - SLIM_SCALING_FACTOR
    low, high = TARGET_FREQUENCY_BAND
    return (
        stdpopsim.DrawMutation(
            time=stdpopsim.GenerationAfter(snap_generations(SPLIT_GENERATIONS)),
            single_site_id=FOCAL_SITE_ID,
            population=ARCHAIC_POPULATION,
        ),
        fitness,
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=check,
            end_time=check,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op=">=",
            allele_frequency=float(low),
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=check,
            end_time=check,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op="<=",
            allele_frequency=float(high),
        ),
        _survival(check),
    )
