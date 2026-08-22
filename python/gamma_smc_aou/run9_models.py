"""Extended events for run9's eight arms.

The model itself is run8's, imported unchanged. What varies is only where the
allele comes from, when selection starts, and what frequency it must have reached
by then.

Three conditioning regimes:

``band_at_pulse``
    The run8 rule for pulse arms whose selection begins at the pulse: the allele
    must land in EAS at 2.0-3.0%. Checked one tick *after* the pulse, as a plain
    time, because the allele only enters when the migration fires and a
    ``GenerationAfter`` on both ends of the window collapses it.

``min_frequency_at_onset``
    The delayed rule: no constraint at the pulse at all, but the allele must have
    drifted to at least 2% by the time selection starts, 400 generations ago.

neither
    De novo arms, where a single copy is the definition of the origin.

Every selected arm additionally conditions on survival to the present.
"""

from __future__ import annotations

from typing import Any

import stdpopsim

from .run7_config import (
    ARCHAIC_POPULATION,
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    FOCAL_SITE_ID,
    SLIM_SCALING_FACTOR,
    SPLIT_GENERATIONS,
    TARGET_FREQUENCY_BAND,
)
from .run8_models import (  # noqa: F401  (re-exported so run9 has one import site)
    build_contig,
    build_stdpopsim_model,
    scoped_slim_patch,
    snap_generations,
)
from .run9_config import Scenario

__all__ = [
    "build_contig",
    "build_events",
    "build_stdpopsim_model",
    "scoped_slim_patch",
    "snap_generations",
]


def _frequency_condition(time, op: str, frequency: float) -> Any:
    return stdpopsim.ConditionOnAlleleFrequency(
        start_time=time,
        end_time=time,
        single_site_id=FOCAL_SITE_ID,
        population=EAS_POPULATION,
        op=op,
        allele_frequency=float(frequency),
    )


def build_events(scenario: Scenario, mode: str) -> tuple[Any, ...]:
    """Extended events for one arm. The neutral mode carries none."""
    if mode == "neutral":
        return ()
    if mode != "selected":
        raise ValueError(f"mode must be 'selected' or 'neutral', got {mode!r}")

    onset = snap_generations(scenario.onset_generations)
    events: list[Any] = []

    if scenario.origin == "de_novo":
        events.append(
            stdpopsim.DrawMutation(
                time=onset, single_site_id=FOCAL_SITE_ID, population=EAS_POPULATION
            )
        )
        survival_from = stdpopsim.GenerationAfter(onset)
    elif scenario.origin == "introgressed_pulse":
        events.append(
            stdpopsim.DrawMutation(
                time=stdpopsim.GenerationAfter(snap_generations(SPLIT_GENERATIONS)),
                single_site_id=FOCAL_SITE_ID,
                population=ARCHAIC_POPULATION,
            )
        )
        survival_from = onset - SLIM_SCALING_FACTOR
    else:
        raise ValueError(f"unknown origin {scenario.origin!r}")

    events.append(
        stdpopsim.ChangeMutationFitness(
            start_time=onset,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            selection_coeff=float(scenario.selection_coefficient),
            dominance_coeff=DOMINANCE_COEFFICIENT,
        )
    )

    if scenario.band_at_pulse:
        # One tick after the migration, so the allele has actually arrived.
        check = onset - SLIM_SCALING_FACTOR
        low, high = TARGET_FREQUENCY_BAND
        events.append(_frequency_condition(check, ">=", low))
        events.append(_frequency_condition(check, "<=", high))

    if scenario.min_frequency_at_onset is not None:
        # Evaluated at the onset itself: the allele has been drifting since the
        # pulse, so there is no migration to straddle here.
        events.append(
            _frequency_condition(onset, ">=", scenario.min_frequency_at_onset)
        )

    events.append(
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=survival_from,
            end_time=0.0,
            single_site_id=FOCAL_SITE_ID,
            population=EAS_POPULATION,
            op=">",
            allele_frequency=0.0,
        )
    )
    return tuple(events)
