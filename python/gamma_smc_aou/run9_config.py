"""Frozen parameters for run9: stronger selection, and a delayed EPAS1-like onset.

run8 established that at s = 0.002 an introgressed sweep is invisible without
carrier conditioning, and that s = 0.001 is dead at every cutoff. run9 pushes the
coefficient up and adds a second question: does it help if selection *starts
recently* rather than at the pulse?

Demography, samples, pairs, rates and decoding are all identical to run8 and are
imported rather than restated, so every difference below is a difference in the
selection model alone.

Arms
----
Eight arms, crossing two origins with four selection regimes:

=============================  ==========  =========  ======  ====================
arm                            origin      onset      s       conditioning
=============================  ==========  =========  ======  ====================
``pulse_s0p003``               2.5% pulse  55 kya     0.003   2-3% at pulse
``pulse_s0p004``               2.5% pulse  55 kya     0.004   2-3% at pulse
``pulse_s0p005``               2.5% pulse  55 kya     0.005   2-3% at pulse
``pulse_delayed_s0p005``       2.5% pulse  **10 kya** 0.005   **>= 2% at onset**
``denovo_s0p003``              de novo     55 kya     0.003   --
``denovo_s0p004``              de novo     55 kya     0.004   --
``denovo_s0p005``              de novo     55 kya     0.005   --
``denovo_delayed_s0p005``      de novo     **10 kya** 0.005   --
=============================  ==========  =========  ======  ====================

Every arm additionally conditions on the allele surviving to the present. The
neutral arm carries no allele and no conditioning at all.

The delayed arms
----------------
The allele still arrives on the 55 kya pulse, then **drifts neutrally for 1,800
generations**, and selection only begins 400 generations ago. The minimum
frequency is enforced at that onset rather than at the pulse, which is what makes
it EPAS1-like: an archaic haplotype that sat at appreciable frequency for a long
time and only recently became advantageous.

This is expensive, and knowingly so. Neutral drift from p = 0.025 over 1,800
generations at a harmonic Ne near 4,400 spans tau ~ 0.20 coalescent units, giving
a diffusion spread of sqrt(p(1-p)tau) ~ 0.07 -- roughly three times the mean -- so
most trajectories drift low or are lost and get rejected. Worse, stdpopsim
checkpoints at ``DrawMutation``, which for a pulse arm sits back at the archaic
split, so each rejection replays the whole simulation rather than a short tail.
A 2% floor was chosen over a cheaper 1% because the floor is the thing that makes
the scenario EPAS1-like rather than a lottery over whatever drift happened to
leave behind.

The de novo delayed arm needs no frequency floor: a de novo allele starts at one
copy by definition, so survival conditioning is the whole of its ascertainment.

Analysis
--------
run8 summarised each arm by its single AF-median replicate, which was its weakest
part. run9 reports the **quantile curve** of the statistic instead -- the p-value
at the 10th, 25th, 50th, 75th and 90th percentile -- and stratifies everything by
final allele frequency. Pooling all eight arms spans a wide AF range, from which
the portable quantity can be estimated: the frequency above which unconditional
detection starts to work.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .run7_config import (
    ADMIXTURE_PROPORTION,
    ARCHAIC_EFFECTIVE_SIZE,
    GENERATION_TIME_YEARS,
    PULSE_GENERATIONS,
    RECOMBINATION_RATE,
    SPLIT_GENERATIONS,
    TARGET_FREQUENCY_BAND,
)
from .run8_config import (  # noqa: F401  (re-exported: run9 changes none of these)
    ALPHA,
    DECODER_SETTINGS,
    FOCAL_POSITION_BP,
    MODES,
    MUTATION_RATE,
    N_DECODE_PAIRS,
    PAIRS_SEED,
    RHO,
    SAMPLE_DIPLOIDS,
    THETA,
)

__all__ = [
    "ALPHA",
    "DECODER_SETTINGS",
    "MODES",
    "N_NEUTRAL",
    "QUANTILES",
    "SCENARIOS",
    "SCENARIOS_BY_ID",
    "Scenario",
    "THRESHOLDS_YEARS",
    "arm_config_record",
    "get_scenario",
    "write_config",
]

SEED_BASE = 20280101

#: 300 rather than run8's 100: the p-value floor is 1/(n+1), and run8's strongest
#: arm sat exactly on 1/101 = 0.0099 with no way to see past it. 300 gives 0.0033.
N_NEUTRAL = 300

#: 75 kya is new. Introgression peaked at 40-50 kya in run8, at the edge of that
#: grid, so the far side of the peak was never observed.
THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    5_000.0,
    10_000.0,
    20_000.0,
    30_000.0,
    40_000.0,
    50_000.0,
    75_000.0,
)

#: Reported in place of run8's single median replicate.
QUANTILES: tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)

#: 10 kya at 25 years per generation.
DELAYED_ONSET_GENERATIONS = 400.0
#: Minimum EAS frequency required when delayed selection begins.
DELAYED_MIN_FREQUENCY = 0.02


@dataclass(frozen=True)
class Scenario:
    scenario_id: str
    label: str
    origin: str  # "de_novo" or "introgressed_pulse"
    onset_generations: float
    selection_coefficient: float
    n_replicates: int
    #: Enforce the 2-3% band at the pulse (the run8 rule).
    band_at_pulse: bool
    #: Enforce a minimum frequency at the onset instead (the delayed rule).
    min_frequency_at_onset: float | None
    rationale: str

    @property
    def onset_years(self) -> float:
        return self.onset_generations * GENERATION_TIME_YEARS

    @property
    def places_allele_in_archaic(self) -> bool:
        return self.origin == "introgressed_pulse"

    @property
    def is_delayed(self) -> bool:
        return self.onset_generations < PULSE_GENERATIONS

    def seed(self, mode: str, replicate_index: int, scenario_index: int) -> int:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        return (
            SEED_BASE
            + MODES.index(mode) * 1_000_000
            + scenario_index * 100_000
            + replicate_index
        )


def _pulse(coefficient: float, index: int) -> Scenario:
    return Scenario(
        scenario_id=f"pulse_s{coefficient:g}".replace(".", "p"),
        label=f"Selected introgressed allele, 55 kya pulse, s = {coefficient:g}",
        origin="introgressed_pulse",
        onset_generations=PULSE_GENERATIONS,
        selection_coefficient=coefficient,
        n_replicates=100,
        band_at_pulse=True,
        min_frequency_at_onset=None,
        rationale=(
            "run8 showed s = 0.002 is invisible unconditionally and s = 0.001 is "
            "dead. Raising s raises final frequency, and since only about AF^2 of "
            "random pairs are carrier-carrier, the unconditional statistic should "
            "scale steeply with it."
        ),
    )


def _denovo(coefficient: float, index: int) -> Scenario:
    return Scenario(
        scenario_id=f"denovo_s{coefficient:g}".replace(".", "p"),
        label=f"De novo sweep at the pulse time, 55 kya, s = {coefficient:g}",
        origin="de_novo",
        onset_generations=PULSE_GENERATIONS,
        selection_coefficient=coefficient,
        n_replicates=100,
        band_at_pulse=False,
        min_frequency_at_onset=None,
        rationale=(
            "Matched to the pulse arm in onset and coefficient so that origin is "
            "the only difference. A de novo allele starts at 1/(2N) rather than "
            "0.025, roughly 190-fold lower, which the softness formula predicts "
            "makes its sweep far harder."
        ),
    )


def _denovo_delayed(coefficient: float) -> Scenario:
    """The delayed pulse arm's counterpart, at a coefficient that can actually move.

    A smoke test at s = 0.005 reached AF 0.010 with no carrier-carrier pairs at
    all, which is the arithmetic rather than bad luck: over 400 generations
    additive selection multiplies the odds by only exp(s t / 2) ~ 2.7, so one copy
    goes nowhere. Two coefficients bracket the useful range instead -- s = 0.01
    should land near 10%, and s = 0.05 has the same ``alpha * t`` as run7's
    LCT-like arm and should reach roughly 0.4-0.5.

    No frequency floor is imposed or needed: a de novo allele starts at one copy,
    so surviving drift is the whole of its ascertainment.
    """
    return Scenario(
        scenario_id=f"denovo_delayed_s{coefficient:g}".replace(".", "p"),
        label=f"De novo sweep from 10 kya, s = {coefficient:g}",
        origin="de_novo",
        onset_generations=DELAYED_ONSET_GENERATIONS,
        selection_coefficient=coefficient,
        n_replicates=100,
        band_at_pulse=False,
        min_frequency_at_onset=None,
        rationale=(
            "Matched in onset to the delayed pulse arm so the contrast isolates "
            "the head start that 2.5% archaic ancestry provides over a single "
            "new copy."
        ),
    )


SCENARIOS: tuple[Scenario, ...] = (
    _pulse(0.003, 0),
    _pulse(0.004, 1),
    _pulse(0.005, 2),
    Scenario(
        scenario_id="pulse_delayed_s0p005",
        label="Introgressed allele, drift then selection from 10 kya, s = 0.005",
        origin="introgressed_pulse",
        onset_generations=DELAYED_ONSET_GENERATIONS,
        selection_coefficient=0.005,
        n_replicates=100,
        band_at_pulse=False,
        min_frequency_at_onset=DELAYED_MIN_FREQUENCY,
        rationale=(
            "EPAS1-like: the allele arrives on the pulse, drifts neutrally for "
            "1,800 generations, and only becomes advantageous 400 generations "
            "ago, by which point it must be at least 2%. Tests whether a recent "
            "onset concentrates coalescence where the neutral baseline is low."
        ),
    ),
    _denovo(0.003, 4),
    _denovo(0.004, 5),
    _denovo(0.005, 6),
    _denovo_delayed(0.01),
    _denovo_delayed(0.05),
)

SCENARIOS_BY_ID = {scenario.scenario_id: scenario for scenario in SCENARIOS}


def get_scenario(scenario_id: str) -> Scenario:
    try:
        return SCENARIOS_BY_ID[scenario_id]
    except KeyError:
        raise ValueError(
            f"unknown scenario {scenario_id!r}; expected one of {sorted(SCENARIOS_BY_ID)}"
        ) from None


def shared_parameter_record() -> dict[str, Any]:
    return {
        "inherits": "run8, unchanged: demography, samples, pairs, rates, decoding",
        "demography": {
            "split_generations": SPLIT_GENERATIONS,
            "pulse_generations": PULSE_GENERATIONS,
            "archaic_effective_size": ARCHAIC_EFFECTIVE_SIZE,
            "admixture_proportion": ADMIXTURE_PROPORTION,
        },
        "simulation": {
            "mutation_rate_per_bp_per_generation": MUTATION_RATE,
            "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
            "generation_time_years": GENERATION_TIME_YEARS,
            "sample_diploids": SAMPLE_DIPLOIDS,
            "sample_haplotypes": 2 * SAMPLE_DIPLOIDS,
        },
        "pairs": {"n_decode_pairs": N_DECODE_PAIRS, "pairs_seed": PAIRS_SEED},
        "decoding": {
            "theta": THETA,
            "rho": RHO,
            "recent_call": DECODER_SETTINGS["recent_call"],
            "implied_constant_ne": THETA / (4.0 * MUTATION_RATE),
        },
        "conditioning": {
            "pulse_band": list(TARGET_FREQUENCY_BAND),
            "delayed_min_frequency_at_onset": DELAYED_MIN_FREQUENCY,
            "delayed_onset_generations": DELAYED_ONSET_GENERATIONS,
            "all_selected_arms_condition_on_survival": True,
            "neutral_conditioning": "none",
        },
        "analysis": {
            "thresholds_years": list(THRESHOLDS_YEARS),
            "quantiles": list(QUANTILES),
            "alpha": ALPHA,
            "n_neutral_replicates": N_NEUTRAL,
            "p_value_floor": 1.0 / (N_NEUTRAL + 1.0),
            "significance_rule": (
                "at most one neutral replicate may match or exceed the observation"
            ),
        },
    }


def arm_config_record() -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run9-config/v1",
        "shared": shared_parameter_record(),
        "scenarios": [
            {
                "scenario_id": s.scenario_id,
                "label": s.label,
                "origin": s.origin,
                "onset_generations": s.onset_generations,
                "onset_years": s.onset_years,
                "selection_coefficient": s.selection_coefficient,
                "band_at_pulse": s.band_at_pulse,
                "min_frequency_at_onset": s.min_frequency_at_onset,
                "n_replicates": s.n_replicates,
                "rationale": s.rationale,
            }
            for s in SCENARIOS
        ],
    }
    json.dumps(record)
    return record


def write_config(destination: str | Path) -> Path:
    directory = Path(destination)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "run9_scenarios.json"
    path.write_text(
        json.dumps(arm_config_record(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
