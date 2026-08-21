"""Frozen parameters for run8: three selection scenarios on one EAS model.

run8 keeps run7's demography unchanged -- the inferred PHLASH EAS history with a
Neanderthal branch grafted on -- and asks three questions with it, each at the
cutoff its own timescale is expected to favour.

Demography (identical to run7, imported rather than restated)
-------------------------------------------------------------
=========================  ==========  =========================================
quantity                   value       source
=========================  ==========  =========================================
human-Neanderthal split    700 kya     Prufer et al. 2014 (550-765 ky)
                           28,000 gen
introgression pulse        55 kya      Sankararaman 2012; Fu 2014 (47-65 ky)
                           2,200 gen
archaic Ne                 3,600       Ragsdale & Gravel 2019
admixture                  2.5%        Prufer 2014; Vernot & Akey 2015
=========================  ==========  =========================================

Scenarios
---------
===========================  ==========  =========  ======  =================
scenario                     origin      onset      s       expected cutoff
===========================  ==========  =========  ======  =================
``recent_denovo``            de novo     5 kya      0.1     ~5 kya
``postintro_denovo``         de novo     55 kya     0.002   ~30 kya
``introgressed_pulse``       2.5% pulse  55 kya     0.002   ~50 kya
``introgressed_pulse_weak``  2.5% pulse  55 kya     0.001   ~50 kya
===========================  ==========  =========  ======  =================

The "expected cutoff" column is a prediction, not a setting: every scenario is
scored on the *same* seven-threshold grid, and which threshold wins is a result.
run7 established that the detectable window tracks when the sweep was
compressing coalescence, so a sweep that finished long ago is invisible at recent
cutoffs and an introgressed allele cannot produce coalescence more recent than
the pulse that delivered it.

``postintro_denovo`` and ``introgressed_pulse`` share an onset *and* a selection
coefficient, so the only thing separating them is where the allele came from: one
new copy on a modern background against 2.5% carried on archaic haplotypes. The
softness formula predicts a large gap, since ``p0`` differs by roughly 190-fold.

Pairs
-----
The empirical scan draws 100,000 random haplotype pairs from a panel of
800-4,000 haplotypes. run8 mirrors that design at simulation scale: 100 diploids
= 200 haplotypes, from which the decoder draws **10,000 random pairs** at a fixed
seed. Tree truth is then evaluated on *exactly* those pairs, read back from the
decoder's own manifest, so truth and decode are comparable pair-for-pair rather
than merely in aggregate.

This is a deliberate break from runs 2-7, which used within-individual pairs
only. Genotype stratification survives the change, restated on haplotype pairs as
carrier-carrier / carrier-noncarrier / noncarrier-noncarrier.

Rates
-----
``MUTATION_RATE`` (1.29e-8) and ``RECOMBINATION_RATE`` (1e-8) apply to the
*simulation* only. Gamma-SMC is given ``theta`` and ``rho`` directly and never
sees these, so the two are decoupled on purpose.

Decoding
--------
``theta = 0.00075`` and ``rho = 0.0006`` are pinned to the Gamma-SMC paper.
``theta = 4 N mu`` implies ``N ~= 14,500`` at this mutation rate, which is far
below EAS's present-day 37,637 but close to its ancient size; since Gamma-SMC
carries a single constant-Ne prior and cannot represent a 10,000-epoch history,
that ancient value is the defensible compromise. The prior biases decoded TMRCAs
young, but the bias is common to the simulated null and the data -- both are
decoded identically -- so it cancels in the calibrated comparison and is not a
threat to validity.

``recent_call = "mean"`` makes ``frac_recent_<T>`` the *proportion of posterior
means below T*, which is the statistic the paper reports. The default "median"
is a different quantity: a Gamma posterior is right-skewed, so its median sits
below its mean and thresholding on the median inflates the statistic at every T.
Runs 3-8 used the median default and are therefore not directly comparable here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .run7_config import (
    ADMIXTURE_PROPORTION,
    ARCHAIC_EFFECTIVE_SIZE,
    ARCHAIC_POPULATION,
    DOMINANCE_COEFFICIENT,
    EAS_POPULATION,
    FOCAL_POSITION_BP,
    FOCAL_SITE_ID,
    GENERATION_TIME_YEARS,
    PROFILE_STRIDE_BP,
    PULSE_GENERATIONS,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    SLIM_BURN_IN,
    SLIM_SCALING_FACTOR,
    SPLIT_GENERATIONS,
    TARGET_FREQUENCY_BAND,
)

__all__ = [
    "DECODER_SETTINGS",
    "MODES",
    "MUTATION_RATE",
    "N_DECODE_PAIRS",
    "N_NEUTRAL",
    "PAIRS_SEED",
    "RHO",
    "SAMPLE_DIPLOIDS",
    "SCENARIOS",
    "SCENARIOS_BY_ID",
    "Scenario",
    "THETA",
    "THRESHOLDS_YEARS",
    "arm_config_record",
    "get_scenario",
    "write_config",
]

MODES = ("selected", "neutral")
SEED_BASE = 20270101
N_NEUTRAL = 100

#: Simulation rates. Gamma-SMC never sees these -- it is given theta/rho.
MUTATION_RATE = 1.29e-8

SAMPLE_DIPLOIDS = 100
N_DECODE_PAIRS = 10_000
PAIRS_SEED = 1729

#: Pinned to the Gamma-SMC paper.
THETA = 0.00075
RHO = 0.0006

DECODER_SETTINGS: dict[str, Any] = {
    "scaled_mutation_rate": THETA,
    "recombination_to_mutation_ratio": RHO / THETA,
    "mutation_rate": MUTATION_RATE,
    "output_at_stride": PROFILE_STRIDE_BP,
    "output_at_hets": False,
    "only_within": False,
    "n_random_pairs": N_DECODE_PAIRS,
    "pairs_seed": PAIRS_SEED,
    "exclude_within": False,
    # The paper's statistic: proportion of posterior *means* below the threshold.
    "recent_call": "mean",
    "cache_size": 1_000,
    "pair_block": 256,
    "exp10": "accurate",
    "backward_alignment": "fixed",
    "threads": 1,
}

#: One grid for every scenario; which threshold wins is a result, not a setting.
THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    5_000.0,
    10_000.0,
    20_000.0,
    30_000.0,
    40_000.0,
    50_000.0,
)

#: Significance: at most one of 100 neutral replicates may exceed the observation.
ALPHA = 0.01


@dataclass(frozen=True)
class Scenario:
    scenario_id: str
    label: str
    origin: str  # "de_novo" or "introgressed_pulse"
    onset_generations: float
    selection_coefficient: float
    expected_cutoff_years: float
    n_replicates: int
    rationale: str

    @property
    def onset_years(self) -> float:
        return self.onset_generations * GENERATION_TIME_YEARS

    @property
    def places_allele_in_archaic(self) -> bool:
        return self.origin == "introgressed_pulse"

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


SCENARIOS: tuple[Scenario, ...] = (
    Scenario(
        scenario_id="recent_denovo",
        label="Recent de novo sweep (LCT-like), 5 kya, s = 0.1",
        origin="de_novo",
        onset_generations=200.0,
        selection_coefficient=0.1,
        expected_cutoff_years=5_000.0,
        n_replicates=100,
        rationale=(
            "A single new copy 200 generations ago, still compressing "
            "coalescence at the moment of sampling, so the signal should sit "
            "below the onset. s = 0.1 rather than an LCT-textbook 0.05 because "
            "run7 measured 0.05 over 200 generations: mean final frequency "
            "0.016. Published LCT estimates span roughly 0.01-0.19."
        ),
    ),
    Scenario(
        scenario_id="postintro_denovo",
        label="De novo sweep at the pulse time, 55 kya, s = 0.002",
        origin="de_novo",
        onset_generations=PULSE_GENERATIONS,
        selection_coefficient=0.002,
        expected_cutoff_years=30_000.0,
        n_replicates=100,
        rationale=(
            "A single new copy on an ordinary modern background, starting when "
            "the pulse lands and under the same coefficient as the introgressed "
            "arm. Holding onset and s fixed leaves origin as the only "
            "difference, so the contrast measures what an archaic starting pool "
            "does and nothing else."
        ),
    ),
    Scenario(
        scenario_id="introgressed_pulse",
        label="Selected introgressed allele, 55 kya pulse, s = 0.002",
        origin="introgressed_pulse",
        onset_generations=PULSE_GENERATIONS,
        selection_coefficient=0.002,
        expected_cutoff_years=50_000.0,
        n_replicates=100,
        rationale=(
            "The allele arrives on archaic haplotypes at 2.5%. Carriers can "
            "only coalesce after the pulse or back in the archaic branch, so "
            "the signal is floored by the pulse and a more recent cutoff cannot "
            "capture it."
        ),
    ),
    Scenario(
        scenario_id="introgressed_pulse_weak",
        label="Selected introgressed allele, 55 kya pulse, s = 0.001",
        origin="introgressed_pulse",
        onset_generations=PULSE_GENERATIONS,
        selection_coefficient=0.001,
        expected_cutoff_years=50_000.0,
        n_replicates=100,
        rationale=(
            "The 'or even lower?' bracket. Expected to end near 7%, which "
            "leaves few carrier-carrier pairs, so it tests where the approach "
            "runs out rather than where it works."
        ),
    ),
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
        "demography": {
            "source": "run7 graft, imported unchanged",
            "split_generations": SPLIT_GENERATIONS,
            "split_years": SPLIT_GENERATIONS * GENERATION_TIME_YEARS,
            "pulse_generations": PULSE_GENERATIONS,
            "pulse_years": PULSE_GENERATIONS * GENERATION_TIME_YEARS,
            "archaic_effective_size": ARCHAIC_EFFECTIVE_SIZE,
            "admixture_proportion": ADMIXTURE_PROPORTION,
            "eas_population": EAS_POPULATION,
            "archaic_population": ARCHAIC_POPULATION,
        },
        "simulation": {
            "sequence_length_bp": SEQUENCE_LENGTH_BP,
            "focal_position_0based": FOCAL_POSITION_BP,
            "focal_site_id": FOCAL_SITE_ID,
            "mutation_rate_per_bp_per_generation": MUTATION_RATE,
            "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
            "generation_time_years": GENERATION_TIME_YEARS,
            "dominance_coefficient": DOMINANCE_COEFFICIENT,
            "sample_diploids": SAMPLE_DIPLOIDS,
            "sample_haplotypes": 2 * SAMPLE_DIPLOIDS,
            "profile_stride_bp": PROFILE_STRIDE_BP,
            "slim_scaling_factor": SLIM_SCALING_FACTOR,
            "slim_burn_in": SLIM_BURN_IN,
        },
        "conditioning": {
            "selected_target_frequency_band_at_pulse": list(TARGET_FREQUENCY_BAND),
            "selected_condition_on_survival": True,
            "neutral_conditioning": "none; the neutral mode carries no focal allele",
        },
        "pairs": {
            "n_decode_pairs": N_DECODE_PAIRS,
            "total_possible_pairs": (2 * SAMPLE_DIPLOIDS) * (2 * SAMPLE_DIPLOIDS - 1) // 2,
            "pairs_seed": PAIRS_SEED,
            "selection": "decoder random draw; truth evaluated on the same manifest",
        },
        "decoding": {
            "theta": THETA,
            "rho": RHO,
            "implied_constant_ne": THETA / (4.0 * MUTATION_RATE),
            "recent_call": DECODER_SETTINGS["recent_call"],
            "primary_statistic": "frac_recent_<T> (proportion of posterior means below T)",
            "secondary_statistic": "mean_p_lt_<T> (mean posterior probability below T)",
            "note": (
                "theta/rho are pinned to the Gamma-SMC paper and are not derived "
                "from the simulation rates; the decoder never sees mu or r."
            ),
        },
        "analysis": {
            "thresholds_years": list(THRESHOLDS_YEARS),
            "alpha": ALPHA,
            "significance_rule": (
                "at most one of the neutral replicates may exceed the observation"
            ),
            "n_neutral_replicates": N_NEUTRAL,
            "scored_at": "focal allele position, per 10 kb stride",
        },
    }


def arm_config_record() -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run8-config/v2",
        "shared": shared_parameter_record(),
        "scenarios": [
            {
                "scenario_id": s.scenario_id,
                "label": s.label,
                "origin": s.origin,
                "onset_generations": s.onset_generations,
                "onset_years": s.onset_years,
                "selection_coefficient": s.selection_coefficient,
                "expected_cutoff_years": s.expected_cutoff_years,
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
    path = directory / "run8_scenarios.json"
    path.write_text(
        json.dumps(arm_config_record(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
