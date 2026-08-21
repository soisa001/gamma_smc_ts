"""Frozen parameters for run7: selection on an introgressed allele, EAS only.

This is the study the graft was built for. One model, one population, no CHB:
the inferred PHLASH EAS history with a Neanderthal branch attached, a focal
allele fixed in that branch, a pulse that delivers it into EAS at ~2.5%, and
selection from the pulse to the present.

The archaic parameters are rounded to their literature values rather than
inherited from another catalogue model:

===================  =============  ==========================================
quantity             value          source
===================  =============  ==========================================
human-Nea split      700 kya        Prufer et al. 2014 (550-765 ky)
introgression pulse  55 kya         Sankararaman 2012; Fu 2014 (47-65 ky)
archaic Ne           3,600          Ragsdale & Gravel 2019
admixture            2.5%           Prufer 2014; Vernot & Akey 2015
===================  =============  ==========================================

At 25 years per generation the split is 28,000 generations and the pulse 2,200.

The ascertainment is explicit and is the point of the design: the focal allele
must actually be segregating at ~2.5% in EAS when selection starts. Fixing it in
the archaic branch and pulsing at 2.5% delivers that in expectation, but drift in
the pulse generation makes the realised frequency vary, so it is *required* by
conditioning rather than assumed. The band is wide enough to be cheap -- with
thousands of genomes the binomial spread around 2.5% is a couple of tenths of a
percent -- and narrow enough that every replicate starts from the same place.

The neutral mode carries no focal allele at all, so the null has nothing
ascertained about it, which is what run6 established the comparison needs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SEQUENCE_LENGTH_BP = 10_000_000
FOCAL_POSITION_BP = 5_000_000
FOCAL_SITE_ID = "run7_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8

GENERATION_TIME_YEARS = 25.0

#: 700 kya at 25 years per generation.
SPLIT_GENERATIONS = 28_000.0
#: 55 kya at 25 years per generation.
PULSE_GENERATIONS = 2_200.0
ARCHAIC_EFFECTIVE_SIZE = 3_600.0
ADMIXTURE_PROPORTION = 0.025

#: A grid rather than one value: at s = 0.01 the allele fixes in essentially
#: every replicate that establishes, which saturates the statistic and hides how
#: power depends on final frequency. Additive selection moves the odds by
#: exp(s t / 2), so over the 2,200 generations since the pulse these land near
#: 0.19, 0.41, 0.86 and 1.00 -- a spread that straddles the f > 2/3 crossover
#: run5 identified.
SELECTION_COEFFICIENTS: tuple[float, ...] = (0.002, 0.003, 0.005, 0.01)
DOMINANCE_COEFFICIENT = 0.5

#: The focal allele must be at ~2.5% in EAS when selection begins.
TARGET_FREQUENCY_BAND = (0.020, 0.030)
#: Replicates are also conditioned on the allele surviving to the present.
#: At the pulse the EAS effective size is only 3,869, so 2.5% is about 39 copies
#: and with h = 0.5 roughly exp(-2 k h s) of them are lost to drift despite
#: selection -- 14% at s = 0.01 and 68% at s = 0.002. Conditioning on survival
#: matches the real analysis, which only ever scans tracts that exist, and the
#: establishment probability is reported separately rather than being buried.
CONDITION_ON_SURVIVAL = True

EAS_POPULATION = "EAS"
ARCHAIC_POPULATION = "Neanderthal"
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
PHLASH_EAS_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"

SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

N_REPLICATES = 80
SAMPLE_DIPLOIDS = 100
PROFILE_STRIDE_BP = 10_000

MODES = ("selected", "neutral")
SEED_BASE = 20261101

THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    4_500.0,
    10_000.0,
    20_000.0,
    30_000.0,
    40_000.0,
    50_000.0,
    100_000.0,
)


@dataclass(frozen=True)
class Run7Arm:
    arm_id: str
    label: str
    description: str

    def seed(self, mode: str, replicate_index: int, s_index: int = 0) -> int:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        return (
            SEED_BASE
            + MODES.index(mode) * 1_000_000
            + s_index * 100_000
            + replicate_index
        )

    def unit_id(self, mode: str, replicate_index: int) -> str:
        return f"{self.arm_id}__{mode}__rep{replicate_index:03d}"


def arm_id_for(selection_coefficient: float) -> str:
    return f"eas_introgressed_s{selection_coefficient:g}".replace(".", "p")


EAS_INTROGRESSED = Run7Arm(
    arm_id="eas_introgressed",
    label="EAS + grafted archaic branch, selection on the introgressed allele",
    description=(
        "The PHLASH EAS median history with a Neanderthal branch splitting at "
        "700 kya and a 2.5% pulse at 55 kya. The focal allele is fixed in the "
        "archaic branch immediately after the split, so what arrives in EAS is "
        "carried on genuine archaic haplotypes rather than on a random set of "
        "modern ones. Selection acts in EAS from the pulse to the present over a "
        "grid of coefficients, the allele is required to be at 2.0-3.0% in "
        "EAS when selection starts, and it must survive to the present. The "
        "neutral mode carries no focal allele at all and is shared across "
        "the grid."
    ),
)

ARMS: tuple[Run7Arm, ...] = (EAS_INTROGRESSED,)


def shared_parameter_record() -> dict[str, Any]:
    return {
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_0based": FOCAL_POSITION_BP,
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        "generation_time_years": GENERATION_TIME_YEARS,
        "split_generations": SPLIT_GENERATIONS,
        "split_years": SPLIT_GENERATIONS * GENERATION_TIME_YEARS,
        "pulse_generations": PULSE_GENERATIONS,
        "pulse_years": PULSE_GENERATIONS * GENERATION_TIME_YEARS,
        "archaic_effective_size": ARCHAIC_EFFECTIVE_SIZE,
        "admixture_proportion": ADMIXTURE_PROPORTION,
        "selection_coefficients": list(SELECTION_COEFFICIENTS),
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "target_frequency_band": list(TARGET_FREQUENCY_BAND),
        "condition_on_survival": CONDITION_ON_SURVIVAL,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "n_replicates_per_mode": N_REPLICATES,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "profile_stride_bp": PROFILE_STRIDE_BP,
        "neutral_contains_focal_allele": False,
        "raw_tmrca_stored": True,
    }


def arm_config_record(arm: Run7Arm) -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run7-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "description": arm.description,
        "shared": shared_parameter_record(),
        "modes": {
            mode: {
                "n_replicates": N_REPLICATES,
                "seeds": [arm.seed(mode, i) for i in range(N_REPLICATES)],
            }
            for mode in MODES
        },
    }
    json.dumps(record)
    return record


def write_arm_configs(destination: str | Path) -> dict[str, Path]:
    directory = Path(destination)
    directory.mkdir(parents=True, exist_ok=True)
    written = {}
    for arm in ARMS:
        path = directory / f"run7_{arm.arm_id}.json"
        path.write_text(
            json.dumps(arm_config_record(arm), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        written[arm.arm_id] = path
    return written


__all__ = [
    "ADMIXTURE_PROPORTION",
    "ARCHAIC_EFFECTIVE_SIZE",
    "ARCHAIC_POPULATION",
    "ARMS",
    "DOMINANCE_COEFFICIENT",
    "EAS_INTROGRESSED",
    "EAS_POPULATION",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "GENERATION_TIME_YEARS",
    "MODES",
    "MUTATION_RATE",
    "N_REPLICATES",
    "PHLASH_EAS_RELATIVE_PATH",
    "PHLASH_EAS_SHA256",
    "PROFILE_STRIDE_BP",
    "PULSE_GENERATIONS",
    "RECOMBINATION_RATE",
    "REQUIRED_SLIM_VERSION",
    "REQUIRED_STDPOPSIM_VERSION",
    "Run7Arm",
    "SAMPLE_DIPLOIDS",
    "SEED_BASE",
    "SELECTION_COEFFICIENTS",
    "CONDITION_ON_SURVIVAL",
    "arm_id_for",
    "SEQUENCE_LENGTH_BP",
    "SLIM_BURN_IN",
    "SLIM_SCALING_FACTOR",
    "SPLIT_GENERATIONS",
    "TARGET_FREQUENCY_BAND",
    "THRESHOLDS_YEARS",
    "arm_config_record",
    "shared_parameter_record",
    "write_arm_configs",
]
