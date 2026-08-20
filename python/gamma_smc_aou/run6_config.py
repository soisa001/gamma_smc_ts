"""Frozen parameters for run6: the genotype-stratified spatial scan.

The CHB side of run6 needs no new simulation -- run5's
``chb_from_introgression`` arm already stored the raw per-pair TMRCA matrices and
the focal genotypes, which is exactly what the scan consumes.  Only the EAS arm
is new here: run5's EAS reference used msprime's structured-coalescent sweep,
which never marks which lineages carry the beneficial allele, so it cannot be
stratified by genotype.

This arm is therefore SLiM, which has two further advantages: it marks the
allele, and it carries the **full PHLASH history** with no epoch restriction, so
the flat-window compromise run5 needed for msprime disappears entirely.

Parameters mirror the CHB arm so the two are comparable: the allele is placed as
standing variation at 2%, selection is ``s = 0.01`` acting to the present, and
the onset is matched to CHB's introgression onset in years (35,960 y, which at
25 years per generation is 1,440 generations).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SEQUENCE_LENGTH_BP = 10_000_000
FOCAL_POSITION_BP = 5_000_000
FOCAL_SITE_ID = "run6_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8

SELECTION_COEFFICIENT = 0.01
DOMINANCE_COEFFICIENT = 0.5

#: Matched to the CHB arm's introgressed frequency at the start of selection.
EAS_PLACEMENT_FREQUENCY = 0.02
#: CHB's onset is 1,240 generations at 29 years; matched in years, not
#: generations, because the two models declare different generation times.
EAS_ONSET_GENERATIONS = 1_440.0
EAS_GENERATION_TIME = 25.0
EAS_POPULATION = "EAS"
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
PHLASH_EAS_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"

SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

N_REPLICATES = 100
SAMPLE_DIPLOIDS = 100
PROFILE_STRIDE_BP = 10_000

MODES = ("selected", "neutral")
SEED_BASE = 20260823

#: The CHB side of the scan reads run5's already-simulated arm.
CHB_SOURCE = {
    "study": "sim_results_run5",
    "arm_id": "chb_from_introgression",
    "generation_time": 29.0,
    "note": "raw TMRCA matrices and focal genotypes were retained by run5",
}


@dataclass(frozen=True)
class Run6Arm:
    arm_id: str
    label: str
    target_population: str
    generation_time: float
    description: str

    def seed(self, mode: str, replicate_index: int) -> int:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        return SEED_BASE + MODES.index(mode) * 1_000_000 + replicate_index

    def unit_id(self, mode: str, replicate_index: int) -> str:
        return f"{self.arm_id}__{mode}__rep{replicate_index:03d}"

    def onset_generations(self) -> float:
        return (
            round(EAS_ONSET_GENERATIONS / SLIM_SCALING_FACTOR) * SLIM_SCALING_FACTOR
        )


EAS_STANDING = Run6Arm(
    arm_id="eas_standing",
    label="EAS standing variation, s = 0.01 (SLiM, full PHLASH history)",
    target_population=EAS_POPULATION,
    generation_time=EAS_GENERATION_TIME,
    description=(
        "The PHLASH EAS median history simulated forwards in SLiM, so the allele "
        "is marked and the genotype stratification is available, and so the full "
        "ten-thousand-epoch history is retained rather than flattened across a "
        "sweep window as msprime required. The allele is placed on a uniformly "
        "random 2% of genomes 1,440 generations ago -- matched in years to the "
        "CHB introgression onset -- and selected at s = 0.01 to the present. The "
        "neutral mode carries no allele at all."
    ),
)

ARMS: tuple[Run6Arm, ...] = (EAS_STANDING,)
ARMS_BY_ID = {arm.arm_id: arm for arm in ARMS}


def get_arm(arm_id: str) -> Run6Arm:
    try:
        return ARMS_BY_ID[arm_id]
    except KeyError:
        raise ValueError(
            f"unknown arm {arm_id!r}; expected one of {sorted(ARMS_BY_ID)}"
        ) from None


def shared_parameter_record() -> dict[str, Any]:
    return {
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_0based": FOCAL_POSITION_BP,
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        "selection_coefficient": SELECTION_COEFFICIENT,
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "eas_placement_frequency": EAS_PLACEMENT_FREQUENCY,
        "eas_onset_generations": EAS_ONSET_GENERATIONS,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "n_replicates_per_mode": N_REPLICATES,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "profile_stride_bp": PROFILE_STRIDE_BP,
        "neutral_contains_focal_allele": False,
        "raw_tmrca_stored": True,
        "chb_source": CHB_SOURCE,
    }


def arm_config_record(arm: Run6Arm) -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run6-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "target_population": arm.target_population,
        "generation_time_years": arm.generation_time,
        "selection_onset_generations": arm.onset_generations(),
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
        path = directory / f"run6_{arm.arm_id}.json"
        path.write_text(
            json.dumps(arm_config_record(arm), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        written[arm.arm_id] = path
    return written


__all__ = [
    "ARMS",
    "ARMS_BY_ID",
    "CHB_SOURCE",
    "DOMINANCE_COEFFICIENT",
    "EAS_GENERATION_TIME",
    "EAS_ONSET_GENERATIONS",
    "EAS_PLACEMENT_FREQUENCY",
    "EAS_POPULATION",
    "EAS_STANDING",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "MODES",
    "MUTATION_RATE",
    "N_REPLICATES",
    "PHLASH_EAS_RELATIVE_PATH",
    "PHLASH_EAS_SHA256",
    "PROFILE_STRIDE_BP",
    "RECOMBINATION_RATE",
    "REQUIRED_SLIM_VERSION",
    "REQUIRED_STDPOPSIM_VERSION",
    "Run6Arm",
    "SAMPLE_DIPLOIDS",
    "SEED_BASE",
    "SELECTION_COEFFICIENT",
    "SEQUENCE_LENGTH_BP",
    "SLIM_BURN_IN",
    "SLIM_SCALING_FACTOR",
    "arm_config_record",
    "get_arm",
    "shared_parameter_record",
    "write_arm_configs",
]
