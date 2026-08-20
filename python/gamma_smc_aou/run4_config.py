"""Frozen parameters for the run4 study: EPAS1-calibrated sweeps.

run2 and run3 both simulated sweeps that were far older, weaker and softer than
the best-characterised real adaptive-introgression locus.  run4 calibrates to
EPAS1 (see ``sim_results_run4/../sim_results_run3/EPAS1_CALIBRATION.md``):

* selection onset **6,000 years ago** rather than 36,000-57,000
* selection coefficient **s = 0.05** rather than 0.01
* a **single founding haplotype** rather than 30-47 independent ones

Those three changes matter for the same reason: a hard sweep's carriers must
coalesce by the sweep onset, so the onset sets a ceiling on carrier TMRCA.  At
6,000 years that ceiling is ~207 generations, which lands inside the original
1-50 ky threshold grid -- the grid run2 found empty because its sweep was 36-57
ky old.  The grid is therefore restored to its original seven cutoffs.

Three arms, each with a plain neutral null that contains **no focal allele at
all**:

1. ``chb_archaic_single_founder`` -- the EPAS1 model. The allele is fixed in
   Neanderthal after the archaic split so it marks archaic ancestry, introgresses
   into CHB, and is then reduced to a single archaic founding haplotype at the
   selection onset. Replicates with no archaic ancestry at the focal base are
   rejected and re-drawn, which is exactly the ascertainment that makes EPAS1 a
   locus worth studying in the first place.
2. ``chb_denovo`` -- positive control. Same demography, same onset, same s, but
   the allele starts as a single de novo copy on a random CHB genome with no
   archaic background. This is the upper bound on detectability.
3. ``eas_denovo`` -- the same de novo hard sweep in the PHLASH EAS history,
   matched to CHB in years.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Shared parameters
# ---------------------------------------------------------------------------

SEQUENCE_LENGTH_BP = 10_000_000
FOCAL_POSITION_BP = 5_000_000
FOCAL_SITE_ID = "run4_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8

#: EPAS1-calibrated. 2*Ne*s ~ 118 with Ne ~ 1,100 implies s ~ 0.053.
SELECTION_COEFFICIENT = 0.05
DOMINANCE_COEFFICIENT = 0.5

#: EPAS1 selection onset, ~9,000 years ago (2,800-10,000 independently).
#: 6,000 sits in the middle of the overlap of the two published intervals.
SELECTION_ONSET_YEARS = 6_000.0

SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

N_REPLICATES = 100
SAMPLE_DIPLOIDS = 100

#: Restored to run2's original grid: an EPAS1-calibrated onset puts the signal
#: back inside it.
TMRCA_THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    4_500.0,
    10_000.0,
    20_000.0,
    30_000.0,
    40_000.0,
    50_000.0,
)
PROFILE_STRIDE_BP = 10_000
SIGNIFICANCE_LEVEL = 0.05

MODES = ("selected", "neutral")
SOURCES = ("tree_truth", "gamma_smc")
SEED_BASE = 20260821

#: Ascertainment: the archaic arm needs archaic ancestry at the focal base.
#: run3 measured that as present in ~37.5% of replicates, so a handful of
#: re-draws per accepted replicate is expected.
MAX_ASCERTAINMENT_ATTEMPTS = 40

# ---------------------------------------------------------------------------
# Catalog contracts
# ---------------------------------------------------------------------------

OOA_MODEL_ID = "OutOfAfricaArchaicAdmixture_5R19"
OOA_GENERATION_TIME = 29.0
NEANDERTHAL_SPLIT_GENERATIONS = 19_275.86217652958
CHB_SPLIT_GENERATIONS = 1_241.379310344827
CHB_ARCHAIC_MIGRATION_END_GENERATIONS = 644.8275862068965
ARCHAIC_POPULATION = "Neanderthal"
CHB_POPULATION = "CHB"

EAS_GENERATION_TIME = 25.0
EAS_POPULATION = "EAS"
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
PHLASH_EAS_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Run4Arm:
    arm_id: str
    label: str
    demography_id: str
    target_population: str
    generation_time: float
    #: "archaic_single_founder" | "denovo"
    origin: str
    #: Population the allele is first introduced into.
    origin_population: str
    #: True when the replicate must be re-drawn if ascertainment fails.
    ascertained: bool
    description: str

    @property
    def archaic(self) -> bool:
        return self.origin == "archaic_single_founder"

    def seed(self, mode: str, replicate_index: int) -> int:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        offset = ARM_ORDER.index(self.arm_id) + 1
        return (
            SEED_BASE
            + offset * 10_000_000
            + MODES.index(mode) * 1_000_000
            + replicate_index
        )

    def unit_id(self, mode: str, replicate_index: int) -> str:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        return f"{self.arm_id}__{mode}__rep{replicate_index:03d}"

    def onset_generations(self) -> float:
        """Selection onset in generations, snapped to the SLiM tick grid."""
        raw = SELECTION_ONSET_YEARS / self.generation_time
        return round(raw / SLIM_SCALING_FACTOR) * SLIM_SCALING_FACTOR


ARM_ORDER = ("chb_archaic_single_founder", "chb_denovo", "eas_denovo")

CHB_ARCHAIC_ARM = Run4Arm(
    arm_id="chb_archaic_single_founder",
    label="CHB archaic, single founder (EPAS1-like)",
    demography_id=OOA_MODEL_ID,
    target_population=CHB_POPULATION,
    generation_time=OOA_GENERATION_TIME,
    origin="archaic_single_founder",
    origin_population=ARCHAIC_POPULATION,
    ascertained=True,
    description=(
        "The allele is fixed in Neanderthal one tick after the archaic split, so "
        "it marks archaic ancestry exactly, and introgresses into CHB by the "
        "catalog's continuous migration. At the selection onset it is reduced to "
        "a single archaic founding haplotype -- every present-day carrier "
        "therefore descends from one introgressed genome, as EPAS1 does. "
        "Replicates with no archaic ancestry at the focal base are re-drawn; "
        "that is the same ascertainment that makes a locus like EPAS1 worth "
        "studying."
    ),
)

CHB_DENOVO_ARM = Run4Arm(
    arm_id="chb_denovo",
    label="CHB de novo hard sweep (positive control)",
    demography_id=OOA_MODEL_ID,
    target_population=CHB_POPULATION,
    generation_time=OOA_GENERATION_TIME,
    origin="denovo",
    origin_population=CHB_POPULATION,
    ascertained=False,
    description=(
        "Same demography, onset and selection coefficient, but the allele begins "
        "as a single de novo copy on a random CHB genome with no archaic "
        "background. This isolates what the statistic can see from a hard sweep "
        "alone, and is the upper bound on detectability."
    ),
)

EAS_DENOVO_ARM = Run4Arm(
    arm_id="eas_denovo",
    label="EAS de novo hard sweep",
    demography_id="PhlashEASMedian",
    target_population=EAS_POPULATION,
    generation_time=EAS_GENERATION_TIME,
    origin="denovo",
    origin_population=EAS_POPULATION,
    ascertained=False,
    description=(
        "The same de novo hard sweep in the PHLASH EAS history, with the onset "
        "matched to the CHB arms in years rather than generations because the "
        "two models declare different generation times (25 versus 29). Replaces "
        "run2/run3's standing-variation EAS arm, which was underpowered by "
        "construction: placing the allele on many independent resident genomes "
        "leaves carriers sharing no recent ancestry, so amplifying them "
        "compresses no coalescent times."
    ),
)

ARMS: tuple[Run4Arm, ...] = (CHB_ARCHAIC_ARM, CHB_DENOVO_ARM, EAS_DENOVO_ARM)
ARMS_BY_ID = {arm.arm_id: arm for arm in ARMS}


def get_arm(arm_id: str) -> Run4Arm:
    try:
        return ARMS_BY_ID[arm_id]
    except KeyError:
        raise ValueError(
            f"unknown arm {arm_id!r}; expected one of {sorted(ARMS_BY_ID)}"
        ) from None


# ---------------------------------------------------------------------------
# Gamma-SMC
# ---------------------------------------------------------------------------

DECODER_SETTINGS: dict[str, Any] = {
    "recombination_to_mutation_ratio": RECOMBINATION_RATE / MUTATION_RATE,
    "mutation_rate": MUTATION_RATE,
    "output_at_stride": PROFILE_STRIDE_BP,
    "output_at_hets": False,
    "only_within": True,
    "recent_call": "median",
    "recent_call_probability": 0.5,
    "cache_size": 1_000,
    "pair_block": 256,
    "exp10": "accurate",
    "backward_alignment": "fixed",
    "threads": 1,
}

PRESENT_NE = {
    "chb_archaic_single_founder": 65_834.77001122756,
    "chb_denovo": 65_834.77001122756,
    "eas_denovo": 37_636.6,
}


def scaled_mutation_rate(arm: Run4Arm) -> float:
    return 4.0 * PRESENT_NE[arm.arm_id] * MUTATION_RATE


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------


def shared_parameter_record() -> dict[str, Any]:
    return {
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_0based": FOCAL_POSITION_BP,
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        "selection_coefficient": SELECTION_COEFFICIENT,
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "selection_onset_years": SELECTION_ONSET_YEARS,
        "calibration_target": "EPAS1",
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "slim_burn_in": SLIM_BURN_IN,
        "n_replicates_per_mode": N_REPLICATES,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "tmrca_thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
        "profile_stride_bp": PROFILE_STRIDE_BP,
        "significance_level": SIGNIFICANCE_LEVEL,
        "seed_base": SEED_BASE,
        "sources": list(SOURCES),
        "neutral_contains_focal_allele": False,
        "present_day_af_conditioning": None,
        "allele_survival_conditioning": None,
        "raw_tmrca_stored": True,
        "statistic": (
            "overall within-individual P(TMRCA < x) at the focal position, "
            "over all sampled diploids"
        ),
        "p_value_rule": "(1 + #{null >= observed}) / (1 + n_null), one-sided upper tail",
        "p_value_floor": 1.0 / (N_REPLICATES + 1),
    }


def arm_config_record(arm: Run4Arm) -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run4-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "demography_id": arm.demography_id,
        "target_population": arm.target_population,
        "origin": arm.origin,
        "origin_population": arm.origin_population,
        "ascertained": arm.ascertained,
        "generation_time_years": arm.generation_time,
        "selection_onset_generations": arm.onset_generations(),
        "selection_onset_years_realized": arm.onset_generations() * arm.generation_time,
        "description": arm.description,
        "scaled_mutation_rate": scaled_mutation_rate(arm),
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
        path = directory / f"run4_{arm.arm_id}.json"
        path.write_text(
            json.dumps(arm_config_record(arm), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        written[arm.arm_id] = path
    return written


__all__ = [
    "ARCHAIC_POPULATION",
    "ARMS",
    "ARMS_BY_ID",
    "ARM_ORDER",
    "CHB_ARCHAIC_ARM",
    "CHB_ARCHAIC_MIGRATION_END_GENERATIONS",
    "CHB_DENOVO_ARM",
    "CHB_POPULATION",
    "CHB_SPLIT_GENERATIONS",
    "DECODER_SETTINGS",
    "DOMINANCE_COEFFICIENT",
    "EAS_DENOVO_ARM",
    "EAS_GENERATION_TIME",
    "EAS_POPULATION",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "MAX_ASCERTAINMENT_ATTEMPTS",
    "MODES",
    "MUTATION_RATE",
    "NEANDERTHAL_SPLIT_GENERATIONS",
    "N_REPLICATES",
    "OOA_GENERATION_TIME",
    "OOA_MODEL_ID",
    "PHLASH_EAS_RELATIVE_PATH",
    "PHLASH_EAS_SHA256",
    "PRESENT_NE",
    "PROFILE_STRIDE_BP",
    "RECOMBINATION_RATE",
    "REQUIRED_SLIM_VERSION",
    "REQUIRED_STDPOPSIM_VERSION",
    "Run4Arm",
    "SAMPLE_DIPLOIDS",
    "SEED_BASE",
    "SELECTION_COEFFICIENT",
    "SELECTION_ONSET_YEARS",
    "SEQUENCE_LENGTH_BP",
    "SIGNIFICANCE_LEVEL",
    "SLIM_BURN_IN",
    "SLIM_SCALING_FACTOR",
    "SOURCES",
    "TMRCA_THRESHOLDS_YEARS",
    "arm_config_record",
    "get_arm",
    "scaled_mutation_rate",
    "shared_parameter_record",
    "write_arm_configs",
]
