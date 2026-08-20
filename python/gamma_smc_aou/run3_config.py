"""Frozen parameters for the run3 selection-detection study.

run3 changes four things relative to run2 (``sim_results_run2/FRAMEWORK.md``):

1. **The CHB arm uses** ``OutOfAfricaArchaicAdmixture_5R19`` **instead of**
   ``AncientEurasia_9K19``.  Archaic ancestry arrives by *continuous migration*
   into CHB over roughly 600 generations rather than a single pulse into an
   ancestral population, which removes run2's "the pulse lands in the Han
   ancestor" problem entirely: CHB is founded at the moment its Neanderthal
   migration begins, so selection can start exactly there.
2. **The focal allele is fixed in the Neanderthal population immediately after
   the archaic split, in BOTH modes.**  It is therefore a marker of archaic
   ancestry, present with or without selection.  The neutral arm now contains
   introgression without selection, so the test isolates selection *given*
   introgression rather than confounding the two.
3. **Raw per-pair TMRCA values are stored**, so any threshold, window or
   summary can be recomputed post hoc without re-simulating.  run2's biggest
   practical failing was that changing a cutoff cost 23 core-hours.
4. **Both tree truth and Gamma-SMC decoding are run**, on the same replicates,
   so the decoder's cost in power can be measured rather than assumed.

The threshold grid is also extended well past the admixture floor, which run2
showed was the reason the original grid saw almost nothing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Shared simulation parameters
# ---------------------------------------------------------------------------

SEQUENCE_LENGTH_BP = 10_000_000
FOCAL_POSITION_BP = 5_000_000
FOCAL_SITE_ID = "run3_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8

SELECTION_COEFFICIENT = 0.01
DOMINANCE_COEFFICIENT = 0.5

#: Frequency the EAS arm places the allele at, matched to the CHB arm.
#:
#: The CHB arm has no analogous knob -- its frequency accumulates from
#: continuous archaic migration -- so this was calibrated from 24 neutral CHB
#: replicates (seeds 9001-9024): mean present-day archaic frequency at the focal
#: base 0.00705, but with **62.5% of replicates carrying no archaic ancestry
#: there at all**, giving a mean of 0.0188 conditional on presence. EAS always
#: places the allele, so it is matched to the conditional value.
EAS_PLACEMENT_FREQUENCY = 0.019
EAS_PLACEMENT_CALIBRATION = {
    "source": "24 neutral CHB replicates, seeds 9001-9024",
    "chb_mean_archaic_frequency_all": 0.00705,
    "chb_mean_archaic_frequency_given_present": 0.0188,
    "chb_fraction_with_no_archaic_ancestry_at_focal_base": 0.625,
    "chb_max_observed": 0.05053,
    "matched_to": "conditional_on_presence",
}

SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

N_REPLICATES = 100
SAMPLE_DIPLOIDS = 100

#: Extended past the admixture floor. run2 stopped at 50 ky and saw almost
#: nothing; raw TMRCA is stored, so this grid is only the default report grid.
TMRCA_THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    4_500.0,
    10_000.0,
    20_000.0,
    30_000.0,
    40_000.0,
    50_000.0,
    75_000.0,
    100_000.0,
    150_000.0,
    200_000.0,
)
PROFILE_STRIDE_BP = 10_000
SIGNIFICANCE_LEVEL = 0.05

MODES = ("selected", "neutral")
SOURCES = ("tree_truth", "gamma_smc")

SEED_BASE = 20260820

# ---------------------------------------------------------------------------
# Catalog contract for OutOfAfricaArchaicAdmixture_5R19
# ---------------------------------------------------------------------------

OOA_MODEL_ID = "OutOfAfricaArchaicAdmixture_5R19"
OOA_GENERATION_TIME = 29.0
#: Forwards in time, Neanderthal branches off the modern human lineage here.
NEANDERTHAL_SPLIT_GENERATIONS = 19_275.86217652958
#: CHB is founded here, and its Neanderthal migration switches on at the same
#: time, so this is both the introgression onset and the selection onset.
CHB_SPLIT_GENERATIONS = 1_241.379310344827
#: Neanderthal migration into CHB switches off here.
CHB_ARCHAIC_MIGRATION_END_GENERATIONS = 644.8275862068965
ARCHAIC_POPULATION = "Neanderthal"
CHB_POPULATION = "CHB"

# ---------------------------------------------------------------------------
# EAS PHLASH artifact
# ---------------------------------------------------------------------------

EAS_GENERATION_TIME = 25.0
EAS_POPULATION = "EAS"
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
PHLASH_EAS_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Run3Arm:
    """One demographic arm of the study."""

    arm_id: str
    label: str
    demography_id: str
    target_population: str
    generation_time: float
    #: True when archaic ancestry is simulated explicitly.
    archaic: bool
    #: Population the focal allele is introduced into.
    origin_population: str
    #: How the SLiM ``add_mut`` replacement chooses carrier genomes.
    placement: str
    description: str

    def seed(self, mode: str, replicate_index: int) -> int:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        arm_offset = 1 if self.archaic else 2
        return (
            SEED_BASE
            + arm_offset * 10_000_000
            + MODES.index(mode) * 1_000_000
            + replicate_index
        )

    def unit_id(self, mode: str, replicate_index: int) -> str:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        return f"{self.arm_id}__{mode}__rep{replicate_index:03d}"

    def years_to_generations(self, years: float) -> float:
        return float(years) / self.generation_time


CHB_ARM = Run3Arm(
    arm_id="chb_ooa_archaic",
    label="CHB (OutOfAfricaArchaicAdmixture_5R19)",
    demography_id=OOA_MODEL_ID,
    target_population=CHB_POPULATION,
    generation_time=OOA_GENERATION_TIME,
    archaic=True,
    origin_population=ARCHAIC_POPULATION,
    placement="fixed_in_archaic_population",
    description=(
        "The focal allele is fixed in the Neanderthal population immediately "
        "after the archaic split, in both modes, so it is an exact marker of "
        "archaic ancestry. It reaches CHB by the catalog's own continuous "
        "Neanderthal migration, which runs from the CHB founding until roughly "
        "645 generations ago. Selection acts in CHB only, from the moment that "
        "migration begins to the present. Because CHB is founded at that same "
        "moment, there is no ancestral-population confound."
    ),
)

EAS_ARM = Run3Arm(
    arm_id="eas_phlash",
    label="EAS (PHLASH median, custom demography)",
    demography_id="PhlashEASMedian",
    target_population=EAS_POPULATION,
    generation_time=EAS_GENERATION_TIME,
    archaic=False,
    origin_population=EAS_POPULATION,
    placement="uniform_random_genomes",
    description=(
        "Unchanged from run2: single-population PHLASH EAS history, no archaic "
        "source. The allele is placed on a uniformly random 2.96% of genomes, "
        "in both modes, at the selection onset. The onset is matched to the CHB "
        "arm in years, not generations, because the two models declare "
        "different generation times (25 vs 29). This arm is NOT an "
        "introgression model and its carriers accumulate no archaic background."
    ),
)

ARMS: tuple[Run3Arm, ...] = (CHB_ARM, EAS_ARM)
ARMS_BY_ID = {arm.arm_id: arm for arm in ARMS}


def get_arm(arm_id: str) -> Run3Arm:
    try:
        return ARMS_BY_ID[arm_id]
    except KeyError:
        raise ValueError(
            f"unknown arm {arm_id!r}; expected one of {sorted(ARMS_BY_ID)}"
        ) from None


# ---------------------------------------------------------------------------
# Gamma-SMC decoding
# ---------------------------------------------------------------------------

#: theta = 4 * Ne * mu. Chosen per arm from the present-day target Ne so the
#: decoder's coalescent prior is not grossly misspecified.
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

#: Present-day effective sizes used to set the decoder's scaled mutation rate.
PRESENT_NE = {
    "chb_ooa_archaic": 65_834.77001122756,
    "eas_phlash": 37_636.6,
}


def scaled_mutation_rate(arm: Run3Arm) -> float:
    """Return theta = 4 * Ne_present * mu for one arm."""
    return 4.0 * PRESENT_NE[arm.arm_id] * MUTATION_RATE


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------


def shared_parameter_record() -> dict[str, Any]:
    return {
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_0based": FOCAL_POSITION_BP,
        "focal_site_id": FOCAL_SITE_ID,
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        "selection_coefficient": SELECTION_COEFFICIENT,
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "slim_scaling_factor": SLIM_SCALING_FACTOR,
        "slim_burn_in": SLIM_BURN_IN,
        "required_slim_version": REQUIRED_SLIM_VERSION,
        "required_stdpopsim_version": REQUIRED_STDPOPSIM_VERSION,
        "n_replicates_per_mode": N_REPLICATES,
        "sample_diploids": SAMPLE_DIPLOIDS,
        "tmrca_thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
        "profile_stride_bp": PROFILE_STRIDE_BP,
        "significance_level": SIGNIFICANCE_LEVEL,
        "seed_base": SEED_BASE,
        "sources": list(SOURCES),
        "allele_present_in_both_modes": True,
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


def arm_config_record(arm: Run3Arm) -> dict[str, Any]:
    record: dict[str, Any] = {
        "schema": "gamma-smc.run3-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "demography_id": arm.demography_id,
        "target_population": arm.target_population,
        "origin_population": arm.origin_population,
        "generation_time_years": arm.generation_time,
        "archaic": arm.archaic,
        "placement": arm.placement,
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
    if arm.archaic:
        record["source"] = {
            "kind": "stdpopsim_catalog",
            "species": "HomSap",
            "model": arm.demography_id,
            "neanderthal_split_generations": NEANDERTHAL_SPLIT_GENERATIONS,
            "chb_split_generations": CHB_SPLIT_GENERATIONS,
            "archaic_migration_end_generations": CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
        }
    else:
        record["source"] = {
            "kind": "tracked_phlash_artifact",
            "path": PHLASH_EAS_RELATIVE_PATH,
            "sha256": PHLASH_EAS_SHA256,
            "curve": "pointwise_median",
            "placement_frequency": EAS_PLACEMENT_FREQUENCY,
            "placement_calibration": EAS_PLACEMENT_CALIBRATION,
        }
    json.dumps(record)
    return record


def write_arm_configs(destination: str | Path) -> dict[str, Path]:
    directory = Path(destination)
    directory.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for arm in ARMS:
        suffix = "chb" if arm.archaic else "eas"
        path = directory / f"run3_{suffix}.json"
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
    "CHB_ARM",
    "CHB_ARCHAIC_MIGRATION_END_GENERATIONS",
    "CHB_POPULATION",
    "CHB_SPLIT_GENERATIONS",
    "DECODER_SETTINGS",
    "DOMINANCE_COEFFICIENT",
    "EAS_ARM",
    "EAS_GENERATION_TIME",
    "EAS_PLACEMENT_CALIBRATION",
    "EAS_PLACEMENT_FREQUENCY",
    "EAS_POPULATION",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
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
    "Run3Arm",
    "SAMPLE_DIPLOIDS",
    "SEED_BASE",
    "SELECTION_COEFFICIENT",
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
