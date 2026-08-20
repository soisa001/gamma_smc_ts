"""Frozen parameters for the run2 selection-detection study.

run2 asks one question: can selection at an introgressed locus be detected from
the within-individual TMRCA distribution alone?

Everything here is a constant on purpose.  The earlier EAS studies swept nine
target allele frequencies, three selection coefficients and three demographic
quantiles, and defined their null by conditioning on the focal allele's present
frequency.  run2 fixes one selection coefficient, one starting frequency, one
onset time and one statistic, and uses an unconditioned null.  See
``sim_results_run2/FRAMEWORK.md``.
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
FOCAL_SITE_ID = "run2_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8
GENERATION_TIME_YEARS = 25.0

# Selection.  s is the homozygote advantage; h = 0.5 makes it additive.
SELECTION_COEFFICIENT = 0.01
DOMINANCE_COEFFICIENT = 0.5

# Standing variation.  This is the AncientEurasia_9K19 Neanderthal -> Loschbour
# pulse proportion, reused verbatim as the starting frequency in both arms.
STANDING_FREQUENCY = 0.0296
PULSE_GENERATIONS = 2_272.0
HAN_SPLIT_GENERATIONS = 2_016.0

# SLiM execution.
SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

# Study size.
N_REPLICATES = 100
SAMPLE_DIPLOIDS = 100

# Analysis.
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

# Seeding.  Every replicate seed is a pure function of (arm, mode, index), so a
# replicate can be re-run in isolation and reproduce byte-for-byte.
SEED_BASE = 20260819


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Run2Arm:
    """One demographic arm of the study."""

    arm_id: str
    label: str
    demography_id: str
    target_population: str
    #: Population the focal allele is placed into at the standing-variation tick.
    origin_population: str
    #: True when the allele rides in on a real archaic introgression pulse.
    introgression: bool
    #: How the SLiM ``add_mut`` replacement chooses carrier genomes.
    placement: str
    description: str

    def seed(self, mode: str, replicate_index: int) -> int:
        """Return the deterministic seed for one replicate."""
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        if not 0 <= replicate_index < 100_000:
            raise ValueError("replicate_index must be in [0, 100000)")
        arm_offset = 1 if self.introgression else 2
        mode_offset = MODES.index(mode)
        return (
            SEED_BASE
            + arm_offset * 10_000_000
            + mode_offset * 1_000_000
            + replicate_index
        )

    def unit_id(self, mode: str, replicate_index: int) -> str:
        if mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
        return f"{self.arm_id}__{mode}__rep{replicate_index:03d}"


CHB_ARM = Run2Arm(
    arm_id="chb_ancient_eurasia",
    label="CHB (Han, AncientEurasia_9K19)",
    demography_id="AncientEurasia_9K19",
    target_population="Han",
    origin_population="Loschbour",
    introgression=True,
    placement="all_introgressing_genomes",
    description=(
        "The catalog's own 2.96% Neanderthal -> Loschbour pulse. The focal "
        "allele is placed on every introgressing genome, which is exactly "
        "equivalent to the allele being fixed in the archaic source at the "
        "pulse: the post-pulse recipient frequency is the pulse proportion, "
        "and carriers are exactly the archaic haplotypes. No rejection "
        "sampling is involved."
    ),
)

EAS_ARM = Run2Arm(
    arm_id="eas_phlash",
    label="EAS (PHLASH median, custom demography)",
    demography_id="PhlashEASMedian",
    target_population="EAS",
    origin_population="EAS",
    introgression=False,
    placement="uniform_random_genomes",
    description=(
        "Single-population PHLASH EAS history. There is no archaic source and "
        "no pulse: at the same tick as the CHB pulse the focal allele is "
        "placed on a uniformly random 2.96% of genomes. This arm is NOT an "
        "introgression model. It matches the CHB arm on starting frequency, "
        "onset time and selection coefficient, and differs in that carriers "
        "sit on ordinary EAS haplotypes rather than divergent archaic ones."
    ),
)

ARMS: tuple[Run2Arm, ...] = (CHB_ARM, EAS_ARM)
ARMS_BY_ID = {arm.arm_id: arm for arm in ARMS}

#: Tracked PHLASH artifact backing the EAS arm, relative to the repository root.
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
PHLASH_EAS_SHA256 = "a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928"


def get_arm(arm_id: str) -> Run2Arm:
    try:
        return ARMS_BY_ID[arm_id]
    except KeyError:
        raise ValueError(
            f"unknown arm {arm_id!r}; expected one of {sorted(ARMS_BY_ID)}"
        ) from None


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------


def shared_parameter_record() -> dict[str, Any]:
    """Return the JSON-serializable block shared by both arms."""

    return {
        "sequence_length_bp": SEQUENCE_LENGTH_BP,
        "focal_position_0based": FOCAL_POSITION_BP,
        "focal_site_id": FOCAL_SITE_ID,
        "mutation_rate_per_bp_per_generation": MUTATION_RATE,
        "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        "generation_time_years": GENERATION_TIME_YEARS,
        "selection_coefficient": SELECTION_COEFFICIENT,
        "dominance_coefficient": DOMINANCE_COEFFICIENT,
        "standing_frequency": STANDING_FREQUENCY,
        "standing_variation_generations_ago_nominal": PULSE_GENERATIONS,
        "han_split_generations_ago_nominal": HAN_SPLIT_GENERATIONS,
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
        "present_day_af_conditioning": None,
        "allele_survival_conditioning": None,
        "null_conditioning": None,
        "statistic": (
            "overall within-individual P(TMRCA < x) at the focal position, "
            "over all sampled diploids"
        ),
        "p_value_rule": "(1 + #{null >= observed}) / (1 + n_null), one-sided upper tail",
        "p_value_floor": 1.0 / (N_REPLICATES + 1),
    }


def arm_config_record(arm: Run2Arm) -> dict[str, Any]:
    """Return the full JSON-serializable configuration for one arm."""

    record: dict[str, Any] = {
        "schema": "gamma-smc.run2-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "demography_id": arm.demography_id,
        "target_population": arm.target_population,
        "origin_population": arm.origin_population,
        "introgression": arm.introgression,
        "standing_variation_placement": arm.placement,
        "description": arm.description,
        "shared": shared_parameter_record(),
        "modes": {
            mode: {
                "n_replicates": N_REPLICATES,
                "seeds": [arm.seed(mode, i) for i in range(N_REPLICATES)],
                "unit_ids": [arm.unit_id(mode, i) for i in range(N_REPLICATES)],
            }
            for mode in MODES
        },
    }
    if arm.introgression:
        record["source"] = {
            "kind": "stdpopsim_catalog",
            "species": "HomSap",
            "model": arm.demography_id,
        }
    else:
        record["source"] = {
            "kind": "tracked_phlash_artifact",
            "path": PHLASH_EAS_RELATIVE_PATH,
            "sha256": PHLASH_EAS_SHA256,
            "curve": "pointwise_median",
        }
    json.dumps(record)
    return record


def write_arm_configs(destination: str | Path) -> dict[str, Path]:
    """Write ``run2_<arm>.json`` for every arm and return the written paths."""

    directory = Path(destination)
    directory.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for arm in ARMS:
        suffix = "chb" if arm.introgression else "eas"
        path = directory / f"run2_{suffix}.json"
        path.write_text(
            json.dumps(arm_config_record(arm), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        written[arm.arm_id] = path
    return written


__all__ = [
    "ARMS",
    "ARMS_BY_ID",
    "CHB_ARM",
    "DOMINANCE_COEFFICIENT",
    "EAS_ARM",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "GENERATION_TIME_YEARS",
    "HAN_SPLIT_GENERATIONS",
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
    "Run2Arm",
    "SAMPLE_DIPLOIDS",
    "SEED_BASE",
    "SELECTION_COEFFICIENT",
    "SEQUENCE_LENGTH_BP",
    "SIGNIFICANCE_LEVEL",
    "SLIM_BURN_IN",
    "SLIM_SCALING_FACTOR",
    "STANDING_FREQUENCY",
    "TMRCA_THRESHOLDS_YEARS",
    "arm_config_record",
    "get_arm",
    "shared_parameter_record",
    "write_arm_configs",
]
