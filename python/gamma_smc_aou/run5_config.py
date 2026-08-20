"""Frozen parameters for the run5 study: selection from introgression onset.

run4 showed that an EPAS1-calibrated sweep (s = 0.05 from 6,000 years ago) is
detectable at 0.73 power, and that what drives detectability is **starting
frequency and sweep age**, not the archaic background as such.  run5 asks the
complementary question: what if selection is weaker (s = 0.01) but runs for the
whole time the allele has been in CHB, starting at introgression rather than
6,000 years ago?

Two arms, each against a plain neutral null carrying no focal allele:

1. ``chb_from_introgression`` -- the allele is fixed in Neanderthal after the
   archaic split, introgresses into CHB, and is selected at s = 0.01 in CHB from
   the moment introgression begins (the CHB founding, 1,240 generations ago)
   through to the present.
2. ``chb_nea_selection`` -- the same, except the allele is *also* under selection
   in the Neanderthal branch itself from just after the human-Neanderthal split.
   A sweep inside the archaic source compresses archaic haplotype diversity, so
   the material that introgresses should descend from fewer, longer haplotypes.
   Whether that raises power is the question the arm exists to answer.

**Ascertainment.** Both arms require the introgressed allele to actually reach an
appreciable frequency in CHB, rather than the near-zero that continuous
low-rate migration usually delivers at any single base (run3 measured 62.5% of
replicates with no archaic ancestry at the focal base at all).  The requirement
is expressed at the end of the archaic migration window and enforced by
stdpopsim's own ``ConditionOnAlleleFrequency``, so SLiM re-draws from its
checkpoint rather than the study rejecting replicates after the fact.
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
FOCAL_SITE_ID = "run5_focal_site"

MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8

#: Weaker than run4's 0.05, but acting for 1,240 generations instead of 205.
SELECTION_COEFFICIENT = 0.01
DOMINANCE_COEFFICIENT = 0.5

#: Minimum CHB frequency required at the end of the archaic migration window.
#: The nominal Neanderthal ancestry proportion in East Asians is ~2%; the
#: model's own expectation from its migration schedule is ~1.2%. Calibrated
#: against the achievable accept rate -- see ASCERTAINMENT_CALIBRATION.
MIN_INTROGRESSED_FREQUENCY = 0.02

#: The archaic arm with Neanderthal selection needs the allele to survive in the
#: source long enough to sweep. Starting it as standing variation in Neanderthal
#: rather than a single copy keeps the accept rate workable: a single copy fixes
#: with probability ~2hs, which would multiply the rejection cost by ~20.
NEA_STARTING_FREQUENCY = 0.05

# ---------------------------------------------------------------------------
# EAS calibration arm (msprime structured coalescent)
# ---------------------------------------------------------------------------

#: An idealized reference, not a model of anything real: a partial sweep from
#: standing variation to a set present-day frequency.
#:
#: msprime rejects demographic events *during* a sweep, and the PHLASH curve has
#: ten thousand of them. Rather than flattening the whole history, the curve is
#: kept as discrete constant-size epochs everywhere older than the sweep window,
#: and held flat (at the harmonic mean, the size governing coalescence over a
#: period of changing size) only across the window itself. That preserves the
#: real deep EAS history, which is what sets the neutral TMRCA background and
#: therefore the null distribution.
EAS_START_FREQUENCY = 0.02
EAS_END_FREQUENCY = 0.75
EAS_SWEEP_DT = 1e-6
EAS_GENERATION_TIME = 25.0
EAS_POPULATION = "EAS"
PHLASH_EAS_RELATIVE_PATH = "no_introgression_EAS_sim/resources/EAS.npz"
#: The sweep window is held flat out to this many generations; the deterministic
#: duration at these parameters is ~998 generations, so this leaves margin.
EAS_FLAT_WINDOW_GENERATIONS = 1_500.0
#: Number of geometric epochs used to represent the PHLASH curve beyond it.
EAS_EPOCH_COUNT = 25
EAS_OLDEST_EPOCH_GENERATIONS = 40_000.0

SLIM_SCALING_FACTOR = 5.0
SLIM_BURN_IN = 0.1
REQUIRED_SLIM_VERSION = "4.2.2"
REQUIRED_STDPOPSIM_VERSION = "0.3.0"

N_REPLICATES = 100
SAMPLE_DIPLOIDS = 100

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
SEED_BASE = 20260822

# ---------------------------------------------------------------------------
# Catalog contract
# ---------------------------------------------------------------------------

OOA_MODEL_ID = "OutOfAfricaArchaicAdmixture_5R19"
OOA_GENERATION_TIME = 29.0
NEANDERTHAL_SPLIT_GENERATIONS = 19_275.86217652958
#: CHB is founded here and its Neanderthal migration switches on at the same
#: time, so this is both the introgression onset and the selection onset.
CHB_SPLIT_GENERATIONS = 1_241.379310344827
#: Neanderthal migration into CHB switches off here; the ascertainment is
#: evaluated at this tick, once all archaic input has arrived.
CHB_ARCHAIC_MIGRATION_END_GENERATIONS = 644.8275862068965
ARCHAIC_POPULATION = "Neanderthal"
CHB_POPULATION = "CHB"

#: Filled in by scripts/calibrate_run5.py.
ASCERTAINMENT_CALIBRATION: dict[str, Any] = {
    "threshold": MIN_INTROGRESSED_FREQUENCY,
    "evaluated_at_generations_ago": CHB_ARCHAIC_MIGRATION_END_GENERATIONS,
    "measured_accept_rate": None,
    "note": "run3 measured 62.5% of replicates with no archaic ancestry at the focal base",
}


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Run5Arm:
    arm_id: str
    label: str
    target_population: str
    generation_time: float
    #: "slim" for the archaic CHB arms, "msprime" for the EAS reference.
    engine: str
    #: True when the allele is also selected inside the Neanderthal branch.
    nea_selection: bool
    description: str

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
        """Selection onset in CHB: the introgression onset, snapped to the grid."""
        return (
            round(CHB_SPLIT_GENERATIONS / SLIM_SCALING_FACTOR) * SLIM_SCALING_FACTOR
        )


ARM_ORDER = ("chb_from_introgression", "chb_nea_selection", "eas_sweep_structured")

CHB_FROM_INTROGRESSION = Run5Arm(
    arm_id="chb_from_introgression",
    label="CHB, selected from introgression onset (s = 0.01)",
    target_population=CHB_POPULATION,
    generation_time=OOA_GENERATION_TIME,
    engine="slim",
    nea_selection=False,
    description=(
        "The allele is fixed in Neanderthal one tick after the archaic split, so "
        "it marks archaic ancestry exactly, and introgresses into CHB by the "
        "catalog's continuous migration. Selection acts in CHB at s = 0.01 from "
        "the CHB founding, which is also the moment archaic migration begins, "
        "through to the present -- 1,240 generations rather than run4's 205. The "
        "introgressed frequency must reach the ascertainment threshold by the end "
        "of the migration window."
    ),
)

CHB_NEA_SELECTION = Run5Arm(
    arm_id="chb_nea_selection",
    label="CHB with Neanderthal selection (s = 0.01 in both branches)",
    target_population=CHB_POPULATION,
    generation_time=OOA_GENERATION_TIME,
    engine="slim",
    nea_selection=True,
    description=(
        "As above, but the allele starts as standing variation in Neanderthal "
        "just after the human-Neanderthal split and is under selection there as "
        "well. A sweep inside the archaic source compresses archaic haplotype "
        "diversity, so the introgressing material should descend from fewer and "
        "longer haplotypes. Whether that raises detection power over the arm "
        "above is what this arm tests."
    ),
)

EAS_SWEEP_STRUCTURED = Run5Arm(
    arm_id="eas_sweep_structured",
    label="EAS structured-coalescent sweep to 75% (calibration reference)",
    target_population=EAS_POPULATION,
    generation_time=EAS_GENERATION_TIME,
    engine="msprime",
    nea_selection=False,
    description=(
        "A deliberately artificial calibration reference, not a model of any "
        "real locus: an msprime structured-coalescent partial sweep from 2% "
        "standing variation to 75% present-day frequency at s = 0.01, matched "
        "to the CHB arms in selection coefficient and roughly in duration "
        "(~1,000 generations). Because msprime rejects demographic events "
        "during a sweep, the population is constant-size at the harmonic mean "
        "of the PHLASH median trajectory over the sweep window. A partial sweep "
        "gives a bimodal TMRCA distribution -- swept pairs share a recent "
        "ancestor, mixed and unswept pairs do not -- which is the point: it "
        "shows what the statistic sees when a known fraction of the sample has "
        "been swept."
    ),
)

ARMS: tuple[Run5Arm, ...] = (
    CHB_FROM_INTROGRESSION,
    CHB_NEA_SELECTION,
    EAS_SWEEP_STRUCTURED,
)
ARMS_BY_ID = {arm.arm_id: arm for arm in ARMS}


def get_arm(arm_id: str) -> Run5Arm:
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
    "chb_from_introgression": 65_834.77001122756,
    "chb_nea_selection": 65_834.77001122756,
    "eas_sweep_structured": 4_884.0,
}


def scaled_mutation_rate(arm: Run5Arm) -> float:
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
        "selection_onset": "introgression onset (CHB founding)",
        "min_introgressed_frequency": MIN_INTROGRESSED_FREQUENCY,
        "nea_starting_frequency": NEA_STARTING_FREQUENCY,
        "ascertainment_calibration": ASCERTAINMENT_CALIBRATION,
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
        "raw_tmrca_stored": True,
        "statistic": (
            "overall within-individual P(TMRCA < x) at the focal position, "
            "over all sampled diploids"
        ),
        "p_value_rule": "(1 + #{null >= observed}) / (1 + n_null), one-sided upper tail",
        "p_value_floor": 1.0 / (N_REPLICATES + 1),
    }


def arm_config_record(arm: Run5Arm) -> dict[str, Any]:
    record = {
        "schema": "gamma-smc.run5-config/v1",
        "arm_id": arm.arm_id,
        "label": arm.label,
        "demography_id": OOA_MODEL_ID,
        "target_population": arm.target_population,
        "generation_time_years": arm.generation_time,
        "nea_selection": arm.nea_selection,
        "selection_onset_generations": arm.onset_generations(),
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
        path = directory / f"run5_{arm.arm_id}.json"
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
    "ASCERTAINMENT_CALIBRATION",
    "CHB_ARCHAIC_MIGRATION_END_GENERATIONS",
    "CHB_FROM_INTROGRESSION",
    "CHB_NEA_SELECTION",
    "CHB_POPULATION",
    "CHB_SPLIT_GENERATIONS",
    "DECODER_SETTINGS",
    "DOMINANCE_COEFFICIENT",
    "EAS_END_FREQUENCY",
    "EAS_GENERATION_TIME",
    "EAS_EPOCH_COUNT",
    "EAS_FLAT_WINDOW_GENERATIONS",
    "EAS_OLDEST_EPOCH_GENERATIONS",
    "EAS_POPULATION",
    "EAS_START_FREQUENCY",
    "EAS_SWEEP_DT",
    "EAS_SWEEP_STRUCTURED",
    "PHLASH_EAS_RELATIVE_PATH",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "MIN_INTROGRESSED_FREQUENCY",
    "MODES",
    "MUTATION_RATE",
    "NEANDERTHAL_SPLIT_GENERATIONS",
    "NEA_STARTING_FREQUENCY",
    "N_REPLICATES",
    "OOA_GENERATION_TIME",
    "OOA_MODEL_ID",
    "PRESENT_NE",
    "PROFILE_STRIDE_BP",
    "RECOMBINATION_RATE",
    "REQUIRED_SLIM_VERSION",
    "REQUIRED_STDPOPSIM_VERSION",
    "Run5Arm",
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
