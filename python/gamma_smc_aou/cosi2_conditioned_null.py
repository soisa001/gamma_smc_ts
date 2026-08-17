"""Conditioned-neutral CoSi2 inputs for the focused EAS/Han comparison.

This module is additive: it does not alter the completed natural-neutral bank.
The two nulls deliberately have different ascertainment because their allele
origins differ.

* EAS uses a neutral Wright--Fisher time-reversal *approximation* conditional
  on the exact selected population count at generation zero.  The resulting
  frozen trajectory, including its stochastic allele age, is replayed by
  CoSi2 with ``COSI_LOAD_TRAJ``.
* Han keeps the fixed Neanderthal birth at generation 2400, uses ``s=0``, and
  accepts a trajectory only when CHB is in [0.011, 0.013] at generation 645
  and is strictly segregating at generation zero.  Accepted trajectories are
  then replayed for the full 10-Mb simulation.

Neither null conditions on sampled carriers or sampled allele frequency.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
import uuid
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import cosi2_framework as framework
from . import cosi2_neutral_simulation as natural_neutral

MODULE_PATH = Path(__file__).resolve()

PLAN_SCHEMA = "gamma-smc.cosi2-conditioned-null-plan/v1"
BRIDGE_SCHEMA = "gamma-smc.cosi2-eas-reverse-bridge/v1"
EAS_PROPOSAL_BLOCK_SCHEMA = "gamma-smc.cosi2-eas-proposal-block/v1"
EAS_BASE_BANK_MANIFEST_SCHEMA = "gamma-smc.cosi2-eas-base-bank-manifest/v2"
EAS_EXTENSION_PLAN_SCHEMA = "gamma-smc.cosi2-eas-extension-plan/v2"
EAS_EXTENSION_BLOCK_BINDING_SCHEMA = "gamma-smc.cosi2-eas-extension-block-binding/v2"
EAS_PROPOSAL_BANK_CONTRACT_ID = "eas-slatkin-global-bank-v2-post-pilot-10m"
EAS_PROPOSAL_KERNEL_ID = "slatkin-reverse-wf-terminal-loss-v1"
EAS_SELECTED_PATH_SCHEMA = "gamma-smc.cosi2-eas-selected-path/v1"
UNIT_SCHEMA = "gamma-smc.cosi2-conditioned-null-unit/v1"
HAN_SCREEN_SCHEMA = "gamma-smc.cosi2-han-conditioned-screen/v1"
HAN_SELECTION_SCHEMA = "gamma-smc.cosi2-han-conditioned-selection/v1"

DEFAULT_WORK_RELATIVE = Path("focused_selection_EAS_sim/cosi2_conditioned_null_work")
PLAN_DIRNAME = "plan"
EAS_RUN_DIRNAME = "eas_units"
HAN_SCREEN_DIRNAME = "han_screen"
HAN_RUN_DIRNAME = "han_units"

EAS_MODEL_ID = "eas_neutral_exact_population_af_reverse_bridge"
HAN_MODEL_ID = "han_introgressed_neutral_g645_gate"

ACCEPTED_PER_MODEL = 100
EAS_SEED_START = 2_026_083_600
HAN_NEW_SEED_START = 2_026_084_600
HAN_CANDIDATE_COUNT = 5_000
HAN_REUSED_CANDIDATE_COUNT = 100
HAN_NEW_CANDIDATE_COUNT = HAN_CANDIDATE_COUNT - HAN_REUSED_CANDIDATE_COUNT

EAS_PRESENT_DIPLOID_NE = 37_637
EAS_PRESENT_CHROMOSOMES = 2 * EAS_PRESENT_DIPLOID_NE
EAS_SELECTED_PRESENT_COUNT = 14_704
EAS_SELECTED_PRESENT_AF = EAS_SELECTED_PRESENT_COUNT / EAS_PRESENT_CHROMOSOMES
EAS_SELECTED_PRESENT_AF_REPORTED = 0.195339692324
EAS_ANCIENT_TERMINAL_GENERATION = 40_000
EAS_ANCIENT_TERMINAL_NE = 13_138
# The oldest admissible one-copy birth row is generation 500,000.  Detecting
# loss immediately older than that row requires one additional reverse draw.
EAS_BRIDGE_MAX_GENERATION = 500_000
EAS_BRIDGE_ABSORPTION_GENERATION = EAS_BRIDGE_MAX_GENERATION + 1
EAS_BRIDGE_SEED_NAMESPACE = 20_260_817
EAS_V1_PROPOSAL_BLOCK_COUNT = 20
EAS_BASE_PROPOSAL_BLOCK_COUNT = EAS_V1_PROPOSAL_BLOCK_COUNT
EAS_BASE_IMPLEMENTATION_SHA256 = (
    "e829711cdd55f822b7db8cb8e56ff16200eb681a3de88a3ba7ac9fcc52b83cc0"
)
EAS_PROPOSALS_PER_BLOCK = 50_000
EAS_BASE_TOTAL_PROPOSALS = EAS_BASE_PROPOSAL_BLOCK_COUNT * EAS_PROPOSALS_PER_BLOCK
EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT = 200
EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT = (
    EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT - EAS_BASE_PROPOSAL_BLOCK_COUNT
)
EAS_V2_TOTAL_PROPOSALS = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT * EAS_PROPOSALS_PER_BLOCK
# Compatibility aliases now name the fixed v2 production target.  Explicit
# base/v2 constants above preserve the append-only boundary in every ledger.
EAS_PROPOSAL_BLOCK_COUNT = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT
EAS_TOTAL_PROPOSALS = EAS_PROPOSAL_BLOCK_COUNT * EAS_PROPOSALS_PER_BLOCK
EAS_PROPOSAL_BLOCK_SEED_START = 2_026_083_500
EAS_CATEGORICAL_RESAMPLE_SEED = 2_026_083_599
# Fail-closed global-bank diagnostics.  These are Monte Carlo quality gates,
# not claims that the importance-corrected approximation is exact.
EAS_MIN_GLOBAL_ESS = 10_000.0
EAS_MIN_WEIGHT_SUM_OVER_MAX = 100.0
EAS_MAX_NORMALIZED_WEIGHT = 0.01
EAS_V1_PILOT_MIN_WEIGHT_SUM_OVER_MAX = 1_000.0
EAS_V1_PILOT_MAX_NORMALIZED_WEIGHT_EXCLUSIVE = 0.01
EAS_MAX_CAP_HIT_UPPER_95 = 0.01
EAS_GENERATION_TIME_YEARS = 25.0

HAN_BIRTH_GENERATION = framework.DEFAULT_HAN_MUTATION_BIRTH_GENERATIONS
HAN_GATE_GENERATION = framework.HAN_MIGRATION_END_GENERATIONS
HAN_GATE_LOWER = 0.011
HAN_GATE_UPPER = 0.013
HAN_PRESENT_LOWER = 1e-12
HAN_PRESENT_UPPER = 1.0 - 1e-12
HAN_PRESENT_TOKEN = "0.000000000001-0.999999999999"
HAN_GENERATION_TIME_YEARS = framework.HAN_GENERATION_TIME_YEARS

SEQUENCE_LENGTH_BP = framework.SEQUENCE_LENGTH_BP
SCREEN_SEQUENCE_LENGTH_BP = 2
SAMPLE_HAPLOIDS = framework.DEFAULT_SAMPLE_HAPLOIDS
MAX_WORKERS = 20


@dataclass(frozen=True)
class ConditionedNullPlan:
    """Immutable high-level simulation contract."""

    accepted_per_model: int = ACCEPTED_PER_MODEL
    eas_seed_start: int = EAS_SEED_START
    eas_present_count: int = EAS_SELECTED_PRESENT_COUNT
    eas_present_chromosomes: int = EAS_PRESENT_CHROMOSOMES
    eas_bridge_max_generation: int = EAS_BRIDGE_MAX_GENERATION
    eas_proposal_bank_contract_id: str = EAS_PROPOSAL_BANK_CONTRACT_ID
    eas_base_proposal_block_count: int = EAS_BASE_PROPOSAL_BLOCK_COUNT
    eas_proposal_block_count: int = EAS_PROPOSAL_BLOCK_COUNT
    eas_extension_proposal_block_count: int = EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT
    eas_proposals_per_block: int = EAS_PROPOSALS_PER_BLOCK
    eas_total_proposals: int = EAS_TOTAL_PROPOSALS
    eas_proposal_block_seed_start: int = EAS_PROPOSAL_BLOCK_SEED_START
    eas_categorical_resample_seed: int = EAS_CATEGORICAL_RESAMPLE_SEED
    eas_min_global_ess: float = EAS_MIN_GLOBAL_ESS
    eas_min_weight_sum_over_max: float = EAS_MIN_WEIGHT_SUM_OVER_MAX
    eas_max_normalized_weight: float = EAS_MAX_NORMALIZED_WEIGHT
    eas_max_cap_hit_upper_95: float = EAS_MAX_CAP_HIT_UPPER_95
    han_candidate_count: int = HAN_CANDIDATE_COUNT
    han_reused_candidate_count: int = HAN_REUSED_CANDIDATE_COUNT
    han_new_seed_start: int = HAN_NEW_SEED_START
    han_gate_lower: float = HAN_GATE_LOWER
    han_gate_upper: float = HAN_GATE_UPPER
    sample_haploids: int = SAMPLE_HAPLOIDS

    def validate(self) -> None:
        if self != ConditionedNullPlan():
            raise ValueError("the conditioned-null production contract is immutable")
        EasProposalBankContract().validate()
        if not math.isclose(
            EAS_SELECTED_PRESENT_AF,
            EAS_SELECTED_PRESENT_AF_REPORTED,
            rel_tol=0.0,
            abs_tol=5e-13,
        ):
            raise RuntimeError(
                "selected EAS endpoint count no longer matches its report"
            )
        eas_seeds = set(range(EAS_SEED_START, EAS_SEED_START + ACCEPTED_PER_MODEL))
        han_seeds = set(
            range(HAN_NEW_SEED_START, HAN_NEW_SEED_START + HAN_NEW_CANDIDATE_COUNT)
        )
        if eas_seeds.intersection(han_seeds):
            raise RuntimeError("conditioned-null seed streams overlap")
        if max(eas_seeds | han_seeds) >= 2**32:
            raise RuntimeError("conditioned-null seeds must fit uint32")

    def to_record(self) -> dict[str, Any]:
        self.validate()
        return asdict(self)


@dataclass(frozen=True)
class EasProposalBankContract:
    """Named, fixed post-pilot proposal-bank contract for inference."""

    contract_id: str = EAS_PROPOSAL_BANK_CONTRACT_ID
    base_block_count: int = EAS_BASE_PROPOSAL_BLOCK_COUNT
    target_block_count: int = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT
    extension_block_count: int = EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT
    proposals_per_block: int = EAS_PROPOSALS_PER_BLOCK
    total_proposals: int = EAS_V2_TOTAL_PROPOSALS
    first_extension_block_index: int = EAS_BASE_PROPOSAL_BLOCK_COUNT
    last_extension_block_index: int = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT - 1
    min_global_effective_sample_size: float = EAS_MIN_GLOBAL_ESS
    min_weight_sum_over_max: float = EAS_MIN_WEIGHT_SUM_OVER_MAX
    max_normalized_weight_inclusive: float = EAS_MAX_NORMALIZED_WEIGHT
    require_zero_cap_hits: bool = True
    fixed_target_selected_before_extension: bool = True
    adaptive_stopping_allowed: bool = False
    quality_gate_rationale: str = (
        "W=sum(weights)/max(weights)=1/max_normalized_weight; W>=100 is the "
        "transparent equivalent of the inclusive max_weight<=0.01 gate. The "
        "v1 pilot used W>=1000; v2 records the post-pilot revision explicitly."
    )

    def validate(self) -> None:
        if self != EasProposalBankContract():
            raise ValueError("the EAS v2 proposal-bank contract is immutable")
        if (
            self.extension_block_count
            != self.target_block_count - self.base_block_count
            or self.total_proposals
            != self.target_block_count * self.proposals_per_block
            or self.first_extension_block_index != self.base_block_count
            or self.last_extension_block_index != self.target_block_count - 1
            or not math.isclose(
                self.min_weight_sum_over_max,
                1.0 / self.max_normalized_weight_inclusive,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or self.adaptive_stopping_allowed
        ):
            raise RuntimeError("the EAS v2 proposal-bank arithmetic changed")

    def to_record(self) -> dict[str, Any]:
        self.validate()
        return asdict(self)


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        temporary.write_text(text, encoding="utf-8", newline="\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n",
    )


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    _atomic_text(
        path,
        frame.to_csv(
            sep="\t",
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        ),
    )


def _positive_workers(value: int) -> int:
    if (
        isinstance(value, bool)
        or int(value) != value
        or not 1 <= int(value) <= MAX_WORKERS
    ):
        raise ValueError(f"workers must be an integer in [1, {MAX_WORKERS}]")
    return int(value)


def eas_integer_ne_by_generation(
    repo_root: str | Path,
    *,
    max_generation: int = EAS_BRIDGE_ABSORPTION_GENERATION,
) -> np.ndarray:
    """Expand PHLASH Ne through the extra generation needed to audit loss."""

    cap = int(max_generation)
    if isinstance(max_generation, bool) or cap != max_generation or cap < 1:
        raise ValueError("max_generation must be a positive integer")
    history = framework.load_eas_history(Path(repo_root).resolve())
    times = history["time_generations"].to_numpy(dtype=float)
    sizes = history["cosi2_integer_ne"].to_numpy(dtype=np.int64)
    generations = np.arange(cap + 1, dtype=np.int64)
    indices = np.searchsorted(times, generations.astype(float), side="right") - 1
    result = sizes[np.clip(indices, 0, len(sizes) - 1)].astype(np.int64)
    result[generations > EAS_ANCIENT_TERMINAL_GENERATION] = EAS_ANCIENT_TERMINAL_NE
    if (
        int(result[0]) != EAS_PRESENT_DIPLOID_NE
        or int(result[EAS_ANCIENT_TERMINAL_GENERATION]) != EAS_ANCIENT_TERMINAL_NE
        or int(result[-1]) != EAS_ANCIENT_TERMINAL_NE
    ):
        raise RuntimeError("expanded EAS integer history changed unexpectedly")
    return result


def _bridge_rng(unit_seed: int, stream: int) -> np.random.Generator:
    if min(unit_seed, stream) < 0:
        raise ValueError("bridge seed and stream must be nonnegative")
    seed_sequence = np.random.SeedSequence(
        [EAS_BRIDGE_SEED_NAMESPACE, int(unit_seed), int(stream)]
    )
    return np.random.Generator(np.random.PCG64(seed_sequence))


def _log_binomial_pmf(
    observed: np.ndarray,
    trials: int,
    probability: np.ndarray,
) -> np.ndarray:
    """Vectorized, boundary-safe log Binomial probability mass."""

    from scipy.special import gammaln, xlog1py, xlogy

    k = np.asarray(observed, dtype=np.int64)
    p = np.asarray(probability, dtype=float)
    n = int(trials)
    if np.any(k < 0) or np.any(k > n) or np.any((p < 0) | (p > 1)):
        raise ValueError("invalid Binomial probability request")
    return (
        gammaln(n + 1.0)
        - gammaln(k + 1.0)
        - gammaln(n - k + 1.0)
        + xlogy(k, p)
        + xlog1py(n - k, -p)
    )


def _weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    probabilities: Sequence[float],
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    requested = np.asarray(probabilities, dtype=float)
    if (
        values.ndim != 1
        or weights.shape != values.shape
        or not len(values)
        or np.any(~np.isfinite(values))
        or np.any(~np.isfinite(weights))
        or np.any(weights < 0)
        or not float(np.sum(weights)) > 0
        or np.any((requested < 0) | (requested > 1))
    ):
        raise ValueError("weighted quantile inputs are invalid")
    order = np.argsort(values, kind="mergesort")
    ordered_values = values[order]
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights)
    cumulative /= cumulative[-1]
    return np.interp(requested, cumulative, ordered_values)


def _reverse_path_log_components(
    integer_ne: Sequence[int],
    derived_counts_ascending: Sequence[int],
) -> dict[str, float]:
    """Recompute the Slatkin components for one present-to-birth path."""

    sizes = np.asarray(integer_ne, dtype=np.int64)
    counts = np.asarray(derived_counts_ascending, dtype=np.int64)
    if counts.ndim != 1 or not len(counts) or len(sizes) < len(counts) + 1:
        raise ValueError("reverse path lacks counts or its older loss generation")
    chromosomes = 2 * sizes[: len(counts)]
    if counts[-1] != 1 or np.any((counts <= 0) | (counts >= chromosomes)):
        raise ValueError("reverse path must be segregating and end in one copy")
    log_forward = 0.0
    log_reverse = 0.0
    for older_generation in range(1, len(counts)):
        younger_count = int(counts[older_generation - 1])
        older_count = int(counts[older_generation])
        younger_chromosomes = int(chromosomes[older_generation - 1])
        older_chromosomes = int(chromosomes[older_generation])
        log_forward += float(
            _log_binomial_pmf(
                np.asarray([younger_count]),
                younger_chromosomes,
                np.asarray([older_count / older_chromosomes]),
            )[0]
        )
        log_reverse += float(
            _log_binomial_pmf(
                np.asarray([older_count]),
                older_chromosomes,
                np.asarray([younger_count / younger_chromosomes]),
            )[0]
        )
    birth = len(counts) - 1
    terminal_ne = int(sizes[birth + 1])
    terminal_chromosomes = 2 * terminal_ne
    birth_chromosomes = int(chromosomes[birth])
    log_reverse += float(
        _log_binomial_pmf(
            np.asarray([0]),
            terminal_chromosomes,
            np.asarray([1.0 / birth_chromosomes]),
        )[0]
    )
    log_mutation = math.log(float(terminal_ne))
    return {
        "log_forward": log_forward,
        "log_reverse": log_reverse,
        "log_mutation": log_mutation,
        "log_weight": log_forward - log_reverse + log_mutation,
    }


def _simulate_reverse_proposal_block(
    integer_ne: Sequence[int],
    *,
    block_index: int,
    block_seed: int,
    proposals_per_block: int = EAS_PROPOSALS_PER_BLOCK,
    capture_indices: Sequence[int] = (),
    present_count: int,
    max_birth_generation: int,
) -> dict[str, Any]:
    """Generic deterministic vectorized reverse-proposal block.

    The reverse process draws through generation ``max_birth_generation + 1``.
    Thus, a one-copy state at the declared maximum birth generation can be
    retained only after its explicit older zero draw.  This private generic
    engine also supports the tiny constant-N oracle test.
    """

    sizes = np.asarray(integer_ne, dtype=np.int64)
    birth_cap = int(max_birth_generation)
    if birth_cap < 1 or sizes.shape != (birth_cap + 2,):
        raise ValueError("integer_ne must cover generations 0..birth_cap+1")
    endpoint_count = int(present_count)
    if not 0 < endpoint_count < 2 * int(sizes[0]):
        raise ValueError("present_count must be strictly segregating")
    block = int(block_index)
    count = int(proposals_per_block)
    if (
        isinstance(block_index, bool)
        or block != block_index
        or block < 0
        or isinstance(proposals_per_block, bool)
        or count != proposals_per_block
        or count < 2
    ):
        raise ValueError("block index/count is invalid")
    requested = tuple(sorted(set(int(value) for value in capture_indices)))
    if any(value < 0 or value >= count for value in requested):
        raise ValueError("capture proposal index lies outside its block")
    seed = int(block_seed)
    if seed < 0:
        raise ValueError("block_seed must be nonnegative")
    rng = _bridge_rng(seed, 0)
    current = np.full(count, endpoint_count, dtype=np.int64)
    active = np.ones(count, dtype=bool)
    log_forward = np.zeros(count, dtype=float)
    log_reverse = np.zeros(count, dtype=float)
    log_mutation = np.zeros(count, dtype=float)
    outcomes = np.full(count, "active", dtype="U12")
    ages = np.full(count, -1, dtype=np.int64)
    captured: dict[int, list[int]] = {index: [endpoint_count] for index in requested}

    for older_generation in range(1, birth_cap + 2):
        indexes = np.flatnonzero(active)
        if not len(indexes):
            break
        younger_ne = int(sizes[older_generation - 1])
        older_ne = int(sizes[older_generation])
        younger_chromosomes = 2 * younger_ne
        older_chromosomes = 2 * older_ne
        younger = current[indexes]
        older = rng.binomial(
            older_chromosomes, younger / float(younger_chromosomes)
        ).astype(np.int64)
        zero = older == 0
        fixation = older == older_chromosomes
        valid_terminal = zero & (younger == 1)
        invalid_terminal = zero & (younger != 1)

        if np.any(valid_terminal):
            terminal_indexes = indexes[valid_terminal]
            log_reverse[terminal_indexes] += older_chromosomes * math.log1p(
                -1.0 / younger_chromosomes
            )
            log_mutation[terminal_indexes] = math.log(float(older_ne))
            outcomes[terminal_indexes] = "valid"
            ages[terminal_indexes] = older_generation - 1
            active[terminal_indexes] = False
        if np.any(invalid_terminal):
            invalid_indexes = indexes[invalid_terminal]
            outcomes[invalid_indexes] = "invalid_loss"
            active[invalid_indexes] = False
        if np.any(fixation):
            fixation_indexes = indexes[fixation]
            outcomes[fixation_indexes] = "fixation"
            active[fixation_indexes] = False

        positive = ~(zero | fixation)
        if np.any(positive):
            positive_indexes = indexes[positive]
            positive_younger = younger[positive]
            positive_older = older[positive]
            forward_step = _log_binomial_pmf(
                positive_younger,
                younger_chromosomes,
                positive_older / float(older_chromosomes),
            )
            reverse_step = _log_binomial_pmf(
                positive_older,
                older_chromosomes,
                positive_younger / float(younger_chromosomes),
            )
            if not np.all(np.isfinite(forward_step - reverse_step)):
                raise RuntimeError("nonfinite EAS proposal transition weight")
            log_forward[positive_indexes] += forward_step
            log_reverse[positive_indexes] += reverse_step
            current[positive_indexes] = positive_older

        for capture_index in requested:
            if active[capture_index]:
                captured[capture_index].append(int(current[capture_index]))

    if np.any(active):
        outcomes[np.flatnonzero(active)] = "cap"
    log_weights = log_forward - log_reverse + log_mutation
    log_weights[outcomes != "valid"] = -np.inf
    for capture_index, path in captured.items():
        if outcomes[capture_index] == "valid":
            expected_rows = int(ages[capture_index]) + 1
            if len(path) != expected_rows or path[-1] != 1:
                raise RuntimeError("captured EAS proposal path is inconsistent")
        else:
            captured[capture_index] = []
    return {
        "schema": EAS_PROPOSAL_BLOCK_SCHEMA,
        "block_index": block,
        "block_seed": seed,
        "proposal_count": count,
        "present_count": endpoint_count,
        "max_birth_generation": birth_cap,
        "outcomes": outcomes,
        "birth_generations": ages,
        "log_forward": log_forward,
        "log_reverse": log_reverse,
        "log_mutation": log_mutation,
        "log_weights": log_weights,
        "captured_paths": {key: tuple(value) for key, value in captured.items()},
    }


def simulate_eas_proposal_block(
    integer_ne: Sequence[int],
    *,
    block_index: int,
    proposals_per_block: int = EAS_PROPOSALS_PER_BLOCK,
    capture_indices: Sequence[int] = (),
) -> dict[str, Any]:
    """Run one fixed EAS pass-A block or deterministic pass-B replay.

    Pass A omits ``capture_indices`` and persists only compact outcome, age,
    and log-probability arrays.  Pass B reruns the exact block and captures
    selected paths.  The endpoint and age cap cannot be overridden here.
    """

    sizes = np.asarray(integer_ne, dtype=np.int64)
    if sizes.shape != (EAS_BRIDGE_ABSORPTION_GENERATION + 1,):
        raise ValueError("integer_ne must cover generations 0..500001")
    block = int(block_index)
    return _simulate_reverse_proposal_block(
        sizes,
        block_index=block,
        block_seed=EAS_PROPOSAL_BLOCK_SEED_START + block,
        proposals_per_block=proposals_per_block,
        capture_indices=capture_indices,
        present_count=EAS_SELECTED_PRESENT_COUNT,
        max_birth_generation=EAS_BRIDGE_MAX_GENERATION,
    )


def validate_eas_proposal_block(
    block: Mapping[str, Any],
    *,
    enforce_production_contract: bool = True,
) -> dict[str, Any]:
    """Validate one compact pass-A block without regenerating it."""

    if block.get("schema") != EAS_PROPOSAL_BLOCK_SCHEMA:
        raise ValueError("EAS proposal-block schema changed")
    index = int(block["block_index"])
    seed = int(block["block_seed"])
    count = int(block["proposal_count"])
    present_count = int(block["present_count"])
    birth_cap = int(block["max_birth_generation"])
    if index < 0 or seed != EAS_PROPOSAL_BLOCK_SEED_START + index:
        raise ValueError("EAS proposal-block index/seed contract changed")
    if present_count != EAS_SELECTED_PRESENT_COUNT:
        raise ValueError("EAS proposal-block endpoint count changed")
    if birth_cap != EAS_BRIDGE_MAX_GENERATION:
        raise ValueError("EAS proposal-block birth cap changed")
    if count < 2 or (
        enforce_production_contract
        and (
            count != EAS_PROPOSALS_PER_BLOCK
            or not 0 <= index < EAS_PROPOSAL_BLOCK_COUNT
        )
    ):
        raise ValueError("EAS proposal-block geometry changed")

    required_arrays = (
        "outcomes",
        "birth_generations",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weights",
    )
    arrays = {name: np.asarray(block[name]) for name in required_arrays}
    if any(value.shape != (count,) for value in arrays.values()):
        raise ValueError("EAS proposal-block array geometry changed")
    outcomes = arrays["outcomes"].astype(str)
    allowed = {"valid", "invalid_loss", "fixation", "cap"}
    if not set(outcomes).issubset(allowed):
        raise ValueError("EAS proposal-block contains an unknown outcome")
    valid = outcomes == "valid"
    ages = arrays["birth_generations"].astype(np.int64)
    if not np.array_equal(np.isfinite(arrays["log_weights"]), valid):
        raise ValueError("EAS proposal finite weights disagree with status")
    if (
        np.any(ages[valid] < 0)
        or np.any(ages[valid] > birth_cap)
        or np.any(ages[~valid] != -1)
    ):
        raise ValueError("EAS proposal birth generations disagree with status")
    expected_weights = (
        arrays["log_forward"] - arrays["log_reverse"] + arrays["log_mutation"]
    )
    if not np.array_equal(arrays["log_weights"][valid], expected_weights[valid]):
        raise ValueError("EAS valid proposal importance weights are inconsistent")
    counts = {label: int(np.count_nonzero(outcomes == label)) for label in allowed}
    if sum(counts.values()) != count:
        raise RuntimeError("EAS proposal outcomes do not sum to the block size")
    valid_ages = ages[valid].astype(float)
    valid_weights = arrays["log_weights"][valid].astype(float)
    summary: dict[str, Any] = {
        "status": "valid",
        "block_index": index,
        "block_seed": seed,
        "proposal_count": count,
        "present_count": present_count,
        "max_birth_generation": birth_cap,
        "absorption_draw_generation": birth_cap + 1,
        "outcome_counts": counts,
        "cap_hit_one_sided_95pct_upper": _tail_interval_upper(counts["cap"], count),
    }
    if len(valid_ages):
        summary["valid_birth_generation_quantiles_q05_q50_q95"] = [
            float(value) for value in np.quantile(valid_ages, (0.05, 0.5, 0.95))
        ]
        summary["valid_log_weight_quantiles_q05_q50_q95"] = [
            float(value) for value in np.quantile(valid_weights, (0.05, 0.5, 0.95))
        ]
    else:
        summary["valid_birth_generation_quantiles_q05_q50_q95"] = []
        summary["valid_log_weight_quantiles_q05_q50_q95"] = []
    return summary


def _eas_proposal_source_binding(repo_root: Path) -> dict[str, Any]:
    import scipy

    return {
        "eas_resource": framework.EAS_RESOURCE_PATH,
        "eas_resource_sha256": sha256_file(repo_root / framework.EAS_RESOURCE_PATH),
        "framework_sha256": sha256_file(repo_root / framework.MODULE_PATH),
        "implementation_sha256": sha256_file(MODULE_PATH),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "python_version": sys.version,
        "bit_generator": "numpy.random.PCG64",
        "seed_sequence_namespace": EAS_BRIDGE_SEED_NAMESPACE,
        "proposal_kernel_id": EAS_PROPOSAL_KERNEL_ID,
    }


def _frozen_eas_base_source_binding(repo_root: Path) -> dict[str, Any]:
    """Return the exact v1 runtime contract with only its old code hash."""

    binding = _eas_proposal_source_binding(repo_root)
    binding["implementation_sha256"] = EAS_BASE_IMPLEMENTATION_SHA256
    # The v1 completion predates this explicit label; deterministic full-block
    # replay below supplies the equivalence attestation instead.
    binding.pop("proposal_kernel_id")
    return binding


def _is_sha256(value: Any) -> bool:
    text = str(value)
    return len(text) == 64 and all(
        character in "0123456789abcdef" for character in text
    )


def evaluate_eas_v2_bank_quality(
    diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate the fixed v2 Monte Carlo gates with exact boundary semantics.

    ``W=sum(weights)/max(weights)`` and the maximum normalized weight are
    reciprocal diagnostics.  Both are retained for readability; equality at
    ``W=100`` / ``max_weight=0.01`` passes the adopted post-pilot contract.
    """

    ess = float(diagnostics["global_effective_sample_size"])
    weight_sum_over_max = float(diagnostics["global_weight_sum_over_max"])
    max_normalized_weight = float(diagnostics["global_max_normalized_weight"])
    cap_hit_count = int(diagnostics["cap_hit_count"])
    if (
        not math.isfinite(ess)
        or ess <= 0
        or not math.isfinite(weight_sum_over_max)
        or weight_sum_over_max <= 0
        or not math.isfinite(max_normalized_weight)
        or max_normalized_weight <= 0
        or cap_hit_count < 0
        or not math.isclose(
            weight_sum_over_max,
            1.0 / max_normalized_weight,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    ):
        raise ValueError("EAS proposal-bank quality diagnostics are inconsistent")
    failures: list[str] = []
    if ess < EAS_MIN_GLOBAL_ESS:
        failures.append("global_effective_sample_size")
    if weight_sum_over_max < EAS_MIN_WEIGHT_SUM_OVER_MAX:
        failures.append("global_weight_sum_over_max")
    if max_normalized_weight > EAS_MAX_NORMALIZED_WEIGHT:
        failures.append("global_max_normalized_weight")
    if cap_hit_count != 0:
        failures.append("cap_hit_count")
    return {
        "contract_id": EAS_PROPOSAL_BANK_CONTRACT_ID,
        "min_global_effective_sample_size": EAS_MIN_GLOBAL_ESS,
        "min_weight_sum_over_max": EAS_MIN_WEIGHT_SUM_OVER_MAX,
        "max_normalized_weight_inclusive": EAS_MAX_NORMALIZED_WEIGHT,
        "require_zero_cap_hits": True,
        "observed_global_effective_sample_size": ess,
        "observed_global_weight_sum_over_max": weight_sum_over_max,
        "observed_global_max_normalized_weight": max_normalized_weight,
        "observed_cap_hit_count": cap_hit_count,
        "failed_quality_gates": failures,
        "quality_gates_passed": not failures,
    }


def _base_bank_failure_record(diagnostics: Mapping[str, Any]) -> dict[str, Any]:
    record = {
        "total_proposals": int(diagnostics["total_proposals"]),
        "valid_proposals": int(diagnostics["valid_proposals"]),
        "global_log_max_weight": float(diagnostics["global_log_max_weight"]),
        "global_weight_sum_over_max": float(diagnostics["global_weight_sum_over_max"]),
        "global_effective_sample_size": float(
            diagnostics["global_effective_sample_size"]
        ),
        "global_max_normalized_weight": float(
            diagnostics["global_max_normalized_weight"]
        ),
        "cap_hit_count": int(diagnostics["cap_hit_count"]),
        "proposal_bank_binding_sha256": str(
            diagnostics["proposal_bank_binding_sha256"]
        ),
    }
    pilot_failures: list[str] = []
    if record["global_effective_sample_size"] < EAS_MIN_GLOBAL_ESS:
        pilot_failures.append("global_effective_sample_size")
    if record["global_weight_sum_over_max"] < EAS_V1_PILOT_MIN_WEIGHT_SUM_OVER_MAX:
        pilot_failures.append("global_weight_sum_over_max")
    if (
        record["global_max_normalized_weight"]
        >= EAS_V1_PILOT_MAX_NORMALIZED_WEIGHT_EXCLUSIVE
    ):
        pilot_failures.append("global_max_normalized_weight")
    if record["cap_hit_count"] != 0:
        pilot_failures.append("cap_hit_count")
    v2_quality = evaluate_eas_v2_bank_quality(record)
    record["pilot_failed_quality_gates"] = pilot_failures
    record["pilot_quality_gates_passed"] = not pilot_failures
    record["v2_failed_quality_gates"] = v2_quality["failed_quality_gates"]
    record["v2_quality_gates_passed"] = v2_quality["quality_gates_passed"]
    return record


def _base_manifest_payload(manifest: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in manifest.items()
        if key != "manifest_payload_sha256"
    }


def _proposal_array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(_canonical_json(list(array.shape)).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def attest_eas_proposal_kernel_equivalence(
    integer_ne: Sequence[int],
    frozen_block: Mapping[str, Any],
) -> dict[str, Any]:
    """Exactly replay a frozen block with the current proposal kernel."""

    validate_eas_proposal_block(frozen_block, enforce_production_contract=False)
    index = int(frozen_block["block_index"])
    count = int(frozen_block["proposal_count"])
    replay = simulate_eas_proposal_block(
        integer_ne,
        block_index=index,
        proposals_per_block=count,
    )
    array_names = (
        "outcomes",
        "birth_generations",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weights",
    )
    array_hashes: dict[str, str] = {}
    for name in array_names:
        frozen_array = np.asarray(frozen_block[name])
        replay_array = np.asarray(replay[name])
        if name == "outcomes":
            equal = np.array_equal(frozen_array, replay_array)
        else:
            equal = np.array_equal(frozen_array, replay_array, equal_nan=True)
        if not equal:
            raise RuntimeError(
                f"current EAS proposal kernel changed frozen block {index}/{name}"
            )
        frozen_hash = _proposal_array_sha256(frozen_array)
        if frozen_hash != _proposal_array_sha256(replay_array):
            raise RuntimeError("equal EAS proposal arrays have different hashes")
        array_hashes[name] = frozen_hash
    record = {
        "status": "exact_match",
        "proposal_kernel_id": EAS_PROPOSAL_KERNEL_ID,
        "attested_block_index": index,
        "attested_block_seed": int(frozen_block["block_seed"]),
        "attested_proposal_count": count,
        "frozen_implementation_sha256": EAS_BASE_IMPLEMENTATION_SHA256,
        "current_implementation_sha256": sha256_file(MODULE_PATH),
        "numpy_version": np.__version__,
        "array_sha256": array_hashes,
    }
    record["attestation_payload_sha256"] = _canonical_sha256(record)
    return record


def _eas_base_bank_contract_record() -> dict[str, Any]:
    return {
        "base_block_count": EAS_BASE_PROPOSAL_BLOCK_COUNT,
        "proposals_per_block": EAS_PROPOSALS_PER_BLOCK,
        "base_total_proposals": EAS_BASE_TOTAL_PROPOSALS,
        "base_first_block_index": 0,
        "base_last_block_index": EAS_BASE_PROPOSAL_BLOCK_COUNT - 1,
        "extension_required": True,
        "quality_gates_revised_post_pilot": True,
        "pilot_quality_gate_contract": {
            "min_global_effective_sample_size": EAS_MIN_GLOBAL_ESS,
            "min_weight_sum_over_max": EAS_V1_PILOT_MIN_WEIGHT_SUM_OVER_MAX,
            "max_normalized_weight_exclusive": (
                EAS_V1_PILOT_MAX_NORMALIZED_WEIGHT_EXCLUSIVE
            ),
            "require_zero_cap_hits": True,
        },
        "v2_quality_gate_contract": {
            "min_global_effective_sample_size": EAS_MIN_GLOBAL_ESS,
            "min_weight_sum_over_max": EAS_MIN_WEIGHT_SUM_OVER_MAX,
            "max_normalized_weight_inclusive": EAS_MAX_NORMALIZED_WEIGHT,
            "require_zero_cap_hits": True,
        },
    }


def build_eas_base_bank_manifest(
    repo_root: str | Path,
    proposal_dir: str | Path,
) -> dict[str, Any]:
    """Revalidate and bind the immutable 20-block bank for v2 extension.

    This is read-only.  It accepts only the frozen v1 implementation/runtime
    tuple and exact block inventories, then independently recomputes the base
    bank diagnostics which motivated extension.
    """

    root = Path(repo_root).resolve()
    proposal_root = Path(proposal_dir).resolve()
    source_binding = _frozen_eas_base_source_binding(root)
    blocks: list[dict[str, Any]] = []
    entries: list[dict[str, Any]] = []
    for index in range(EAS_BASE_PROPOSAL_BLOCK_COUNT):
        bundle = proposal_root / f"block_{index:02d}"
        completion_path = bundle / "completion.json"
        expected_inventory = {
            "completion.json",
            "proposal_block.npz",
            "summary.json",
        }
        if (
            bundle.is_symlink()
            or not bundle.is_dir()
            or {path.name for path in bundle.iterdir()} != expected_inventory
            or not completion_path.is_file()
        ):
            raise ValueError(f"frozen EAS base block inventory changed: {index}")
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
        if (
            completion.get("schema") != EAS_PROPOSAL_BLOCK_SCHEMA
            or completion.get("status") != "complete"
            or int(completion.get("block_index", -1)) != index
            or int(completion.get("block_seed", -1))
            != EAS_PROPOSAL_BLOCK_SEED_START + index
            or int(completion.get("proposal_count", -1)) != EAS_PROPOSALS_PER_BLOCK
            or int(completion.get("present_count", -1)) != EAS_SELECTED_PRESENT_COUNT
            or int(completion.get("max_birth_generation", -1))
            != EAS_BRIDGE_MAX_GENERATION
            or int(completion.get("absorption_draw_generation", -1))
            != EAS_BRIDGE_ABSORPTION_GENERATION
            or completion.get("source_binding") != source_binding
        ):
            raise ValueError(f"frozen EAS base block contract changed: {index}")
        outputs = completion.get("outputs")
        if not isinstance(outputs, dict) or set(outputs) != {
            "proposal_block.npz",
            "summary.json",
        }:
            raise ValueError(f"frozen EAS base output manifest changed: {index}")
        for name, output in outputs.items():
            path = bundle / name
            if (
                path.is_symlink()
                or not path.is_file()
                or path.stat().st_size != int(output.get("size_bytes", -1))
                or sha256_file(path) != output.get("sha256")
            ):
                raise ValueError(f"frozen EAS base checksum failed: {index}/{name}")
        summary = json.loads((bundle / "summary.json").read_text(encoding="utf-8"))
        if completion.get("summary_sha256") != _canonical_sha256(summary):
            raise ValueError(f"frozen EAS base summary binding changed: {index}")
        with np.load(bundle / "proposal_block.npz", allow_pickle=False) as archive:
            arrays = {name: np.asarray(archive[name]) for name in archive.files}
        if set(arrays) != {
            "outcomes",
            "birth_generations",
            "log_forward",
            "log_reverse",
            "log_mutation",
            "log_weights",
        }:
            raise ValueError(f"frozen EAS base array inventory changed: {index}")
        block: dict[str, Any] = {
            "schema": EAS_PROPOSAL_BLOCK_SCHEMA,
            "block_index": index,
            "block_seed": EAS_PROPOSAL_BLOCK_SEED_START + index,
            "proposal_count": EAS_PROPOSALS_PER_BLOCK,
            "present_count": EAS_SELECTED_PRESENT_COUNT,
            "max_birth_generation": EAS_BRIDGE_MAX_GENERATION,
            **arrays,
            "captured_paths": {},
            "artifact_dir": str(bundle),
            "proposal_block_sha256": outputs["proposal_block.npz"]["sha256"],
            "completion_sha256": sha256_file(completion_path),
            "source_binding_kind": "frozen_v1_base_manifest",
        }
        observed_summary = validate_eas_proposal_block(
            block, enforce_production_contract=False
        )
        if observed_summary != summary:
            raise ValueError(f"frozen EAS base summary is stale: {index}")
        blocks.append(block)
        entries.append(
            {
                "block_index": index,
                "block_seed": EAS_PROPOSAL_BLOCK_SEED_START + index,
                "proposal_count": EAS_PROPOSALS_PER_BLOCK,
                "completion_size_bytes": completion_path.stat().st_size,
                "completion_sha256": sha256_file(completion_path),
                "proposal_block_size_bytes": int(
                    outputs["proposal_block.npz"]["size_bytes"]
                ),
                "proposal_block_sha256": str(outputs["proposal_block.npz"]["sha256"]),
                "summary_size_bytes": int(outputs["summary.json"]["size_bytes"]),
                "summary_file_sha256": str(outputs["summary.json"]["sha256"]),
                "summary_payload_sha256": str(completion["summary_sha256"]),
            }
        )
    _, diagnostics = select_eas_global_proposals(
        blocks,
        accepted_count=ACCEPTED_PER_MODEL,
        enforce_production_contract=False,
    )
    failure = _base_bank_failure_record(diagnostics)
    if failure["pilot_quality_gates_passed"]:
        raise RuntimeError("frozen million-particle bank no longer motivates extension")
    kernel_attestation = attest_eas_proposal_kernel_equivalence(
        eas_integer_ne_by_generation(root), blocks[0]
    )
    manifest: dict[str, Any] = {
        "schema": EAS_BASE_BANK_MANIFEST_SCHEMA,
        "status": "complete",
        "contract": _eas_base_bank_contract_record(),
        "source_binding": source_binding,
        "source_binding_sha256": _canonical_sha256(source_binding),
        "blocks": entries,
        "base_failure_diagnostics": failure,
        "kernel_equivalence_attestation": kernel_attestation,
    }
    manifest["manifest_payload_sha256"] = _canonical_sha256(
        _base_manifest_payload(manifest)
    )
    validate_eas_base_bank_manifest(root, manifest)
    return manifest


def validate_eas_base_bank_manifest(
    repo_root: str | Path,
    manifest: Mapping[str, Any],
    *,
    proposal_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Validate a frozen base manifest and optionally all referenced bytes."""

    root = Path(repo_root).resolve()
    expected_contract = _eas_base_bank_contract_record()
    expected_keys = {
        "schema",
        "status",
        "contract",
        "source_binding",
        "source_binding_sha256",
        "blocks",
        "base_failure_diagnostics",
        "kernel_equivalence_attestation",
        "manifest_payload_sha256",
    }
    if (
        set(manifest) != expected_keys
        or manifest.get("schema") != EAS_BASE_BANK_MANIFEST_SCHEMA
        or manifest.get("status") != "complete"
        or manifest.get("contract") != expected_contract
        or manifest.get("source_binding") != _frozen_eas_base_source_binding(root)
        or manifest.get("source_binding_sha256")
        != _canonical_sha256(manifest["source_binding"])
        or manifest.get("manifest_payload_sha256")
        != _canonical_sha256(_base_manifest_payload(manifest))
    ):
        raise ValueError("frozen EAS base-bank manifest contract changed")
    entries = manifest.get("blocks")
    if not isinstance(entries, list) or len(entries) != EAS_BASE_PROPOSAL_BLOCK_COUNT:
        raise ValueError("frozen EAS base-bank block inventory changed")
    for index, entry in enumerate(entries):
        if (
            not isinstance(entry, dict)
            or int(entry.get("block_index", -1)) != index
            or int(entry.get("block_seed", -1)) != EAS_PROPOSAL_BLOCK_SEED_START + index
            or int(entry.get("proposal_count", -1)) != EAS_PROPOSALS_PER_BLOCK
            or any(
                not _is_sha256(entry.get(key))
                for key in (
                    "completion_sha256",
                    "proposal_block_sha256",
                    "summary_file_sha256",
                    "summary_payload_sha256",
                )
            )
            or any(
                int(entry.get(key, 0)) <= 0
                for key in (
                    "completion_size_bytes",
                    "proposal_block_size_bytes",
                    "summary_size_bytes",
                )
            )
        ):
            raise ValueError(f"frozen EAS base manifest entry changed: {index}")
    failure = manifest.get("base_failure_diagnostics")
    if (
        not isinstance(failure, dict)
        or _base_bank_failure_record(failure) != failure
        or int(failure.get("total_proposals", -1)) != EAS_BASE_TOTAL_PROPOSALS
        or bool(failure.get("pilot_quality_gates_passed", True))
        or not failure.get("pilot_failed_quality_gates")
        or bool(failure.get("v2_quality_gates_passed", True))
        or not failure.get("v2_failed_quality_gates")
        or not _is_sha256(failure.get("proposal_bank_binding_sha256"))
        or not math.isclose(
            float(failure.get("global_weight_sum_over_max", math.nan)),
            1.0 / float(failure.get("global_max_normalized_weight", math.nan)),
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    ):
        raise ValueError("frozen EAS base failure diagnostics changed")
    attestation = manifest.get("kernel_equivalence_attestation")
    if (
        not isinstance(attestation, dict)
        or attestation.get("status") != "exact_match"
        or attestation.get("proposal_kernel_id") != EAS_PROPOSAL_KERNEL_ID
        or int(attestation.get("attested_block_index", -1)) != 0
        or int(attestation.get("attested_block_seed", -1))
        != EAS_PROPOSAL_BLOCK_SEED_START
        or int(attestation.get("attested_proposal_count", -1))
        != EAS_PROPOSALS_PER_BLOCK
        or attestation.get("frozen_implementation_sha256")
        != EAS_BASE_IMPLEMENTATION_SHA256
        or attestation.get("current_implementation_sha256") != sha256_file(MODULE_PATH)
        or attestation.get("numpy_version") != np.__version__
        or set(attestation.get("array_sha256", {}))
        != {
            "outcomes",
            "birth_generations",
            "log_forward",
            "log_reverse",
            "log_mutation",
            "log_weights",
        }
        or any(
            not _is_sha256(value)
            for value in attestation.get("array_sha256", {}).values()
        )
        or attestation.get("attestation_payload_sha256")
        != _canonical_sha256(
            {
                key: value
                for key, value in attestation.items()
                if key != "attestation_payload_sha256"
            }
        )
    ):
        raise ValueError("frozen/current EAS proposal-kernel attestation changed")
    if proposal_dir is not None:
        # Byte verification is intentionally cheaper than manifest creation:
        # the exact 50k-particle kernel replay is a one-time build attestation,
        # not work repeated by every restart/status phase.
        proposal_root = Path(proposal_dir).resolve()
        observed_blocks = [
            verify_eas_proposal_block(
                root,
                proposal_root / f"block_{index:02d}",
                expected_block_index=index,
                expected_proposals=EAS_PROPOSALS_PER_BLOCK,
                enforce_production_contract=False,
                base_manifest=manifest,
            )
            for index in range(EAS_BASE_PROPOSAL_BLOCK_COUNT)
        ]
        attested_arrays = attestation["array_sha256"]
        for name, expected_hash in attested_arrays.items():
            if _proposal_array_sha256(np.asarray(observed_blocks[0][name])) != (
                expected_hash
            ):
                raise ValueError(
                    "frozen EAS base block no longer matches its kernel attestation"
                )
    return {
        "status": "valid",
        "manifest_payload_sha256": str(manifest["manifest_payload_sha256"]),
        "base_block_count": EAS_BASE_PROPOSAL_BLOCK_COUNT,
        "base_total_proposals": EAS_BASE_TOTAL_PROPOSALS,
        "pilot_failed_quality_gates": list(failure["pilot_failed_quality_gates"]),
        "v2_failed_quality_gates": list(failure["v2_failed_quality_gates"]),
    }


def _make_eas_extension_plan(
    repo_root: Path,
    base_manifest: Mapping[str, Any],
    *,
    target_block_count: int,
) -> dict[str, Any]:
    validate_eas_base_bank_manifest(repo_root, base_manifest)
    contract = EasProposalBankContract()
    contract.validate()
    if int(target_block_count) != contract.target_block_count:
        raise ValueError(
            "EAS extension target must equal the fixed v2 200-block contract"
        )
    ledger: list[dict[str, Any]] = []
    for index in range(contract.target_block_count):
        is_base = index < contract.base_block_count
        base_entry = base_manifest["blocks"][index] if is_base else None
        ledger.append(
            {
                "block_index": index,
                "block_seed": EAS_PROPOSAL_BLOCK_SEED_START + index,
                "proposal_count": contract.proposals_per_block,
                "block_dirname": f"block_{index:02d}",
                "source_kind": (
                    "frozen_v1_base_manifest"
                    if is_base
                    else "current_v2_extension_implementation"
                ),
                "required_action": (
                    "reuse_and_verify_exact_bytes" if is_base else "generate_or_verify"
                ),
                "immutable_existing": is_base,
                "expected_completion_sha256": (
                    base_entry["completion_sha256"] if is_base else ""
                ),
                "expected_proposal_block_sha256": (
                    base_entry["proposal_block_sha256"] if is_base else ""
                ),
            }
        )
    plan: dict[str, Any] = {
        "schema": EAS_EXTENSION_PLAN_SCHEMA,
        "status": "complete",
        "contract": contract.to_record(),
        "base_manifest_payload_sha256": base_manifest["manifest_payload_sha256"],
        "base_proposal_bank_binding_sha256": base_manifest["base_failure_diagnostics"][
            "proposal_bank_binding_sha256"
        ],
        "base_kernel_attestation_sha256": base_manifest[
            "kernel_equivalence_attestation"
        ]["attestation_payload_sha256"],
        "extension_source_binding": _eas_proposal_source_binding(repo_root),
        "extension_source_binding_sha256": _canonical_sha256(
            _eas_proposal_source_binding(repo_root)
        ),
        "block_ledger": ledger,
        "block_ledger_sha256": _canonical_sha256(ledger),
        "post_pilot_decision": {
            "fixed_total_proposals": EAS_V2_TOTAL_PROPOSALS,
            "fixed_target_block_count": EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT,
            "adaptive_stopping_allowed": False,
            "pilot_min_weight_sum_over_max": (EAS_V1_PILOT_MIN_WEIGHT_SUM_OVER_MAX),
            "v2_min_weight_sum_over_max": EAS_MIN_WEIGHT_SUM_OVER_MAX,
            "gate_revision_rationale": contract.quality_gate_rationale,
        },
        "invalidation": {
            "requires_new_global_selection": True,
            "reuse_v1_selected_particles": False,
            "reuse_v1_selected_path_bundles": False,
            "invalidate_downstream_eas_decode_and_results": True,
            "reason": (
                "the v1 million-particle bank failed its pilot diagnostics; "
                "all inference must derive from the complete fixed v2 bank"
            ),
        },
    }
    plan["plan_payload_sha256"] = _canonical_sha256(
        {key: value for key, value in plan.items() if key != "plan_payload_sha256"}
    )
    return plan


def build_eas_extension_plan(
    repo_root: str | Path,
    base_manifest: Mapping[str, Any],
    *,
    target_block_count: int = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT,
) -> dict[str, Any]:
    """Build the deterministic append-only v2 10-million-particle plan."""

    plan = _make_eas_extension_plan(
        Path(repo_root).resolve(),
        base_manifest,
        target_block_count=int(target_block_count),
    )
    validate_eas_extension_plan(
        repo_root,
        plan,
        base_manifest=base_manifest,
        expected_target_block_count=int(target_block_count),
    )
    return plan


def validate_eas_extension_plan(
    repo_root: str | Path,
    plan: Mapping[str, Any],
    *,
    base_manifest: Mapping[str, Any],
    expected_target_block_count: int = EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT,
) -> dict[str, Any]:
    """Validate the exact named v2 target; alternate/adaptive targets fail."""

    expected = _make_eas_extension_plan(
        Path(repo_root).resolve(),
        base_manifest,
        target_block_count=int(expected_target_block_count),
    )
    if dict(plan) != expected:
        raise ValueError("EAS v2 extension plan differs from its fixed contract")
    return {
        "status": "valid",
        "contract_id": EAS_PROPOSAL_BANK_CONTRACT_ID,
        "base_block_count": EAS_BASE_PROPOSAL_BLOCK_COUNT,
        "target_block_count": EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT,
        "extension_block_count": EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT,
        "total_proposals": EAS_V2_TOTAL_PROPOSALS,
        "adaptive_stopping_allowed": False,
        "plan_payload_sha256": str(plan["plan_payload_sha256"]),
    }


def _eas_extension_block_binding(
    base_manifest: Mapping[str, Any],
    extension_plan: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": EAS_EXTENSION_BLOCK_BINDING_SCHEMA,
        "contract_id": EAS_PROPOSAL_BANK_CONTRACT_ID,
        "base_manifest_payload_sha256": str(base_manifest["manifest_payload_sha256"]),
        "extension_plan_payload_sha256": str(extension_plan["plan_payload_sha256"]),
        "extension_source_binding_sha256": str(
            extension_plan["extension_source_binding_sha256"]
        ),
    }


def write_eas_proposal_block(
    repo_root: str | Path,
    output_dir: str | Path,
    *,
    block_index: int,
    proposals_per_block: int = EAS_PROPOSALS_PER_BLOCK,
    enforce_production_contract: bool = True,
    base_manifest: Mapping[str, Any] | None = None,
    extension_plan: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Atomically write or checksum-verify one restartable pass-A block."""

    root = Path(repo_root).resolve()
    destination = Path(output_dir).resolve()
    index = int(block_index)
    proposal_count = int(proposals_per_block)
    plan_validation: dict[str, Any] | None = None
    if enforce_production_contract:
        if base_manifest is None or extension_plan is None:
            raise ValueError("production EAS blocks require the v2 manifest and plan")
        plan_validation = validate_eas_extension_plan(
            root,
            extension_plan,
            base_manifest=base_manifest,
        )
        if not 0 <= index < EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT:
            raise ValueError("production EAS block index must be in 0..199")
        if proposal_count != EAS_PROPOSALS_PER_BLOCK:
            raise ValueError("production EAS blocks require exactly 50000 proposals")
    if destination.exists() and (destination / "completion.json").is_file():
        return verify_eas_proposal_block(
            root,
            destination,
            expected_block_index=index,
            expected_proposals=proposal_count,
            enforce_production_contract=enforce_production_contract,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
        )
    if destination.exists() and any(destination.iterdir()):
        raise ValueError(
            "EAS proposal-block destination is nonempty without completion"
        )
    if enforce_production_contract and index < EAS_BASE_PROPOSAL_BLOCK_COUNT:
        raise ValueError(
            "frozen EAS base blocks cannot be regenerated; supply their manifest"
        )

    block = simulate_eas_proposal_block(
        eas_integer_ne_by_generation(root),
        block_index=index,
        proposals_per_block=proposal_count,
    )
    summary = validate_eas_proposal_block(
        block, enforce_production_contract=enforce_production_contract
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.rmdir()
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.stage.", dir=destination.parent)
    )
    try:
        archive = stage / "proposal_block.npz"
        np.savez_compressed(
            archive,
            outcomes=np.asarray(block["outcomes"]),
            birth_generations=np.asarray(block["birth_generations"]),
            log_forward=np.asarray(block["log_forward"]),
            log_reverse=np.asarray(block["log_reverse"]),
            log_mutation=np.asarray(block["log_mutation"]),
            log_weights=np.asarray(block["log_weights"]),
        )
        _atomic_json(stage / "summary.json", summary)
        outputs = {
            name: {
                "size_bytes": (stage / name).stat().st_size,
                "sha256": sha256_file(stage / name),
            }
            for name in ("proposal_block.npz", "summary.json")
        }
        completion = {
            "schema": EAS_PROPOSAL_BLOCK_SCHEMA,
            "status": "complete",
            "created_utc": _utc_now(),
            "block_index": int(block["block_index"]),
            "block_seed": int(block["block_seed"]),
            "proposal_count": int(block["proposal_count"]),
            "present_count": int(block["present_count"]),
            "present_chromosomes": EAS_PRESENT_CHROMOSOMES,
            "max_birth_generation": EAS_BRIDGE_MAX_GENERATION,
            "absorption_draw_generation": EAS_BRIDGE_ABSORPTION_GENERATION,
            "source_binding": _eas_proposal_source_binding(root),
            "summary_sha256": _canonical_sha256(summary),
            "outputs": outputs,
        }
        if enforce_production_contract:
            if (
                plan_validation is None
                or base_manifest is None
                or extension_plan is None
            ):
                raise RuntimeError("validated EAS extension binding is absent")
            completion["extension_bank_binding"] = _eas_extension_block_binding(
                base_manifest, extension_plan
            )
        _atomic_json(stage / "completion.json", completion)
        os.replace(stage, destination)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return verify_eas_proposal_block(
        root,
        destination,
        expected_block_index=index,
        expected_proposals=proposal_count,
        enforce_production_contract=enforce_production_contract,
        base_manifest=base_manifest,
        extension_plan=extension_plan,
    )


def verify_eas_proposal_block(
    repo_root: str | Path,
    block_dir: str | Path,
    *,
    expected_block_index: int,
    expected_proposals: int = EAS_PROPOSALS_PER_BLOCK,
    enforce_production_contract: bool = True,
    base_manifest: Mapping[str, Any] | None = None,
    extension_plan: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Load one pass-A block only after source and checksum validation."""

    root = Path(repo_root).resolve()
    bundle = Path(block_dir).resolve()
    plan_validation: dict[str, Any] | None = None
    if enforce_production_contract:
        if base_manifest is None or extension_plan is None:
            raise ValueError("production EAS blocks require the v2 manifest and plan")
        plan_validation = validate_eas_extension_plan(
            root,
            extension_plan,
            base_manifest=base_manifest,
        )
    completion_path = bundle / "completion.json"
    if bundle.is_symlink() or not bundle.is_dir() or not completion_path.is_file():
        raise ValueError(f"EAS proposal-block bundle is absent: {bundle}")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if (
        completion.get("schema") != EAS_PROPOSAL_BLOCK_SCHEMA
        or completion.get("status") != "complete"
        or int(completion.get("block_index", -1)) != int(expected_block_index)
        or int(completion.get("proposal_count", -1)) != int(expected_proposals)
        or int(completion.get("block_seed", -1))
        != EAS_PROPOSAL_BLOCK_SEED_START + int(expected_block_index)
        or int(completion.get("present_count", -1)) != EAS_SELECTED_PRESENT_COUNT
        or int(completion.get("max_birth_generation", -1)) != EAS_BRIDGE_MAX_GENERATION
        or int(completion.get("absorption_draw_generation", -1))
        != EAS_BRIDGE_ABSORPTION_GENERATION
    ):
        raise ValueError("EAS proposal-block completion contract changed")
    index = int(expected_block_index)
    observed_source_binding = completion.get("source_binding")
    current_source_binding = _eas_proposal_source_binding(root)
    manifest_entry: Mapping[str, Any] | None = None
    source_binding_kind = "current_extension_implementation"
    if enforce_production_contract and index < EAS_BASE_PROPOSAL_BLOCK_COUNT:
        if base_manifest is None:
            raise ValueError("frozen EAS base block requires its v2 manifest")
        manifest_validation = validate_eas_base_bank_manifest(root, base_manifest)
        if observed_source_binding != _frozen_eas_base_source_binding(root):
            raise ValueError("frozen EAS base block source/runtime binding changed")
        manifest_entry = base_manifest["blocks"][index]
        source_binding_kind = "frozen_v1_base_manifest"
    elif observed_source_binding != current_source_binding:
        if base_manifest is None or index >= EAS_BASE_PROPOSAL_BLOCK_COUNT:
            raise ValueError("EAS proposal-block source/runtime binding is stale")
        manifest_validation = validate_eas_base_bank_manifest(root, base_manifest)
        if observed_source_binding != _frozen_eas_base_source_binding(root):
            raise ValueError("frozen EAS base block source/runtime binding changed")
        manifest_entry = base_manifest["blocks"][index]
        source_binding_kind = "frozen_v1_base_manifest"
    else:
        manifest_validation = None
    if enforce_production_contract and index >= EAS_BASE_PROPOSAL_BLOCK_COUNT:
        if base_manifest is None or extension_plan is None:
            raise RuntimeError("validated EAS extension contract is absent")
        if completion.get("extension_bank_binding") != _eas_extension_block_binding(
            base_manifest, extension_plan
        ):
            raise ValueError("EAS extension block is not bound to the fixed v2 plan")
    outputs = completion.get("outputs")
    if not isinstance(outputs, dict) or set(outputs) != {
        "proposal_block.npz",
        "summary.json",
    }:
        raise ValueError("EAS proposal-block output inventory changed")
    observed_inventory = {path.name for path in bundle.iterdir()}
    if observed_inventory != {
        "completion.json",
        "proposal_block.npz",
        "summary.json",
    }:
        raise ValueError("EAS proposal-block directory inventory changed")
    for name, record in outputs.items():
        path = bundle / name
        if (
            path.is_symlink()
            or not path.is_file()
            or path.stat().st_size != int(record.get("size_bytes", -1))
            or sha256_file(path) != record.get("sha256")
        ):
            raise ValueError(f"EAS proposal-block checksum failed: {name}")
    summary = json.loads((bundle / "summary.json").read_text(encoding="utf-8"))
    if completion.get("summary_sha256") != _canonical_sha256(summary):
        raise ValueError("EAS proposal-block summary binding changed")
    if manifest_entry is not None and (
        completion_path.stat().st_size != int(manifest_entry["completion_size_bytes"])
        or sha256_file(completion_path) != manifest_entry["completion_sha256"]
        or outputs["proposal_block.npz"]["size_bytes"]
        != int(manifest_entry["proposal_block_size_bytes"])
        or outputs["proposal_block.npz"]["sha256"]
        != manifest_entry["proposal_block_sha256"]
        or outputs["summary.json"]["size_bytes"]
        != int(manifest_entry["summary_size_bytes"])
        or outputs["summary.json"]["sha256"] != manifest_entry["summary_file_sha256"]
        or completion["summary_sha256"] != manifest_entry["summary_payload_sha256"]
    ):
        raise ValueError("frozen EAS base block differs from its v2 manifest")
    with np.load(bundle / "proposal_block.npz", allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]) for name in archive.files}
    expected_arrays = {
        "outcomes",
        "birth_generations",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weights",
    }
    if set(arrays) != expected_arrays:
        raise ValueError("EAS proposal-block array inventory changed")
    block: dict[str, Any] = {
        "schema": EAS_PROPOSAL_BLOCK_SCHEMA,
        "block_index": int(completion["block_index"]),
        "block_seed": int(completion["block_seed"]),
        "proposal_count": int(completion["proposal_count"]),
        "present_count": int(completion["present_count"]),
        "max_birth_generation": int(completion["max_birth_generation"]),
        **arrays,
        "captured_paths": {},
        "artifact_dir": str(bundle),
        "proposal_block_sha256": outputs["proposal_block.npz"]["sha256"],
        "completion_sha256": sha256_file(completion_path),
        "source_binding_kind": source_binding_kind,
        "base_manifest_payload_sha256": (
            str(base_manifest["manifest_payload_sha256"])
            if enforce_production_contract and base_manifest is not None
            else (
                manifest_validation["manifest_payload_sha256"]
                if manifest_validation is not None
                else ""
            )
        ),
        "extension_plan_payload_sha256": (
            plan_validation["plan_payload_sha256"]
            if plan_validation is not None
            else ""
        ),
    }
    observed_summary = validate_eas_proposal_block(
        block, enforce_production_contract=enforce_production_contract
    )
    if observed_summary != summary:
        raise ValueError("EAS proposal-block summary is stale")
    return block


def load_eas_extended_proposal_bank(
    repo_root: str | Path,
    proposal_dir: str | Path,
    *,
    base_manifest: Mapping[str, Any] | None,
    extension_plan: Mapping[str, Any] | None,
    enforce_production_contract: bool = True,
) -> list[dict[str, Any]]:
    """Load the exact fixed bank, rejecting partial or mixed-provenance roots."""

    root = Path(repo_root).resolve()
    proposal_root = Path(proposal_dir).resolve()
    if not enforce_production_contract:
        raise ValueError("the extended-bank loader is production-contract only")
    if base_manifest is None or extension_plan is None:
        raise ValueError("the extended EAS bank requires its v2 manifest and plan")
    expected_names = {
        f"block_{index:02d}" for index in range(EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT)
    }
    children = list(proposal_root.iterdir()) if proposal_root.is_dir() else []
    if (
        proposal_root.is_symlink()
        or {path.name for path in children} != expected_names
        or any(path.is_symlink() or not path.is_dir() for path in children)
    ):
        raise ValueError(
            "EAS proposal root must contain exactly fixed blocks 0 through 199"
        )
    validate_eas_base_bank_manifest(
        root,
        base_manifest,
        proposal_dir=proposal_root,
    )
    plan_validation = validate_eas_extension_plan(
        root,
        extension_plan,
        base_manifest=base_manifest,
    )
    blocks = [
        verify_eas_proposal_block(
            root,
            proposal_root / f"block_{index:02d}",
            expected_block_index=index,
            expected_proposals=EAS_PROPOSALS_PER_BLOCK,
            enforce_production_contract=True,
            base_manifest=base_manifest,
            extension_plan=extension_plan,
        )
        for index in range(EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT)
    ]
    manifest_hash = str(base_manifest["manifest_payload_sha256"])
    plan_hash = str(plan_validation["plan_payload_sha256"])
    for index, block in enumerate(blocks):
        expected_kind = (
            "frozen_v1_base_manifest"
            if index < EAS_BASE_PROPOSAL_BLOCK_COUNT
            else "current_extension_implementation"
        )
        if (
            block.get("source_binding_kind") != expected_kind
            or block.get("base_manifest_payload_sha256") != manifest_hash
            or block.get("extension_plan_payload_sha256") != plan_hash
        ):
            raise RuntimeError("loaded EAS proposal block lost v2 provenance")
    return blocks


def _categorical_resample_indices(
    normalized_weights: np.ndarray,
    *,
    count: int,
    seed: int,
) -> np.ndarray:
    weights = np.asarray(normalized_weights, dtype=float)
    if (
        weights.ndim != 1
        or not len(weights)
        or np.any(~np.isfinite(weights))
        or np.any(weights < 0)
        or not math.isclose(float(weights.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)
    ):
        raise ValueError("categorical resampling weights are invalid")
    draws = int(count)
    if isinstance(count, bool) or draws != count or draws < 1:
        raise ValueError("categorical resampling count must be positive")
    rng = _bridge_rng(int(seed), 2)
    return rng.choice(len(weights), size=draws, replace=True, p=weights).astype(
        np.int64
    )


def select_eas_global_proposals(
    blocks: Sequence[Mapping[str, Any]],
    *,
    accepted_count: int = ACCEPTED_PER_MODEL,
    enforce_production_contract: bool = True,
    repo_root: str | Path | None = None,
    base_manifest: Mapping[str, Any] | None = None,
    extension_plan: Mapping[str, Any] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Select occurrence units by iid categorical draws from the global bank.

    Sampling is with replacement.  Repeated source particles are legitimate;
    every occurrence receives a distinct CoSi2 genealogy seed and is retained
    in the ledger.  This preserves the iid conditional-resampling semantics
    required by the downstream ``+1`` empirical p-value calculation.
    """

    manifest_validation: dict[str, Any] | None = None
    plan_validation: dict[str, Any] | None = None
    if enforce_production_contract:
        if repo_root is None or base_manifest is None or extension_plan is None:
            raise ValueError(
                "production EAS selection requires repo root, base manifest, and v2 plan"
            )
        root = Path(repo_root).resolve()
        manifest_validation = validate_eas_base_bank_manifest(root, base_manifest)
        plan_validation = validate_eas_extension_plan(
            root,
            extension_plan,
            base_manifest=base_manifest,
        )
    if not blocks:
        raise ValueError("EAS proposal block collection is empty")
    ordered = sorted(blocks, key=lambda value: int(value["block_index"]))
    block_indexes = [int(value["block_index"]) for value in ordered]
    if block_indexes != list(range(len(ordered))):
        raise ValueError("EAS proposal blocks are not consecutive from zero")
    proposal_counts = [int(value["proposal_count"]) for value in ordered]
    if len(set(proposal_counts)) != 1:
        raise ValueError("EAS proposal block sizes differ")
    if enforce_production_contract and (
        len(ordered) != EAS_PROPOSAL_BLOCK_COUNT
        or proposal_counts[0] != EAS_PROPOSALS_PER_BLOCK
        or int(accepted_count) != ACCEPTED_PER_MODEL
    ):
        raise ValueError(
            "EAS global bank differs from the 200x50000 production contract"
        )

    compact_rows: list[pd.DataFrame] = []
    block_bindings: list[dict[str, Any]] = []
    for block in ordered:
        validate_eas_proposal_block(
            block, enforce_production_contract=enforce_production_contract
        )
        binding = {
            "block_index": int(block["block_index"]),
            "block_seed": int(block["block_seed"]),
            "proposal_count": int(block["proposal_count"]),
            "proposal_block_sha256": str(block.get("proposal_block_sha256", "")),
            "completion_sha256": str(block.get("completion_sha256", "")),
        }
        block_bindings.append(binding)
        if enforce_production_contract:
            if manifest_validation is None or plan_validation is None:
                raise RuntimeError("validated EAS bank contract is absent")
            index = int(block["block_index"])
            expected_kind = (
                "frozen_v1_base_manifest"
                if index < EAS_BASE_PROPOSAL_BLOCK_COUNT
                else "current_extension_implementation"
            )
            if (
                block.get("source_binding_kind") != expected_kind
                or block.get("base_manifest_payload_sha256")
                != manifest_validation["manifest_payload_sha256"]
                or block.get("extension_plan_payload_sha256")
                != plan_validation["plan_payload_sha256"]
            ):
                raise RuntimeError(
                    f"EAS proposal block {index} lacks fixed-v2 provenance"
                )
            for key in ("proposal_block_sha256", "completion_sha256"):
                value = binding[key]
                if len(value) != 64 or any(
                    character not in "0123456789abcdef" for character in value
                ):
                    raise RuntimeError(
                        f"EAS proposal block {binding['block_index']} lacks "
                        f"hash-bound {key}"
                    )
        count = int(block["proposal_count"])
        required_arrays = (
            "outcomes",
            "birth_generations",
            "log_forward",
            "log_reverse",
            "log_mutation",
            "log_weights",
        )
        arrays = {name: np.asarray(block[name]) for name in required_arrays}
        if any(value.shape != (count,) for value in arrays.values()):
            raise ValueError("EAS proposal block array geometry changed")
        valid = arrays["outcomes"].astype(str) == "valid"
        if not np.array_equal(np.isfinite(arrays["log_weights"]), valid):
            raise ValueError("EAS proposal finite weights disagree with status")
        indexes = np.flatnonzero(valid)
        compact_rows.append(
            pd.DataFrame(
                {
                    "block_index": int(block["block_index"]),
                    "block_seed": int(block["block_seed"]),
                    "proposal_index": indexes,
                    "global_proposal_index": (
                        int(block["block_index"]) * count + indexes
                    ),
                    "birth_generation": arrays["birth_generations"][indexes],
                    "log_forward": arrays["log_forward"][indexes],
                    "log_reverse": arrays["log_reverse"][indexes],
                    "log_mutation": arrays["log_mutation"][indexes],
                    "log_weight": arrays["log_weights"][indexes],
                    "source_proposal_block_sha256": str(
                        block.get("proposal_block_sha256", "")
                    ),
                    "source_proposal_completion_sha256": str(
                        block.get("completion_sha256", "")
                    ),
                }
            )
        )
    valid_frame = pd.concat(compact_rows, ignore_index=True).sort_values(
        "global_proposal_index", kind="mergesort"
    )
    if valid_frame.empty:
        raise RuntimeError("EAS global proposal bank has no valid weighted paths")
    log_weights = valid_frame["log_weight"].to_numpy(dtype=float)
    log_max = float(np.max(log_weights))
    scaled = np.exp(log_weights - log_max)
    weight_sum_over_max = float(scaled.sum())
    normalized = scaled / weight_sum_over_max
    ess = float(1.0 / np.sum(normalized**2))
    max_normalized = float(np.max(normalized))
    total_proposals = int(sum(proposal_counts))
    cap_count = sum(
        int(np.count_nonzero(np.asarray(block["outcomes"]).astype(str) == "cap"))
        for block in ordered
    )
    quarter_records: list[dict[str, Any]] = []
    all_global = valid_frame["global_proposal_index"].to_numpy(dtype=np.int64)
    all_ages = valid_frame["birth_generation"].to_numpy(dtype=float)
    quarter_size = total_proposals // 4
    for quarter in range(4):
        lower = quarter * quarter_size
        upper = total_proposals if quarter == 3 else (quarter + 1) * quarter_size
        mask = (all_global >= lower) & (all_global < upper)
        if not np.any(mask):
            raise RuntimeError(
                f"EAS proposal-bank quarter {quarter + 1} has no valid paths"
            )
        quarter_log = log_weights[mask]
        quarter_scaled = np.exp(quarter_log - np.max(quarter_log))
        quarter_weights = quarter_scaled / quarter_scaled.sum()
        quantiles = _weighted_quantile(
            all_ages[mask], quarter_weights, (0.05, 0.5, 0.95)
        )
        quarter_records.append(
            {
                "quarter": quarter + 1,
                "proposal_lower_inclusive": lower,
                "proposal_upper_exclusive": upper,
                "valid_proposals": int(mask.sum()),
                "weighted_age_q05": float(quantiles[0]),
                "weighted_age_median": float(quantiles[1]),
                "weighted_age_q95": float(quantiles[2]),
            }
        )
    global_age_quantiles = _weighted_quantile(
        all_ages, normalized, (0.05, 0.25, 0.5, 0.75, 0.95)
    )
    diagnostics = {
        "status": "valid",
        "proposal_blocks": len(ordered),
        "proposals_per_block": proposal_counts[0],
        "total_proposals": total_proposals,
        "valid_proposals": len(valid_frame),
        "global_log_max_weight": log_max,
        "global_weight_sum_over_max": weight_sum_over_max,
        "global_effective_sample_size": ess,
        "global_max_normalized_weight": max_normalized,
        "cap_hit_count": cap_count,
        "cap_hit_fraction": cap_count / total_proposals,
        "cap_hit_one_sided_95pct_upper": _tail_interval_upper(
            cap_count, total_proposals
        ),
        "weighted_age_quantiles_q05_q25_q50_q75_q95": [
            float(value) for value in global_age_quantiles
        ],
        "quarter_bank_stability": quarter_records,
        "proposal_block_bindings": block_bindings,
        "proposal_bank_binding_sha256": _canonical_sha256(block_bindings),
        "categorical_resample_seed": EAS_CATEGORICAL_RESAMPLE_SEED,
        "resampling_with_replacement": True,
        "proposal_bank_contract_id": EAS_PROPOSAL_BANK_CONTRACT_ID,
        "base_manifest_payload_sha256": (
            manifest_validation["manifest_payload_sha256"]
            if manifest_validation is not None
            else ""
        ),
        "extension_plan_payload_sha256": (
            plan_validation["plan_payload_sha256"]
            if plan_validation is not None
            else ""
        ),
        "fixed_target_evaluated_without_adaptive_stopping": bool(
            enforce_production_contract
        ),
        "importance_method": (
            "Slatkin_forward_reverse_transition_ratio_plus_N_birth_mutation_opportunity"
        ),
        "time_reversal_is_importance_corrected_approximation": True,
    }
    quality = evaluate_eas_v2_bank_quality(diagnostics)
    diagnostics["v2_quality_gate"] = quality
    if enforce_production_contract and not quality["quality_gates_passed"]:
        raise RuntimeError(
            "EAS global proposal bank fails ESS/W/max-weight/cap gates: "
            f"ESS={ess:.6g}, W={weight_sum_over_max:.6g}, "
            f"max_weight={max_normalized:.6g}, cap={cap_count}; "
            f"failed={quality['failed_quality_gates']}"
        )
    selected_rows = _categorical_resample_indices(
        normalized,
        count=int(accepted_count),
        seed=EAS_CATEGORICAL_RESAMPLE_SEED,
    )
    selected = valid_frame.iloc[selected_rows].copy().reset_index(drop=True)
    selected.insert(0, "accepted_index", np.arange(len(selected), dtype=np.int64))
    selected.insert(
        1,
        "unit_id",
        [f"{EAS_MODEL_ID}_r{index:03d}" for index in range(len(selected))],
    )
    selected.insert(2, "cosi2_seed", EAS_SEED_START + np.arange(len(selected)))
    selected["normalized_global_weight"] = normalized[selected_rows]
    selected["source_particle_occurrence_index"] = selected.groupby(
        "global_proposal_index", sort=False
    ).cumcount()
    multiplicities = selected["global_proposal_index"].value_counts()
    selected["source_particle_draw_multiplicity"] = selected[
        "global_proposal_index"
    ].map(multiplicities)
    selected["selection_method"] = (
        "seeded_iid_categorical_with_replacement_from_global_weighted_bank"
    )
    diagnostics.update(
        {
            "selected_count": len(selected),
            "selected_unique_global_proposals": int(
                selected["global_proposal_index"].nunique()
            ),
            "selected_duplicate_occurrences": int(
                len(selected) - selected["global_proposal_index"].nunique()
            ),
            "selected_particles_with_multiplicity_gt1": int(
                np.count_nonzero(multiplicities.to_numpy(dtype=np.int64) > 1)
            ),
            "selected_max_particle_multiplicity": int(multiplicities.max()),
        }
    )
    if enforce_production_contract:
        for column in (
            "source_proposal_block_sha256",
            "source_proposal_completion_sha256",
        ):
            values = set(selected[column].astype(str))
            if any(
                len(value) != 64
                or any(character not in "0123456789abcdef" for character in value)
                for value in values
            ):
                raise RuntimeError(f"EAS selected rows lack hash-bound {column}")
    return selected, diagnostics


def replay_eas_global_selection(
    integer_ne: Sequence[int],
    blocks: Sequence[Mapping[str, Any]],
    selected: pd.DataFrame,
) -> dict[str, tuple[int, ...]]:
    """Pass B: rerun selected blocks and reproduce every chosen full path."""

    required = {
        "unit_id",
        "block_index",
        "proposal_index",
        "birth_generation",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weight",
    }
    if (
        not required.issubset(selected.columns)
        or selected["unit_id"].duplicated().any()
    ):
        raise ValueError("EAS selected proposal table is malformed")
    indexed_blocks = {int(block["block_index"]): block for block in blocks}
    paths: dict[str, tuple[int, ...]] = {}
    for block_index, group in selected.groupby("block_index", sort=True):
        index = int(block_index)
        if index not in indexed_blocks:
            raise ValueError(f"selected EAS block is absent: {index}")
        proposal_indices = tuple(group["proposal_index"].astype(int))
        replay = simulate_eas_proposal_block(
            integer_ne,
            block_index=index,
            proposals_per_block=int(indexed_blocks[index]["proposal_count"]),
            capture_indices=proposal_indices,
        )
        original = indexed_blocks[index]
        for key in (
            "outcomes",
            "birth_generations",
            "log_forward",
            "log_reverse",
            "log_mutation",
            "log_weights",
        ):
            original_array = np.asarray(original[key])
            replay_array = np.asarray(replay[key])
            if key == "outcomes":
                equal = np.array_equal(original_array, replay_array)
            else:
                equal = np.array_equal(original_array, replay_array, equal_nan=True)
            if not equal:
                raise RuntimeError(f"EAS block {index} replay changed {key}")
        for row in group.to_dict(orient="records"):
            proposal = int(row["proposal_index"])
            comparisons = {
                "birth_generation": replay["birth_generations"][proposal],
                "log_forward": replay["log_forward"][proposal],
                "log_reverse": replay["log_reverse"][proposal],
                "log_mutation": replay["log_mutation"][proposal],
                "log_weight": replay["log_weights"][proposal],
            }
            for label, observed in comparisons.items():
                expected = row[label]
                if label == "birth_generation":
                    matches = int(observed) == int(expected)
                else:
                    matches = float(observed) == float(expected)
                if not matches:
                    raise RuntimeError(
                        f"selected EAS proposal changed {label}: {row['unit_id']}"
                    )
            path = tuple(replay["captured_paths"][proposal])
            if len(path) != int(row["birth_generation"]) + 1 or path[-1] != 1:
                raise RuntimeError("selected EAS replay path is incomplete")
            recomputed = _reverse_path_log_components(integer_ne, path)
            for label, observed in recomputed.items():
                if not math.isclose(
                    float(observed),
                    float(row[label]),
                    rel_tol=0.0,
                    abs_tol=1e-9,
                ):
                    raise RuntimeError(
                        f"selected EAS full path changed {label}: {row['unit_id']}"
                    )
            paths[str(row["unit_id"])] = path
    if set(paths) != set(selected["unit_id"].astype(str)):
        raise RuntimeError("EAS selected-path replay inventory is incomplete")
    return paths


def eas_selected_path_frames(
    integer_ne: Sequence[int],
    derived_counts_ascending: Sequence[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert one selected global-bank particle to audit and CoSi2 tables."""

    sizes = np.asarray(integer_ne, dtype=np.int64)
    if sizes.shape != (EAS_BRIDGE_ABSORPTION_GENERATION + 1,):
        raise ValueError("integer_ne must cover generations 0..500001")
    counts = np.asarray(derived_counts_ascending, dtype=np.int64)
    if counts.ndim != 1 or not len(counts):
        raise ValueError("selected EAS path must be a nonempty vector")
    birth = len(counts) - 1
    if birth > EAS_BRIDGE_MAX_GENERATION:
        raise ValueError("selected EAS path exceeds the maximum birth generation")
    ne = sizes[: len(counts)]
    chromosomes = 2 * ne
    if (
        int(counts[0]) != EAS_SELECTED_PRESENT_COUNT
        or int(counts[-1]) != 1
        or np.any((counts <= 0) | (counts >= chromosomes))
    ):
        raise ValueError("selected EAS path violates endpoint/birth/segregation")
    audit = pd.DataFrame(
        {
            "generation": np.arange(len(counts), dtype=np.int64),
            "diploid_ne": ne,
            "derived_count": counts,
            "chromosomes": chromosomes,
            "frequency": counts / chromosomes,
        }
    )
    trajectory = pd.DataFrame(
        {
            "sim": 1,
            "gen": audit["generation"].to_numpy(dtype=np.int64)[::-1],
            "selfreq_1": audit["frequency"].to_numpy(dtype=float)[::-1],
        }
    )
    validate_eas_selected_path_frames(audit, trajectory)
    return audit, trajectory


def validate_eas_selected_path_frames(
    audit: pd.DataFrame,
    trajectory: pd.DataFrame,
    *,
    expected_birth_generation: int | None = None,
) -> dict[str, Any]:
    """Fail closed on a globally selected EAS path and load trajectory."""

    required_audit = {
        "generation",
        "diploid_ne",
        "derived_count",
        "chromosomes",
        "frequency",
    }
    required_trajectory = {"sim", "gen", "selfreq_1"}
    if not required_audit.issubset(audit.columns) or not required_trajectory.issubset(
        trajectory.columns
    ):
        raise ValueError("EAS selected-path tables lack required columns")
    generations = pd.to_numeric(audit["generation"], errors="raise").to_numpy(
        dtype=np.int64
    )
    if not len(generations) or not np.array_equal(
        generations, np.arange(len(generations), dtype=np.int64)
    ):
        raise ValueError("EAS selected-path generations are not exact")
    birth = int(generations[-1])
    if birth > EAS_BRIDGE_MAX_GENERATION or (
        expected_birth_generation is not None
        and birth != int(expected_birth_generation)
    ):
        raise ValueError("EAS selected-path birth generation changed")
    ne = pd.to_numeric(audit["diploid_ne"], errors="raise").to_numpy(dtype=np.int64)
    counts = pd.to_numeric(audit["derived_count"], errors="raise").to_numpy(
        dtype=np.int64
    )
    chromosomes = pd.to_numeric(audit["chromosomes"], errors="raise").to_numpy(
        dtype=np.int64
    )
    frequencies = pd.to_numeric(audit["frequency"], errors="raise").to_numpy(
        dtype=float
    )
    if not np.array_equal(chromosomes, 2 * ne):
        raise ValueError("EAS selected-path chromosomes disagree with Ne")
    if not np.allclose(frequencies, counts / chromosomes, rtol=0.0, atol=5e-16):
        raise ValueError("EAS selected-path frequencies disagree with counts")
    if (
        int(ne[0]) != EAS_PRESENT_DIPLOID_NE
        or int(counts[0]) != EAS_SELECTED_PRESENT_COUNT
        or int(counts[-1]) != 1
        or np.any((counts <= 0) | (counts >= chromosomes))
    ):
        raise ValueError("EAS selected path violates endpoint/birth/segregation")
    if set(pd.to_numeric(trajectory["sim"], errors="raise")) != {1}:
        raise ValueError("EAS load trajectory must contain simulation 1 only")
    observed_generations = pd.to_numeric(trajectory["gen"], errors="raise").to_numpy(
        dtype=np.int64
    )
    observed_frequencies = pd.to_numeric(
        trajectory["selfreq_1"], errors="raise"
    ).to_numpy(dtype=float)
    if not np.array_equal(observed_generations, generations[::-1]) or not np.allclose(
        observed_frequencies, frequencies[::-1], rtol=0.0, atol=5e-16
    ):
        raise ValueError("EAS load trajectory differs from the count audit")
    return {
        "status": "valid",
        "bridge_engine": (
            "neutral_wright_fisher_reverse_proposal_with_slatkin_transition_ratio_"
            "importance_weighting_and_global_categorical_resampling"
        ),
        "birth_generation": birth,
        "birth_age_years": birth * EAS_GENERATION_TIME_YEARS,
        "present_count": int(counts[0]),
        "present_chromosomes": int(chromosomes[0]),
        "present_population_af": float(counts[0] / chromosomes[0]),
        "birth_count": int(counts[-1]),
        "max_birth_generation": EAS_BRIDGE_MAX_GENERATION,
        "absorption_draw_generation": EAS_BRIDGE_ABSORPTION_GENERATION,
        "cap_reached": False,
        "importance_corrected_reverse_proposal": True,
        "native_cosi2_age_inference": False,
        "time_reversal_is_approximation": True,
    }


def _eas_selected_particle_record(selected_row: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "accepted_index",
        "unit_id",
        "cosi2_seed",
        "block_index",
        "block_seed",
        "proposal_index",
        "global_proposal_index",
        "birth_generation",
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weight",
        "normalized_global_weight",
        "source_particle_occurrence_index",
        "source_particle_draw_multiplicity",
        "selection_method",
        "source_proposal_block_sha256",
        "source_proposal_completion_sha256",
    }
    if not required.issubset(selected_row):
        raise ValueError("EAS selected-particle row lacks required fields")
    integer_fields = (
        "accepted_index",
        "cosi2_seed",
        "block_index",
        "block_seed",
        "proposal_index",
        "global_proposal_index",
        "birth_generation",
        "source_particle_occurrence_index",
        "source_particle_draw_multiplicity",
    )
    float_fields = (
        "log_forward",
        "log_reverse",
        "log_mutation",
        "log_weight",
        "normalized_global_weight",
    )
    record: dict[str, Any] = {key: int(selected_row[key]) for key in integer_fields}
    record.update({key: float(selected_row[key]) for key in float_fields})
    for key in (
        "unit_id",
        "selection_method",
        "source_proposal_block_sha256",
        "source_proposal_completion_sha256",
    ):
        record[key] = str(selected_row[key])
    index = record["accepted_index"]
    if (
        not 0 <= index < ACCEPTED_PER_MODEL
        or record["unit_id"] != f"{EAS_MODEL_ID}_r{index:03d}"
        or record["cosi2_seed"] != EAS_SEED_START + index
        or record["block_seed"] != EAS_PROPOSAL_BLOCK_SEED_START + record["block_index"]
        or record["global_proposal_index"]
        != record["block_index"] * EAS_PROPOSALS_PER_BLOCK + record["proposal_index"]
        or not 0 <= record["proposal_index"] < EAS_PROPOSALS_PER_BLOCK
        or not 0 <= record["birth_generation"] <= EAS_BRIDGE_MAX_GENERATION
        or record["source_particle_occurrence_index"] < 0
        or record["source_particle_draw_multiplicity"]
        <= record["source_particle_occurrence_index"]
    ):
        raise ValueError("EAS selected-particle index/seed/occurrence contract changed")
    if record["selection_method"] != (
        "seeded_iid_categorical_with_replacement_from_global_weighted_bank"
    ):
        raise ValueError("EAS selected-particle resampling method changed")
    if not all(math.isfinite(record[key]) for key in float_fields) or not (
        0.0 < record["normalized_global_weight"] <= 1.0
    ):
        raise ValueError("EAS selected-particle weights are invalid")
    if not math.isclose(
        record["log_weight"],
        record["log_forward"] - record["log_reverse"] + record["log_mutation"],
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        raise ValueError("EAS selected-particle log weight is inconsistent")
    for key in (
        "source_proposal_block_sha256",
        "source_proposal_completion_sha256",
    ):
        value = record[key]
        if len(value) != 64 or any(
            character not in "0123456789abcdef" for character in value
        ):
            raise ValueError(f"EAS selected-particle {key} is invalid")
    return record


def write_eas_selected_path_bundle(
    repo_root: str | Path,
    output_dir: str | Path,
    *,
    selected_row: Mapping[str, Any],
    derived_counts_ascending: Sequence[int],
    global_diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    """Atomically persist one selected occurrence's trajectory and config."""

    root = Path(repo_root).resolve()
    destination = Path(output_dir).resolve()
    record = _eas_selected_particle_record(selected_row)
    if destination.exists() and (destination / "completion.json").is_file():
        return verify_eas_selected_path_bundle(
            root,
            destination,
            expected_selected_row=record,
            expected_counts=derived_counts_ascending,
            expected_global_diagnostics=global_diagnostics,
        )
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("EAS selected-path destination is nonempty without completion")
    audit, trajectory = eas_selected_path_frames(
        eas_integer_ne_by_generation(root), derived_counts_ascending
    )
    validation = validate_eas_selected_path_frames(
        audit,
        trajectory,
        expected_birth_generation=record["birth_generation"],
    )
    diagnostics = json.loads(_canonical_json(dict(global_diagnostics)))
    selection = {"selected_particle": record, "global_diagnostics": diagnostics}
    parameters = render_eas_replay_parameter_text(
        root, birth_generation=record["birth_generation"]
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.rmdir()
    stage = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.stage.", dir=destination.parent)
    )
    try:
        _atomic_frame(stage / "bridge_counts.tsv", audit)
        _atomic_frame(stage / "trajectory.tsv", trajectory)
        _atomic_text(stage / "parameters.par", parameters)
        _atomic_json(stage / "selection.json", selection)
        outputs = {
            name: {
                "size_bytes": (stage / name).stat().st_size,
                "sha256": sha256_file(stage / name),
            }
            for name in (
                "bridge_counts.tsv",
                "trajectory.tsv",
                "parameters.par",
                "selection.json",
            )
        }
        completion = {
            "schema": EAS_SELECTED_PATH_SCHEMA,
            "status": "complete",
            "created_utc": _utc_now(),
            "unit_id": record["unit_id"],
            "cosi2_seed": record["cosi2_seed"],
            "birth_generation": record["birth_generation"],
            "selected_particle_sha256": _canonical_sha256(record),
            "global_diagnostics_sha256": _canonical_sha256(diagnostics),
            "validation": validation,
            "source_binding": _eas_proposal_source_binding(root),
            "outputs": outputs,
        }
        _atomic_json(stage / "completion.json", completion)
        os.replace(stage, destination)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return verify_eas_selected_path_bundle(
        root,
        destination,
        expected_selected_row=record,
        expected_counts=derived_counts_ascending,
        expected_global_diagnostics=diagnostics,
    )


def verify_eas_selected_path_bundle(
    repo_root: str | Path,
    bundle_dir: str | Path,
    *,
    expected_selected_row: Mapping[str, Any],
    expected_counts: Sequence[int] | None = None,
    expected_global_diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify one occurrence bundle, including exact directory inventory."""

    root = Path(repo_root).resolve()
    bundle = Path(bundle_dir).resolve()
    record = _eas_selected_particle_record(expected_selected_row)
    diagnostics = json.loads(_canonical_json(dict(expected_global_diagnostics)))
    completion_path = bundle / "completion.json"
    if bundle.is_symlink() or not bundle.is_dir() or not completion_path.is_file():
        raise ValueError(f"EAS selected-path bundle is absent: {bundle}")
    expected_inventory = {
        "bridge_counts.tsv",
        "trajectory.tsv",
        "parameters.par",
        "selection.json",
        "completion.json",
    }
    observed_inventory = {path.name for path in bundle.iterdir()}
    if observed_inventory != expected_inventory:
        raise ValueError("EAS selected-path directory inventory changed")
    completion = json.loads(completion_path.read_text(encoding="utf-8"))
    if (
        completion.get("schema") != EAS_SELECTED_PATH_SCHEMA
        or completion.get("status") != "complete"
        or completion.get("unit_id") != record["unit_id"]
        or int(completion.get("cosi2_seed", -1)) != record["cosi2_seed"]
        or int(completion.get("birth_generation", -1)) != record["birth_generation"]
        or completion.get("selected_particle_sha256") != _canonical_sha256(record)
        or completion.get("global_diagnostics_sha256") != _canonical_sha256(diagnostics)
        or completion.get("source_binding") != _eas_proposal_source_binding(root)
    ):
        raise ValueError("EAS selected-path completion contract changed")
    outputs = completion.get("outputs")
    expected_outputs = expected_inventory - {"completion.json"}
    if not isinstance(outputs, dict) or set(outputs) != expected_outputs:
        raise ValueError("EAS selected-path output inventory changed")
    for name, output in outputs.items():
        path = bundle / name
        if (
            path.is_symlink()
            or not path.is_file()
            or path.stat().st_size != int(output.get("size_bytes", -1))
            or sha256_file(path) != output.get("sha256")
        ):
            raise ValueError(f"EAS selected-path checksum failed: {name}")
    selection = json.loads((bundle / "selection.json").read_text(encoding="utf-8"))
    if selection != {"selected_particle": record, "global_diagnostics": diagnostics}:
        raise ValueError("EAS selected-path selection record changed")
    audit = pd.read_csv(bundle / "bridge_counts.tsv", sep="\t")
    trajectory = pd.read_csv(bundle / "trajectory.tsv", sep="\t")
    validation = validate_eas_selected_path_frames(
        audit,
        trajectory,
        expected_birth_generation=record["birth_generation"],
    )
    if completion.get("validation") != validation:
        raise ValueError("EAS selected-path validation record changed")
    if expected_counts is not None:
        regenerated_audit, regenerated_trajectory = eas_selected_path_frames(
            eas_integer_ne_by_generation(root), expected_counts
        )
        pd.testing.assert_frame_equal(
            audit,
            regenerated_audit,
            check_dtype=False,
            rtol=0.0,
            atol=5e-16,
        )
        pd.testing.assert_frame_equal(
            trajectory,
            regenerated_trajectory,
            check_dtype=False,
            rtol=0.0,
            atol=5e-16,
        )
    expected_parameters = render_eas_replay_parameter_text(
        root, birth_generation=record["birth_generation"]
    )
    if (bundle / "parameters.par").read_text(encoding="utf-8") != expected_parameters:
        raise ValueError("EAS selected-path CoSi2 parameters changed")
    return {
        "schema": EAS_SELECTED_PATH_SCHEMA,
        "status": "complete",
        "unit_id": record["unit_id"],
        "cosi2_seed": record["cosi2_seed"],
        "birth_generation": record["birth_generation"],
        "trajectory_path": str((bundle / "trajectory.tsv").resolve()),
        "trajectory_sha256": outputs["trajectory.tsv"]["sha256"],
        "parameters_path": str((bundle / "parameters.par").resolve()),
        "parameters_sha256": outputs["parameters.par"]["sha256"],
        "completion_sha256": sha256_file(completion_path),
    }


def _source_framework_cells(repo_root: Path) -> pd.DataFrame:
    plan = framework.FrameworkPlan(selection_coefficients=(0.01,))
    plan.validate()
    return framework.build_model_cells(framework.load_eas_history(repo_root), plan)


def _replace_sweep_and_seed(
    base_text: str,
    *,
    sweep_line: str,
    header_lines: Sequence[str],
) -> str:
    result: list[str] = []
    sweep_count = 0
    seed_count = 0
    for line in base_text.splitlines():
        if line.startswith("# Generated CoSi2"):
            result.extend(header_lines)
        elif line.startswith("pop_event sweep_mult_standing "):
            result.append(sweep_line)
            sweep_count += 1
        elif line.startswith("random_seed "):
            result.append("random_seed 0")
            seed_count += 1
        else:
            result.append(line)
    if sweep_count != 1 or seed_count != 1:
        raise RuntimeError("base CoSi2 parameter contract changed")
    return "\n".join(result) + "\n"


def render_eas_replay_parameter_text(
    repo_root: str | Path,
    *,
    birth_generation: int,
) -> str:
    """Render the full 10-Mb s=0 EAS config for one loaded bridge."""

    root = Path(repo_root).resolve()
    birth = int(birth_generation)
    if (
        isinstance(birth_generation, bool)
        or birth != birth_generation
        or not (1 <= birth <= EAS_BRIDGE_MAX_GENERATION)
    ):
        raise ValueError("birth_generation must be within the bridge cap")
    source_plan = framework.FrameworkPlan(selection_coefficients=(0.01,))
    cells = _source_framework_cells(root)
    source = cells[cells["sample_population"].astype(str) == "EAS"].iloc[0].to_dict()
    source.update(
        {
            "cell_id": EAS_MODEL_ID,
            "selection_coefficient": 0.0,
            "mutation_birth_generations_ago": birth,
            "selection_onset_generations_ago": birth,
            "present_af_lower": EAS_SELECTED_PRESENT_AF,
            "present_af_upper": EAS_SELECTED_PRESENT_AF,
        }
    )
    base = framework.render_eas_parameter_file(
        framework.rle_eas_history(framework.load_eas_history(root)),
        source,
        source_plan,
    )
    sweep = (
        f'pop_event sweep_mult_standing "{EAS_MODEL_ID}" 1 {birth} 0 '
        f"{framework.FOCAL_RELATIVE_POSITION:g} "
        f"{EAS_SELECTED_PRESENT_AF:.17g} 1 {birth}"
    )
    text = _replace_sweep_and_seed(
        base,
        sweep_line=sweep,
        header_lines=(
            "# Generated CoSi2 EAS conditioned-neutral replay; do not hand-edit.",
            "# s=0; trajectory is supplied with COSI_LOAD_TRAJ.",
            "# Endpoint is the exact accepted selected population count 14704/75274.",
            "# Allele age comes from the frozen reverse-time WF approximation.",
        ),
    )
    lint_conditioned_parameter_text(
        text,
        model_id=EAS_MODEL_ID,
        birth_generation=birth,
        screen=False,
    )
    return text


def render_han_parameter_text(
    repo_root: str | Path,
    *,
    screen: bool,
) -> str:
    """Render the Han s=0 screen or full replay config."""

    root = Path(repo_root).resolve()
    source_plan = framework.FrameworkPlan(
        selection_coefficients=(0.01,),
        sequence_length_bp=(
            SCREEN_SEQUENCE_LENGTH_BP if screen else SEQUENCE_LENGTH_BP
        ),
        sample_haploids=(2 if screen else SAMPLE_HAPLOIDS),
    )
    cells = framework.build_model_cells(framework.load_eas_history(root), source_plan)
    source = cells[cells["sample_population"].astype(str) == "Han"].iloc[0].to_dict()
    source.update(
        {
            "cell_id": HAN_MODEL_ID,
            "selection_coefficient": 0.0,
            "present_af_lower": HAN_PRESENT_LOWER,
            "present_af_upper": HAN_PRESENT_UPPER,
        }
    )
    ledger = framework.build_han_event_ledger(source_plan)
    base = framework.render_han_parameter_file(ledger, source, source_plan)
    sweep = (
        f'pop_event sweep_mult_standing "{HAN_MODEL_ID}" 4 '
        f"{HAN_BIRTH_GENERATION} 0 {framework.FOCAL_RELATIVE_POSITION:g} "
        f"{HAN_PRESENT_TOKEN} 3 {HAN_GATE_GENERATION}"
    )
    text = _replace_sweep_and_seed(
        base,
        sweep_line=sweep,
        header_lines=(
            "# Generated CoSi2 Han conditioned-neutral model; do not hand-edit.",
            "# Neanderthal one-copy birth; s=0 everywhere.",
            "# Native endpoint requires strict CHB segregation only.",
            "# CHB AF 0.011-0.013 at generation 645 is validated post hoc.",
        ),
    )
    lint_conditioned_parameter_text(
        text,
        model_id=HAN_MODEL_ID,
        birth_generation=HAN_BIRTH_GENERATION,
        screen=screen,
    )
    return text


def lint_conditioned_parameter_text(
    text: str,
    *,
    model_id: str,
    birth_generation: int,
    screen: bool,
) -> dict[str, Any]:
    """Lint the exact conditioned-null sweep subset and geometry."""

    directives = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    sweeps = [
        line for line in directives if line.startswith("pop_event sweep_mult_standing ")
    ]
    if len(sweeps) != 1:
        raise ValueError("conditioned config must contain one sweep directive")
    tokens = [
        token
        for token in next(csv.reader([sweeps[0]], delimiter=" ", quotechar='"'))
        if token
    ]
    if tokens[:3] != ["pop_event", "sweep_mult_standing", model_id]:
        raise ValueError("conditioned config model ID changed")
    if int(tokens[4]) != int(birth_generation) or float(tokens[5]) != 0.0:
        raise ValueError("conditioned config birth or s=0 contract changed")
    if model_id == EAS_MODEL_ID:
        expected = {
            "birth_population": 1,
            "endpoint": f"{EAS_SELECTED_PRESENT_AF:.17g}",
            "selection_population": 1,
            "selection_onset": int(birth_generation),
            "length": SEQUENCE_LENGTH_BP,
            "sample": SAMPLE_HAPLOIDS,
        }
    elif model_id == HAN_MODEL_ID:
        expected = {
            "birth_population": 4,
            "endpoint": HAN_PRESENT_TOKEN,
            "selection_population": 3,
            "selection_onset": HAN_GATE_GENERATION,
            "length": SCREEN_SEQUENCE_LENGTH_BP if screen else SEQUENCE_LENGTH_BP,
            "sample": 2 if screen else SAMPLE_HAPLOIDS,
        }
    else:
        raise ValueError(f"unsupported conditioned-null model: {model_id}")
    if (
        int(tokens[3]) != expected["birth_population"]
        or tokens[7] != expected["endpoint"]
        or int(tokens[8]) != expected["selection_population"]
        or int(tokens[9]) != expected["selection_onset"]
    ):
        raise ValueError(
            "conditioned config sweep population/endpoint contract changed"
        )
    lengths = [line for line in directives if line.startswith("length ")]
    samples = [line for line in directives if line.startswith("sample_size ")]
    if lengths != [f"length {expected['length']}"]:
        raise ValueError("conditioned config sequence length changed")
    expected_population = 1 if model_id == EAS_MODEL_ID else 3
    matching_sample = [
        line
        for line in samples
        if line.split()[1:] == [str(expected_population), str(expected["sample"])]
    ]
    if len(matching_sample) != 1 or "random_seed 0" not in directives:
        raise ValueError("conditioned config sample/seed contract changed")
    return {
        "status": "valid",
        "model_id": model_id,
        "selection_coefficient": 0.0,
        "birth_generation": int(birth_generation),
        "endpoint_token": expected["endpoint"],
        "sequence_length_bp": expected["length"],
        "sample_haploids": expected["sample"],
        "sample_af_conditioned": False,
    }


def validate_han_trajectory(
    trajectory_path: str | Path,
    *,
    require_gate: bool = True,
) -> dict[str, Any]:
    """Validate birth, g645 CHB introgression gate, and present segregation."""

    path = Path(trajectory_path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Han trajectory is absent: {path}")
    frame = pd.read_csv(path, sep="\t")
    required = {
        "sim",
        "gen",
        *(f"selfreq_{identifier}" for identifier in range(1, 6)),
        *(f"popsize_{identifier}" for identifier in range(1, 6)),
    }
    if not required.issubset(frame.columns):
        raise ValueError("Han trajectory lacks required population columns")
    if set(frame["sim"].astype(str)) not in ({"1"}, {"1.0"}):
        raise ValueError("Han trajectory must contain simulation 1 only")
    generations = pd.to_numeric(frame["gen"], errors="raise").to_numpy(dtype=float)
    expected = np.arange(HAN_BIRTH_GENERATION, -1, -1, dtype=float)
    if not np.array_equal(generations, expected):
        raise ValueError("Han trajectory does not contain every birth-to-present row")
    frequency_columns = [f"selfreq_{identifier}" for identifier in range(1, 6)]
    size_columns = [f"popsize_{identifier}" for identifier in range(1, 6)]
    frequencies = frame[frequency_columns].apply(pd.to_numeric, errors="raise")
    sizes = frame[size_columns].apply(pd.to_numeric, errors="raise")
    if (
        not np.isfinite(frequencies.to_numpy(dtype=float)).all()
        or ((frequencies < 0) | (frequencies > 1)).any(axis=None)
        or not np.isfinite(sizes.to_numpy(dtype=float)).all()
        or (sizes <= 0).any(axis=None)
    ):
        raise ValueError("Han trajectory has invalid frequency or size")
    birth_frequency = float(frequencies.iloc[0]["selfreq_4"])
    birth_ne = float(sizes.iloc[0]["popsize_4"])
    if not math.isclose(
        birth_frequency,
        1.0 / (2.0 * birth_ne),
        rel_tol=1e-9,
        abs_tol=1e-15,
    ):
        raise ValueError("Han allele does not begin from one Neanderthal copy")
    if any(
        float(frequencies.iloc[0][column]) != 0.0
        for column in frequency_columns
        if column != "selfreq_4"
    ):
        raise ValueError("Han allele is present outside Neanderthal at birth")
    endpoint = float(frequencies.iloc[-1]["selfreq_3"])
    if not 0.0 < endpoint < 1.0:
        raise ValueError("Han allele is not strictly segregating in CHB at present")
    boundary = frame.loc[frame["gen"] == HAN_GATE_GENERATION]
    if len(boundary) != 1:
        raise ValueError("Han trajectory lacks exactly one generation-645 row")
    boundary_af = float(boundary.iloc[0]["selfreq_3"])
    gate_passed = HAN_GATE_LOWER <= boundary_af <= HAN_GATE_UPPER
    if require_gate and not gate_passed:
        raise ValueError("Han trajectory fails the generation-645 CHB gate")
    return {
        "status": "valid",
        "model_id": HAN_MODEL_ID,
        "birth_generation": HAN_BIRTH_GENERATION,
        "birth_population": "Neanderthal",
        "birth_population_frequency": birth_frequency,
        "han_generation_645_chb_af": boundary_af,
        "han_generation_645_gate_lower": HAN_GATE_LOWER,
        "han_generation_645_gate_upper": HAN_GATE_UPPER,
        "han_generation_645_gate_passed": gate_passed,
        "present_chb_af": endpoint,
        "present_strictly_segregating": True,
        "final_af_conditioning": "strict_segregation_only",
        "sample_af_conditioned": False,
    }


def build_eas_unit_ledger() -> pd.DataFrame:
    """Return the fixed 100-unit EAS bridge/replay ledger."""

    rows = [
        {
            "plan_order": index,
            "unit_id": f"{EAS_MODEL_ID}_r{index:03d}",
            "model_id": EAS_MODEL_ID,
            "replicate_index": index,
            "seed": EAS_SEED_START + index,
            "selection_coefficient": 0.0,
            "conditioning": "exact_present_population_count_reverse_bridge",
            "present_population_count": EAS_SELECTED_PRESENT_COUNT,
            "present_population_chromosomes": EAS_PRESENT_CHROMOSOMES,
            "present_population_af": EAS_SELECTED_PRESENT_AF,
            "sample_af_conditioned": False,
            "bridge_cap_generation": EAS_BRIDGE_MAX_GENERATION,
        }
        for index in range(ACCEPTED_PER_MODEL)
    ]
    return pd.DataFrame(rows)


def build_han_candidate_ledger(
    *,
    existing_run_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Return all 5,000 predeclared Han candidates in deterministic order.

    The first 100 rows represent the complete old survival bank in its original
    order.  They are not prefiltered.  Reuse is allowed only after the caller
    binds and verifies each source completion and output checksum.
    """

    rows: list[dict[str, Any]] = []
    existing_root = Path(existing_run_dir).resolve() if existing_run_dir else None
    for index in range(HAN_REUSED_CANDIDATE_COUNT):
        unit_id = f"{natural_neutral.HAN_MODEL_ID}_r{index:03d}"
        unit_dir = existing_root / "units" / unit_id if existing_root else None
        rows.append(
            {
                "candidate_order": index,
                "candidate_id": f"han_candidate_{index:04d}",
                "seed": natural_neutral.HAN_SEED_START + index,
                "source_kind": "existing_complete_natural_neutral_candidate",
                "source_unit_id": unit_id,
                "source_unit_dir": str(unit_dir) if unit_dir else None,
                "selection_coefficient": 0.0,
                "birth_generation": HAN_BIRTH_GENERATION,
                "gate_generation": HAN_GATE_GENERATION,
                "gate_lower": HAN_GATE_LOWER,
                "gate_upper": HAN_GATE_UPPER,
                "present_condition": "strictly_segregating",
            }
        )
    for new_index in range(HAN_NEW_CANDIDATE_COUNT):
        order = HAN_REUSED_CANDIDATE_COUNT + new_index
        rows.append(
            {
                "candidate_order": order,
                "candidate_id": f"han_candidate_{order:04d}",
                "seed": HAN_NEW_SEED_START + new_index,
                "source_kind": "new_screen_then_replay",
                "source_unit_id": None,
                "source_unit_dir": None,
                "selection_coefficient": 0.0,
                "birth_generation": HAN_BIRTH_GENERATION,
                "gate_generation": HAN_GATE_GENERATION,
                "gate_lower": HAN_GATE_LOWER,
                "gate_upper": HAN_GATE_UPPER,
                "present_condition": "strictly_segregating",
            }
        )
    frame = pd.DataFrame(rows)
    if (
        len(frame) != HAN_CANDIDATE_COUNT
        or frame["candidate_id"].duplicated().any()
        or frame["seed"].duplicated().any()
        or not np.array_equal(frame["candidate_order"], np.arange(HAN_CANDIDATE_COUNT))
    ):
        raise RuntimeError("Han candidate ledger is not deterministic and unique")
    return frame


def select_han_acceptances(
    candidate_ledger: pd.DataFrame,
    screen_results: pd.DataFrame,
) -> pd.DataFrame:
    """Freeze the first 100 joint passes from the complete candidate order."""

    required_results = {
        "candidate_id",
        "trajectory_sha256",
        "present_chb_af",
        "han_generation_645_chb_af",
        "gate_passed",
        "strictly_segregating",
        "screen_completion_sha256",
    }
    if not required_results.issubset(screen_results.columns):
        raise ValueError("Han screen results lack required provenance columns")
    if screen_results["candidate_id"].duplicated().any():
        raise ValueError("Han screen results contain duplicate candidates")
    merged = candidate_ledger.merge(
        screen_results,
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )

    def strict_flags(values: pd.Series, label: str) -> pd.Series:
        normalized: list[bool] = []
        for value in values.tolist():
            if value is None or (isinstance(value, float) and math.isnan(value)):
                normalized.append(False)
            elif isinstance(value, (bool, np.bool_)):
                normalized.append(bool(value))
            elif isinstance(value, str) and value.casefold() in {"true", "false"}:
                normalized.append(value.casefold() == "true")
            else:
                raise ValueError(f"Han {label} contains a non-boolean value: {value!r}")
        return pd.Series(normalized, index=values.index, dtype=bool)

    gate_passed = strict_flags(merged["gate_passed"], "gate_passed")
    segregating = strict_flags(merged["strictly_segregating"], "strictly_segregating")
    passed = merged[gate_passed & segregating].sort_values(
        "candidate_order", kind="mergesort"
    )
    if len(passed) < ACCEPTED_PER_MODEL:
        raise RuntimeError(
            f"Han candidate bank has only {len(passed)} joint passes; need 100"
        )
    accepted = passed.iloc[:ACCEPTED_PER_MODEL].copy()
    accepted.insert(0, "accepted_index", np.arange(ACCEPTED_PER_MODEL))
    accepted.insert(
        1,
        "unit_id",
        [f"{HAN_MODEL_ID}_r{index:03d}" for index in range(ACCEPTED_PER_MODEL)],
    )
    accepted["selection_rule"] = (
        "first_100_candidate_order_joint_g645_gate_and_present_segregation"
    )
    return accepted.reset_index(drop=True)


def _tail_interval_upper(cap_count: int, attempts: int, alpha: float = 0.05) -> float:
    if attempts <= 0 or not 0 <= cap_count <= attempts:
        raise ValueError("tail audit counts are invalid")
    if cap_count == attempts:
        return 1.0
    # One-sided Clopper--Pearson upper bound.
    from scipy.stats import beta

    return float(beta.ppf(1.0 - alpha, cap_count + 1, attempts - cap_count))


def build_plan_tables(
    *,
    existing_han_run_dir: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    ConditionedNullPlan().validate()
    return {
        "eas_units": build_eas_unit_ledger(),
        "han_candidates": build_han_candidate_ledger(
            existing_run_dir=existing_han_run_dir
        ),
    }


__all__ = [
    "ACCEPTED_PER_MODEL",
    "ConditionedNullPlan",
    "EasProposalBankContract",
    "EAS_BASE_PROPOSAL_BLOCK_COUNT",
    "EAS_BASE_TOTAL_PROPOSALS",
    "EAS_BRIDGE_ABSORPTION_GENERATION",
    "EAS_BRIDGE_MAX_GENERATION",
    "EAS_CATEGORICAL_RESAMPLE_SEED",
    "EAS_MAX_NORMALIZED_WEIGHT",
    "EAS_MIN_GLOBAL_ESS",
    "EAS_MIN_WEIGHT_SUM_OVER_MAX",
    "EAS_MODEL_ID",
    "EAS_PROPOSAL_BLOCK_COUNT",
    "EAS_PROPOSAL_BANK_CONTRACT_ID",
    "EAS_PROPOSAL_BLOCK_SEED_START",
    "EAS_PROPOSAL_KERNEL_ID",
    "EAS_PROPOSALS_PER_BLOCK",
    "EAS_SEED_START",
    "EAS_SELECTED_PRESENT_AF",
    "EAS_SELECTED_PRESENT_COUNT",
    "EAS_TOTAL_PROPOSALS",
    "EAS_V2_EXTENSION_PROPOSAL_BLOCK_COUNT",
    "EAS_V2_TARGET_PROPOSAL_BLOCK_COUNT",
    "EAS_V2_TOTAL_PROPOSALS",
    "HAN_CANDIDATE_COUNT",
    "HAN_BIRTH_GENERATION",
    "HAN_GATE_GENERATION",
    "HAN_GATE_LOWER",
    "HAN_GATE_UPPER",
    "HAN_MODEL_ID",
    "attest_eas_proposal_kernel_equivalence",
    "build_eas_base_bank_manifest",
    "build_eas_extension_plan",
    "build_eas_unit_ledger",
    "build_han_candidate_ledger",
    "build_plan_tables",
    "eas_integer_ne_by_generation",
    "eas_selected_path_frames",
    "evaluate_eas_v2_bank_quality",
    "lint_conditioned_parameter_text",
    "load_eas_extended_proposal_bank",
    "replay_eas_global_selection",
    "render_eas_replay_parameter_text",
    "render_han_parameter_text",
    "select_eas_global_proposals",
    "select_han_acceptances",
    "simulate_eas_proposal_block",
    "validate_eas_proposal_block",
    "validate_eas_base_bank_manifest",
    "validate_eas_extension_plan",
    "validate_eas_selected_path_frames",
    "validate_han_trajectory",
    "verify_eas_proposal_block",
    "verify_eas_selected_path_bundle",
    "write_eas_proposal_block",
    "write_eas_selected_path_bundle",
]
