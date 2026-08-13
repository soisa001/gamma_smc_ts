"""Demographic and selected-site building blocks for the EAS sweep study.

This module deliberately stops below orchestration and simulation execution.  It
validates the tracked PHLASH artifact, constructs demographic models, and emits
stdpopsim 0.3.0 extended events plus JSON-serializable records describing those
objects.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import msprime
import numpy as np
import stdpopsim


PHLASH_MVN_SCHEMA = "phlash.aou.log-ne-mvn/v1"
STDPOPSIM_VERSION = "0.3.0"

TMRCA_THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
GENERATION_TIME_YEARS = 25.0
TMRCA_THRESHOLDS_GENERATIONS = tuple(
    years / GENERATION_TIME_YEARS for years in TMRCA_THRESHOLDS_YEARS
)
SEQUENCE_LENGTH_BP = 10_000_000
FOCAL_POSITION_BP = 5_000_000
MUTATION_RATE = 1.25e-8
RECOMBINATION_RATE = 1.0e-8
SELECTION_COEFFICIENTS = (0.01, 0.005, 0.001)
DOMINANCE_COEFFICIENT = 0.5
DE_NOVO_TARGET_FREQUENCY = 0.50
DE_NOVO_DEFAULT_LOWER_FREQUENCY = 0.45
DE_NOVO_DEFAULT_UPPER_FREQUENCY = 0.55
INTROGRESSION_TARGET_FREQUENCIES = tuple(i / 10 for i in range(1, 10))

EAS_POPULATION = "EAS"
FOCAL_SITE_ID = "eas_selected_site"
ANCIENT_EURASIA_MODEL_ID = "AncientEurasia_9K19"
ARCHAIC_POPULATION = "Neanderthal"
INTROGRESSION_RECIPIENT_POPULATION = "Loschbour"
INTROGRESSION_TARGET_POPULATION = "Han"
INTROGRESSION_PULSE_GENERATIONS = 2_272.0
HAN_SPLIT_GENERATIONS = 2_016.0
HUMAN_NEANDERTHAL_SPLIT_GENERATIONS = 27_840.0
PRE_PULSE_SOURCE_CHECK_GENERATIONS = INTROGRESSION_PULSE_GENERATIONS + 1.0
DEFAULT_ARCHAIC_ORIGIN_GENERATIONS = 10_000.0
# A fixed source allele is supported by passing 1.0, but is an expensive
# rejection target.  The executable default only requires the allele to be
# established in the small Neanderthal source before the pulse.
DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY = 0.10

_QUANTILES = (
    ("q025", 0.025),
    ("median", 0.5),
    ("q975", 0.975),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise_sha256(value: str) -> str:
    value = str(value).strip().lower()
    if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise ValueError("expected_sha256 must be a 64-character hexadecimal digest")
    return value


def _scalar_string(array: np.ndarray, *, field: str) -> str:
    array = np.asarray(array)
    if array.shape != () and array.size != 1:
        raise ValueError(f"{field} must be a scalar string")
    value = array.item() if array.shape == () else array.reshape(-1)[0]
    if not isinstance(value, (str, np.str_)):
        raise ValueError(f"{field} must be a scalar string")
    return str(value)


def _float_list(values: np.ndarray) -> list[float]:
    return [float(value) for value in np.asarray(values).reshape(-1)]


@dataclass(frozen=True)
class PhlashEasArtifact:
    """Validated EAS PHLASH bootstrap artifact and its pointwise quantiles."""

    path: Path
    expected_sha256: str
    actual_sha256: str
    schema: str
    population: str
    time_generations: np.ndarray
    mean_log_ne: np.ndarray
    covariance_factor: np.ndarray
    bootstrap_ne: np.ndarray
    jitter: float
    q025_ne: np.ndarray
    median_ne: np.ndarray
    q975_ne: np.ndarray

    def quantile_trajectories(self) -> dict[str, np.ndarray]:
        return {
            "q025": self.q025_ne.copy(),
            "median": self.median_ne.copy(),
            "q975": self.q975_ne.copy(),
        }

    def provenance_record(self) -> dict[str, Any]:
        """Return a JSON-serializable description of the validated input."""

        return {
            "schema": "gamma-smc.eas-phlash-source/v1",
            "artifact_schema": self.schema,
            "population": self.population,
            "source_path": str(self.path),
            "source_sha256_expected": self.expected_sha256,
            "source_sha256_actual": self.actual_sha256,
            "source_sha256_verified": self.expected_sha256 == self.actual_sha256,
            "time_units": "generations",
            "time_generations": _float_list(self.time_generations),
            "number_of_time_points": int(self.time_generations.size),
            "number_of_bootstrap_fits": int(self.bootstrap_ne.shape[0]),
            "covariance_factor_shape": [int(x) for x in self.covariance_factor.shape],
            "jitter_log_ne_sd": float(self.jitter),
            "pointwise_quantiles": [probability for _, probability in _QUANTILES],
        }


def load_phlash_eas_npz(
    path: str | Path,
    *,
    expected_sha256: str,
) -> PhlashEasArtifact:
    """Load and validate a tracked ``phlash.aou.log-ne-mvn/v1`` EAS artifact.

    The caller must supply the tracked file digest.  This makes artifact identity
    an explicit run input rather than an assumption based on a local filename.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    expected = _normalise_sha256(expected_sha256)
    actual = _sha256(source)
    if actual != expected:
        raise ValueError(
            f"PHLASH EAS artifact SHA-256 mismatch: expected {expected}, got {actual}"
        )

    required = {
        "schema",
        "population",
        "time",
        "mean_log_ne",
        "covariance_factor",
        "bootstrap_ne",
    }
    try:
        with np.load(source, allow_pickle=False) as data:
            missing = required.difference(data.files)
            if missing:
                raise ValueError(
                    f"PHLASH EAS artifact is missing fields: {sorted(missing)}"
                )
            schema = _scalar_string(np.asarray(data["schema"]), field="schema")
            population = _scalar_string(
                np.asarray(data["population"]), field="population"
            )
            time = np.asarray(data["time"], dtype=np.float64)
            mean_log_ne = np.asarray(data["mean_log_ne"], dtype=np.float64)
            covariance_factor = np.asarray(data["covariance_factor"], dtype=np.float64)
            bootstrap_ne = np.asarray(data["bootstrap_ne"], dtype=np.float64)
            jitter_array = (
                np.asarray(data["jitter"], dtype=np.float64)
                if "jitter" in data
                else np.asarray(0.0)
            )
    except ValueError:
        raise
    except Exception as error:
        raise ValueError(
            f"could not read PHLASH EAS artifact {source}: {error}"
        ) from error

    if schema != PHLASH_MVN_SCHEMA:
        raise ValueError(
            f"PHLASH artifact schema is {schema!r}; expected {PHLASH_MVN_SCHEMA!r}"
        )
    if population != EAS_POPULATION:
        raise ValueError(
            f"PHLASH artifact population is {population!r}; expected {EAS_POPULATION!r}"
        )
    if time.ndim != 1 or time.size < 2:
        raise ValueError("time must be a one-dimensional array of length at least two")
    if not np.isfinite(time).all() or np.any(time <= 0):
        raise ValueError("time must contain only positive finite generations")
    if np.any(np.diff(time) <= 0):
        raise ValueError("time must be strictly increasing")
    if not math.isclose(float(time[0]), 100.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(
            "the tracked PHLASH EAS time grid must begin at 100 generations"
        )
    if mean_log_ne.shape != time.shape or not np.isfinite(mean_log_ne).all():
        raise ValueError("mean_log_ne must be finite and have the same shape as time")
    mean_ne = np.exp(mean_log_ne)
    if not np.isfinite(mean_ne).all() or np.any(mean_ne <= 0):
        raise ValueError("exp(mean_log_ne) must be positive and finite")
    if (
        covariance_factor.ndim != 2
        or covariance_factor.shape[0] < 1
        or covariance_factor.shape[1] != time.size
        or not np.isfinite(covariance_factor).all()
    ):
        raise ValueError(
            "covariance_factor must be finite with shape (positive rank, time points)"
        )
    if (
        bootstrap_ne.ndim != 2
        or bootstrap_ne.shape[0] < 2
        or bootstrap_ne.shape[1] != time.size
        or not np.isfinite(bootstrap_ne).all()
        or np.any(bootstrap_ne <= 0)
    ):
        raise ValueError(
            "bootstrap_ne must be positive and finite with shape "
            "(at least two fits, time points)"
        )
    if jitter_array.shape != () and jitter_array.size != 1:
        raise ValueError("jitter must be scalar")
    jitter = float(jitter_array.item())
    if not math.isfinite(jitter) or jitter < 0:
        raise ValueError("jitter must be nonnegative and finite")

    quantiles = np.quantile(
        bootstrap_ne,
        [probability for _, probability in _QUANTILES],
        axis=0,
        method="linear",
    )
    return PhlashEasArtifact(
        path=source,
        expected_sha256=expected,
        actual_sha256=actual,
        schema=schema,
        population=population,
        time_generations=time,
        mean_log_ne=mean_log_ne,
        covariance_factor=covariance_factor,
        bootstrap_ne=bootstrap_ne,
        jitter=jitter,
        q025_ne=quantiles[0],
        median_ne=quantiles[1],
        q975_ne=quantiles[2],
    )


@dataclass(frozen=True)
class EasDemographyModel:
    """One pointwise PHLASH quantile as msprime and stdpopsim models."""

    label: str
    quantile: float
    time_generations: np.ndarray
    ne: np.ndarray
    msprime_demography: msprime.Demography
    stdpopsim_model: stdpopsim.DemographicModel
    record: dict[str, Any]

    def provenance_record(self) -> dict[str, Any]:
        return copy.deepcopy(self.record)


def _require_stdpopsim_030() -> None:
    if stdpopsim.__version__ != STDPOPSIM_VERSION:
        raise RuntimeError(
            f"this study requires stdpopsim {STDPOPSIM_VERSION}; "
            f"found {stdpopsim.__version__}"
        )


def _piecewise_msprime_demography(
    time_generations: np.ndarray,
    ne: np.ndarray,
    *,
    population: str,
) -> msprime.Demography:
    demography = msprime.Demography()
    demography.add_population(name=population, initial_size=float(ne[0]))
    for time, size in zip(time_generations[1:], ne[1:]):
        demography.add_population_parameters_change(
            time=float(time),
            initial_size=float(size),
            population=population,
        )
    demography.sort_events()
    demography.debug()
    return demography


def build_eas_demography_models(
    artifact: PhlashEasArtifact,
) -> dict[str, EasDemographyModel]:
    """Construct q2.5, median, and q97.5 single-population demographic models.

    PHLASH's first supported point is 100 generations.  Each returned trajectory
    therefore contains an explicit time-zero value copied from the 100-generation
    estimate, and its msprime model contains the corresponding (redundant but
    provenance-visible) change event at 100 generations.
    """

    _require_stdpopsim_030()
    models: dict[str, EasDemographyModel] = {}
    trajectories = artifact.quantile_trajectories()
    for label, probability in _QUANTILES:
        source_ne = trajectories[label]
        time = np.concatenate(([0.0], artifact.time_generations.astype(np.float64)))
        ne = np.concatenate(([float(source_ne[0])], source_ne.astype(np.float64)))
        msprime_model = _piecewise_msprime_demography(
            time, ne, population=EAS_POPULATION
        )
        stdpopsim_model = stdpopsim.DemographicModel(
            id=f"PhlashEAS{label.capitalize()}",
            description=f"PHLASH EAS pointwise {probability:g} quantile",
            long_description=(
                "Single-population piecewise-constant EAS history constructed "
                "from pointwise quantiles of tracked PHLASH bootstrap fits."
            ),
            generation_time=GENERATION_TIME_YEARS,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            model=msprime_model,
        )
        record = {
            "schema": "gamma-smc.eas-demography/v1",
            "population": EAS_POPULATION,
            "trajectory_label": label,
            "pointwise_bootstrap_quantile": probability,
            "time_units": "generations",
            "time_generations": _float_list(time),
            "effective_population_size": _float_list(ne),
            "generation_time_years": GENERATION_TIME_YEARS,
            "mutation_rate_per_bp_per_generation": MUTATION_RATE,
            "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
            "presentward_extrapolation": {
                "method": "constant_first_supported_ne",
                "from_generations_ago": float(artifact.time_generations[0]),
                "to_generations_ago": 0.0,
                "effective_population_size": float(source_ne[0]),
            },
            "ancient_extrapolation": {
                "method": "constant_last_supported_ne",
                "from_generations_ago": float(artifact.time_generations[-1]),
                "effective_population_size": float(source_ne[-1]),
            },
            "source": artifact.provenance_record(),
        }
        json.dumps(record)
        models[label] = EasDemographyModel(
            label=label,
            quantile=probability,
            time_generations=time,
            ne=ne,
            msprime_demography=msprime_model,
            stdpopsim_model=stdpopsim_model,
            record=record,
        )
    return models


@dataclass(frozen=True)
class OriginAgeEstimate:
    """Approximate de novo origin time for one demographic trajectory and s."""

    trajectory_label: str
    selection_coefficient: float
    dominance_coefficient: float
    target_frequency: float
    age_generations: float
    age_years: float
    ne_at_origin: float
    initial_frequency: float
    slim_scaling_factor: float
    fixed_point_residual_generations: float
    grid_step_generations: float
    grid_max_generations: float
    boundary_limited: bool

    def to_record(self) -> dict[str, Any]:
        record = {
            "schema": "gamma-smc.de-novo-origin-age/v2",
            "trajectory_label": self.trajectory_label,
            "selection_coefficient": self.selection_coefficient,
            "dominance_coefficient": self.dominance_coefficient,
            "target_frequency": self.target_frequency,
            "target_frequency_is_placeholder": True,
            "age_generations": self.age_generations,
            "age_years": self.age_years,
            "ne_at_origin": self.ne_at_origin,
            "de_novo_initial_frequency": self.initial_frequency,
            "slim_scaling_factor": self.slim_scaling_factor,
            "de_novo_initial_frequency_semantics": (
                "one copy in the Q-scaled diploid population, approximated as "
                "Q/(2*Ne) before integral population-size rounding"
            ),
            "fixed_point_residual_generations": self.fixed_point_residual_generations,
            "grid_step_generations": self.grid_step_generations,
            "grid_max_generations": self.grid_max_generations,
            "boundary_limited": self.boundary_limited,
            "method": "deterministic_logit_fixed_frequency_grid_approximation",
            "equation": ("t ~= [logit(f) - logit(Q/(2*Ne(t)))] / (h*s)"),
            "is_exact": False,
            "caveat": (
                "Placement approximation only; it ignores drift and stochastic "
                "trajectory conditioning. Realized frequency must be checked in SLiM."
            ),
        }
        json.dumps(record)
        return record


def _logit(value: np.ndarray | float) -> np.ndarray | float:
    return np.log(np.asarray(value) / (1.0 - np.asarray(value)))


def estimate_de_novo_origin_age(
    time_generations: Sequence[float],
    ne: Sequence[float],
    *,
    selection_coefficient: float,
    trajectory_label: str,
    target_frequency: float = DE_NOVO_TARGET_FREQUENCY,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    slim_scaling_factor: float = 1.0,
    grid_step_generations: float = 1.0,
    grid_max_generations: float | None = None,
) -> OriginAgeEstimate:
    """Estimate an origin time by a deterministic logistic fixed-point grid.

    Under stdpopsim's SLiM fitness convention, heterozygotes have fitness
    ``1 + h*s``.  For a rare additive allele, the deterministic logit growth
    rate is therefore approximately ``h*s`` rather than ``s``.  stdpopsim's
    ``Q`` rescaling changes a diploid population from ``Ne`` to approximately
    ``Ne/Q`` while :class:`stdpopsim.DrawMutation` still introduces one copy.
    The corresponding starting frequency is therefore approximated as
    ``Q/(2*Ne)``.  In original-generation units, the reciprocal scaling of
    elapsed ticks and selection strength leaves the logit-growth denominator
    as ``h*s``.
    """

    times = np.asarray(time_generations, dtype=np.float64)
    sizes = np.asarray(ne, dtype=np.float64)
    if times.ndim != 1 or sizes.shape != times.shape or times.size < 2:
        raise ValueError("time_generations and ne must be equal-length vectors")
    if times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("time_generations must begin at zero and increase")
    if not np.isfinite(sizes).all() or np.any(sizes <= 0.5):
        raise ValueError("Ne must be finite and greater than 0.5")
    if not math.isfinite(selection_coefficient) or selection_coefficient <= 0:
        raise ValueError("selection_coefficient must be positive and finite")
    if not math.isfinite(dominance_coefficient) or dominance_coefficient <= 0:
        raise ValueError("dominance_coefficient must be positive and finite")
    if not math.isfinite(slim_scaling_factor) or slim_scaling_factor <= 0:
        raise ValueError("slim_scaling_factor must be positive and finite")
    if not 0 < target_frequency < 1:
        raise ValueError("target_frequency must be in (0, 1)")
    if not math.isfinite(grid_step_generations) or grid_step_generations <= 0:
        raise ValueError("grid_step_generations must be positive and finite")
    if grid_max_generations is None:
        grid_max_generations = float(times[-1])
    if not math.isfinite(grid_max_generations) or grid_max_generations <= 0:
        raise ValueError("grid_max_generations must be positive and finite")

    ages = np.arange(
        grid_step_generations,
        grid_max_generations + 0.5 * grid_step_generations,
        grid_step_generations,
        dtype=np.float64,
    )
    if ages.size == 0:
        raise ValueError("origin-age search grid is empty")
    epoch = np.searchsorted(times, ages, side="right") - 1
    epoch = np.clip(epoch, 0, sizes.size - 1)
    ne_at_age = sizes[epoch]
    initial_frequency = float(slim_scaling_factor) / (2.0 * ne_at_age)
    if np.any(initial_frequency >= 1.0):
        raise ValueError(
            "slim_scaling_factor leaves fewer than one diploid copy at an "
            "origin-age grid point"
        )
    duration = (_logit(target_frequency) - _logit(initial_frequency)) / (
        dominance_coefficient * selection_coefficient
    )
    signed_residual = np.asarray(duration - ages, dtype=np.float64)
    best_index = int(np.argmin(np.abs(signed_residual)))
    best_age = float(ages[best_index])
    boundary_limited = best_index in (0, ages.size - 1)
    return OriginAgeEstimate(
        trajectory_label=str(trajectory_label),
        selection_coefficient=float(selection_coefficient),
        dominance_coefficient=float(dominance_coefficient),
        target_frequency=float(target_frequency),
        age_generations=best_age,
        age_years=best_age * GENERATION_TIME_YEARS,
        ne_at_origin=float(ne_at_age[best_index]),
        initial_frequency=float(initial_frequency[best_index]),
        slim_scaling_factor=float(slim_scaling_factor),
        fixed_point_residual_generations=float(signed_residual[best_index]),
        grid_step_generations=float(grid_step_generations),
        grid_max_generations=float(grid_max_generations),
        boundary_limited=boundary_limited,
    )


def build_de_novo_origin_age_grid(
    demographies: Mapping[str, EasDemographyModel],
    *,
    selection_coefficients: Sequence[float] = SELECTION_COEFFICIENTS,
    target_frequency: float = DE_NOVO_TARGET_FREQUENCY,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    slim_scaling_factor: float = 1.0,
    grid_step_generations: float = 1.0,
) -> tuple[OriginAgeEstimate, ...]:
    """Return Q-aware deterministic origin ages for every s and quantile."""

    missing = {label for label, _ in _QUANTILES}.difference(demographies)
    if missing:
        raise ValueError(f"demographies are missing trajectories: {sorted(missing)}")
    estimates = []
    for selection_coefficient in selection_coefficients:
        for label, _ in _QUANTILES:
            model = demographies[label]
            estimates.append(
                estimate_de_novo_origin_age(
                    model.time_generations,
                    model.ne,
                    selection_coefficient=float(selection_coefficient),
                    trajectory_label=label,
                    target_frequency=target_frequency,
                    dominance_coefficient=dominance_coefficient,
                    slim_scaling_factor=slim_scaling_factor,
                    grid_step_generations=grid_step_generations,
                )
            )
    return tuple(estimates)


def build_focal_contig(
    *,
    single_site_id: str = FOCAL_SITE_ID,
    sequence_length: int = SEQUENCE_LENGTH_BP,
    focal_position: int = FOCAL_POSITION_BP,
) -> stdpopsim.Contig:
    """Build the 10-Mb contig and reserve the focal base as a single site."""

    _require_stdpopsim_030()
    if sequence_length < 2 or not 0 <= focal_position < sequence_length:
        raise ValueError("focal_position must fall inside sequence_length")
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig(
        length=int(sequence_length),
        mutation_rate=MUTATION_RATE,
        recombination_rate=RECOMBINATION_RATE,
    )
    contig.add_single_site(
        id=single_site_id,
        coordinate=int(focal_position),
        description="EAS focal selected mutation",
        long_description=(
            "A single mutation drawn by a stdpopsim extended event; the focal "
            "base is excluded from the ordinary neutral mutation overlay."
        ),
    )
    return contig


def _validate_selection_parameters(
    selection_coefficient: float,
    dominance_coefficient: float,
) -> None:
    if not math.isfinite(selection_coefficient) or selection_coefficient <= 0:
        raise ValueError("selection_coefficient must be positive and finite")
    if not math.isfinite(dominance_coefficient) or not 0 <= dominance_coefficient <= 1:
        raise ValueError("dominance_coefficient must be finite and in [0, 1]")


def _validate_frequency_interval(
    lower: float,
    upper: float | None,
    *,
    target: float | None = None,
) -> None:
    if not math.isfinite(lower) or not 0 < lower <= 1:
        raise ValueError("lower allele-frequency bound must be in (0, 1]")
    if upper is not None:
        if not math.isfinite(upper) or not 0 <= upper < 1:
            raise ValueError("upper allele-frequency bound must be in [0, 1)")
        if upper < lower:
            raise ValueError("upper allele-frequency bound must be >= lower bound")
    if target is not None and not (lower <= target <= (upper or 1.0)):
        raise ValueError("target_frequency must fall inside the conditioning interval")


def build_no_introgression_events(
    *,
    origin_age_generations: float,
    selection_coefficient: float,
    target_frequency: float = DE_NOVO_TARGET_FREQUENCY,
    lower_frequency: float = DE_NOVO_DEFAULT_LOWER_FREQUENCY,
    upper_frequency: float = DE_NOVO_DEFAULT_UPPER_FREQUENCY,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    population: str = EAS_POPULATION,
    single_site_id: str = FOCAL_SITE_ID,
) -> tuple[stdpopsim.ExtendedEvent, ...]:
    """Construct a de novo sweep with present-day lower and upper AF bounds."""

    _require_stdpopsim_030()
    _validate_selection_parameters(selection_coefficient, dominance_coefficient)
    if not math.isfinite(origin_age_generations) or origin_age_generations <= 1:
        raise ValueError("origin_age_generations must be finite and greater than one")
    if not 0 < target_frequency < 1:
        raise ValueError("target_frequency must be in (0, 1)")
    _validate_frequency_interval(
        lower_frequency, upper_frequency, target=target_frequency
    )

    return (
        stdpopsim.DrawMutation(
            time=float(origin_age_generations),
            single_site_id=single_site_id,
            population=population,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(float(origin_age_generations)),
            end_time=0.0,
            single_site_id=single_site_id,
            population=population,
            op=">",
            allele_frequency=0.0,
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=float(origin_age_generations),
            end_time=0.0,
            single_site_id=single_site_id,
            population=population,
            selection_coeff=float(selection_coefficient),
            dominance_coeff=float(dominance_coefficient),
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=0.0,
            end_time=0.0,
            single_site_id=single_site_id,
            population=population,
            op=">=",
            allele_frequency=float(lower_frequency),
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=0.0,
            end_time=0.0,
            single_site_id=single_site_id,
            population=population,
            op="<=",
            allele_frequency=float(upper_frequency),
        ),
    )


def _serialize_time(value: float | stdpopsim.GenerationAfter) -> dict[str, Any]:
    return {
        "generations_ago": float(value),
        "semantics": (
            "generation_after_in_forward_time"
            if isinstance(value, stdpopsim.GenerationAfter)
            else "exact_generation"
        ),
    }


def serialize_extended_events(
    events: Iterable[stdpopsim.ExtendedEvent],
) -> list[dict[str, Any]]:
    """Serialize the stdpopsim event attributes without losing time markers."""

    records = []
    for event in events:
        record: dict[str, Any] = {
            "event_type": type(event).__name__,
            "single_site_id": event.single_site_id,
        }
        if hasattr(event, "population"):
            record["population"] = event.population
        if hasattr(event, "time"):
            record["time"] = _serialize_time(event.time)
        if hasattr(event, "start_time"):
            record["start_time"] = _serialize_time(event.start_time)
        if hasattr(event, "end_time"):
            record["end_time"] = _serialize_time(event.end_time)
        if isinstance(event, stdpopsim.ChangeMutationFitness):
            record["selection_coefficient"] = float(event.selection_coeff)
            record["dominance_coefficient"] = float(event.dominance_coeff)
        if isinstance(event, stdpopsim.ConditionOnAlleleFrequency):
            record["operator"] = event.op
            record["allele_frequency"] = float(event.allele_frequency)
        records.append(record)
    json.dumps(records)
    return records


@dataclass(frozen=True)
class SingleSiteSweepSpec:
    """Runtime objects and a separate JSON-safe specification record."""

    scenario: str
    contig: stdpopsim.Contig
    extended_events: tuple[stdpopsim.ExtendedEvent, ...]
    record: dict[str, Any]
    demographic_model: stdpopsim.DemographicModel | None = None

    def provenance_record(self) -> dict[str, Any]:
        return copy.deepcopy(self.record)


def _contig_record(contig: stdpopsim.Contig, single_site_id: str) -> dict[str, Any]:
    index = [i for i, dfe in enumerate(contig.dfe_list) if dfe.id == single_site_id]
    if len(index) != 1:
        raise ValueError(f"single site {single_site_id!r} is not unique on the contig")
    interval = np.asarray(contig.interval_list[index[0]], dtype=np.int64)
    if interval.shape != (1, 2) or interval[0, 1] - interval[0, 0] != 1:
        raise ValueError(
            f"single site {single_site_id!r} does not cover exactly one base"
        )
    return {
        "sequence_length_bp": int(contig.length),
        "mutation_rate_per_bp_per_generation": float(contig.mutation_rate),
        "recombination_rate_per_bp_per_generation": float(
            contig.recombination_map.mean_rate
        ),
        "single_site_id": single_site_id,
        "focal_interval_zero_based_half_open": [int(x) for x in interval[0]],
        "focal_mutation_draw_count": 1,
        "ordinary_mutation_overlay_at_focal_base": False,
        "recurrence_at_focal_base": False,
        "single_site_api": "stdpopsim.Contig.add_single_site",
    }


def build_no_introgression_sweep_spec(
    *,
    origin_age_generations: float,
    selection_coefficient: float,
    target_frequency: float = DE_NOVO_TARGET_FREQUENCY,
    lower_frequency: float = DE_NOVO_DEFAULT_LOWER_FREQUENCY,
    upper_frequency: float = DE_NOVO_DEFAULT_UPPER_FREQUENCY,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    single_site_id: str = FOCAL_SITE_ID,
) -> SingleSiteSweepSpec:
    contig = build_focal_contig(single_site_id=single_site_id)
    events = build_no_introgression_events(
        origin_age_generations=origin_age_generations,
        selection_coefficient=selection_coefficient,
        target_frequency=target_frequency,
        lower_frequency=lower_frequency,
        upper_frequency=upper_frequency,
        dominance_coefficient=dominance_coefficient,
        single_site_id=single_site_id,
    )
    record = {
        "schema": "gamma-smc.no-introgression-sweep-spec/v1",
        "scenario": "no_introgression_de_novo",
        "population": EAS_POPULATION,
        "origin_age_generations": float(origin_age_generations),
        "origin_age_years": float(origin_age_generations) * GENERATION_TIME_YEARS,
        "selection_coefficient": float(selection_coefficient),
        "dominance_coefficient": float(dominance_coefficient),
        "target_frequency": float(target_frequency),
        "target_frequency_is_placeholder": True,
        "present_frequency_interval": {
            "lower_inclusive": float(lower_frequency),
            "upper_inclusive": float(upper_frequency),
        },
        "generation_time_years": GENERATION_TIME_YEARS,
        "tmrca_thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
        "tmrca_thresholds_generations": list(TMRCA_THRESHOLDS_GENERATIONS),
        "contig": _contig_record(contig, single_site_id),
        "extended_events": serialize_extended_events(events),
    }
    json.dumps(record)
    return SingleSiteSweepSpec(
        scenario="no_introgression_de_novo",
        contig=contig,
        extended_events=events,
        record=record,
    )


def _population_name(model: stdpopsim.DemographicModel, population: int | str) -> str:
    if isinstance(population, str):
        return population
    return model.model.populations[int(population)].name


def _find_mass_migration(
    model: stdpopsim.DemographicModel,
    *,
    time: float,
    source: str,
    destination: str,
    proportion: float,
) -> None:
    for event in model.model.events:
        if not isinstance(event, msprime.MassMigration):
            continue
        if (
            math.isclose(float(event.time), time)
            and _population_name(model, event.source) == source
            and _population_name(model, event.dest) == destination
            and math.isclose(float(event.proportion), proportion)
        ):
            return
    raise RuntimeError(
        f"{model.id} is missing expected mass migration at {time:g}: "
        f"{source} -> {destination} ({proportion:g})"
    )


def load_ancient_eurasia_model() -> stdpopsim.DemographicModel:
    """Load and contract-check stdpopsim 0.3.0 ``AncientEurasia_9K19``."""

    _require_stdpopsim_030()
    catalog_model = stdpopsim.get_species("HomSap").get_demographic_model(
        ANCIENT_EURASIA_MODEL_ID
    )
    model = copy.deepcopy(catalog_model)
    if model.generation_time != GENERATION_TIME_YEARS:
        raise RuntimeError(
            f"{ANCIENT_EURASIA_MODEL_ID} generation time changed from 25 years"
        )
    required_populations = {
        ARCHAIC_POPULATION,
        INTROGRESSION_RECIPIENT_POPULATION,
        INTROGRESSION_TARGET_POPULATION,
    }
    names = {population.name for population in model.populations}
    if not required_populations.issubset(names):
        raise RuntimeError(
            f"{ANCIENT_EURASIA_MODEL_ID} is missing populations: "
            f"{sorted(required_populations.difference(names))}"
        )
    _find_mass_migration(
        model,
        time=INTROGRESSION_PULSE_GENERATIONS,
        source=INTROGRESSION_RECIPIENT_POPULATION,
        destination=ARCHAIC_POPULATION,
        proportion=0.0296,
    )
    _find_mass_migration(
        model,
        time=HAN_SPLIT_GENERATIONS,
        source=INTROGRESSION_TARGET_POPULATION,
        destination=INTROGRESSION_RECIPIENT_POPULATION,
        proportion=1.0,
    )
    _find_mass_migration(
        model,
        time=HUMAN_NEANDERTHAL_SPLIT_GENERATIONS,
        source=ARCHAIC_POPULATION,
        destination=INTROGRESSION_RECIPIENT_POPULATION,
        proportion=1.0,
    )
    # The generic study contig supplies the explicit rates. Keep the copied
    # model's metadata consistent with those simulation inputs as well.
    model.mutation_rate = MUTATION_RATE
    model.recombination_rate = RECOMBINATION_RATE
    return model


def build_introgression_events(
    *,
    selection_coefficient: float,
    minimum_han_frequency: float,
    maximum_han_frequency: float | None = None,
    mutation_age_generations: float = DEFAULT_ARCHAIC_ORIGIN_GENERATIONS,
    source_minimum_frequency: float = DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY,
    realized_han_split_generations: float = HAN_SPLIT_GENERATIONS,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    single_site_id: str = FOCAL_SITE_ID,
) -> tuple[stdpopsim.ExtendedEvent, ...]:
    """Construct an archaic-specific standing-variation sweep into Han.

    The mutation is drawn on the Neanderthal branch after the human-Neanderthal
    split and before the 2.96% pulse into Loschbour.  Selection is absent in the
    source, starts in Loschbour at the pulse, and continues in Han after the
    explicitly supplied realized split tick. The default is the unscaled
    2,016-generation catalog split.
    """

    _require_stdpopsim_030()
    _validate_selection_parameters(selection_coefficient, dominance_coefficient)
    _validate_frequency_interval(minimum_han_frequency, maximum_han_frequency)
    if not math.isfinite(source_minimum_frequency) or not (
        0 < source_minimum_frequency <= 1
    ):
        raise ValueError("source_minimum_frequency must be in (0, 1]")
    if not math.isfinite(mutation_age_generations) or not (
        PRE_PULSE_SOURCE_CHECK_GENERATIONS
        < mutation_age_generations
        < HUMAN_NEANDERTHAL_SPLIT_GENERATIONS
    ):
        raise ValueError(
            "mutation_age_generations must be after the human-Neanderthal split "
            "and before the introgression pulse"
        )
    if not math.isfinite(realized_han_split_generations) or not (
        0 < realized_han_split_generations < INTROGRESSION_PULSE_GENERATIONS
    ):
        raise ValueError(
            "realized_han_split_generations must be positive and more recent "
            "than the introgression pulse"
        )

    events: list[stdpopsim.ExtendedEvent] = [
        stdpopsim.DrawMutation(
            time=float(mutation_age_generations),
            single_site_id=single_site_id,
            population=ARCHAIC_POPULATION,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=stdpopsim.GenerationAfter(float(mutation_age_generations)),
            end_time=PRE_PULSE_SOURCE_CHECK_GENERATIONS,
            single_site_id=single_site_id,
            population=ARCHAIC_POPULATION,
            op=">",
            allele_frequency=0.0,
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=PRE_PULSE_SOURCE_CHECK_GENERATIONS,
            end_time=PRE_PULSE_SOURCE_CHECK_GENERATIONS,
            single_site_id=single_site_id,
            population=ARCHAIC_POPULATION,
            op=">=",
            allele_frequency=float(source_minimum_frequency),
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=INTROGRESSION_PULSE_GENERATIONS,
            end_time=float(realized_han_split_generations),
            single_site_id=single_site_id,
            population=INTROGRESSION_RECIPIENT_POPULATION,
            selection_coeff=float(selection_coefficient),
            dominance_coeff=float(dominance_coefficient),
        ),
        stdpopsim.ChangeMutationFitness(
            start_time=float(realized_han_split_generations),
            end_time=0.0,
            single_site_id=single_site_id,
            population=INTROGRESSION_TARGET_POPULATION,
            selection_coeff=float(selection_coefficient),
            dominance_coeff=float(dominance_coefficient),
        ),
        stdpopsim.ConditionOnAlleleFrequency(
            start_time=0.0,
            end_time=0.0,
            single_site_id=single_site_id,
            population=INTROGRESSION_TARGET_POPULATION,
            op=">=",
            allele_frequency=float(minimum_han_frequency),
        ),
    ]
    if maximum_han_frequency is not None:
        events.append(
            stdpopsim.ConditionOnAlleleFrequency(
                start_time=0.0,
                end_time=0.0,
                single_site_id=single_site_id,
                population=INTROGRESSION_TARGET_POPULATION,
                op="<=",
                allele_frequency=float(maximum_han_frequency),
            )
        )
    return tuple(events)


def build_introgression_sweep_spec(
    *,
    selection_coefficient: float,
    minimum_han_frequency: float,
    maximum_han_frequency: float | None = None,
    mutation_age_generations: float = DEFAULT_ARCHAIC_ORIGIN_GENERATIONS,
    source_minimum_frequency: float = DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY,
    realized_han_split_generations: float = HAN_SPLIT_GENERATIONS,
    dominance_coefficient: float = DOMINANCE_COEFFICIENT,
    single_site_id: str = FOCAL_SITE_ID,
) -> SingleSiteSweepSpec:
    model = load_ancient_eurasia_model()
    contig = build_focal_contig(single_site_id=single_site_id)
    events = build_introgression_events(
        selection_coefficient=selection_coefficient,
        minimum_han_frequency=minimum_han_frequency,
        maximum_han_frequency=maximum_han_frequency,
        mutation_age_generations=mutation_age_generations,
        source_minimum_frequency=source_minimum_frequency,
        realized_han_split_generations=realized_han_split_generations,
        dominance_coefficient=dominance_coefficient,
        single_site_id=single_site_id,
    )
    record = {
        "schema": "gamma-smc.introgression-sweep-spec/v1",
        "scenario": "neanderthal_introgression_standing_variation",
        "stdpopsim_version": stdpopsim.__version__,
        "demographic_model": ANCIENT_EURASIA_MODEL_ID,
        "generation_time_years": GENERATION_TIME_YEARS,
        "mutation_origin": {
            "population": ARCHAIC_POPULATION,
            "age_generations": float(mutation_age_generations),
            "age_years": float(mutation_age_generations) * GENERATION_TIME_YEARS,
            "after_human_neanderthal_split": True,
            "before_introgression_pulse": True,
            "archaic_specific_no_ils_by_construction": True,
        },
        "source_frequency_condition": {
            "population": ARCHAIC_POPULATION,
            "time_generations": PRE_PULSE_SOURCE_CHECK_GENERATIONS,
            "timing": "one_generation_before_introgression_pulse",
            "minimum_inclusive": float(source_minimum_frequency),
        },
        "introgression": {
            "recipient_population": INTROGRESSION_RECIPIENT_POPULATION,
            "pulse_generations_ago": INTROGRESSION_PULSE_GENERATIONS,
            "pulse_proportion": 0.0296,
        },
        "selection": {
            "selection_coefficient": float(selection_coefficient),
            "dominance_coefficient": float(dominance_coefficient),
            "mode": "standing_variation_at_introgression_onset",
            "requested_han_split_generations_ago": HAN_SPLIT_GENERATIONS,
            "realized_han_split_generations_ago": float(realized_han_split_generations),
            "split_time_semantics": (
                "requested catalog time and caller-supplied realized simulation "
                "tick are recorded separately"
            ),
            "loschbour_interval_generations_ago": [
                INTROGRESSION_PULSE_GENERATIONS,
                float(realized_han_split_generations),
            ],
            "han_interval_generations_ago": [
                float(realized_han_split_generations),
                0.0,
            ],
        },
        "present_han_frequency": {
            "minimum_inclusive": float(minimum_han_frequency),
            "maximum_inclusive": (
                float(maximum_han_frequency)
                if maximum_han_frequency is not None
                else None
            ),
        },
        "tmrca_thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
        "tmrca_thresholds_generations": list(TMRCA_THRESHOLDS_GENERATIONS),
        "contig": _contig_record(contig, single_site_id),
        "extended_events": serialize_extended_events(events),
    }
    json.dumps(record)
    return SingleSiteSweepSpec(
        scenario="neanderthal_introgression_standing_variation",
        contig=contig,
        extended_events=events,
        record=record,
        demographic_model=model,
    )


__all__ = [
    "ANCIENT_EURASIA_MODEL_ID",
    "ARCHAIC_POPULATION",
    "DE_NOVO_DEFAULT_LOWER_FREQUENCY",
    "DE_NOVO_DEFAULT_UPPER_FREQUENCY",
    "DE_NOVO_TARGET_FREQUENCY",
    "DEFAULT_ARCHAIC_ORIGIN_GENERATIONS",
    "DEFAULT_ARCHAIC_SOURCE_MINIMUM_FREQUENCY",
    "DOMINANCE_COEFFICIENT",
    "EAS_POPULATION",
    "EasDemographyModel",
    "FOCAL_POSITION_BP",
    "FOCAL_SITE_ID",
    "GENERATION_TIME_YEARS",
    "HAN_SPLIT_GENERATIONS",
    "HUMAN_NEANDERTHAL_SPLIT_GENERATIONS",
    "INTROGRESSION_PULSE_GENERATIONS",
    "INTROGRESSION_RECIPIENT_POPULATION",
    "INTROGRESSION_TARGET_FREQUENCIES",
    "INTROGRESSION_TARGET_POPULATION",
    "MUTATION_RATE",
    "OriginAgeEstimate",
    "PHLASH_MVN_SCHEMA",
    "PRE_PULSE_SOURCE_CHECK_GENERATIONS",
    "PhlashEasArtifact",
    "RECOMBINATION_RATE",
    "SELECTION_COEFFICIENTS",
    "SEQUENCE_LENGTH_BP",
    "STDPOPSIM_VERSION",
    "SingleSiteSweepSpec",
    "TMRCA_THRESHOLDS_GENERATIONS",
    "TMRCA_THRESHOLDS_YEARS",
    "build_de_novo_origin_age_grid",
    "build_eas_demography_models",
    "build_focal_contig",
    "build_introgression_events",
    "build_introgression_sweep_spec",
    "build_no_introgression_events",
    "build_no_introgression_sweep_spec",
    "estimate_de_novo_origin_age",
    "load_ancient_eurasia_model",
    "load_phlash_eas_npz",
    "serialize_extended_events",
]
