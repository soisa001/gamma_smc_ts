"""Immutable plan for the focused selected-versus-neutral EAS campaign.

This module deliberately owns the campaign namespace instead of reusing the
legacy study directories.  Its first phase defines stable execution units and
their analysis contract.  Simulation, decoding, and analysis phases are added
behind that frozen plan; they must never reinterpret an existing unit ID.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import zstandard

from .defaults import DEFAULT_CACHE_SIZE
from .eas_sweep_models import (
    FOCAL_POSITION_BP,
    GENERATION_TIME_YEARS,
    MUTATION_RATE,
    RECOMBINATION_RATE,
    SEQUENCE_LENGTH_BP,
    TMRCA_THRESHOLDS_GENERATIONS,
    TMRCA_THRESHOLDS_YEARS,
)
from .eas_sweep_study import BASE_SEED

SCHEMA_VERSION = "gamma-smc.focused-selection-campaign/v1"
DEFAULT_CAMPAIGN_DIR = "focused_selection_EAS_sim"
DEFAULT_SELECTION_COEFFICIENTS = (0.01, 0.005)
DEFAULT_SELECTED_REPLICATES = 10
DEFAULT_NEUTRAL_REPLICATES = 100
DEFAULT_AF_TARGETS = (0.10, 0.20, 0.30)
DEFAULT_AF_HALF_WIDTH = 0.025
DEFAULT_SAMPLE_DIPLOIDS = 100
DEFAULT_SELECTED_POOL_DIPLOIDS = 500
DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE = 4
DEFAULT_SLIM_SCALING_FACTOR = 5.0
DEFAULT_SELECTED_MAX_DRAWS = 50
DEFAULT_SELECTED_DRAW_TIMEOUT_MINUTES = 5.0
DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_MINUTES = 45.0
DEFAULT_HAN_SELECTION_CALIBRATION = (
    "focused_selection_EAS_sim/calibration/han_selection_end_frozen.json"
)
DECODE_IMPLEMENTATION_PATHS = (
    "python/gamma_smc_aou/focused_selection_campaign.py",
    "python/gamma_smc_aou/focused_selection_decode.py",
    "python/gamma_smc_aou/decoder.py",
    "python/gamma_smc_aou/tree_sequence.py",
    "python/gamma_smc_aou/eas_sweep_models.py",
)
DECODE_AUXILIARY_OUTPUTS = {
    "overall_summary": "overall.summary.tsv",
    "run_manifest": "overall.summary.tsv.run.json",
    "decoder_stdout": "decoder.stdout.zst",
    "decoder_stderr": "decoder.stderr.zst",
}
DEMOGRAPHIES = (
    {
        "demography_id": "eas_phlash_median",
        "demography_kind": "no_introgression",
        "population": "EAS",
        "source_model": "PHLASH EAS pointwise median trajectory",
        "selection_origin": "de_novo",
    },
    {
        "demography_id": "ancient_eurasia_han_introgression",
        "demography_kind": "introgression",
        "population": "Han",
        "source_model": "stdpopsim AncientEurasia_9K19",
        "selection_origin": "archaic_specific_introgressed_standing_variation",
    },
)
EXECUTION_COLUMNS = (
    "unit_id",
    "demography_id",
    "demography_kind",
    "population",
    "source_model",
    "selection_origin",
    "simulation_class",
    "selection_coefficient",
    "target_allele_frequency",
    "population_af_lower",
    "population_af_upper",
    "exact_sample_alt_count",
    "sample_diploids",
    "replicate_index",
    "seed",
    "sequence_length_bp",
    "focal_position_bp",
)


@dataclass(frozen=True)
class FocusedCampaignPlan:
    """Result-affecting dimensions of the focused campaign."""

    repo_root: Path
    campaign_dir: Path
    selected_replicates: int = DEFAULT_SELECTED_REPLICATES
    neutral_replicates: int = DEFAULT_NEUTRAL_REPLICATES
    base_seed: int = BASE_SEED

    def validate(self) -> None:
        repo_root = self.repo_root.resolve()
        campaign_dir = self.campaign_dir.resolve()
        try:
            campaign_dir.relative_to(repo_root)
        except ValueError as error:
            raise ValueError("campaign_dir must be inside repo_root") from error
        if any("onedrive" in part.casefold() for part in campaign_dir.parts):
            raise ValueError("focused campaign output must not be inside OneDrive")
        for name in ("selected_replicates", "neutral_replicates"):
            value = getattr(self, name)
            if isinstance(value, bool) or int(value) != value or int(value) < 1:
                raise ValueError(f"{name} must be a positive integer")
        if int(self.base_seed) < 1:
            raise ValueError("base_seed must be positive")


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    compression: str | Mapping[str, Any] | None = None
    if path.suffix == ".gz":
        compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            compression=compression,
            float_format="%.12g",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _decode_source_hashes(repo_root: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in DECODE_IMPLEMENTATION_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise ValueError(f"decode implementation source is absent: {path}")
        result[relative] = _sha256_file(path)
    return result


def _portable_decoder_path(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _decode_settings(present_ne: float, *, threads: int) -> dict[str, Any]:
    if not np.isfinite(float(present_ne)) or float(present_ne) <= 0:
        raise ValueError("decode present-day Ne must be finite and positive")
    if int(threads) < 1:
        raise ValueError("focused decode threads must be positive")
    return {
        "input_format": "trees_via_streamed_vcf",
        "output_at_stride": 10_000,
        "output_at_hets": False,
        "pair_selector": "explicit_within_individual_pairs",
        "cache_size": int(DEFAULT_CACHE_SIZE),
        "pair_block": 256,
        "exp10": "accurate",
        "backward_alignment": "fixed",
        "recent_call": "median",
        "recent_call_probability": 0.5,
        "scaled_mutation_rate": 4 * float(present_ne) * MUTATION_RATE,
        "recombination_to_mutation_ratio": RECOMBINATION_RATE / MUTATION_RATE,
        "unscaled_mutation_rate": MUTATION_RATE,
        "generation_time_years": GENERATION_TIME_YEARS,
        "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
        "threads": int(threads),
    }


def _compact_decoder_logs(summary_path: Path, run: Mapping[str, Any]) -> dict[str, Any]:
    """Compress verbose decoder streams and replace them with checksum records."""

    compact = dict(run)
    for stream in ("stdout", "stderr"):
        content = str(compact.pop(stream, "") or "").encode("utf-8")
        path = summary_path.with_name(f"decoder.{stream}.zst")
        temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
        try:
            temporary.write_bytes(zstandard.ZstdCompressor(level=6).compress(content))
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        compact[f"{stream}_artifact"] = {
            "path": path.name,
            "uncompressed_size_bytes": len(content),
            "compressed_size_bytes": path.stat().st_size,
            "sha256": _sha256_file(path),
        }
    run_path = summary_path.with_suffix(summary_path.suffix + ".run.json")
    _atomic_json(run_path, compact)
    return compact


def _stable_seed(base_seed: int, unit_id: str) -> int:
    token = f"{base_seed}:{unit_id}".encode()
    offset = int.from_bytes(hashlib.sha256(token).digest()[:8], "big")
    return int(1 + offset % (2**31 - 2))


def _frequency_slug(value: float) -> str:
    return f"af{round(100 * value):02d}"


def _selection_slug(value: float) -> str:
    return f"s{value:.3f}".replace(".", "p")


def build_execution_units(plan: FocusedCampaignPlan) -> pd.DataFrame:
    """Build stable units; defaults yield 720 simulations across the grid."""
    plan.validate()
    rows: list[dict[str, Any]] = []
    for demography in DEMOGRAPHIES:
        for target_af in DEFAULT_AF_TARGETS:
            classes = [
                ("selected", int(plan.selected_replicates), selection)
                for selection in DEFAULT_SELECTION_COEFFICIENTS
            ]
            # The neutral distribution has no selection coefficient.  Keep one
            # independent null bank per demography-by-AF cell and reuse that
            # same bank for both prespecified selected-s comparisons; copying
            # it under two s labels would be pseudoreplication.
            classes.append(("neutral", int(plan.neutral_replicates), 0.0))
            for simulation_class, replicates, selection in classes:
                for replicate_index in range(1, replicates + 1):
                    class_slug = (
                        f"selected_{_selection_slug(selection)}"
                        if simulation_class == "selected"
                        else "neutral"
                    )
                    unit_id = (
                        f"{demography['demography_id']}__"
                        f"{_frequency_slug(target_af)}__{class_slug}__"
                        f"rep{replicate_index:03d}"
                    )
                    rows.append(
                        {
                            "unit_id": unit_id,
                            **demography,
                            "simulation_class": simulation_class,
                            "selection_coefficient": selection,
                            "target_allele_frequency": target_af,
                            "population_af_lower": target_af - DEFAULT_AF_HALF_WIDTH,
                            "population_af_upper": target_af + DEFAULT_AF_HALF_WIDTH,
                            "exact_sample_alt_count": round(
                                2 * DEFAULT_SAMPLE_DIPLOIDS * target_af
                            ),
                            "sample_diploids": DEFAULT_SAMPLE_DIPLOIDS,
                            "replicate_index": replicate_index,
                            "seed": _stable_seed(plan.base_seed, unit_id),
                            "sequence_length_bp": SEQUENCE_LENGTH_BP,
                            "focal_position_bp": FOCAL_POSITION_BP,
                        }
                    )
    frame = pd.DataFrame.from_records(rows, columns=EXECUTION_COLUMNS)
    if frame["unit_id"].duplicated().any() or frame["seed"].duplicated().any():
        raise RuntimeError("focused campaign unit IDs and seeds must be unique")
    return frame


def _tsv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")


def build_study_design(
    plan: FocusedCampaignPlan, units: pd.DataFrame
) -> dict[str, Any]:
    """Return the JSON-safe scientific and execution contract."""
    selected = units[units["simulation_class"] == "selected"]
    neutral = units[units["simulation_class"] == "neutral"]
    return {
        "schema": SCHEMA_VERSION,
        "phase": "plan",
        "campaign_root": ".",
        "scientific_grid": {
            "demographies": [dict(record) for record in DEMOGRAPHIES],
            "selection_coefficients": list(DEFAULT_SELECTION_COEFFICIENTS),
            "target_allele_frequencies": list(DEFAULT_AF_TARGETS),
            "selected_replicates_per_cell": int(plan.selected_replicates),
            "neutral_replicates_per_cell": int(plan.neutral_replicates),
            "cells": len(DEMOGRAPHIES) * len(DEFAULT_AF_TARGETS),
            "selected_units": len(selected),
            "neutral_units": len(neutral),
            "total_units": len(units),
        },
        "sequence": {
            "length_bp": SEQUENCE_LENGTH_BP,
            "focal_position_0based": FOCAL_POSITION_BP,
            "mutation_rate_per_bp_per_generation": MUTATION_RATE,
            "recombination_rate_per_bp_per_generation": RECOMBINATION_RATE,
        },
        "time": {
            "generation_time_years": GENERATION_TIME_YEARS,
            "thresholds_years": list(TMRCA_THRESHOLDS_YEARS),
            "thresholds_generations": list(TMRCA_THRESHOLDS_GENERATIONS),
        },
        "rng": {
            "base_seed": int(plan.base_seed),
            "derivation": "sha256(base_seed:unit_id)",
            "integer_mapping": "1 + first_8_digest_bytes mod (2^31 - 2)",
        },
        "simulation_contract": {
            "selected_engine": "SLiM through stdpopsim 0.3 model contracts",
            "slim_scaling_factor": DEFAULT_SLIM_SCALING_FACTOR,
            "slim_burn_in": 0.1,
            "selected_acceptance_execution": {
                "planned_accepted_replicates_per_s_by_af_by_demography_cell": int(
                    plan.selected_replicates
                ),
                "immutable_unit_semantics": (
                    "each planned selected unit retries a deterministic disjoint draw "
                    "stream until one trajectory passes focal identity, candidate-pool "
                    "AF band, exact-k sampling, and genotype QC"
                ),
                "maximum_draws_per_unit": DEFAULT_SELECTED_MAX_DRAWS,
                "draw_timeout_minutes": DEFAULT_SELECTED_DRAW_TIMEOUT_MINUTES,
                "cumulative_timeout_minutes_per_unit": (
                    DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_MINUTES
                ),
                "minimum_calibration_band_and_panel_hit_rate": 0.20,
                "exhaustion_probability_at_minimum_hit_rate": (
                    (1.0 - 0.20) ** DEFAULT_SELECTED_MAX_DRAWS
                ),
                "replacement_units": False,
            },
            "neutral_engine": (
                "msprime DTWF for 200 generations then SmcPrimeApproxCoalescent "
                "under matched demography"
            ),
            "neutral_recent_dtwf_duration_generations": 200,
            "neutral_coalescent_approximation": "SMC-prime after 200 generations",
            "neutral_independence": (
                "four independently seeded ancestry candidates per fixed batch; "
                "simulation units have disjoint seed streams"
            ),
            "neutral_af_semantics": (
                "a genuine conditioned neutral SNP whose 500-diploid candidate "
                "pool is inside the same AF band as selected trajectories, followed "
                "by the same exact-k 100-diploid panel sampler; genealogy-independent "
                "pseudo-genotype labels are forbidden in the primary null"
            ),
            "population_af_half_width": DEFAULT_AF_HALF_WIDTH,
            "population_af_intervals_are_nonoverlapping": True,
            "sample_diploids": DEFAULT_SAMPLE_DIPLOIDS,
            "candidate_pool_diploids": DEFAULT_SELECTED_POOL_DIPLOIDS,
            "selected_candidate_pool_diploids": DEFAULT_SELECTED_POOL_DIPLOIDS,
            "neutral_candidate_pool_diploids": DEFAULT_SELECTED_POOL_DIPLOIDS,
            "candidate_pool_frequency_gate_scope": "selected_and_neutral",
            "candidate_pool_frequency_gate_order": "before_exact_sample_panel",
            "selected_terminal_population_af_conditioning_in_slim": False,
            "neutral_ancestry_batch_size": DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE,
            "exact_sample_af_semantics": (
                "exact alternate-copy count k=2*N*target_af; k is 20, 40, or "
                "60 among 200 sampled haplotypes"
            ),
            "minimum_hom_ref_diploids": 2,
            "minimum_hom_alt_diploids": 2,
            "minimum_heterozygous_diploids": 1,
            "neutral_ascertainment": (
                "simulate fixed batches of four ancestries with 200 recent "
                "generations of DTWF followed by SMC-prime; for EAS, condition "
                "branches by mutation-time length in the fixed 5-Mb "
                "marginal tree of the requested 10-Mb sequence; for Han, scan a "
                "20-Mb ancestry by Loschbour-to-Neanderthal pulse migrations, "
                "ascend to the archaic branch crossing 2,400 generations, weight "
                "by eligible genomic span, and crop a 10-Mb window; in both cases "
                "the 500-diploid pool must pass the AF band and admit the same "
                "uniform exact-k 100-diploid panel sampler"
            ),
            "introgression_source_condition": (
                "single archaic origin present as nonzero standing variation at "
                "the pulse; neutral branches may feed one or several migrant lineages"
            ),
            "introgression_mutation_age_generations": 2_400.0,
            "introgression_selection_episode": {
                "starts_generations_ago": 2_272.0,
                "cessation_rule": (
                    "fixed s-by-AF endpoint from the disjoint terminal-AF calibration; "
                    "neutral thereafter"
                ),
                "frozen_calibration_default_path": (DEFAULT_HAN_SELECTION_CALIBRATION),
                "frozen_artifact_and_selected_endpoint_table_sha256_required": True,
                "pulse_proportion": 0.0296,
                "dominance_coefficient": 0.5,
                "selection_end_is_target_and_s_specific": True,
                "post_pulse_recipient_nonloss_conditioning": True,
            },
        },
        "decode_contract": {
            "gamma_raw_invocations_per_unit": 1,
            "postprocessing": (
                "derive overall, hom_ref, heterozygous, and hom_alt summaries "
                "from the same within-individual pair posterior"
            ),
            "raw_center_profile_required": True,
            "output_stride_bp": 10_000,
            "flow_field_cache_bp": int(DEFAULT_CACHE_SIZE),
            "raw_posterior_retained_by_default": True,
            "auxiliary_outputs_checksummed": sorted(DECODE_AUXILIARY_OUTPUTS),
            "spatial_profiles_aggregated_across_units": True,
        },
        "analysis_contract": {
            "within_selected_sim": "hom_alt versus hom_ref",
            "between_sim": "selected versus matched neutral",
            "tmrca_statistic": "P(TMRCA < x)",
            "roc_auc": "Mann-Whitney probability selected statistic exceeds neutral",
            "p_value": "one-sided upper-tail Monte Carlo with +1 correction",
            "neutral_leave_one_out_calibration": True,
            "raw_profile_key": [
                "unit_id",
                "demography_id",
                "simulation_class",
                "target_allele_frequency",
                "source",
                "genotype_class",
                "threshold_years",
            ],
        },
        "hashes": {
            "execution_units_tsv_sha256": hashlib.sha256(_tsv_bytes(units)).hexdigest()
        },
    }


def _read_existing_plan(campaign_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    table_path = campaign_dir / "execution_units.tsv"
    design_path = campaign_dir / "study_design.json"
    if table_path.is_file() != design_path.is_file():
        raise ValueError("focused campaign plan is incomplete")
    if not table_path.is_file():
        raise FileNotFoundError
    try:
        frame = pd.read_csv(table_path, sep="\t")
        design = json.loads(design_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise ValueError("focused campaign plan is unreadable") from error
    return frame, design


def write_plan(plan: FocusedCampaignPlan) -> tuple[Path, Path]:
    """Write once; subsequent calls checksum-validate without rewriting."""
    plan.validate()
    units = build_execution_units(plan)
    design = build_study_design(plan, units)
    campaign_dir = plan.campaign_dir.resolve()
    table_path = campaign_dir / "execution_units.tsv"
    design_path = campaign_dir / "study_design.json"
    try:
        recorded_units, recorded_design = _read_existing_plan(campaign_dir)
    except FileNotFoundError:
        campaign_dir.mkdir(parents=True, exist_ok=True)
        work_dir = campaign_dir / "work"
        results_dir = campaign_dir / "results"
        work_dir.mkdir(exist_ok=True)
        results_dir.mkdir(exist_ok=True)
        temporary_table = table_path.with_suffix(".tsv.tmp")
        temporary_table.write_bytes(_tsv_bytes(units))
        temporary_table.replace(table_path)
        temporary_design = design_path.with_suffix(".json.tmp")
        temporary_design.write_text(
            json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary_design.replace(design_path)
    else:
        if _canonical_sha256(recorded_design) != _canonical_sha256(design):
            raise ValueError(
                "existing focused campaign design is incompatible; choose a new directory"
            )
        try:
            pd.testing.assert_frame_equal(
                recorded_units,
                units,
                check_dtype=False,
                check_exact=False,
                rtol=1e-12,
                atol=1e-12,
            )
        except AssertionError as error:
            raise ValueError(
                "existing focused execution table is incompatible"
            ) from error
        # Hash the LF-normalized on-disk bytes so a Git checkout using CRLF and
        # a WSL checkout validate the same immutable plan without normalizing
        # float spellings through a pandas read/write round trip.
        actual = hashlib.sha256(
            table_path.read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest()
        if actual != recorded_design["hashes"]["execution_units_tsv_sha256"]:
            raise ValueError("existing focused execution table checksum is invalid")
    return table_path, design_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plan the focused EAS selected-versus-neutral campaign"
    )
    parser.add_argument(
        "phase",
        choices=("plan", "simulate", "decode", "analyze", "report", "status", "all"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument(
        "--selected-replicates", type=int, default=DEFAULT_SELECTED_REPLICATES
    )
    parser.add_argument(
        "--neutral-replicates", type=int, default=DEFAULT_NEUTRAL_REPLICATES
    )
    parser.add_argument("--base-seed", type=int, default=BASE_SEED)
    parser.add_argument(
        "--simulation-class", choices=("selected", "neutral", "all"), default="all"
    )
    parser.add_argument("--demography", action="append", default=[])
    parser.add_argument("--af", action="append", type=float, default=[])
    parser.add_argument(
        "--selection-coefficient", action="append", type=float, default=[]
    )
    parser.add_argument("--unit", action="append", default=[])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--slim-bin", type=Path)
    parser.add_argument("--decoder-bin", type=Path)
    parser.add_argument("--neutral-max-attempts", type=int, default=100)
    parser.add_argument(
        "--selected-max-draws", type=int, default=DEFAULT_SELECTED_MAX_DRAWS
    )
    parser.add_argument(
        "--selected-draw-timeout-minutes",
        type=float,
        default=DEFAULT_SELECTED_DRAW_TIMEOUT_MINUTES,
    )
    parser.add_argument(
        "--selected-cumulative-timeout-minutes",
        type=float,
        default=DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_MINUTES,
    )
    parser.add_argument(
        "--slim-scaling-factor", type=float, default=DEFAULT_SLIM_SCALING_FACTOR
    )
    parser.add_argument("--han-selection-calibration", type=Path)
    parser.add_argument("--empirical-tsv", type=Path)
    parser.add_argument("--remove-raw-after-success", action="store_true")
    return parser


def _filtered_units(units: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    selected = units.copy()
    if args.simulation_class != "all":
        selected = selected[selected["simulation_class"] == args.simulation_class]
    if args.demography:
        selected = selected[selected["demography_id"].isin(args.demography)]
    if args.af:
        selected = selected[
            np.logical_or.reduce(
                [
                    np.isclose(selected["target_allele_frequency"], value)
                    for value in args.af
                ]
            )
        ]
    if args.selection_coefficient:
        selected = selected[
            np.logical_or.reduce(
                [
                    np.isclose(selected["selection_coefficient"], value)
                    for value in args.selection_coefficient
                ]
            )
        ]
    if args.unit:
        selected = selected[
            selected["unit_id"].map(
                lambda value: any(token in str(value) for token in args.unit)
            )
        ]
    return selected.reset_index(drop=True)


def _simulation_status(campaign_dir: Path, units: pd.DataFrame) -> pd.DataFrame:
    from .focused_selection_simulation import (
        SCHEMA_VERSION as SIMULATION_SCHEMA_VERSION,
    )
    from .focused_selection_simulation import (
        unit_lock_is_held,
    )

    output_paths = {
        "tree": Path("simulation.trees"),
        "pair_table": Path("sample_manifest.tsv"),
        "overall_pairs": Path("pairs/overall.pairs.tsv"),
        "truth_profiles": Path("truth_profiles.tsv.gz"),
        "truth_class_summaries": Path("truth_class_summaries.tsv"),
    }
    rows = []
    for unit in units.to_dict(orient="records"):
        directory = campaign_dir / "work" / unit["unit_id"]
        completion = directory / "simulation_complete.json"
        failure = directory / "simulation_failed.json"
        error = ""
        failure_category = ""
        if completion.is_file():
            try:
                payload = json.loads(completion.read_text(encoding="utf-8"))
                contract = payload.get("contract")
                if (
                    payload.get("schema") != SIMULATION_SCHEMA_VERSION
                    or payload.get("status") != "complete"
                    or not isinstance(contract, Mapping)
                    or payload.get("contract_sha256") != _canonical_sha256(contract)
                ):
                    raise ValueError("completion schema/status/contract is invalid")
                recorded_unit = contract.get("unit")
                planned_unit = {key: unit[key] for key in EXECUTION_COLUMNS}
                if not isinstance(recorded_unit, Mapping) or _canonical_sha256(
                    recorded_unit
                ) != _canonical_sha256(planned_unit):
                    raise ValueError(
                        "completion unit does not match the immutable plan"
                    )
                outputs = payload.get("outputs")
                if not isinstance(outputs, Mapping):
                    raise TypeError("completion output manifest is invalid")
                for label, relative in output_paths.items():
                    record = outputs.get(label)
                    path = directory / relative
                    if (
                        not isinstance(record, Mapping)
                        or str(record.get("path", "")).replace("\\", "/")
                        != relative.as_posix()
                        or not path.is_file()
                        or _sha256_file(path) != record.get("sha256")
                    ):
                        raise ValueError(f"completion output is invalid: {label}")
            except (OSError, ValueError, json.JSONDecodeError) as caught:
                status = "invalid_completion"
                error = str(caught)
            else:
                status = "complete"
        elif unit_lock_is_held(directory / ".simulation.lock"):
            status = "running"
        elif failure.is_file():
            status = "failed_or_exhausted"
            try:
                failure_payload = json.loads(failure.read_text(encoding="utf-8"))
                failure_error = str(failure_payload.get("error", ""))
                failure_type = str(failure_payload.get("error_type", ""))
                if failure_type == "SimulationExhausted":
                    if "cumulative selected-simulation budget" in failure_error:
                        failure_category = "selected_cumulative_timeout_exhausted"
                    elif "produced no exact sampled AF" in failure_error:
                        failure_category = "selected_biological_draws_exhausted"
                    elif "neutral branch" in failure_error:
                        failure_category = "neutral_ascertainment_attempts_exhausted"
                    else:
                        failure_category = "bounded_search_exhausted"
                elif failure_type == "SimulationDrawTimeout":
                    failure_category = "selected_draw_timeout"
                else:
                    failure_category = "execution_failure"
                error = failure_error
            except (OSError, json.JSONDecodeError, TypeError, ValueError):
                failure_category = "unreadable_failure_record"
        elif directory.exists():
            status = "partial"
        else:
            status = "not_started"
        rows.append(
            {
                **unit,
                "status": status,
                "failure_category": failure_category,
                "status_error": error,
            }
        )
    frame = pd.DataFrame(rows)
    results = campaign_dir / "results"
    results.mkdir(parents=True, exist_ok=True)
    frame.to_csv(results / "simulation_status.tsv", sep="\t", index=False)
    return frame


def _decode_task_unlocked(payload: Mapping[str, Any]) -> dict[str, Any]:
    from .decoder import run_within_decoder
    from .focused_selection_decode import postprocess_raw_posteriors, sha256_file

    unit = dict(payload["unit"])
    directory = Path(payload["campaign_dir"]) / "work" / unit["unit_id"]
    simulation_completion = directory / "simulation_complete.json"
    if not simulation_completion.is_file():
        return {"unit_id": unit["unit_id"], "status": "not_simulated", "error": ""}
    decoded = directory / "decoded"
    completion_path = decoded / "completion.json"
    auxiliary_outputs = {
        label: decoded / filename
        for label, filename in DECODE_AUXILIARY_OUTPUTS.items()
    }
    if completion_path.is_file():
        try:
            artifacts = postprocess_raw_posteriors(
                decoded / "posterior.zst",
                directory / "sample_manifest.tsv",
                decoded,
                unit_record=unit,
                focal_position_0based=FOCAL_POSITION_BP,
                remove_raw_after_success=bool(payload["remove_raw_after_success"]),
                auxiliary_output_paths=auxiliary_outputs,
            )
            return {
                "unit_id": unit["unit_id"],
                "status": "cached",
                "error": "",
                "completion": str(artifacts.completion),
            }
        except ValueError as error:
            return {"unit_id": unit["unit_id"], "status": "failed", "error": str(error)}
    decoded.mkdir(parents=True, exist_ok=True)
    raw = decoded / "posterior.zst"
    summary = decoded / "overall.summary.tsv"
    run = run_within_decoder(
        payload["decoder_bin"],
        directory / "simulation.trees",
        summary,
        scaled_mutation_rate=4 * float(unit["present_ne"]) * MUTATION_RATE,
        recombination_to_mutation_ratio=RECOMBINATION_RATE / MUTATION_RATE,
        mutation_rate=MUTATION_RATE,
        threshold_years=TMRCA_THRESHOLDS_YEARS,
        generation_time=GENERATION_TIME_YEARS,
        input_format="trees",
        raw_output=raw,
        output_at_stride=10_000,
        output_at_hets=False,
        only_within=False,
        pairs_file=directory / "pairs" / "overall.pairs.tsv",
        threads=int(payload["threads"]),
        cache_size=DEFAULT_CACHE_SIZE,
        pair_block=256,
        exp10="accurate",
        backward_alignment="fixed",
    )
    compact_run = _compact_decoder_logs(summary, run)
    artifacts = postprocess_raw_posteriors(
        raw,
        directory / "sample_manifest.tsv",
        decoded,
        unit_record=unit,
        focal_position_0based=FOCAL_POSITION_BP,
        remove_raw_after_success=bool(payload["remove_raw_after_success"]),
        auxiliary_output_paths=auxiliary_outputs,
    )
    return {
        "unit_id": unit["unit_id"],
        "status": "complete",
        "error": "",
        "completion": str(artifacts.completion),
        "decoder_sha256": sha256_file(payload["decoder_bin"]),
        "decode_seconds": float(compact_run["decode_seconds"]),
    }


def _decode_task(payload: Mapping[str, Any]) -> dict[str, Any]:
    from .focused_selection_simulation import exclusive_unit_lock

    unit = dict(payload["unit"])
    directory = Path(payload["campaign_dir"]) / "work" / unit["unit_id"]
    with exclusive_unit_lock(
        directory / "decoded" / ".decode.lock",
        unit_id=str(unit["unit_id"]),
        phase="decode",
    ):
        return _decode_task_unlocked(payload)


def _decode_units(
    units: pd.DataFrame,
    campaign_dir: Path,
    decoder_bin: Path,
    *,
    workers: int,
    threads: int,
    remove_raw_after_success: bool,
) -> pd.DataFrame:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from .focused_selection_simulation import load_cell_demography

    decoder_bin = decoder_bin.resolve()
    if not decoder_bin.is_file():
        raise ValueError(f"Gamma-SMC decoder is absent: {decoder_bin}")
    decoder_sha256 = _sha256_file(decoder_bin)
    repo_root = campaign_dir.parent.resolve()
    decode_sources = _decode_source_hashes(repo_root)
    decoder_path = _portable_decoder_path(decoder_bin, repo_root)

    records = []
    for row in units.to_dict(orient="records"):
        unit_dir = campaign_dir / "work" / str(row["unit_id"])
        completion_path = unit_dir / "simulation_complete.json"
        if not completion_path.is_file():
            records.append({**row, "present_ne": np.nan})
            continue
        try:
            simulation_completion = json.loads(
                completion_path.read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(
                f"simulation completion is unreadable: {row['unit_id']}"
            ) from error
        if simulation_completion.get("status") != "complete" or (
            simulation_completion.get("contract_sha256")
            != _canonical_sha256(simulation_completion.get("contract"))
        ):
            raise ValueError(
                f"simulation completion contract is invalid: {row['unit_id']}"
            )
        required_inputs = {
            "tree": unit_dir / "simulation.trees",
            "pair_table": unit_dir / "sample_manifest.tsv",
            "overall_pairs": unit_dir / "pairs" / "overall.pairs.tsv",
        }
        for label, path in required_inputs.items():
            expected = (
                simulation_completion.get("outputs", {}).get(label, {}).get("sha256")
            )
            if not path.is_file() or _sha256_file(path) != expected:
                raise ValueError(
                    f"simulation input checksum failed for {row['unit_id']}: {label}"
                )
        _, _, present_ne = load_cell_demography(row, campaign_dir.parent)
        records.append(
            {
                **row,
                "present_ne": present_ne,
                "decode_provenance": {
                    "decoder_sha256": decoder_sha256,
                    "decoder_path": decoder_path,
                    "simulation_completion_sha256": _sha256_file(completion_path),
                    "tree_sha256": _sha256_file(required_inputs["tree"]),
                    "pair_table_sha256": _sha256_file(required_inputs["pair_table"]),
                    "overall_pairs_sha256": _sha256_file(
                        required_inputs["overall_pairs"]
                    ),
                    "implementation_sources": decode_sources,
                    "settings": _decode_settings(present_ne, threads=threads),
                },
            }
        )
    tasks = [
        {
            "unit": row,
            "campaign_dir": str(campaign_dir),
            "decoder_bin": str(decoder_bin),
            "threads": int(threads),
            "remove_raw_after_success": remove_raw_after_success,
        }
        for row in records
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=int(workers)) as pool:
        futures = {pool.submit(_decode_task, task): task for task in tasks}
        for future in as_completed(futures):
            try:
                rows.append(future.result())
            except Exception as error:  # noqa: BLE001 - worker status captures failures
                task = futures[future]
                rows.append(
                    {
                        "unit_id": str(task["unit"]["unit_id"]),
                        "status": "failed",
                        "error": str(error),
                    }
                )
    frame = pd.DataFrame(rows).sort_values("unit_id")
    frame.to_csv(campaign_dir / "results" / "decode_status.tsv", sep="\t", index=False)
    return frame


def _validated_output_path(
    root: Path,
    outputs: Mapping[str, Any],
    label: str,
    relative: Path,
) -> Path:
    record = outputs.get(label)
    path = root / relative
    if (
        not isinstance(record, Mapping)
        or str(record.get("path", "")).replace("\\", "/") != relative.as_posix()
        or not path.is_file()
        or _sha256_file(path) != record.get("sha256")
    ):
        raise ValueError(f"completion output is invalid: {label}")
    if "size_bytes" in record and path.stat().st_size != int(record["size_bytes"]):
        raise ValueError(f"completion output size is invalid: {label}")
    return path


def _update_spatial_accumulator(
    accumulator: dict[tuple[Any, ...], dict[str, Any]],
    frame: pd.DataFrame,
    unit: Mapping[str, Any],
) -> None:
    required = {
        "unit_id",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "source",
        "genotype_class",
        "threshold_years",
        "position_0based",
        "mean_p_tmrca_lt_threshold",
        "mean_tmrca_generations",
    }
    missing = sorted(required.difference(frame.columns))
    if missing or frame.empty:
        detail = ", ".join(missing) if missing else "empty table"
        raise ValueError(f"spatial profile is invalid for {unit['unit_id']}: {detail}")
    for column in ("unit_id", "demography_id", "simulation_class"):
        if set(frame[column].astype(str)) != {str(unit[column])}:
            raise ValueError(f"spatial profile {column} disagrees with the plan")
    for column in ("selection_coefficient", "target_allele_frequency"):
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or not np.allclose(
            values, float(unit[column]), rtol=0, atol=1e-12
        ):
            raise ValueError(f"spatial profile {column} disagrees with the plan")
    if set(frame["genotype_class"].astype(str)) != {
        "overall",
        "hom_ref",
        "heterozygous",
        "hom_alt",
    }:
        raise ValueError("spatial profile lacks a focal genotype class")
    if set(frame["threshold_years"].astype(float)) != set(TMRCA_THRESHOLDS_YEARS):
        raise ValueError("spatial profile lacks a prespecified threshold")
    key_columns = ["genotype_class", "threshold_years", "position_0based"]
    if frame.duplicated(key_columns).any():
        raise ValueError(f"spatial profile contains duplicate rows: {unit['unit_id']}")
    source_values = set(frame["source"].astype(str))
    if len(source_values) != 1:
        raise ValueError(f"spatial profile source changes within {unit['unit_id']}")
    source = next(iter(source_values))
    for (genotype_class, threshold), group in frame.groupby(
        ["genotype_class", "threshold_years"], sort=True
    ):
        ordered = group.sort_values("position_0based")
        raw_positions = pd.to_numeric(
            ordered["position_0based"], errors="coerce"
        ).to_numpy(dtype=float)
        positions = raw_positions.astype(np.int64)
        probability = pd.to_numeric(
            ordered["mean_p_tmrca_lt_threshold"], errors="coerce"
        ).to_numpy(dtype=float)
        tmrca = pd.to_numeric(
            ordered["mean_tmrca_generations"], errors="coerce"
        ).to_numpy(dtype=float)
        if (
            not np.isfinite(raw_positions).all()
            or not np.array_equal(raw_positions, positions)
            or np.any(positions < 0)
            or np.any(np.diff(positions) <= 0)
            or not np.isfinite(probability).all()
            or np.any((probability < 0) | (probability > 1))
            or not np.isfinite(tmrca).all()
            or np.any(tmrca <= 0)
        ):
            raise ValueError(f"spatial profile values are invalid: {unit['unit_id']}")
        key = (
            str(source),
            str(unit["demography_id"]),
            str(unit["simulation_class"]),
            float(unit["selection_coefficient"]),
            float(unit["target_allele_frequency"]),
            str(genotype_class),
            float(threshold),
        )
        state = accumulator.get(key)
        if state is None:
            accumulator[key] = {
                "positions": positions,
                "probability_sum": probability.copy(),
                "probability_sumsq": probability**2,
                "tmrca_sum": tmrca.copy(),
                "tmrca_sumsq": tmrca**2,
                "n_units": 1,
            }
        else:
            if not np.array_equal(state["positions"], positions):
                raise ValueError(f"spatial output grid changes within cell: {key}")
            state["probability_sum"] += probability
            state["probability_sumsq"] += probability**2
            state["tmrca_sum"] += tmrca
            state["tmrca_sumsq"] += tmrca**2
            state["n_units"] += 1


def _spatial_accumulator_frame(
    accumulator: Mapping[tuple[Any, ...], Mapping[str, Any]],
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    columns = [
        "source",
        "demography_id",
        "simulation_class",
        "selection_coefficient",
        "target_allele_frequency",
        "genotype_class",
        "threshold_years",
    ]
    for key, state in sorted(accumulator.items()):
        n_units = int(state["n_units"])
        probability_mean = state["probability_sum"] / n_units
        tmrca_mean = state["tmrca_sum"] / n_units
        if n_units > 1:
            probability_variance = np.maximum(
                (state["probability_sumsq"] - n_units * probability_mean**2)
                / (n_units - 1),
                0,
            )
            tmrca_variance = np.maximum(
                (state["tmrca_sumsq"] - n_units * tmrca_mean**2) / (n_units - 1),
                0,
            )
        else:
            probability_variance = np.zeros_like(probability_mean)
            tmrca_variance = np.zeros_like(tmrca_mean)
        prefix = dict(zip(columns, key, strict=True))
        frames.append(
            pd.DataFrame(
                {
                    **{label: value for label, value in prefix.items()},
                    "position_0based": state["positions"],
                    "n_units": n_units,
                    "mean_p_tmrca_lt_threshold": probability_mean,
                    "sd_p_tmrca_lt_threshold": np.sqrt(probability_variance),
                    "sem_p_tmrca_lt_threshold": np.sqrt(probability_variance / n_units),
                    "mean_tmrca_generations": tmrca_mean,
                    "sd_tmrca_generations": np.sqrt(tmrca_variance),
                    "sem_tmrca_generations": np.sqrt(tmrca_variance / n_units),
                }
            )
        )
    if not frames:
        raise ValueError("no validated spatial profiles are available")
    return pd.concat(frames, ignore_index=True)


def _aggregate_compact_decode(
    campaign_dir: Path, units: pd.DataFrame
) -> tuple[Path, Path, Path, Path, Path]:
    from .focused_selection_decode import POSTPROCESSOR_SOURCE
    from .focused_selection_decode import SCHEMA_VERSION as DECODE_SCHEMA_VERSION
    from .focused_selection_simulation import (
        SCHEMA_VERSION as SIMULATION_SCHEMA_VERSION,
    )
    from .focused_selection_simulation import (
        _simulate_unit_unlocked,
        load_cell_demography,
    )

    repo_root = campaign_dir.parent.resolve()
    current_sources = _decode_source_hashes(repo_root)
    class_frames = []
    pair_frames = []
    truth_frames = []
    gamma_spatial: dict[tuple[Any, ...], dict[str, Any]] = {}
    truth_spatial: dict[tuple[Any, ...], dict[str, Any]] = {}
    missing: list[str] = []
    for unit in units.to_dict(orient="records"):
        unit_id = str(unit["unit_id"])
        unit_root = campaign_dir / "work" / unit_id
        directory = unit_root / "decoded"
        completion_path = directory / "completion.json"
        simulation_completion_path = unit_root / "simulation_complete.json"
        if not completion_path.is_file() or not simulation_completion_path.is_file():
            missing.append(str(unit_id))
            continue
        try:
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            simulation_completion = json.loads(
                simulation_completion_path.read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as error:
            raise ValueError(f"completion record is unreadable: {unit_id}") from error
        simulation_contract = simulation_completion.get("contract")
        if (
            simulation_completion.get("schema") != SIMULATION_SCHEMA_VERSION
            or simulation_completion.get("status") != "complete"
            or not isinstance(simulation_contract, Mapping)
            or simulation_completion.get("contract_sha256")
            != _canonical_sha256(simulation_contract)
            or _canonical_sha256(simulation_contract.get("unit"))
            != _canonical_sha256({key: unit[key] for key in EXECUTION_COLUMNS})
        ):
            raise ValueError(f"simulation completion contract is invalid: {unit_id}")
        simulation_implementation = simulation_contract.get("implementation")
        if not isinstance(simulation_implementation, Mapping):
            raise TypeError(f"simulation implementation contract is invalid: {unit_id}")
        slim_path: Path | None = None
        if str(unit["simulation_class"]) == "selected":
            slim_record = simulation_implementation.get("slim")
            if not isinstance(slim_record, Mapping) or not str(
                slim_record.get("path", "")
            ):
                raise ValueError(f"simulation SLiM contract is invalid: {unit_id}")
            recorded_slim_path = Path(str(slim_record["path"]))
            slim_path = (
                recorded_slim_path
                if recorded_slim_path.is_absolute()
                else repo_root / recorded_slim_path
            ).resolve()
        validated_simulation = _simulate_unit_unlocked(
            unit,
            repo_root,
            campaign_dir,
            slim_path=slim_path,
        )
        if not validated_simulation.cache_hit:
            raise ValueError(
                f"simulation was not validated as an immutable cache hit: {unit_id}"
            )
        simulation_outputs = simulation_completion.get("outputs")
        if not isinstance(simulation_outputs, Mapping):
            raise TypeError(f"simulation output manifest is invalid: {unit_id}")
        tree_path = _validated_output_path(
            unit_root, simulation_outputs, "tree", Path("simulation.trees")
        )
        pair_table_path = _validated_output_path(
            unit_root, simulation_outputs, "pair_table", Path("sample_manifest.tsv")
        )
        overall_pairs_path = _validated_output_path(
            unit_root,
            simulation_outputs,
            "overall_pairs",
            Path("pairs/overall.pairs.tsv"),
        )
        truth_profile_path = _validated_output_path(
            unit_root,
            simulation_outputs,
            "truth_profiles",
            Path("truth_profiles.tsv.gz"),
        )
        truth_path = _validated_output_path(
            unit_root,
            simulation_outputs,
            "truth_class_summaries",
            Path("truth_class_summaries.tsv"),
        )

        decode_contract = completion.get("contract")
        if (
            completion.get("schema") != DECODE_SCHEMA_VERSION
            or completion.get("status") != "complete"
            or not isinstance(decode_contract, Mapping)
            or completion.get("contract_sha256") != _canonical_sha256(decode_contract)
        ):
            raise ValueError(f"decode completion contract is invalid: {unit_id}")
        static = decode_contract.get("static")
        if not isinstance(static, Mapping):
            raise TypeError(f"decode static contract is invalid: {unit_id}")
        recorded_unit = static.get("unit_record")
        if not isinstance(recorded_unit, Mapping):
            raise TypeError(f"decode unit contract is invalid: {unit_id}")
        provenance = recorded_unit.get("decode_provenance")
        if not isinstance(provenance, Mapping):
            raise TypeError(f"decode provenance is absent: {unit_id}")
        decoder_record = Path(str(provenance.get("decoder_path", "")))
        decoder_path = (
            decoder_record
            if decoder_record.is_absolute()
            else repo_root / decoder_record
        ).resolve()
        if not decoder_path.is_file() or _sha256_file(decoder_path) != provenance.get(
            "decoder_sha256"
        ):
            raise ValueError(f"current decoder binary differs from decode: {unit_id}")
        if provenance.get("implementation_sources") != current_sources:
            raise ValueError(f"current decode sources differ from decode: {unit_id}")
        settings = provenance.get("settings")
        if not isinstance(settings, Mapping):
            raise TypeError(f"decode settings are absent: {unit_id}")
        try:
            threads = int(settings["threads"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f"decode thread provenance is invalid: {unit_id}"
            ) from error
        _, _, present_ne = load_cell_demography(unit, repo_root)
        expected_provenance = {
            "decoder_sha256": _sha256_file(decoder_path),
            "decoder_path": _portable_decoder_path(decoder_path, repo_root),
            "simulation_completion_sha256": _sha256_file(simulation_completion_path),
            "tree_sha256": _sha256_file(tree_path),
            "pair_table_sha256": _sha256_file(pair_table_path),
            "overall_pairs_sha256": _sha256_file(overall_pairs_path),
            "implementation_sources": current_sources,
            "settings": _decode_settings(present_ne, threads=threads),
        }
        expected_unit = {
            **unit,
            "present_ne": present_ne,
            "decode_provenance": expected_provenance,
        }
        if _canonical_sha256(recorded_unit) != _canonical_sha256(expected_unit):
            raise ValueError(
                f"decode provenance differs from current inputs: {unit_id}"
            )
        pair_contract = static.get("pair_table")
        if (
            not isinstance(pair_contract, Mapping)
            or Path(str(pair_contract.get("path", ""))).resolve()
            != pair_table_path.resolve()
            or pair_contract.get("sha256") != _sha256_file(pair_table_path)
        ):
            raise ValueError(f"decode pair-table contract is invalid: {unit_id}")
        expected_parameters = {
            "unscaled_mutation_rate_per_bp_per_generation": MUTATION_RATE,
            "generation_time_years": GENERATION_TIME_YEARS,
            "thresholds_years": [float(value) for value in TMRCA_THRESHOLDS_YEARS],
            "focal_position_0based": float(FOCAL_POSITION_BP),
        }
        if _canonical_sha256(static.get("parameters")) != _canonical_sha256(
            expected_parameters
        ):
            raise ValueError(f"decode postprocess parameters are invalid: {unit_id}")
        implementation = static.get("implementation")
        if not isinstance(implementation, Mapping) or implementation.get(
            "postprocessor_sha256"
        ) != _sha256_file(POSTPROCESSOR_SOURCE):
            raise ValueError(f"current postprocessor differs from decode: {unit_id}")
        expected_auxiliary = {
            label: str((directory / filename).resolve())
            for label, filename in DECODE_AUXILIARY_OUTPUTS.items()
        }
        if static.get("auxiliary_output_paths") != expected_auxiliary:
            raise ValueError(f"decode auxiliary contract is invalid: {unit_id}")
        raw_inputs = decode_contract.get("raw_inputs")
        if not isinstance(raw_inputs, Mapping):
            raise TypeError(f"decode raw-input manifest is invalid: {unit_id}")
        raw_path = directory / "posterior.zst"
        raw_meta_path = directory / "posterior.zst.meta"
        if (
            Path(str(raw_inputs.get("posterior_path", ""))).resolve()
            != raw_path.resolve()
            or Path(str(raw_inputs.get("metadata_path", ""))).resolve()
            != raw_meta_path.resolve()
        ):
            raise ValueError(f"decode raw-input paths are invalid: {unit_id}")
        if raw_path.is_file() != raw_meta_path.is_file():
            raise ValueError(f"only one raw decode input remains: {unit_id}")
        if raw_path.is_file() and (
            _sha256_file(raw_path) != raw_inputs.get("posterior_sha256")
            or _sha256_file(raw_meta_path) != raw_inputs.get("metadata_sha256")
        ):
            raise ValueError(f"decode raw-input checksum failed: {unit_id}")
        decode_paths = {
            "spatial_profiles": directory / "spatial_profiles.tsv.gz",
            "class_summaries": directory / "class_summaries.tsv",
            "pair_summaries": directory / "pair_summaries.tsv.gz",
            **{
                label: directory / filename
                for label, filename in DECODE_AUXILIARY_OUTPUTS.items()
            },
        }
        decode_outputs = completion.get("outputs")
        if not isinstance(decode_outputs, Mapping) or set(decode_outputs) != set(
            decode_paths
        ):
            raise ValueError(f"decode output manifest is invalid: {unit_id}")
        for label, path in decode_paths.items():
            record = decode_outputs.get(label)
            if (
                not isinstance(record, Mapping)
                or record.get("path") != path.name
                or not path.is_file()
                or _sha256_file(path) != record.get("sha256")
                or path.stat().st_size != int(record.get("size_bytes", -1))
            ):
                raise ValueError(f"decode output checksum failed: {unit_id}: {label}")
        class_frames.append(pd.read_csv(decode_paths["class_summaries"], sep="\t"))
        pair_frames.append(pd.read_csv(decode_paths["pair_summaries"], sep="\t"))
        truth_frames.append(pd.read_csv(truth_path, sep="\t"))
        _update_spatial_accumulator(
            gamma_spatial,
            pd.read_csv(decode_paths["spatial_profiles"], sep="\t"),
            unit,
        )
        _update_spatial_accumulator(
            truth_spatial, pd.read_csv(truth_profile_path, sep="\t"), unit
        )
    if missing:
        raise ValueError(
            f"{len(missing)} planned units lack validated simulation/decode completion; "
            f"first: {', '.join(missing[:5])}"
        )
    if not class_frames:
        raise ValueError("no validated decoded units are available for analysis")
    results = campaign_dir / "results"
    classes_path = results / "combined_class_summaries.tsv.gz"
    pairs_path = results / "combined_pair_summaries.tsv.gz"
    truth_path = results / "combined_truth_class_summaries.tsv.gz"
    spatial_path = results / "aggregated_spatial_profiles.tsv.gz"
    truth_spatial_path = results / "aggregated_truth_spatial_profiles.tsv.gz"
    _atomic_frame(classes_path, pd.concat(class_frames, ignore_index=True))
    _atomic_frame(pairs_path, pd.concat(pair_frames, ignore_index=True))
    _atomic_frame(truth_path, pd.concat(truth_frames, ignore_index=True))
    _atomic_frame(spatial_path, _spatial_accumulator_frame(gamma_spatial))
    _atomic_frame(truth_spatial_path, _spatial_accumulator_frame(truth_spatial))
    return classes_path, pairs_path, truth_path, spatial_path, truth_spatial_path


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    repo_root = args.repo_root.resolve()
    campaign_dir = (args.campaign_dir or repo_root / DEFAULT_CAMPAIGN_DIR).resolve()
    plan = FocusedCampaignPlan(
        repo_root=repo_root,
        campaign_dir=campaign_dir,
        selected_replicates=args.selected_replicates,
        neutral_replicates=args.neutral_replicates,
        base_seed=args.base_seed,
    )
    table_path, _ = write_plan(plan)
    units = pd.read_csv(table_path, sep="\t")
    selected_units = _filtered_units(units, args)
    if args.phase == "plan":
        return 0
    if args.phase == "status":
        status = _simulation_status(campaign_dir, units)
        print(status.groupby("status").size().to_string())
        return 0
    if args.phase in {"simulate", "all"}:
        from .focused_selection_simulation import simulate_units

        if selected_units.empty:
            raise ValueError("unit filters selected no simulations")
        status = simulate_units(
            selected_units,
            repo_root,
            campaign_dir,
            slim_path=args.slim_bin,
            workers=args.workers,
            neutral_max_attempts=args.neutral_max_attempts,
            selected_max_draws=args.selected_max_draws,
            selected_draw_timeout_seconds=(args.selected_draw_timeout_minutes * 60.0),
            selected_cumulative_timeout_seconds=(
                args.selected_cumulative_timeout_minutes * 60.0
            ),
            slim_scaling_factor=args.slim_scaling_factor,
            han_selection_calibration_path=(
                args.han_selection_calibration
                or campaign_dir / "calibration" / "han_selection_end_frozen.json"
            ),
        )
        status.to_csv(
            campaign_dir / "results" / "last_simulation_invocation.tsv",
            sep="\t",
            index=False,
        )
        _simulation_status(campaign_dir, units)
        if (status["status"] == "failed").any():
            return 1
    if args.phase in {"decode", "all"}:
        if args.decoder_bin is None:
            raise ValueError("decode requires --decoder-bin")
        if selected_units.empty:
            raise ValueError("unit filters selected no decodes")
        decode_status = _decode_units(
            selected_units if args.phase == "decode" else units,
            campaign_dir,
            args.decoder_bin.resolve(),
            workers=args.workers,
            threads=args.threads,
            remove_raw_after_success=args.remove_raw_after_success,
        )
        if not decode_status["status"].isin({"complete", "cached"}).all():
            return 1
    if args.phase in {"analyze", "all"}:
        from .focused_selection_analysis import (
            run_focused_analysis,
            run_truth_gamma_comparison,
        )

        classes, pairs, truth, spatial, truth_spatial = _aggregate_compact_decode(
            campaign_dir, units
        )
        run_focused_analysis(
            classes,
            campaign_dir / "results" / "analysis",
            pair_summaries_path=pairs,
            spatial_summaries_path=spatial,
            empirical_path=args.empirical_tsv,
        )
        run_focused_analysis(
            truth,
            campaign_dir / "results" / "analysis_tree_truth",
            spatial_summaries_path=truth_spatial,
            minimum_selected_per_cell=DEFAULT_SELECTED_REPLICATES,
            minimum_neutral_per_cell=DEFAULT_NEUTRAL_REPLICATES,
        )
        run_truth_gamma_comparison(
            truth,
            classes,
            campaign_dir / "results" / "truth_gamma_comparison",
        )
    if args.phase in {"report", "analyze", "all"}:
        from .focused_selection_report import generate_run_report

        print(generate_run_report(campaign_dir))
    return 0


__all__ = [
    "DEFAULT_AF_TARGETS",
    "DEFAULT_CAMPAIGN_DIR",
    "DEFAULT_HAN_SELECTION_CALIBRATION",
    "DEFAULT_NEUTRAL_ANCESTRY_BATCH_SIZE",
    "DEFAULT_NEUTRAL_REPLICATES",
    "DEFAULT_SELECTED_CUMULATIVE_TIMEOUT_MINUTES",
    "DEFAULT_SELECTED_DRAW_TIMEOUT_MINUTES",
    "DEFAULT_SELECTED_MAX_DRAWS",
    "DEFAULT_SELECTED_POOL_DIPLOIDS",
    "DEFAULT_SELECTED_REPLICATES",
    "DEFAULT_SELECTION_COEFFICIENTS",
    "DEFAULT_SLIM_SCALING_FACTOR",
    "DEMOGRAPHIES",
    "FocusedCampaignPlan",
    "build_execution_units",
    "build_parser",
    "build_study_design",
    "main",
    "write_plan",
]
