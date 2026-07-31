from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import bitmatrix as bitmatrix_reader


AUTOSOMES = tuple(range(1, 23))
POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "MID", "SAS")
SOFT_COLUMN = "mean_p_tmrca_lt_threshold"
TMRCA_COLUMN = "mean_tmrca_generations"
ANCESTRY_LABEL_COLUMN = "ancestry_pred_other"
ANCESTRY_ID_COLUMNS = ("research_id", "sample_id", "s", "person_id")
QC_ID_COLUMNS = ("s", "research_id", "sample_id", "sample_id.s", "person_id")
RELATEDNESS_ID_COLUMNS = (
    "sample_id.s",
    "sample_id",
    "s",
    "research_id",
    "person_id",
)


def threshold_suffix(years: float) -> str:
    return str(int(years)) if float(years).is_integer() else str(years)


def called_fraction_column(threshold_years: float) -> str:
    return f"frac_recent_{threshold_suffix(threshold_years)}"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def contract_sha256(contract: dict) -> str:
    encoded = json.dumps(
        contract, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_controlled_tsv(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame.columns = [str(column).lstrip("\ufeff").strip() for column in frame.columns]
    if len(set(frame.columns)) != len(frame.columns):
        raise ValueError(f"table has duplicate columns after normalization: {path}")
    return frame


def _id_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(
        f"{label} table has none of the supported sample-ID columns "
        f"{list(candidates)}; observed {list(frame.columns)}"
    )


def _clean_ids(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().loc[lambda series: series.ne("")]


def build_workbench_sample_list(
    *,
    ancestry_path: str | Path,
    qc_exclusions_path: str | Path,
    relatedness_exclusions_path: str | Path,
    bcf_samples_path: str | Path,
    population: str,
    output_path: str | Path,
    audit_path: str | Path,
) -> dict:
    population = population.upper()
    if population not in POPULATIONS:
        raise ValueError(f"unsupported population: {population}")
    ancestry = _read_controlled_tsv(ancestry_path)
    qc = _read_controlled_tsv(qc_exclusions_path)
    relatedness = _read_controlled_tsv(relatedness_exclusions_path)
    if ANCESTRY_LABEL_COLUMN not in ancestry.columns:
        raise ValueError(
            f"ancestry table is missing required column {ANCESTRY_LABEL_COLUMN}"
        )
    ancestry_id_column = _id_column(ancestry, ANCESTRY_ID_COLUMNS, "ancestry")
    qc_id_column = _id_column(qc, QC_ID_COLUMNS, "QC exclusion")
    relatedness_id_column = _id_column(
        relatedness, RELATEDNESS_ID_COLUMNS, "relatedness exclusion"
    )

    ancestry_pairs = ancestry[[ancestry_id_column, ANCESTRY_LABEL_COLUMN]].copy()
    ancestry_pairs.columns = ["sample_id", "population"]
    ancestry_pairs["sample_id"] = ancestry_pairs["sample_id"].astype(str).str.strip()
    ancestry_pairs["population"] = (
        ancestry_pairs["population"].astype(str).str.strip().str.upper()
    )
    ancestry_pairs = ancestry_pairs.loc[ancestry_pairs["sample_id"].ne("")]
    conflicts = (
        ancestry_pairs.groupby("sample_id", sort=False)["population"].nunique().gt(1)
    )
    if conflicts.any():
        examples = conflicts.index[conflicts].tolist()[:5]
        raise ValueError(f"ancestry table assigns conflicting labels to: {examples}")
    ancestry_pairs = ancestry_pairs.drop_duplicates("sample_id", keep="first")
    ancestry_by_id = dict(
        zip(ancestry_pairs["sample_id"], ancestry_pairs["population"])
    )

    qc_ids = set(_clean_ids(qc[qc_id_column]))
    relatedness_ids = set(_clean_ids(relatedness[relatedness_id_column]))
    excluded_ids = qc_ids | relatedness_ids

    bcf_samples_file = Path(bcf_samples_path)
    if not bcf_samples_file.is_file():
        raise ValueError(f"BCF sample list is absent: {bcf_samples_file}")
    bcf_samples = [
        line.strip()
        for line in bcf_samples_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not bcf_samples:
        raise ValueError("BCF contains no samples")
    if len(set(bcf_samples)) != len(bcf_samples):
        raise ValueError("BCF sample list contains duplicate sample IDs")

    selected = [
        sample_id
        for sample_id in bcf_samples
        if ancestry_by_id.get(sample_id) == population
        and sample_id not in excluded_ids
    ]
    if not selected:
        raise ValueError(
            f"no {population} samples remain after ancestry, BCF, QC, and relatedness filters"
        )

    target_ancestry = {
        sample_id for sample_id, label in ancestry_by_id.items() if label == population
    }
    bcf_sample_set = set(bcf_samples)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    temporary.write_text("\n".join(selected) + "\n", encoding="utf-8")
    os.replace(temporary, output)

    source_paths = {
        "ancestry": Path(ancestry_path),
        "qc_exclusions": Path(qc_exclusions_path),
        "relatedness_exclusions": Path(relatedness_exclusions_path),
        "bcf_samples": bcf_samples_file,
    }
    audit = {
        "schema_version": 1,
        "population": population,
        "ancestry_label_column": ANCESTRY_LABEL_COLUMN,
        "id_columns": {
            "ancestry": ancestry_id_column,
            "qc_exclusions": qc_id_column,
            "relatedness_exclusions": relatedness_id_column,
        },
        "sources": {
            name: {
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for name, path in source_paths.items()
        },
        "counts": {
            "ancestry_rows": int(len(ancestry)),
            "ancestry_unique_ids": int(len(ancestry_by_id)),
            "target_ancestry_ids": int(len(target_ancestry)),
            "qc_exclusion_ids": int(len(qc_ids)),
            "relatedness_exclusion_ids": int(len(relatedness_ids)),
            "union_exclusion_ids": int(len(excluded_ids)),
            "bcf_samples": int(len(bcf_samples)),
            "bcf_samples_absent_from_ancestry": int(
                len(bcf_sample_set.difference(ancestry_by_id))
            ),
            "target_ancestry_absent_from_bcf": int(
                len(target_ancestry.difference(bcf_sample_set))
            ),
            "target_bcf_samples_before_exclusions": int(
                len(target_ancestry.intersection(bcf_sample_set))
            ),
            "target_bcf_samples_excluded_by_qc": int(
                len(target_ancestry & bcf_sample_set & qc_ids)
            ),
            "target_bcf_samples_excluded_by_relatedness": int(
                len(target_ancestry & bcf_sample_set & relatedness_ids)
            ),
            "selected_samples": int(len(selected)),
        },
        "sample_list": {
            "path": str(output.resolve()),
            "sha256": sha256_file(output),
            "n_samples": int(len(selected)),
            "order": "source BCF header order",
        },
    }
    audit_destination = Path(audit_path)
    audit_destination.parent.mkdir(parents=True, exist_ok=True)
    audit_temporary = audit_destination.with_name(
        audit_destination.name + f".tmp.{os.getpid()}"
    )
    audit_temporary.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    os.replace(audit_temporary, audit_destination)
    return audit


def build_workbench_callable_mask(
    *,
    hardmask_path: str | Path,
    contig: str,
    sequence_length: int,
    output_path: str | Path,
    audit_path: str | Path,
) -> dict:
    if not contig or sequence_length < 1:
        raise ValueError("contig and positive sequence length are required")
    hardmask = Path(hardmask_path)
    if not hardmask.is_file() or hardmask.stat().st_size == 0:
        raise ValueError(f"hard-mask BED is absent or empty: {hardmask}")
    aliases = {contig}
    aliases.add(contig[3:] if contig.lower().startswith("chr") else f"chr{contig}")
    intervals = []
    matching_records = 0
    clipped_records = 0
    with hardmask.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("track ") or stripped.startswith("browser "):
                continue
            fields = stripped.split()
            if fields[0] not in aliases:
                continue
            matching_records += 1
            if len(fields) < 3:
                raise ValueError(
                    f"hard-mask BED line {line_number} has fewer than three fields"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"hard-mask BED line {line_number} has non-integer coordinates"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(f"hard-mask BED line {line_number} has invalid interval")
            clipped_start = max(0, start)
            clipped_end = min(sequence_length, end)
            if clipped_start != start or clipped_end != end:
                clipped_records += 1
            if clipped_start < clipped_end:
                intervals.append((clipped_start, clipped_end))

    if matching_records == 0:
        raise ValueError(
            f"hard-mask BED has no records for {contig} (accepted aliases: {sorted(aliases)})"
        )

    intervals.sort()
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    callable_intervals = []
    cursor = 0
    for start, end in merged:
        if cursor < start:
            callable_intervals.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < sequence_length:
        callable_intervals.append((cursor, sequence_length))
    if not callable_intervals:
        raise ValueError(f"hard mask excludes the entire {contig} sequence")

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    temporary.write_text(
        "".join(f"{contig}\t{start}\t{end}\n" for start, end in callable_intervals),
        encoding="utf-8",
    )
    os.replace(temporary, output)
    excluded_bases = int(sum(end - start for start, end in merged))
    callable_bases = int(sequence_length - excluded_bases)
    audit = {
        "schema_version": 1,
        "source_semantics": "excluded_intervals",
        "decoder_mask_semantics": "included_intervals",
        "contig": contig,
        "accepted_contig_aliases": sorted(aliases),
        "sequence_length": int(sequence_length),
        "source": {
            "path": str(hardmask.resolve()),
            "sha256": sha256_file(hardmask),
            "size_bytes": hardmask.stat().st_size,
        },
        "counts": {
            "matching_source_records": int(matching_records),
            "clipped_source_records": int(clipped_records),
            "merged_excluded_intervals": int(len(merged)),
            "callable_intervals": int(len(callable_intervals)),
            "excluded_bases": excluded_bases,
            "callable_bases": callable_bases,
            "callable_fraction": float(callable_bases / sequence_length),
        },
        "callable_mask": {
            "path": str(output.resolve()),
            "sha256": sha256_file(output),
            "size_bytes": output.stat().st_size,
        },
    }
    audit_destination = Path(audit_path)
    audit_destination.parent.mkdir(parents=True, exist_ok=True)
    audit_temporary = audit_destination.with_name(
        audit_destination.name + f".tmp.{os.getpid()}"
    )
    audit_temporary.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    os.replace(audit_temporary, audit_destination)
    return audit


def build_workbench_contract(
    *,
    population: str,
    chromosome: int,
    input_uri: str,
    input_fingerprint: str,
    local_input: str | Path,
    index_uri: str,
    index_fingerprint: str,
    local_index: str | Path,
    output_summary: str | Path,
    pairs_manifest: str | Path,
    sample_list: str | Path,
    sample_audit: str | Path,
    ancestry_uri: str,
    ancestry_fingerprint: str,
    local_ancestry: str | Path,
    qc_exclusions_uri: str,
    qc_exclusions_fingerprint: str,
    local_qc_exclusions: str | Path,
    relatedness_exclusions_uri: str,
    relatedness_exclusions_fingerprint: str,
    local_relatedness_exclusions: str | Path,
    mask_uri: str | None,
    mask_fingerprint: str | None,
    local_mask_source: str | Path | None,
    local_mask: str | Path | None,
    mask_audit: str | Path | None,
    theta: float,
    rho_over_theta: float,
    mutation_rate: float,
    generation_time: float,
    threshold_years: float,
    recent_call: str,
    stride: int,
    cache_size: int,
    threads: int,
    pair_block: int,
    exp10: str,
    backward_alignment: str,
    code_commit: str,
    bitmatrix: str | Path | None = None,
    candidate_regions: str | Path | None = None,
    candidate_manifest: str | Path | None = None,
    n_random_pairs: int = 0,
    pairs_seed: int = 1729,
    exclude_within: bool = False,
    signal_fraction: float = 0.05,
    merge_gap: int = 20_000,
    profile_half_width: int = 500_000,
    variant_half_width: int = 100_000,
    minimum_genotype_pairs: int = 20,
) -> dict:
    population = population.upper()
    if population not in POPULATIONS:
        raise ValueError(f"unsupported population: {population}")
    if chromosome not in AUTOSOMES:
        raise ValueError(f"chromosome must be 1-22, observed {chromosome}")
    if recent_call not in {"mean", "median"}:
        raise ValueError(f"unsupported recent-call rule: {recent_call}")
    if stride < 1 or cache_size < 1 or threads < 1 or pair_block < 1:
        raise ValueError("stride, cache size, threads, and pair block must be positive")
    if n_random_pairs < 0 or pairs_seed < 0:
        raise ValueError("random-pair count and seed must be nonnegative")
    if not 0 <= signal_fraction <= 1:
        raise ValueError("candidate signal fraction must be in [0, 1]")
    if min(merge_gap, profile_half_width, variant_half_width) < 0:
        raise ValueError("candidate gap/profile widths must be nonnegative")
    if minimum_genotype_pairs < 1:
        raise ValueError("minimum genotype-pair count must be positive")
    candidate_fields = (bitmatrix, candidate_regions, candidate_manifest)
    if any(value is not None for value in candidate_fields) and not all(
        value is not None for value in candidate_fields
    ):
        raise ValueError(
            "bit matrix, candidate regions, and candidate manifest must be supplied together"
        )
    mask_fields = (
        mask_uri,
        mask_fingerprint,
        local_mask_source,
        local_mask,
        mask_audit,
    )
    if any(value is not None for value in mask_fields) and not all(
        value is not None for value in mask_fields
    ):
        raise ValueError("mask URI, fingerprint, and local path must be supplied together")
    return {
        "schema_version": 1,
        "population": population,
        "chromosome": int(chromosome),
        "input": {
            "uri": input_uri,
            "cloud_fingerprint": input_fingerprint,
            "local_path": str(Path(local_input).resolve()),
            "index": {
                "uri": index_uri,
                "cloud_fingerprint": index_fingerprint,
                "local_path": str(Path(local_index).resolve()),
            },
        },
        "mask": (
            {
                "uri": mask_uri,
                "cloud_fingerprint": mask_fingerprint,
                "source_local_path": str(Path(local_mask_source).resolve()),
                "local_path": str(Path(local_mask).resolve()),
                "audit": str(Path(mask_audit).resolve()),
                "source_semantics": "excluded_intervals_complemented",
            }
            if mask_uri is not None
            else None
        ),
        "output_summary": str(Path(output_summary).resolve()),
        "pairs_manifest": str(Path(pairs_manifest).resolve()),
        "analysis_outputs": (
            {
                "bitmatrix": str(Path(bitmatrix).resolve()),
                "candidate_regions": str(Path(candidate_regions).resolve()),
                "candidate_manifest": str(Path(candidate_manifest).resolve()),
            }
            if bitmatrix is not None
            else None
        ),
        "sample_selection": {
            "population_label_column": ANCESTRY_LABEL_COLUMN,
            "sample_list": str(Path(sample_list).resolve()),
            "audit": str(Path(sample_audit).resolve()),
            "sources": {
                "ancestry": {
                    "uri": ancestry_uri,
                    "cloud_fingerprint": ancestry_fingerprint,
                    "local_path": str(Path(local_ancestry).resolve()),
                },
                "qc_exclusions": {
                    "uri": qc_exclusions_uri,
                    "cloud_fingerprint": qc_exclusions_fingerprint,
                    "local_path": str(Path(local_qc_exclusions).resolve()),
                },
                "relatedness_exclusions": {
                    "uri": relatedness_exclusions_uri,
                    "cloud_fingerprint": relatedness_exclusions_fingerprint,
                    "local_path": str(Path(local_relatedness_exclusions).resolve()),
                },
            },
        },
        "decoder": {
            "theta": float(theta),
            "rho_over_theta": float(rho_over_theta),
            "mutation_rate": float(mutation_rate),
            "generation_time": float(generation_time),
            "threshold_years": float(threshold_years),
            "recent_call": recent_call,
            "stride_bp": int(stride),
            "cache_size_bp": int(cache_size),
            "threads": int(threads),
            "pair_block": int(pair_block),
            "exp10": exp10,
            "backward_alignment": backward_alignment,
            "output_at_hets": False,
            "pair_selection": (
                "random_haplotype_pairs" if n_random_pairs > 0 else "within_individual"
            ),
            "n_random_pairs": int(n_random_pairs),
            "pairs_seed": int(pairs_seed),
            "exclude_within": bool(exclude_within),
        },
        "candidate_screen": {
            "called_fraction_strictly_greater_than": float(signal_fraction),
            "merge_gap_bp": int(merge_gap),
            "profile_half_width_bp": int(profile_half_width),
            "variant_half_width_bp": int(variant_half_width),
            "minimum_hom_genotype_pairs": int(minimum_genotype_pairs),
            "representative_variant_statistic": "raw_tmrca_r_squared",
        },
        "code_commit": code_commit,
    }


def validate_summary(
    summary_path: str | Path,
    *,
    threshold_years: float,
) -> tuple[pd.DataFrame, dict]:
    path = Path(summary_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"summary is absent or empty: {path}")
    frame = pd.read_csv(path, sep="\t")
    called_column = called_fraction_column(threshold_years)
    required = {
        "position_0based",
        "position_1based",
        "n_pairs",
        SOFT_COLUMN,
        TMRCA_COLUMN,
        called_column,
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"summary is missing columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("summary contains no positions")
    positions = frame["position_0based"].to_numpy(dtype=float)
    if not np.all(np.isfinite(positions)) or np.any(positions < 0):
        raise ValueError("summary positions must be finite and nonnegative")
    if np.any(np.diff(positions) <= 0):
        raise ValueError("summary positions must be unique and strictly increasing")
    one_based = frame["position_1based"].to_numpy(dtype=float)
    if not np.allclose(one_based, positions + 1, rtol=0, atol=1e-8):
        raise ValueError("position_1based is not position_0based + 1")
    n_pairs = frame["n_pairs"].to_numpy(dtype=float)
    if np.any(~np.isfinite(n_pairs)) or np.any(n_pairs < 0):
        raise ValueError("n_pairs must be finite and nonnegative")
    if not np.any(n_pairs > 0):
        raise ValueError("no output position has a usable pair")
    for column in (SOFT_COLUMN, called_column):
        values = frame[column].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if len(finite) == 0 or np.any((finite < 0) | (finite > 1)):
            raise ValueError(f"{column} must contain probabilities in [0, 1]")
    tmrca = frame[TMRCA_COLUMN].to_numpy(dtype=float)
    finite_tmrca = tmrca[np.isfinite(tmrca)]
    if len(finite_tmrca) == 0 or np.any(finite_tmrca < 0):
        raise ValueError("mean TMRCA must contain finite nonnegative values")
    metrics = {
        "n_output_positions": int(len(frame)),
        "first_position_0based": float(positions[0]),
        "last_position_0based": float(positions[-1]),
        "n_pairs_min": int(np.nanmin(n_pairs)),
        "n_pairs_median": float(np.nanmedian(n_pairs)),
        "n_pairs_max": int(np.nanmax(n_pairs)),
    }
    return frame, metrics


def _command_value(command: list[str], flag: str) -> str:
    try:
        index = command.index(flag)
    except ValueError as error:
        raise ValueError(f"run command is missing {flag}") from error
    if index + 1 >= len(command):
        raise ValueError(f"run command has no value after {flag}")
    return command[index + 1]


def _assert_float(command: list[str], flag: str, expected: float) -> None:
    observed = float(_command_value(command, flag))
    if not np.isclose(observed, expected, rtol=0, atol=max(1e-15, abs(expected) * 1e-12)):
        raise ValueError(f"run command {flag}={observed}; expected {expected}")


def validate_run_json(run_json_path: str | Path, contract: dict) -> dict:
    path = Path(run_json_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"run metadata is absent or empty: {path}")
    run = json.loads(path.read_text(encoding="utf-8"))
    command = run.get("command")
    if not isinstance(command, list) or not all(isinstance(item, str) for item in command):
        raise ValueError("run metadata command must be a string list")
    decoder = contract["decoder"]
    expected_input = contract["input"]["local_path"]
    if str(Path(_command_value(command, "--input")).resolve()) != expected_input:
        raise ValueError("run command input does not match the contract")
    expected_output = contract["output_summary"]
    if str(Path(_command_value(command, "--recent_summary")).resolve()) != expected_output:
        raise ValueError("run command output does not match the contract")
    expected_manifest = contract["pairs_manifest"]
    if str(Path(_command_value(command, "--pairs_manifest")).resolve()) != expected_manifest:
        raise ValueError("run command pair manifest does not match the contract")
    expected_samples = contract["sample_selection"]["sample_list"]
    if str(Path(_command_value(command, "--samples")).resolve()) != expected_samples:
        raise ValueError("run command sample list does not match the contract")
    expected_mask = contract["mask"]
    if expected_mask is None:
        if "--mask" in command:
            raise ValueError("run command unexpectedly supplied a mask")
    elif str(Path(_command_value(command, "--mask")).resolve()) != str(
        Path(expected_mask["local_path"]).resolve()
    ):
        raise ValueError("run command mask does not match the contract")
    _assert_float(command, "--scaled_mutation_rate", decoder["theta"])
    _assert_float(
        command,
        "--recombination_to_mutation_ratio",
        decoder["rho_over_theta"],
    )
    _assert_float(command, "--unscaled_mutation_rate", decoder["mutation_rate"])
    _assert_float(command, "--generation_time", decoder["generation_time"])
    thresholds = _command_value(command, "--recent_threshold_years").split(",")
    if len(thresholds) != 1 or not np.isclose(
        float(thresholds[0]), decoder["threshold_years"], rtol=0, atol=1e-12
    ):
        raise ValueError("run command threshold does not match the contract")
    exact_values = {
        "--recent_call": str(decoder["recent_call"]),
        "--output_at_stride": str(decoder["stride_bp"]),
        "--cache_size": str(decoder["cache_size_bp"]),
        "--threads": str(decoder["threads"]),
        "--pair_block": str(decoder["pair_block"]),
        "--exp10": str(decoder["exp10"]),
        "--backward_alignment": str(decoder["backward_alignment"]),
    }
    for flag, expected in exact_values.items():
        if _command_value(command, flag) != expected:
            raise ValueError(f"run command {flag} does not match the contract")
    if "--output_at_hets=false" not in command:
        raise ValueError("run command did not disable heterozygous-site output")
    if decoder["pair_selection"] == "within_individual":
        if "--only_within" not in command:
            raise ValueError("run command did not select within-individual pairs")
        forbidden_selectors = {"--n_random_pairs", "--pairs_file", "--exclude_within"}
        observed_forbidden = forbidden_selectors.intersection(command)
        if observed_forbidden:
            raise ValueError(
                f"run command mixed within-individual pairs with {sorted(observed_forbidden)}"
            )
    else:
        if "--only_within" in command or "--pairs_file" in command:
            raise ValueError("random scan command mixed incompatible pair selectors")
        if _command_value(command, "--n_random_pairs") != str(
            decoder["n_random_pairs"]
        ):
            raise ValueError("run command random-pair count does not match the contract")
        if _command_value(command, "--pairs_seed") != str(decoder["pairs_seed"]):
            raise ValueError("run command pair seed does not match the contract")
        if decoder["exclude_within"] != ("--exclude_within" in command):
            raise ValueError("run command within-pair exclusion does not match the contract")
    outputs = contract.get("analysis_outputs")
    if outputs is not None:
        if str(Path(_command_value(command, "--recent_bitmatrix")).resolve()) != outputs[
            "bitmatrix"
        ]:
            raise ValueError("run command bit matrix does not match the contract")
    if int(run.get("stride_bp", -1)) != decoder["stride_bp"]:
        raise ValueError("run metadata stride does not match the contract")
    if int(run.get("cache_size_bp", -1)) != decoder["cache_size_bp"]:
        raise ValueError("run metadata cache size does not match the contract")
    run_manifest = run.get("pairs_manifest")
    if run_manifest is None or str(Path(run_manifest).resolve()) != expected_manifest:
        raise ValueError("run metadata pair manifest does not match the contract")
    return run


def validate_pairs_manifest(
    path: str | Path,
    *,
    expected_samples: list[str] | None = None,
    pair_mode: str = "within_individual",
    expected_count: int | None = None,
) -> dict:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"pair manifest is absent or empty: {path}")
    pairs = set()
    pair_index = 0
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) < 2:
            raise ValueError(f"pair manifest line {line_number} has fewer than two fields")
        try:
            haplotype_i, haplotype_j = int(fields[0]), int(fields[1])
        except ValueError as error:
            raise ValueError(
                f"pair manifest line {line_number} has non-integer indices"
            ) from error
        pair = (haplotype_i, haplotype_j)
        if haplotype_i < 0 or haplotype_j <= haplotype_i:
            raise ValueError(f"pair manifest line {line_number} is not an ordered pair")
        if pair_mode == "within_individual" and (
            haplotype_i % 2 != 0 or haplotype_j != haplotype_i + 1
        ):
            raise ValueError(
                f"pair manifest line {line_number} is not a within-individual pair"
            )
        if pair in pairs:
            raise ValueError(f"pair manifest repeats pair {pair}")
        if expected_samples is not None:
            if haplotype_j >= 2 * len(expected_samples):
                raise ValueError("pair manifest references outside the sample list")
            if pair_mode == "within_individual":
                if pair_index >= len(expected_samples):
                    raise ValueError("pair manifest has more pairs than the sample list")
                expected_pair = (2 * pair_index, 2 * pair_index + 1)
                if pair != expected_pair:
                    raise ValueError(
                        "pair manifest order/indices do not match the BCF-ordered sample list"
                    )
            if len(fields) < 4:
                raise ValueError("pair manifest is missing haplotype labels")
            expected_i = (
                f"{expected_samples[haplotype_i // 2]}.{haplotype_i % 2}"
            )
            expected_j = (
                f"{expected_samples[haplotype_j // 2]}.{haplotype_j % 2}"
            )
            if fields[2] != expected_i or fields[3] != expected_j:
                raise ValueError(
                    "pair manifest haplotype labels do not match the sample list"
                )
        pairs.add(pair)
        pair_index += 1
    if not pairs:
        raise ValueError("pair manifest contains no pairs")
    if expected_count is not None and len(pairs) != expected_count:
        raise ValueError(
            f"pair manifest has {len(pairs)} pairs; expected {expected_count}"
        )
    return {"n_pairs": len(pairs), "pair_mode": pair_mode}


def read_sample_list(path: str | Path) -> list[str]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"sample list is absent or empty: {path}")
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    if any(not line.strip() for line in raw_lines):
        raise ValueError("sample list contains an empty line")
    samples = [line.strip() for line in raw_lines]
    if len(set(samples)) != len(samples):
        raise ValueError("sample list contains duplicate IDs")
    return samples


def validate_sample_list(path: str | Path) -> dict:
    samples = read_sample_list(path)
    return {"n_samples": len(samples), "sha256": sha256_file(path)}


def validate_sample_selection(contract: dict, *, verify_sources: bool) -> dict:
    selection = contract["sample_selection"]
    sample_path = Path(selection["sample_list"])
    audit_path = Path(selection["audit"])
    sample_metrics = validate_sample_list(sample_path)
    if not audit_path.is_file() or audit_path.stat().st_size == 0:
        raise ValueError(f"sample-selection audit is absent or empty: {audit_path}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("population") != contract["population"]:
        raise ValueError("sample-selection audit population does not match the contract")
    if audit.get("ancestry_label_column") != ANCESTRY_LABEL_COLUMN:
        raise ValueError("sample-selection audit used the wrong ancestry label column")
    audit_sample = audit.get("sample_list", {})
    if str(Path(audit_sample.get("path", "")).resolve()) != str(sample_path.resolve()):
        raise ValueError("sample-selection audit points to a different sample list")
    if audit_sample.get("sha256") != sample_metrics["sha256"]:
        raise ValueError("sample-selection audit sample-list hash does not match")
    if int(audit_sample.get("n_samples", -1)) != sample_metrics["n_samples"]:
        raise ValueError("sample-selection audit sample count does not match")
    if verify_sources:
        for name, source in selection["sources"].items():
            local_path = Path(source["local_path"])
            if not local_path.is_file():
                raise ValueError(f"sample-selection source is absent: {local_path}")
            audit_source = audit.get("sources", {}).get(name, {})
            if audit_source.get("sha256") != sha256_file(local_path):
                raise ValueError(f"sample-selection source hash changed: {name}")
    return {
        **sample_metrics,
        "audit_sha256": sha256_file(audit_path),
        "audit_size_bytes": audit_path.stat().st_size,
    }


def validate_callable_mask(contract: dict, *, required: bool) -> dict | None:
    mask = contract["mask"]
    if mask is None:
        return None
    source_path = Path(mask["source_local_path"])
    callable_path = Path(mask["local_path"])
    audit_path = Path(mask["audit"])
    if not callable_path.exists() and not audit_path.exists() and not required:
        return None
    for label, path in (
        ("hard-mask source", source_path),
        ("callable mask", callable_path),
        ("callable-mask audit", audit_path),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{label} is absent or empty: {path}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("source_semantics") != "excluded_intervals":
        raise ValueError("mask audit did not treat the hard mask as exclusions")
    if audit.get("decoder_mask_semantics") != "included_intervals":
        raise ValueError("mask audit did not create a positive decoder mask")
    audit_contig = str(audit.get("contig", ""))
    normalized_contig = (
        audit_contig[3:] if audit_contig.lower().startswith("chr") else audit_contig
    )
    if normalized_contig != str(contract["chromosome"]):
        raise ValueError("callable-mask audit contig does not match the contract")
    if audit.get("source", {}).get("sha256") != sha256_file(source_path):
        raise ValueError("hard-mask source hash does not match its audit")
    callable_record = audit.get("callable_mask", {})
    if str(Path(callable_record.get("path", "")).resolve()) != str(
        callable_path.resolve()
    ):
        raise ValueError("mask audit points to a different callable mask")
    callable_sha256 = sha256_file(callable_path)
    if callable_record.get("sha256") != callable_sha256:
        raise ValueError("callable-mask hash does not match its audit")
    counts = audit.get("counts", {})
    if int(counts.get("callable_bases", 0)) < 1:
        raise ValueError("callable-mask audit reports no callable bases")
    return {
        "callable_sha256": callable_sha256,
        "callable_size_bytes": callable_path.stat().st_size,
        "audit_sha256": sha256_file(audit_path),
        "audit_size_bytes": audit_path.stat().st_size,
        "contig": audit.get("contig"),
        "sequence_length": int(audit.get("sequence_length", -1)),
        "excluded_bases": int(counts.get("excluded_bases", -1)),
        "callable_bases": int(counts.get("callable_bases", -1)),
        "callable_fraction": float(counts.get("callable_fraction", -1)),
    }


def _manifest_pairs(path: str | Path) -> np.ndarray:
    rows = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        rows.append((int(fields[0]), int(fields[1])))
    return np.asarray(rows, dtype=np.int64)


def validate_workbench_bitmatrix(
    contract: dict,
    *,
    summary_frame: pd.DataFrame,
    required: bool,
) -> dict | None:
    outputs = contract.get("analysis_outputs")
    if outputs is None:
        return None
    path = Path(outputs["bitmatrix"])
    metadata_path = bitmatrix_reader.meta_path(path)
    if not path.exists() and not required:
        return None
    for label, candidate in (("bit matrix", path), ("bit-matrix metadata", metadata_path)):
        if not candidate.is_file() or candidate.stat().st_size == 0:
            raise ValueError(f"{label} is absent or empty: {candidate}")
    metadata = bitmatrix_reader.read_meta(path)
    decoder = contract["decoder"]
    pairs = _manifest_pairs(contract["pairs_manifest"])
    metadata_pairs = np.asarray(metadata.get("pairs", []), dtype=np.int64)
    if metadata_pairs.shape != pairs.shape or not np.array_equal(metadata_pairs, pairs):
        raise ValueError("bit-matrix pair order does not match the pair manifest")
    if int(metadata.get("n_pairs", -1)) != len(pairs):
        raise ValueError("bit-matrix pair count does not match the pair manifest")
    if metadata.get("recent_call") != decoder["recent_call"]:
        raise ValueError("bit-matrix call rule does not match the contract")
    thresholds = metadata.get("thresholds_years", [])
    if len(thresholds) != 1 or not np.isclose(
        float(thresholds[0]), decoder["threshold_years"], rtol=0, atol=1e-12
    ):
        raise ValueError("bit-matrix threshold does not match the contract")
    expected_two_ne = decoder["theta"] / (2.0 * decoder["mutation_rate"])
    if not np.isclose(
        float(metadata.get("two_ne_generations", -1)),
        expected_two_ne,
        rtol=1e-9,
        atol=1e-9,
    ):
        raise ValueError("bit-matrix time scale does not match the contract")
    observed_positions = bitmatrix_reader.output_positions(metadata)
    expected_positions = summary_frame["position_0based"].to_numpy(dtype=np.int64)
    if not np.array_equal(observed_positions, expected_positions):
        raise ValueError("bit-matrix positions do not match the summary")
    sample_names = read_sample_list(contract["sample_selection"]["sample_list"])
    labels = metadata.get("sample_names", {})
    expected_labels = {
        str(2 * index): f"{sample}.0"
        for index, sample in enumerate(sample_names)
    }
    expected_labels.update({
        str(2 * index + 1): f"{sample}.1"
        for index, sample in enumerate(sample_names)
    })
    if labels != expected_labels:
        raise ValueError("bit-matrix sample order does not match the sample list")
    return {
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "metadata_sha256": sha256_file(metadata_path),
        "metadata_size_bytes": metadata_path.stat().st_size,
        "n_pairs": int(len(pairs)),
        "n_positions": int(len(observed_positions)),
    }


def validate_candidate_analysis(contract: dict, *, required: bool) -> dict | None:
    outputs = contract.get("analysis_outputs")
    if outputs is None:
        return None
    regions_path = Path(outputs["candidate_regions"])
    manifest_path = Path(outputs["candidate_manifest"])
    if not regions_path.exists() and not manifest_path.exists() and not required:
        return None
    for label, path in (("candidate regions", regions_path), ("candidate manifest", manifest_path)):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{label} is absent or empty: {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("population") != contract["population"]:
        raise ValueError("candidate manifest population does not match the contract")
    if int(manifest.get("chromosome", -1)) != contract["chromosome"]:
        raise ValueError("candidate manifest chromosome does not match the contract")
    if manifest.get("regions_sha256") != sha256_file(regions_path):
        raise ValueError("candidate-region table does not match its analysis manifest")
    pair_digest = sha256_file(contract["pairs_manifest"])
    if manifest.get("pair_contract", {}).get("pairs_manifest_sha256") != pair_digest:
        raise ValueError("candidate analysis did not use the decoded pair manifest")
    artifact_records = manifest.get("artifacts")
    if not isinstance(artifact_records, list):
        raise ValueError("candidate manifest artifact list is malformed")
    root = manifest_path.parent.resolve()
    present = 0
    for record in artifact_records:
        relative = Path(str(record.get("path", "")))
        artifact = (root / relative).resolve()
        try:
            artifact.relative_to(root)
        except ValueError as error:
            raise ValueError("candidate artifact escapes its output directory") from error
        if not artifact.exists() and not required:
            continue
        if not artifact.is_file():
            raise ValueError(f"candidate artifact is absent: {artifact}")
        present += 1
        if artifact.stat().st_size != int(record.get("size_bytes", -1)):
            raise ValueError(f"candidate artifact size changed: {artifact}")
        if sha256_file(artifact) != record.get("sha256"):
            raise ValueError(f"candidate artifact hash changed: {artifact}")
    if required and present != len(artifact_records):
        raise ValueError("candidate analysis is missing one or more artifacts")
    return {
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_size_bytes": manifest_path.stat().st_size,
        "regions_sha256": sha256_file(regions_path),
        "regions_size_bytes": regions_path.stat().st_size,
        "n_regions": int(manifest.get("n_regions", -1)),
        "n_regions_with_representative_variant": int(
            manifest.get("n_regions_with_representative_variant", -1)
        ),
        "n_artifacts": int(len(artifact_records)),
    }


def write_workbench_completion(
    *,
    summary_path: str | Path,
    run_json_path: str | Path,
    completion_path: str | Path,
    contract: dict,
) -> dict:
    summary_path = Path(summary_path)
    run_json_path = Path(run_json_path)
    pairs_manifest_path = Path(contract["pairs_manifest"])
    sample_selection = validate_sample_selection(contract, verify_sources=True)
    mask_metrics = validate_callable_mask(contract, required=True)
    summary_frame, summary_metrics = validate_summary(
        summary_path,
        threshold_years=contract["decoder"]["threshold_years"],
    )
    run = validate_run_json(run_json_path, contract)
    expected_samples = read_sample_list(contract["sample_selection"]["sample_list"])
    decoder = contract["decoder"]
    pair_metrics = validate_pairs_manifest(
        pairs_manifest_path,
        expected_samples=expected_samples,
        pair_mode=decoder["pair_selection"],
        expected_count=(
            decoder["n_random_pairs"]
            if decoder["pair_selection"] == "random_haplotype_pairs"
            else len(expected_samples)
        ),
    )
    if int(run.get("n_output_positions", -1)) != summary_metrics["n_output_positions"]:
        raise ValueError("run metadata output-position count does not match the summary")
    if int(run.get("n_pairs_recorded", -1)) != pair_metrics["n_pairs"]:
        raise ValueError("run metadata pair count does not match the pair manifest")
    if (
        decoder["pair_selection"] == "within_individual"
        and pair_metrics["n_pairs"] != sample_selection["n_samples"]
    ):
        raise ValueError("within-individual pair count does not match the sample list")
    if summary_metrics["n_pairs_max"] > pair_metrics["n_pairs"]:
        raise ValueError("summary uses more pairs than the pair manifest contains")
    bitmatrix_metrics = validate_workbench_bitmatrix(
        contract, summary_frame=summary_frame, required=True
    )
    candidate_metrics = validate_candidate_analysis(contract, required=True)
    completion = {
        "schema_version": 1,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "contract_sha256": contract_sha256(contract),
        "contract": contract,
        "summary": {
            **summary_metrics,
            "sha256": sha256_file(summary_path),
            "size_bytes": summary_path.stat().st_size,
        },
        "run_json": {
            "sha256": sha256_file(run_json_path),
            "size_bytes": run_json_path.stat().st_size,
            "decode_seconds": float(run["decode_seconds"]),
        },
        "pairs_manifest": {
            **pair_metrics,
            "sha256": sha256_file(pairs_manifest_path),
            "size_bytes": pairs_manifest_path.stat().st_size,
        },
        "sample_selection": {
            **sample_selection,
            "sample_list_size_bytes": Path(
                contract["sample_selection"]["sample_list"]
            ).stat().st_size,
        },
        "mask": mask_metrics,
        "bitmatrix": bitmatrix_metrics,
        "candidate_analysis": candidate_metrics,
    }
    destination = Path(completion_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, destination)
    return completion


def validate_workbench_completion(
    *,
    summary_path: str | Path,
    run_json_path: str | Path,
    completion_path: str | Path,
    contract: dict,
) -> dict:
    destination = Path(completion_path)
    if not destination.is_file():
        raise ValueError(f"completion record is absent: {destination}")
    completion = json.loads(destination.read_text(encoding="utf-8"))
    expected_digest = contract_sha256(contract)
    if completion.get("contract_sha256") != expected_digest:
        raise ValueError("completion contract does not match requested inputs/settings")
    if completion.get("contract") != contract:
        raise ValueError("completion contract payload does not match requested contract")
    summary_frame, _ = validate_summary(
        summary_path,
        threshold_years=contract["decoder"]["threshold_years"],
    )
    validate_run_json(run_json_path, contract)
    pairs_manifest_path = Path(contract["pairs_manifest"])
    expected_samples = read_sample_list(contract["sample_selection"]["sample_list"])
    decoder = contract["decoder"]
    pair_metrics = validate_pairs_manifest(
        pairs_manifest_path,
        expected_samples=expected_samples,
        pair_mode=decoder["pair_selection"],
        expected_count=(
            decoder["n_random_pairs"]
            if decoder["pair_selection"] == "random_haplotype_pairs"
            else len(expected_samples)
        ),
    )
    sample_selection = validate_sample_selection(contract, verify_sources=False)
    mask_metrics = validate_callable_mask(contract, required=False)
    bitmatrix_metrics = validate_workbench_bitmatrix(
        contract, summary_frame=summary_frame, required=False
    )
    candidate_metrics = validate_candidate_analysis(contract, required=False)
    if pair_metrics["n_pairs"] != completion["pairs_manifest"]["n_pairs"]:
        raise ValueError("pair-manifest count no longer matches its completion record")
    if sample_selection["n_samples"] != completion["sample_selection"]["n_samples"]:
        raise ValueError("sample-list count no longer matches its completion record")
    if sample_selection["sha256"] != completion["sample_selection"]["sha256"]:
        raise ValueError("sample-list hash no longer matches its completion record")
    if (
        sample_selection["audit_sha256"]
        != completion["sample_selection"]["audit_sha256"]
    ):
        raise ValueError("sample-selection audit hash no longer matches its completion record")
    if mask_metrics is not None and mask_metrics != completion.get("mask"):
        raise ValueError("callable mask no longer matches its completion record")
    if bitmatrix_metrics is not None and bitmatrix_metrics != completion.get("bitmatrix"):
        raise ValueError("bit matrix no longer matches its completion record")
    if (
        candidate_metrics is not None
        and candidate_metrics != completion.get("candidate_analysis")
    ):
        raise ValueError("candidate analysis no longer matches its completion record")
    if sha256_file(summary_path) != completion["summary"]["sha256"]:
        raise ValueError("summary hash no longer matches its completion record")
    if sha256_file(run_json_path) != completion["run_json"]["sha256"]:
        raise ValueError("run-metadata hash no longer matches its completion record")
    if sha256_file(pairs_manifest_path) != completion["pairs_manifest"]["sha256"]:
        raise ValueError("pair-manifest hash no longer matches its completion record")
    return completion


def _plot_chromosome(
    frame: pd.DataFrame,
    *,
    population: str,
    chromosome: int,
    threshold_years: float,
    output_stem: Path,
) -> None:
    called_column = called_fraction_column(threshold_years)
    position_mb = frame["position_0based"].to_numpy(dtype=float) / 1e6
    figure, axes = plt.subplots(3, 1, figsize=(16, 11), sharex=True, constrained_layout=True)
    panels = (
        (called_column, f"Fraction: posterior mean TMRCA < {threshold_years:g} years", "#5b21b6"),
        (SOFT_COLUMN, f"Mean posterior P(TMRCA < {threshold_years:g} years)", "#0369a1"),
        (TMRCA_COLUMN, "Mean posterior TMRCA (generations)", "#374151"),
    )
    for axis, (column, ylabel, color) in zip(axes, panels):
        axis.plot(position_mb, frame[column], color=color, linewidth=0.8, rasterized=True)
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.16, linewidth=0.6)
    finite_tmrca = frame[TMRCA_COLUMN].to_numpy(dtype=float)
    if np.all(finite_tmrca[np.isfinite(finite_tmrca)] > 0):
        axes[2].set_yscale("log")
    axes[2].set_xlabel(f"Chromosome {chromosome} position (Mb)")
    figure.suptitle(
        f"{population} Gamma-SMC chromosome {chromosome} scan\n"
        "Descriptive posterior summaries; significance requires matched null calibration",
        fontsize=17,
    )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(Path(f"{output_stem}.png"), dpi=180, bbox_inches="tight")
    figure.savefig(Path(f"{output_stem}.pdf"), bbox_inches="tight")
    plt.close(figure)


def plot_workbench_population(
    summary_paths: dict[int, str | Path],
    *,
    population: str,
    output_dir: str | Path,
    threshold_years: float = 4500,
    whole_genome: bool = False,
    top_n: int = 100,
) -> dict:
    population = population.upper()
    if population not in POPULATIONS:
        raise ValueError(f"unsupported population: {population}")
    chromosomes = sorted(summary_paths)
    if not chromosomes:
        raise ValueError("no chromosome summaries were supplied")
    if any(chromosome not in AUTOSOMES for chromosome in chromosomes):
        raise ValueError("chromosome summaries must be autosomes 1-22")
    if whole_genome and chromosomes != list(AUTOSOMES):
        raise ValueError("whole-genome plotting requires all autosomes 1-22")

    output_dir = Path(output_dir)
    chromosome_plot_dir = output_dir / "chromosomes"
    chromosome_plot_dir.mkdir(parents=True, exist_ok=True)
    called_column = called_fraction_column(threshold_years)
    frames: dict[int, pd.DataFrame] = {}
    chromosome_rows = []
    input_records = []
    for chromosome in chromosomes:
        path = Path(summary_paths[chromosome]).resolve()
        frame, validation = validate_summary(path, threshold_years=threshold_years)
        frames[chromosome] = frame
        _plot_chromosome(
            frame,
            population=population,
            chromosome=chromosome,
            threshold_years=threshold_years,
            output_stem=chromosome_plot_dir / f"chr{chromosome}.gamma_smc",
        )
        hard_index = int(frame[called_column].astype(float).idxmax())
        soft_index = int(frame[SOFT_COLUMN].astype(float).idxmax())
        tmrca_index = int(frame[TMRCA_COLUMN].astype(float).idxmin())
        chromosome_rows.append({
            "population": population,
            "chromosome": chromosome,
            **validation,
            "maximum_called_fraction": float(frame.loc[hard_index, called_column]),
            "maximum_called_fraction_position_1based": int(frame.loc[hard_index, "position_1based"]),
            "maximum_soft_probability": float(frame.loc[soft_index, SOFT_COLUMN]),
            "maximum_soft_probability_position_1based": int(frame.loc[soft_index, "position_1based"]),
            "minimum_mean_tmrca_generations": float(frame.loc[tmrca_index, TMRCA_COLUMN]),
            "minimum_mean_tmrca_position_1based": int(frame.loc[tmrca_index, "position_1based"]),
        })
        input_records.append({
            "chromosome": chromosome,
            "path": str(path),
            "sha256": sha256_file(path),
            "rows": int(len(frame)),
        })

    chromosome_summary = pd.DataFrame(chromosome_rows)
    chromosome_summary.to_csv(
        output_dir / "chromosome_scan_summary.tsv", sep="\t", index=False
    )
    result = {
        "population": population,
        "chromosomes": chromosomes,
        "threshold_years": float(threshold_years),
        "called_fraction_column": called_column,
        "whole_genome_complete": bool(whole_genome),
        "inputs": input_records,
    }

    if whole_genome:
        genome_parts = []
        offset = 0.0
        ticks = []
        boundaries = []
        for chromosome in chromosomes:
            frame = frames[chromosome].copy()
            sequence_length = float(frame["position_0based"].max()) + 1.0
            frame.insert(0, "population", population)
            frame.insert(1, "chromosome", chromosome)
            frame["genome_position_0based"] = frame["position_0based"] + offset
            genome_parts.append(frame)
            ticks.append(offset + sequence_length / 2)
            offset += sequence_length
            boundaries.append(offset)
        genome = pd.concat(genome_parts, ignore_index=True)
        figure, axes = plt.subplots(
            3, 1, figsize=(22, 12), sharex=True, constrained_layout=True
        )
        panels = (
            (called_column, f"Fraction: posterior mean TMRCA < {threshold_years:g} years"),
            (SOFT_COLUMN, f"Mean posterior P(TMRCA < {threshold_years:g} years)"),
            (TMRCA_COLUMN, "Mean posterior TMRCA (generations)"),
        )
        colors = ("#2563eb", "#7c3aed")
        for axis, (column, ylabel) in zip(axes, panels):
            for index, chromosome in enumerate(chromosomes):
                subset = genome.loc[genome["chromosome"].eq(chromosome)]
                axis.scatter(
                    subset["genome_position_0based"] / 1e9,
                    subset[column],
                    s=1.2,
                    alpha=0.75,
                    color=colors[index % 2],
                    rasterized=True,
                )
            for boundary in boundaries[:-1]:
                axis.axvline(boundary / 1e9, color="0.8", linewidth=0.45)
            axis.set_ylabel(ylabel)
            axis.grid(axis="y", alpha=0.16, linewidth=0.6)
        finite_tmrca = genome[TMRCA_COLUMN].to_numpy(dtype=float)
        if np.all(finite_tmrca[np.isfinite(finite_tmrca)] > 0):
            axes[2].set_yscale("log")
        axes[2].set_xticks(np.asarray(ticks) / 1e9, [str(value) for value in chromosomes])
        axes[2].set_xlabel("Chromosome")
        figure.suptitle(
            f"{population} whole-genome Gamma-SMC scan\n"
            "Descriptive posterior summaries; significance requires matched null calibration",
            fontsize=18,
        )
        figure.savefig(
            output_dir / f"{population}.whole_genome.gamma_smc.png",
            dpi=180,
            bbox_inches="tight",
        )
        figure.savefig(
            output_dir / f"{population}.whole_genome.gamma_smc.pdf",
            bbox_inches="tight",
        )
        plt.close(figure)

        hard_top = genome.nlargest(top_n, called_column).copy()
        hard_top.insert(3, "ranking_statistic", called_column)
        hard_top.insert(4, "ranking_value", hard_top[called_column])
        soft_top = genome.nlargest(top_n, SOFT_COLUMN).copy()
        soft_top.insert(3, "ranking_statistic", SOFT_COLUMN)
        soft_top.insert(4, "ranking_value", soft_top[SOFT_COLUMN])
        top = pd.concat([hard_top, soft_top], ignore_index=True)
        top.to_csv(output_dir / "whole_genome_top_windows.tsv", sep="\t", index=False)
        result.update({
            "whole_genome_rows": int(len(genome)),
            "whole_genome_span_bp": int(offset),
            "top_windows_per_statistic": int(top_n),
        })

    (output_dir / "plot_manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result
