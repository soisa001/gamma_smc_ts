from __future__ import annotations

import gzip
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
from matplotlib.ticker import PercentFormatter

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
PLOT_MANIFEST_SCHEMA = "gamma_smc_aou.workbench-population-plots/v2"
REPORT_MANIFEST_SCHEMA = "gamma_smc_aou.workbench-run-report/v3"

CANDIDATE_REGION_COLUMNS = [
    "population",
    "chromosome",
    "region_id",
    "start_0based",
    "end_0based_exclusive",
    "n_signal_windows",
    "peak_position_0based",
    "peak_position_1based",
    "peak_fraction_recent",
    "peak_mean_tmrca_generations",
    "peak_mean_p_tmrca_lt_threshold",
]
RANKED_GENE_LIST_COLUMNS = [
    "population",
    "hit_id",
    "chromosome",
    "merged_start_0based",
    "merged_end_0based_exclusive",
    "merged_start_1based",
    "merged_end_1based_inclusive",
    "n_consecutive_1mb_bins",
    "n_top_windows",
    "ranking_statistic",
    "peak_ranking_value",
    "peak_position_0based",
    "peak_position_1based",
    "peak_genome_position_0based",
    "peak_mean_tmrca_generations",
    "peak_mean_p_tmrca_lt_threshold",
    "probable_gene",
    "probable_gene_id",
    "probable_gene_relation",
    "probable_gene_distance_bp",
    "probable_gene_start_0based",
    "probable_gene_end_0based_exclusive",
    "candidate_genes",
    "candidate_gene_coordinates_grch38",
    "candidate_gene_distances_to_peak_bp",
]


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


def _atomic_text(path: str | Path, text: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_frame(path: str | Path, frame: pd.DataFrame) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_figure(figure, path: str | Path, **savefig_kwargs) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.stem}.tmp.{os.getpid()}{destination.suffix}"
    )
    try:
        figure.savefig(temporary, **savefig_kwargs)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _artifact_record(path: Path, root: Path) -> dict:
    return {
        "path": path.resolve().relative_to(root.resolve()).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _reuse_artifact_manifest(
    manifest_path: Path,
    *,
    schema: str,
    contract: dict,
) -> dict | None:
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != schema or manifest.get("contract") != contract:
            return None
        root = manifest_path.parent.resolve()
        artifacts = manifest.get("artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            return None
        for record in artifacts:
            relative = Path(str(record["path"]))
            artifact = (root / relative).resolve()
            artifact.relative_to(root)
            if not artifact.is_file():
                return None
            if artifact.stat().st_size != int(record["size_bytes"]):
                return None
            if sha256_file(artifact) != record["sha256"]:
                return None
        result = dict(manifest["result"])
        result["reused"] = True
        return result
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _write_artifact_manifest(
    manifest_path: Path,
    *,
    schema: str,
    contract: dict,
    result: dict,
    artifacts: list[Path],
) -> None:
    manifest = {
        "schema": schema,
        "contract": contract,
        "result": result,
        "artifacts": [
            _artifact_record(path, manifest_path.parent) for path in artifacts
        ],
    }
    _atomic_text(manifest_path, json.dumps(manifest, indent=2) + "\n")


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
        if ancestry_by_id.get(sample_id) == population and sample_id not in excluded_ids
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
                raise ValueError(
                    f"hard-mask BED line {line_number} has invalid interval"
                )
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
    signal_fraction: float = 0.02,
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
        raise ValueError(
            "mask URI, fingerprint, and local path must be supplied together"
        )
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
    if not np.isclose(
        observed, expected, rtol=0, atol=max(1e-15, abs(expected) * 1e-12)
    ):
        raise ValueError(f"run command {flag}={observed}; expected {expected}")


def validate_run_json(run_json_path: str | Path, contract: dict) -> dict:
    path = Path(run_json_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"run metadata is absent or empty: {path}")
    run = json.loads(path.read_text(encoding="utf-8"))
    command = run.get("command")
    if not isinstance(command, list) or not all(
        isinstance(item, str) for item in command
    ):
        raise ValueError("run metadata command must be a string list")
    decoder = contract["decoder"]
    expected_input = contract["input"]["local_path"]
    if str(Path(_command_value(command, "--input")).resolve()) != expected_input:
        raise ValueError("run command input does not match the contract")
    expected_output = contract["output_summary"]
    if (
        str(Path(_command_value(command, "--recent_summary")).resolve())
        != expected_output
    ):
        raise ValueError("run command output does not match the contract")
    expected_manifest = contract["pairs_manifest"]
    if (
        str(Path(_command_value(command, "--pairs_manifest")).resolve())
        != expected_manifest
    ):
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
            raise ValueError(
                "run command random-pair count does not match the contract"
            )
        if _command_value(command, "--pairs_seed") != str(decoder["pairs_seed"]):
            raise ValueError("run command pair seed does not match the contract")
        if decoder["exclude_within"] != ("--exclude_within" in command):
            raise ValueError(
                "run command within-pair exclusion does not match the contract"
            )
    outputs = contract.get("analysis_outputs")
    if outputs is not None:
        if (
            str(Path(_command_value(command, "--recent_bitmatrix")).resolve())
            != outputs["bitmatrix"]
        ):
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
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) < 2:
            raise ValueError(
                f"pair manifest line {line_number} has fewer than two fields"
            )
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
                    raise ValueError(
                        "pair manifest has more pairs than the sample list"
                    )
                expected_pair = (2 * pair_index, 2 * pair_index + 1)
                if pair != expected_pair:
                    raise ValueError(
                        "pair manifest order/indices do not match the BCF-ordered sample list"
                    )
            if len(fields) < 4:
                raise ValueError("pair manifest is missing haplotype labels")
            expected_i = f"{expected_samples[haplotype_i // 2]}.{haplotype_i % 2}"
            expected_j = f"{expected_samples[haplotype_j // 2]}.{haplotype_j % 2}"
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
        raise ValueError(
            "sample-selection audit population does not match the contract"
        )
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
    for label, candidate in (
        ("bit matrix", path),
        ("bit-matrix metadata", metadata_path),
    ):
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
    # The native decoder represents scaled_mutation_rate as a C++ float and
    # derives 2Ne from that actual value. Permit that single-precision input
    # rounding while continuing to reject a materially different time scale.
    if not np.isclose(
        float(metadata.get("two_ne_generations", -1)),
        expected_two_ne,
        rtol=np.finfo(np.float32).eps,
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
        str(2 * index): f"{sample}.0" for index, sample in enumerate(sample_names)
    }
    expected_labels.update(
        {str(2 * index + 1): f"{sample}.1" for index, sample in enumerate(sample_names)}
    )
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
    for label, path in (
        ("candidate regions", regions_path),
        ("candidate manifest", manifest_path),
    ):
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
            raise ValueError(
                "candidate artifact escapes its output directory"
            ) from error
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
        raise ValueError(
            "run metadata output-position count does not match the summary"
        )
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
            "sample_list_size_bytes": Path(contract["sample_selection"]["sample_list"])
            .stat()
            .st_size,
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
        raise ValueError(
            "completion contract payload does not match requested contract"
        )
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
        raise ValueError(
            "sample-selection audit hash no longer matches its completion record"
        )
    if mask_metrics is not None and mask_metrics != completion.get("mask"):
        raise ValueError("callable mask no longer matches its completion record")
    if bitmatrix_metrics is not None and bitmatrix_metrics != completion.get(
        "bitmatrix"
    ):
        raise ValueError("bit matrix no longer matches its completion record")
    if candidate_metrics is not None and candidate_metrics != completion.get(
        "candidate_analysis"
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
    figure, axes = plt.subplots(
        3, 1, figsize=(16, 11), sharex=True, constrained_layout=True
    )
    panels = (
        (
            called_column,
            f"Fraction: posterior mean TMRCA < {threshold_years:g} years",
            "#5b21b6",
        ),
        (
            SOFT_COLUMN,
            f"Mean posterior P(TMRCA < {threshold_years:g} years)",
            "#0369a1",
        ),
        (TMRCA_COLUMN, "Mean posterior TMRCA (generations)", "#374151"),
    )
    for axis, (column, ylabel, color) in zip(axes, panels):
        axis.plot(
            position_mb, frame[column], color=color, linewidth=0.8, rasterized=True
        )
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
    _atomic_figure(figure, Path(f"{output_stem}.png"), dpi=180, bbox_inches="tight")
    _atomic_figure(figure, Path(f"{output_stem}.pdf"), bbox_inches="tight")
    plt.close(figure)


def _summary_stride(frame: pd.DataFrame) -> int:
    positions = frame["position_0based"].to_numpy(dtype=float)
    differences = np.diff(positions)
    positive = differences[differences > 0]
    return max(1, int(round(float(np.median(positive))))) if len(positive) else 1


def _sequence_span(frame: pd.DataFrame) -> float:
    positions = frame["position_0based"].to_numpy(dtype=float)
    return float(positions.max()) + _summary_stride(frame)


def candidate_regions_from_summary(
    frame: pd.DataFrame,
    *,
    population: str,
    chromosome: int,
    sequence_length: int,
    threshold_years: float = 4500,
    minimum_fraction: float = 0.02,
    merge_gap: int = 20_000,
    stride: int | None = None,
) -> pd.DataFrame:
    """Derive merged strict-threshold regions from a validated scan summary."""
    if not 0 <= minimum_fraction <= 1:
        raise ValueError("minimum fraction must be in [0, 1]")
    if sequence_length < 1 or merge_gap < 0:
        raise ValueError("sequence length must be positive and merge gap nonnegative")
    stride = _summary_stride(frame) if stride is None else int(stride)
    if stride < 1:
        raise ValueError("stride must be positive")

    called_column = called_fraction_column(threshold_years)
    signal = frame.loc[frame[called_column].astype(float).gt(minimum_fraction)].copy()
    groups: list[list[int]] = []
    current: list[int] = []
    current_end = -1
    for index, row in signal.sort_values("position_0based").iterrows():
        start = int(row["position_0based"])
        end = min(sequence_length, start + stride)
        if current and start > current_end + merge_gap:
            groups.append(current)
            current = []
        current.append(int(index))
        current_end = max(current_end, end)
    if current:
        groups.append(current)

    rows = []
    for region_number, indices in enumerate(groups, 1):
        subset = frame.loc[indices].copy().sort_values(
            [called_column, "mean_tmrca_generations", "position_0based"],
            ascending=[False, True, True],
        )
        peak = subset.iloc[0]
        positions = frame.loc[indices, "position_0based"].astype(int)
        rows.append(
            {
                "population": population.upper(),
                "chromosome": int(chromosome),
                "region_id": f"chr{chromosome}_{region_number:04d}",
                "start_0based": int(positions.min()),
                "end_0based_exclusive": int(
                    min(sequence_length, positions.max() + stride)
                ),
                "n_signal_windows": int(len(indices)),
                "peak_position_0based": int(peak["position_0based"]),
                "peak_position_1based": int(peak["position_1based"]),
                "peak_fraction_recent": float(peak[called_column]),
                "peak_mean_tmrca_generations": float(
                    peak["mean_tmrca_generations"]
                ),
                "peak_mean_p_tmrca_lt_threshold": float(
                    peak["mean_p_tmrca_lt_threshold"]
                ),
            }
        )
    return pd.DataFrame(rows, columns=CANDIDATE_REGION_COLUMNS)


def _assemble_genome(
    frames: dict[int, pd.DataFrame],
    *,
    population: str,
    chromosome_spans: dict[int, float] | None = None,
) -> tuple[pd.DataFrame, list[float], list[float], float]:
    genome_parts = []
    offset = 0.0
    ticks = []
    boundaries = []
    for chromosome in sorted(frames):
        frame = frames[chromosome].copy()
        sequence_length = (
            float(chromosome_spans[chromosome])
            if chromosome_spans is not None
            else _sequence_span(frame)
        )
        frame.insert(0, "population", population)
        frame.insert(1, "chromosome", chromosome)
        frame["genome_position_0based"] = frame["position_0based"] + offset
        genome_parts.append(frame)
        ticks.append(offset + sequence_length / 2)
        offset += sequence_length
        boundaries.append(offset)
    return pd.concat(genome_parts, ignore_index=True), ticks, boundaries, offset


def _draw_recent_genome_axis(
    axis,
    genome: pd.DataFrame,
    *,
    chromosomes: list[int],
    called_column: str,
    ticks: list[float],
    boundaries: list[float],
    signal_fraction: float,
    ymax: float,
    show_chromosomes: bool,
) -> None:
    colors = ("#111827", "#2563eb")
    for index, chromosome in enumerate(chromosomes):
        subset = genome.loc[genome["chromosome"].eq(chromosome)]
        axis.scatter(
            subset["genome_position_0based"] / 1e9,
            subset[called_column],
            s=1.2,
            alpha=0.8,
            color=colors[index % 2],
            linewidths=0,
            rasterized=True,
        )
    for boundary in boundaries[:-1]:
        axis.axvline(boundary / 1e9, color="0.78", linewidth=0.45, linestyle=":")
    axis.axhline(signal_fraction, color="#059669", linewidth=0.8, linestyle="--")
    axis.set_ylim(0, ymax)
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    axis.grid(axis="y", alpha=0.18, linewidth=0.55)
    axis.margins(x=0)
    if show_chromosomes:
        axis.set_xticks(
            np.asarray(ticks) / 1e9,
            [str(chromosome) for chromosome in chromosomes],
        )
        axis.set_xlabel("Chromosome")
    else:
        axis.tick_params(axis="x", labelbottom=False)


def _plot_population_recent_genome(
    genome: pd.DataFrame,
    *,
    population: str,
    chromosomes: list[int],
    called_column: str,
    threshold_years: float,
    ticks: list[float],
    boundaries: list[float],
    signal_fraction: float,
    output_dir: Path,
) -> list[Path]:
    finite = genome[called_column].to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    ymax = min(1.0, max(signal_fraction, float(finite.max())) * 1.08)
    figure, axis = plt.subplots(figsize=(22, 5.2), constrained_layout=True)
    _draw_recent_genome_axis(
        axis,
        genome,
        chromosomes=chromosomes,
        called_column=called_column,
        ticks=ticks,
        boundaries=boundaries,
        signal_fraction=signal_fraction,
        ymax=ymax,
        show_chromosomes=True,
    )
    axis.set_ylabel(f"% pairs called coalesced < {threshold_years:g}y")
    axis.set_title(
        f"{population} genome-wide Gamma-SMC scan: posterior mean TMRCA < "
        f"{threshold_years:g} years"
    )
    outputs = [
        output_dir / f"{population}.whole_genome.gamma_smc.png",
        output_dir / f"{population}.whole_genome.gamma_smc.pdf",
    ]
    _atomic_figure(figure, outputs[0], dpi=180, bbox_inches="tight")
    _atomic_figure(figure, outputs[1], bbox_inches="tight")
    plt.close(figure)
    return outputs


def _plot_population_diagnostics(
    genome: pd.DataFrame,
    *,
    population: str,
    chromosomes: list[int],
    called_column: str,
    threshold_years: float,
    ticks: list[float],
    boundaries: list[float],
    output_dir: Path,
) -> list[Path]:
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
                linewidths=0,
                rasterized=True,
            )
        for boundary in boundaries[:-1]:
            axis.axvline(boundary / 1e9, color="0.8", linewidth=0.45)
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.16, linewidth=0.6)
    finite_tmrca = genome[TMRCA_COLUMN].to_numpy(dtype=float)
    finite_tmrca = finite_tmrca[np.isfinite(finite_tmrca)]
    if len(finite_tmrca) and np.all(finite_tmrca > 0):
        axes[2].set_yscale("log")
    axes[2].set_xticks(np.asarray(ticks) / 1e9, [str(value) for value in chromosomes])
    axes[2].set_xlabel("Chromosome")
    figure.suptitle(
        f"{population} whole-genome Gamma-SMC diagnostics\n"
        "Descriptive posterior summaries; significance requires matched null calibration",
        fontsize=18,
    )
    outputs = [
        output_dir / f"{population}.whole_genome.diagnostics.png",
        output_dir / f"{population}.whole_genome.diagnostics.pdf",
    ]
    _atomic_figure(figure, outputs[0], dpi=180, bbox_inches="tight")
    _atomic_figure(figure, outputs[1], bbox_inches="tight")
    plt.close(figure)
    return outputs


def plot_workbench_population(
    summary_paths: dict[int, str | Path],
    *,
    population: str,
    output_dir: str | Path,
    threshold_years: float = 4500,
    whole_genome: bool = False,
    top_n: int = 100,
    signal_fraction: float = 0.02,
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
    if top_n < 1:
        raise ValueError("top_n must be positive")
    if not 0 <= signal_fraction <= 1:
        raise ValueError("signal_fraction must be in [0, 1]")

    output_dir = Path(output_dir)
    chromosome_plot_dir = output_dir / "chromosomes"
    chromosome_plot_dir.mkdir(parents=True, exist_ok=True)
    called_column = called_fraction_column(threshold_years)
    input_records = []
    for chromosome in chromosomes:
        path = Path(summary_paths[chromosome]).resolve()
        if not path.is_file():
            raise ValueError(f"summary is absent: {path}")
        input_records.append(
            {
                "chromosome": chromosome,
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    contract = {
        "population": population,
        "chromosomes": chromosomes,
        "threshold_years": float(threshold_years),
        "called_fraction_column": called_column,
        "whole_genome": bool(whole_genome),
        "top_n": int(top_n),
        "signal_fraction": float(signal_fraction),
        "inputs": input_records,
    }
    manifest_path = output_dir / "plot_manifest.json"
    reused = _reuse_artifact_manifest(
        manifest_path, schema=PLOT_MANIFEST_SCHEMA, contract=contract
    )
    if reused is not None:
        return reused

    frames: dict[int, pd.DataFrame] = {}
    chromosome_rows = []
    artifacts: list[Path] = []
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
        artifacts.extend(
            [
                chromosome_plot_dir / f"chr{chromosome}.gamma_smc.png",
                chromosome_plot_dir / f"chr{chromosome}.gamma_smc.pdf",
            ]
        )
        hard_index = int(frame[called_column].astype(float).idxmax())
        soft_index = int(frame[SOFT_COLUMN].astype(float).idxmax())
        tmrca_index = int(frame[TMRCA_COLUMN].astype(float).idxmin())
        chromosome_rows.append(
            {
                "population": population,
                "chromosome": chromosome,
                **validation,
                "maximum_called_fraction": float(frame.loc[hard_index, called_column]),
                "maximum_called_fraction_position_1based": int(
                    frame.loc[hard_index, "position_1based"]
                ),
                "maximum_soft_probability": float(frame.loc[soft_index, SOFT_COLUMN]),
                "maximum_soft_probability_position_1based": int(
                    frame.loc[soft_index, "position_1based"]
                ),
                "minimum_mean_tmrca_generations": float(
                    frame.loc[tmrca_index, TMRCA_COLUMN]
                ),
                "minimum_mean_tmrca_position_1based": int(
                    frame.loc[tmrca_index, "position_1based"]
                ),
            }
        )

    chromosome_summary = pd.DataFrame(chromosome_rows)
    chromosome_summary_path = output_dir / "chromosome_scan_summary.tsv"
    _atomic_frame(chromosome_summary_path, chromosome_summary)
    artifacts.append(chromosome_summary_path)
    result = {
        "population": population,
        "chromosomes": chromosomes,
        "threshold_years": float(threshold_years),
        "called_fraction_column": called_column,
        "whole_genome_complete": bool(whole_genome),
        "inputs": input_records,
        "reused": False,
    }

    if whole_genome:
        genome, ticks, boundaries, offset = _assemble_genome(
            frames, population=population
        )
        artifacts.extend(
            _plot_population_recent_genome(
                genome,
                population=population,
                chromosomes=chromosomes,
                called_column=called_column,
                threshold_years=threshold_years,
                ticks=ticks,
                boundaries=boundaries,
                signal_fraction=signal_fraction,
                output_dir=output_dir,
            )
        )
        artifacts.extend(
            _plot_population_diagnostics(
                genome,
                population=population,
                chromosomes=chromosomes,
                called_column=called_column,
                threshold_years=threshold_years,
                ticks=ticks,
                boundaries=boundaries,
                output_dir=output_dir,
            )
        )

        hard_top = genome.nlargest(top_n, called_column).copy()
        hard_top.insert(3, "ranking_statistic", called_column)
        hard_top.insert(4, "ranking_value", hard_top[called_column])
        soft_top = genome.nlargest(top_n, SOFT_COLUMN).copy()
        soft_top.insert(3, "ranking_statistic", SOFT_COLUMN)
        soft_top.insert(4, "ranking_value", soft_top[SOFT_COLUMN])
        top = pd.concat([hard_top, soft_top], ignore_index=True)
        top_path = output_dir / "whole_genome_top_windows.tsv"
        _atomic_frame(top_path, top)
        artifacts.append(top_path)
        result.update(
            {
                "whole_genome_rows": int(len(genome)),
                "whole_genome_span_bp": int(offset),
                "top_windows_per_statistic": int(top_n),
            }
        )

    _write_artifact_manifest(
        manifest_path,
        schema=PLOT_MANIFEST_SCHEMA,
        contract=contract,
        result=result,
        artifacts=artifacts,
    )
    return result


def _plot_combined_recent_genome(
    genomes: dict[str, pd.DataFrame],
    *,
    chromosomes: list[int],
    called_column: str,
    threshold_years: float,
    ticks: list[float],
    boundaries: list[float],
    signal_fraction: float,
    output_stem: Path,
    fixed_ymax: float | None = None,
    hit_labels: pd.DataFrame | None = None,
    hit_label_min_fraction: float = 0.0,
) -> list[Path]:
    finite = np.concatenate(
        [genome[called_column].to_numpy(dtype=float) for genome in genomes.values()]
    )
    finite = finite[np.isfinite(finite)]
    ymax = (
        float(fixed_ymax)
        if fixed_ymax is not None
        else min(1.0, max(signal_fraction, float(finite.max())) * 1.08)
    )
    figure, raw_axes = plt.subplots(
        len(genomes),
        1,
        figsize=(22, max(4.5, 2.4 * len(genomes))),
        sharex=True,
        sharey=True,
        constrained_layout=True,
        squeeze=False,
    )
    axes = raw_axes[:, 0]
    for index, (population, genome) in enumerate(genomes.items()):
        axis = axes[index]
        _draw_recent_genome_axis(
            axis,
            genome,
            chromosomes=chromosomes,
            called_column=called_column,
            ticks=ticks,
            boundaries=boundaries,
            signal_fraction=signal_fraction,
            ymax=ymax,
            show_chromosomes=index == len(axes) - 1,
        )
        if fixed_ymax is not None:
            clipped = genome.loc[genome[called_column].astype(float).ge(ymax)]
            if not clipped.empty:
                axis.scatter(
                    clipped["genome_position_0based"] / 1e9,
                    np.full(len(clipped), ymax),
                    marker="v",
                    s=9,
                    color="#9f1239",
                    linewidths=0,
                    rasterized=True,
                    zorder=4,
                )
        if hit_labels is not None and not hit_labels.empty:
            population_hits = hit_labels.loc[
                hit_labels["population"].astype(str).eq(population)
                & hit_labels["peak_ranking_value"].astype(float).gt(
                    hit_label_min_fraction
                )
            ].sort_values(["chromosome", "peak_genome_position_0based"])
            for hit in population_hits.itertuples(index=False):
                label_y = min(float(hit.peak_ranking_value), ymax * 0.94)
                axis.annotate(
                    str(hit.probable_gene),
                    xy=(float(hit.peak_genome_position_0based) / 1e9, label_y),
                    xytext=(6, 0),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation=0,
                    fontsize=8.5,
                    fontweight="bold",
                    color="#1f2937",
                    annotation_clip=True,
                )
        axis.text(
            0.003,
            0.92,
            population,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
        )
    figure.supylabel(f"% pairs called coalesced < {threshold_years:g}y")
    subtitle = "Posterior-mean TMRCA call; dashed line is the candidate-region screen"
    if fixed_ymax is not None:
        subtitle = (
            f"0-{fixed_ymax:.0%} detail; triangles mark values at or above the "
            f"axis ceiling; horizontal labels mark peaks >{hit_label_min_fraction:.0%} "
            "and their nearest protein-coding genes"
        )
    figure.suptitle(
        "Gamma-SMC genome-wide recent-coalescence scan by population\n" + subtitle
    )
    outputs = [
        Path(f"{output_stem}.png"),
        Path(f"{output_stem}.pdf"),
    ]
    _atomic_figure(figure, outputs[0], dpi=180, bbox_inches="tight")
    _atomic_figure(figure, outputs[1], bbox_inches="tight")
    plt.close(figure)
    return outputs


def _gtf_attributes(value: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in value.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, raw_value = item.split(None, 1)
        attributes[key] = raw_value.strip().strip('"')
    return attributes


def _load_protein_coding_genes(path: str | Path) -> pd.DataFrame:
    """Read GRCh38 GENCODE GTF genes into 0-based half-open coordinates."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"gene annotation is absent or empty: {path}")
    rows: list[dict] = []
    handle = (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open("rt", encoding="utf-8")
    )
    with handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chromosome_label = fields[0].removeprefix("chr")
            if not chromosome_label.isdigit():
                continue
            chromosome = int(chromosome_label)
            if chromosome not in AUTOSOMES:
                continue
            attributes = _gtf_attributes(fields[8])
            gene_type = attributes.get("gene_type", attributes.get("gene_biotype", ""))
            if gene_type != "protein_coding":
                continue
            gene_name = attributes.get("gene_name")
            gene_id = attributes.get("gene_id")
            if not gene_name or not gene_id:
                continue
            rows.append(
                {
                    "chromosome": chromosome,
                    "gene_start_0based": int(fields[3]) - 1,
                    "gene_end_0based_exclusive": int(fields[4]),
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "gene_type": gene_type,
                }
            )
    genes = pd.DataFrame(rows)
    if genes.empty:
        raise ValueError("gene annotation contains no autosomal protein-coding genes")
    if genes.duplicated(["gene_id", "chromosome"]).any():
        raise ValueError("gene annotation contains duplicate autosomal gene records")
    return genes.sort_values(
        ["chromosome", "gene_start_0based", "gene_end_0based_exclusive", "gene_name"]
    ).reset_index(drop=True)


def _point_gene_distance(position: int, start: pd.Series, end: pd.Series) -> pd.Series:
    return pd.Series(
        np.where(
            position < start,
            start - position,
            np.where(position >= end, position - end + 1, 0),
        ),
        index=start.index,
        dtype=int,
    )


def _annotate_ranked_hit(
    *,
    chromosome: int,
    region_start: int,
    region_end: int,
    peak_position: int,
    genes: pd.DataFrame,
    context_flank: int,
) -> dict:
    chromosome_genes = genes.loc[genes["chromosome"].eq(chromosome)].copy()
    if chromosome_genes.empty:
        raise ValueError(f"gene annotation contains no protein-coding chr{chromosome} genes")
    chromosome_genes["distance_to_peak_bp"] = _point_gene_distance(
        peak_position,
        chromosome_genes["gene_start_0based"],
        chromosome_genes["gene_end_0based_exclusive"],
    )
    chromosome_genes["overlaps_hit"] = (
        chromosome_genes["gene_start_0based"].lt(region_end)
        & chromosome_genes["gene_end_0based_exclusive"].gt(region_start)
    )
    context_start = max(0, region_start - context_flank)
    context_end = region_end + context_flank
    candidates = chromosome_genes.loc[
        chromosome_genes["gene_start_0based"].lt(context_end)
        & chromosome_genes["gene_end_0based_exclusive"].gt(context_start)
    ].copy()
    if candidates.empty:
        candidates = chromosome_genes.nsmallest(5, "distance_to_peak_bp").copy()
    candidates.sort_values(
        ["distance_to_peak_bp", "gene_start_0based", "gene_name"], inplace=True
    )
    probable = candidates.iloc[0]
    if int(probable["distance_to_peak_bp"]) == 0:
        relation = "peak_overlap"
    elif bool(probable["overlaps_hit"]):
        relation = "merged_hit_overlap"
    elif int(probable["gene_start_0based"]) < context_end and int(
        probable["gene_end_0based_exclusive"]
    ) > context_start:
        relation = "within_context_flank"
    else:
        relation = "nearest_outside_context"
    return {
        "probable_gene": str(probable["gene_name"]),
        "probable_gene_id": str(probable["gene_id"]),
        "probable_gene_relation": relation,
        "probable_gene_distance_bp": int(probable["distance_to_peak_bp"]),
        "probable_gene_start_0based": int(probable["gene_start_0based"]),
        "probable_gene_end_0based_exclusive": int(
            probable["gene_end_0based_exclusive"]
        ),
        "candidate_genes": ";".join(candidates["gene_name"].astype(str)),
        "candidate_gene_coordinates_grch38": ";".join(
            f"{row.gene_name}=chr{chromosome}:{int(row.gene_start_0based) + 1}-"
            f"{int(row.gene_end_0based_exclusive)}"
            for row in candidates.itertuples(index=False)
        ),
        "candidate_gene_distances_to_peak_bp": ";".join(
            f"{row.gene_name}={int(row.distance_to_peak_bp)}"
            for row in candidates.itertuples(index=False)
        ),
    }


def _build_ranked_gene_list(
    genomes: dict[str, pd.DataFrame],
    *,
    called_column: str,
    genes: pd.DataFrame,
    top_n: int,
    hit_bin_size: int,
    context_flank: int,
) -> pd.DataFrame:
    """Group each population's top scan windows through consecutive 1-Mb bins."""
    rows: list[dict] = []
    for population, genome in genomes.items():
        ranked = genome.loc[genome[called_column].astype(float).gt(0)].sort_values(
            [called_column, TMRCA_COLUMN, "position_0based"],
            ascending=[False, True, True],
        ).head(top_n).copy()
        if ranked.empty:
            continue
        ranked["hit_bin_start_0based"] = (
            ranked["position_0based"].astype(int) // hit_bin_size
        ) * hit_bin_size
        hit_number = 0
        for chromosome, chromosome_ranked in ranked.groupby("chromosome", sort=True):
            chromosome = int(chromosome)
            scan_chromosome = genome.loc[genome["chromosome"].eq(chromosome)]
            chromosome_end = int(np.ceil(_sequence_span(scan_chromosome)))
            bins = sorted(chromosome_ranked["hit_bin_start_0based"].unique())
            groups: list[list[int]] = []
            current: list[int] = []
            for bin_start in bins:
                bin_start = int(bin_start)
                if current and bin_start > current[-1] + hit_bin_size:
                    groups.append(current)
                    current = []
                current.append(bin_start)
            if current:
                groups.append(current)
            for group in groups:
                hit_number += 1
                region_start = int(group[0])
                region_end = min(chromosome_end, int(group[-1]) + hit_bin_size)
                hit_windows = chromosome_ranked.loc[
                    chromosome_ranked["hit_bin_start_0based"].isin(group)
                ].sort_values(
                    [called_column, TMRCA_COLUMN, "position_0based"],
                    ascending=[False, True, True],
                )
                peak = hit_windows.iloc[0]
                peak_position = int(peak["position_0based"])
                rows.append(
                    {
                        "population": population,
                        "hit_id": f"{population}_hit_{hit_number:03d}",
                        "chromosome": chromosome,
                        "merged_start_0based": region_start,
                        "merged_end_0based_exclusive": region_end,
                        "merged_start_1based": region_start + 1,
                        "merged_end_1based_inclusive": region_end,
                        "n_consecutive_1mb_bins": len(group),
                        "n_top_windows": len(hit_windows),
                        "ranking_statistic": called_column,
                        "peak_ranking_value": float(peak[called_column]),
                        "peak_position_0based": peak_position,
                        "peak_position_1based": int(peak["position_1based"]),
                        "peak_genome_position_0based": float(
                            peak["genome_position_0based"]
                        ),
                        "peak_mean_tmrca_generations": float(peak[TMRCA_COLUMN]),
                        "peak_mean_p_tmrca_lt_threshold": float(peak[SOFT_COLUMN]),
                        **_annotate_ranked_hit(
                            chromosome=chromosome,
                            region_start=region_start,
                            region_end=region_end,
                            peak_position=peak_position,
                            genes=genes,
                            context_flank=context_flank,
                        ),
                    }
                )
    if not rows:
        return pd.DataFrame(columns=RANKED_GENE_LIST_COLUMNS)
    return (
        pd.DataFrame(rows, columns=RANKED_GENE_LIST_COLUMNS)
        .sort_values(["population", "chromosome", "merged_start_0based"])
        .reset_index(drop=True)
    )


def summarize_workbench_run(
    results_root: str | Path,
    *,
    populations: list[str] | tuple[str, ...],
    chromosomes: list[int] | tuple[int, ...],
    output_dir: str | Path,
    threshold_years: float = 4500,
    signal_fraction: float = 0.02,
    merge_gap: int = 20_000,
    whole_genome: bool = False,
    top_n: int = 100,
    gene_annotation: str | Path | None = None,
    hit_bin_size: int = 1_000_000,
    gene_context_flank: int = 500_000,
    zoom_ymax: float = 0.04,
    hit_label_min_fraction: float = 0.02,
) -> dict:
    results_root = Path(results_root).resolve()
    output_dir = Path(output_dir).resolve()
    populations = [population.upper() for population in populations]
    chromosomes = sorted(int(chromosome) for chromosome in chromosomes)
    if not populations or len(populations) != len(set(populations)):
        raise ValueError("populations must be nonempty and unique")
    if any(population not in POPULATIONS for population in populations):
        raise ValueError("unsupported population in run report")
    if not chromosomes or len(chromosomes) != len(set(chromosomes)):
        raise ValueError("chromosomes must be nonempty and unique")
    if any(chromosome not in AUTOSOMES for chromosome in chromosomes):
        raise ValueError("run report chromosomes must be autosomes 1-22")
    if whole_genome and chromosomes != list(AUTOSOMES):
        raise ValueError("whole-genome report requires all autosomes 1-22")
    if not 0 <= signal_fraction <= 1:
        raise ValueError("signal_fraction must be in [0, 1]")
    if top_n < 1 or hit_bin_size < 1 or gene_context_flank < 0 or merge_gap < 0:
        raise ValueError("top_n/bin size must be positive and gene flank nonnegative")
    if not 0 < zoom_ymax <= 1:
        raise ValueError("zoom_ymax must be in (0, 1]")
    if not 0 <= hit_label_min_fraction <= 1:
        raise ValueError("hit_label_min_fraction must be in [0, 1]")

    gene_annotation_record = None
    if gene_annotation is not None:
        gene_annotation = Path(gene_annotation).resolve()
        if not gene_annotation.is_file() or gene_annotation.stat().st_size == 0:
            raise ValueError(f"gene annotation is absent or empty: {gene_annotation}")
        gene_annotation_record = {
            "path": str(gene_annotation),
            "size_bytes": gene_annotation.stat().st_size,
            "sha256": sha256_file(gene_annotation),
        }

    input_records = []
    summary_paths: dict[str, dict[int, Path]] = {}
    for population in populations:
        summary_paths[population] = {}
        for chromosome in chromosomes:
            root = results_root / population / "chromosomes"
            summary = root / f"chr{chromosome}.gamma_smc.tsv"
            if not summary.is_file():
                raise ValueError(f"summary input is absent: {summary}")
            input_records.append(
                {
                    "population": population,
                    "chromosome": chromosome,
                    "kind": "summary",
                    "path": str(summary.resolve()),
                    "size_bytes": summary.stat().st_size,
                    "sha256": sha256_file(summary),
                }
            )
            summary_paths[population][chromosome] = summary

    called_column = called_fraction_column(threshold_years)
    contract = {
        "populations": populations,
        "chromosomes": chromosomes,
        "threshold_years": float(threshold_years),
        "called_fraction_column": called_column,
        "signal_fraction": float(signal_fraction),
        "merge_gap": int(merge_gap),
        "whole_genome": bool(whole_genome),
        "top_n": int(top_n),
        "gene_annotation": gene_annotation_record,
        "hit_bin_size": int(hit_bin_size),
        "gene_context_flank": int(gene_context_flank),
        "zoom_ymax": float(zoom_ymax),
        "hit_label_min_fraction": float(hit_label_min_fraction),
        "inputs": input_records,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "run_report_manifest.json"
    reused = _reuse_artifact_manifest(
        manifest_path, schema=REPORT_MANIFEST_SCHEMA, contract=contract
    )
    if reused is not None:
        return reused

    frames_by_population: dict[str, dict[int, pd.DataFrame]] = {}
    region_frames = []
    population_rows = []
    for population in populations:
        frames_by_population[population] = {}
        population_regions = []
        for chromosome in chromosomes:
            summary, _ = validate_summary(
                summary_paths[population][chromosome],
                threshold_years=threshold_years,
            )
            frames_by_population[population][chromosome] = summary
            regions = candidate_regions_from_summary(
                summary,
                population=population,
                chromosome=chromosome,
                sequence_length=int(_sequence_span(summary)),
                threshold_years=threshold_years,
                minimum_fraction=signal_fraction,
                merge_gap=merge_gap,
            )
            population_regions.append(regions)
            region_frames.append(regions)
        combined_regions = pd.concat(population_regions, ignore_index=True)
        if combined_regions.empty:
            peak_chromosome = peak_position = np.nan
            peak_fraction = np.nan
        else:
            peak = combined_regions.loc[
                combined_regions["peak_fraction_recent"].astype(float).idxmax()
            ]
            peak_chromosome = int(peak["chromosome"])
            peak_position = int(peak["peak_position_1based"])
            peak_fraction = float(peak["peak_fraction_recent"])
        population_rows.append(
            {
                "population": population,
                "chromosomes_complete": len(chromosomes),
                "regions_found": int(len(combined_regions)),
                "chromosomes_with_regions": int(
                    combined_regions["chromosome"].nunique()
                ),
                "signal_windows": int(combined_regions["n_signal_windows"].sum()),
                "region_span_bp": int(
                    (
                        combined_regions["end_0based_exclusive"]
                        - combined_regions["start_0based"]
                    ).sum()
                ),
                "maximum_peak_fraction_recent": peak_fraction,
                "maximum_peak_chromosome": peak_chromosome,
                "maximum_peak_position_1based": peak_position,
            }
        )

    population_summary = pd.DataFrame(population_rows)
    all_regions = pd.concat(region_frames, ignore_index=True)
    population_summary_path = output_dir / "regions_by_population.tsv"
    all_regions_path = output_dir / "all_candidate_regions.tsv"
    _atomic_frame(population_summary_path, population_summary)
    _atomic_frame(all_regions_path, all_regions)

    chromosome_spans = {
        chromosome: max(
            _sequence_span(frames_by_population[population][chromosome])
            for population in populations
        )
        for chromosome in chromosomes
    }
    genomes: dict[str, pd.DataFrame] = {}
    ticks = boundaries = None
    for population in populations:
        genome, population_ticks, population_boundaries, _ = _assemble_genome(
            frames_by_population[population],
            population=population,
            chromosome_spans=chromosome_spans,
        )
        genomes[population] = genome
        ticks = population_ticks
        boundaries = population_boundaries
    scope = "whole_genome" if whole_genome else "requested_chromosomes"
    figure_paths = _plot_combined_recent_genome(
        genomes,
        chromosomes=chromosomes,
        called_column=called_column,
        threshold_years=threshold_years,
        ticks=ticks or [],
        boundaries=boundaries or [],
        signal_fraction=signal_fraction,
        output_stem=output_dir / f"combined.{scope}.gamma_smc",
    )
    gene_list_path = None
    gene_list = pd.DataFrame()
    if gene_annotation is not None:
        genes = _load_protein_coding_genes(gene_annotation)
        gene_list = _build_ranked_gene_list(
            genomes,
            called_column=called_column,
            genes=genes,
            top_n=top_n,
            hit_bin_size=hit_bin_size,
            context_flank=gene_context_flank,
        )
        gene_list_path = output_dir / "gene_list.tsv"
        _atomic_frame(gene_list_path, gene_list)
    zoom_percent = f"{zoom_ymax * 100:g}".replace(".", "p")
    zoom_stem = output_dir / f"combined.{scope}.gamma_smc.zoom{zoom_percent}pct"
    zoom_paths = _plot_combined_recent_genome(
        genomes,
        chromosomes=chromosomes,
        called_column=called_column,
        threshold_years=threshold_years,
        ticks=ticks or [],
        boundaries=boundaries or [],
        signal_fraction=signal_fraction,
        output_stem=zoom_stem,
        fixed_ymax=zoom_ymax,
        hit_labels=gene_list,
        hit_label_min_fraction=hit_label_min_fraction,
    )
    counts = {row["population"]: int(row["regions_found"]) for row in population_rows}
    result = {
        "populations": populations,
        "chromosomes": chromosomes,
        "whole_genome_complete": bool(whole_genome),
        "threshold_years": float(threshold_years),
        "signal_fraction": float(signal_fraction),
        "merge_gap": int(merge_gap),
        "population_region_counts": counts,
        "total_regions": int(len(all_regions)),
        "ranked_gene_hits": int(len(gene_list)),
        "plotted_gene_hits": int(
            gene_list["peak_ranking_value"].astype(float).gt(
                hit_label_min_fraction
            ).sum()
        ) if not gene_list.empty else 0,
        "hit_label_min_fraction": float(hit_label_min_fraction),
        "reused": False,
    }
    report_artifacts = [
        population_summary_path,
        all_regions_path,
        *figure_paths,
        *zoom_paths,
    ]
    if gene_list_path is not None:
        report_artifacts.append(gene_list_path)
    _write_artifact_manifest(
        manifest_path,
        schema=REPORT_MANIFEST_SCHEMA,
        contract=contract,
        result=result,
        artifacts=report_artifacts,
    )
    return result
