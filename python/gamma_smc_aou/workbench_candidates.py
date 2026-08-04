from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zstandard
from scipy.stats import t

from .workbench import (
    CANDIDATE_REGION_COLUMNS,
    POPULATIONS,
    candidate_regions_from_summary,
    called_fraction_column,
    sha256_file,
    validate_summary,
)


REGION_COLUMNS = CANDIDATE_REGION_COLUMNS

VARIANT_COLUMNS = [
    "population",
    "chromosome",
    "region_id",
    "chrom",
    "position_0based",
    "position_1based",
    "id",
    "ref",
    "alt",
    "alt_index",
    "variant_type",
    "qual",
    "filter",
    "alternate_allele_frequency",
    "tmrca_position_0based",
    "n_hom_ref_pairs",
    "n_hom_alt_pairs",
    "mean_hom_ref_tmrca_generations",
    "mean_hom_alt_tmrca_generations",
    "hom_alt_minus_hom_ref_generations",
    "r2_tmrca",
    "r2_log1p_tmrca",
]

PROFILE_COLUMNS = [
    "population",
    "chromosome",
    "region_id",
    "representative_position_1based",
    "representative_ref",
    "representative_alt",
    "position_0based",
    "position_1based",
    "pair_class",
    "n_pairs",
    "mean_tmrca_generations",
    "median_tmrca_generations",
    "standard_error_generations",
    "ci95_lower_generations",
    "ci95_upper_generations",
]

DISTRIBUTION_COLUMNS = [
    "population",
    "chromosome",
    "region_id",
    "position_0based",
    "position_1based",
    "n_pairs",
    "tmrca_p05_generations",
    "tmrca_p25_generations",
    "tmrca_median_generations",
    "tmrca_p75_generations",
    "tmrca_p95_generations",
]


def _atomic_text(path: str | Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _atomic_frame(path: str | Path, frame: pd.DataFrame, *, gzip: bool = False) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    compression = {"method": "gzip", "compresslevel": 6, "mtime": 0} if gzip else None
    frame.to_csv(temporary, sep="\t", index=False, compression=compression)
    os.replace(temporary, path)


def build_candidate_regions(
    summary_path: str | Path,
    *,
    population: str,
    chromosome: int,
    sequence_length: int,
    output_path: str | Path,
    positions_path: str | Path,
    threshold_years: float = 4500,
    minimum_fraction: float = 0.02,
    merge_gap: int = 20_000,
    stride: int = 10_000,
    profile_half_width: int = 500_000,
) -> pd.DataFrame:
    """Find >threshold scan windows and merge gaps up to ``merge_gap`` bases."""
    population = population.upper()
    if population not in POPULATIONS:
        raise ValueError(f"unsupported population: {population}")
    if chromosome < 1 or chromosome > 22 or sequence_length < 1:
        raise ValueError("an autosome and positive sequence length are required")
    if not 0 <= minimum_fraction <= 1:
        raise ValueError("minimum fraction must be in [0, 1]")
    if min(merge_gap, stride, profile_half_width) < 0 or stride < 1:
        raise ValueError("gap/profile widths must be nonnegative and stride positive")

    frame, _ = validate_summary(summary_path, threshold_years=threshold_years)
    regions = candidate_regions_from_summary(
        frame,
        population=population,
        chromosome=chromosome,
        sequence_length=sequence_length,
        threshold_years=threshold_years,
        minimum_fraction=minimum_fraction,
        merge_gap=merge_gap,
        stride=stride,
    )
    _atomic_frame(output_path, regions)
    requested: set[int] = set()
    for row in regions.itertuples(index=False):
        start = max(0, int(row.peak_position_0based) - profile_half_width)
        end = min(sequence_length - 1, int(row.peak_position_0based) + profile_half_width)
        requested.update(
            frame.loc[
                frame["position_0based"].between(start, end, inclusive="both"),
                "position_0based",
            ].astype(int)
        )
    _atomic_text(positions_path, "".join(f"{value}\n" for value in sorted(requested)))
    return regions


def read_pair_manifest(path: str | Path) -> np.ndarray:
    pairs = []
    for line_number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) < 2:
            raise ValueError(f"pair manifest line {line_number} is malformed")
        pairs.append((int(fields[0]), int(fields[1])))
    if not pairs:
        raise ValueError("pair manifest contains no pairs")
    result = np.asarray(pairs, dtype=np.int64)
    if np.any(result[:, 0] < 0) or np.any(result[:, 1] <= result[:, 0]):
        raise ValueError("pair manifest must contain ordered distinct haplotype pairs")
    if len(np.unique(result, axis=0)) != len(result):
        raise ValueError("pair manifest contains duplicate pairs")
    return result


def _read_exact(reader, size: int) -> bytes:
    chunks = []
    remaining = size
    while remaining:
        chunk = reader.read(remaining)
        if not chunk:
            raise EOFError(f"posterior stream ended {remaining} bytes early")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def convert_raw_to_ordered_tmrca(
    raw_path: str | Path,
    *,
    pairs_manifest_path: str | Path,
    mutation_rate: float,
    output_path: str | Path,
    metadata_path: str | Path,
    working_path: str | Path,
) -> tuple[np.memmap, dict]:
    """Convert chunk-major alpha/beta to position-major TMRCA in manifest order."""
    raw_path = Path(raw_path)
    raw_meta_path = raw_path.with_name(raw_path.name + ".meta")
    if not raw_path.is_file() or not raw_meta_path.is_file():
        raise ValueError("raw posterior and metadata are required")
    raw_meta = json.loads(raw_meta_path.read_text(encoding="utf-8"))
    manifest_pairs = read_pair_manifest(pairs_manifest_path)
    raw_pairs = np.asarray(raw_meta.get("pairs", []), dtype=np.int64)
    if raw_pairs.shape != manifest_pairs.shape or not np.array_equal(
        raw_pairs, manifest_pairs
    ):
        raise ValueError("raw posterior pair order does not exactly match the input manifest")
    positions = np.asarray(raw_meta.get("output_positions", []), dtype=np.int64)
    if positions.ndim != 1 or not len(positions) or np.any(np.diff(positions) <= 0):
        raise ValueError("raw posterior positions are absent or not strictly ordered")
    n_pairs = int(raw_meta["num_pairs"])
    chunk_size = int(raw_meta["chunk_size"])
    n_positions = int(raw_meta["sequence_length"])
    if n_pairs != len(manifest_pairs) or n_positions != len(positions):
        raise ValueError("raw posterior dimensions do not match its metadata")
    if mutation_rate <= 0:
        raise ValueError("mutation rate must be positive")
    two_ne = float(raw_meta["scaled_mutation_rate"]) / (2.0 * mutation_rate)

    working_path = Path(working_path)
    working_path.parent.mkdir(parents=True, exist_ok=True)
    tmrca = np.memmap(
        working_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_positions, n_pairs),
    )
    bytes_per_chunk = 2 * n_positions * chunk_size * np.dtype(np.float32).itemsize
    with raw_path.open("rb") as compressed:
        with zstandard.ZstdDecompressor().stream_reader(compressed) as reader:
            for chunk_index, pair_start in enumerate(range(0, n_pairs, chunk_size)):
                raw = _read_exact(reader, bytes_per_chunk)
                values = np.frombuffer(raw, dtype=np.float32).reshape(
                    2, n_positions, chunk_size
                )
                pair_end = min(n_pairs, pair_start + chunk_size)
                real = pair_end - pair_start
                alpha = values[0, :, :real]
                beta = values[1, :, :real]
                valid = (
                    np.isfinite(alpha)
                    & np.isfinite(beta)
                    & (alpha > 0)
                    & (beta > 0)
                )
                means = np.full(alpha.shape, np.nan, dtype=np.float32)
                np.divide(alpha, beta, out=means, where=valid)
                means *= np.float32(two_ne)
                tmrca[:, pair_start:pair_end] = means
            if reader.read(1):
                raise ValueError("raw posterior stream has trailing decompressed bytes")
    tmrca.flush()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + f".tmp.{os.getpid()}")
    compressor = zstandard.ZstdCompressor(level=3, write_checksum=True)
    with temporary.open("wb") as destination:
        with compressor.stream_writer(destination, closefd=False) as writer:
            for start in range(0, n_positions, 64):
                block = np.ascontiguousarray(tmrca[start : start + 64])
                writer.write(memoryview(block).cast("B"))
    os.replace(temporary, output_path)

    metadata = {
        "format": "gamma_smc_candidate_tmrca_v1",
        "dtype": "float32",
        "layout": "position_major_then_pair",
        "shape": [n_positions, n_pairs],
        "positions_0based": positions.tolist(),
        "pair_order": "column k is data row k of pairs_manifest",
        "pairs_manifest": str(Path(pairs_manifest_path).resolve()),
        "pairs_manifest_sha256": sha256_file(pairs_manifest_path),
        "raw_metadata_sha256": sha256_file(raw_meta_path),
        "scaled_mutation_rate": float(raw_meta["scaled_mutation_rate"]),
        "unscaled_mutation_rate": float(mutation_rate),
        "two_ne_generations": float(two_ne),
        "array_sha256": sha256_file(output_path),
        "array_size_bytes": output_path.stat().st_size,
    }
    _atomic_text(metadata_path, json.dumps(metadata, indent=2) + "\n")
    return tmrca, metadata


def _sample_names_from_raw_meta(raw_path: Path) -> list[str]:
    metadata = json.loads(
        raw_path.with_name(raw_path.name + ".meta").read_text(encoding="utf-8")
    )
    labels = metadata.get("sample_names", {})
    if not isinstance(labels, dict) or len(labels) % 2:
        raise ValueError("raw metadata sample-name mapping is malformed")
    samples = []
    for haplotype in range(0, len(labels), 2):
        first = str(labels.get(str(haplotype), ""))
        second = str(labels.get(str(haplotype + 1), ""))
        if not first.endswith(".0") or not second.endswith(".1"):
            raise ValueError("raw metadata haplotype labels are malformed")
        if first[:-2] != second[:-2]:
            raise ValueError("raw metadata diploid haplotype labels disagree")
        samples.append(first[:-2])
    return samples


def _parse_genotypes(genotypes: list[str]) -> np.ndarray:
    alleles = np.full(2 * len(genotypes), -1, dtype=np.int16)
    for sample_index, genotype in enumerate(genotypes):
        genotype = genotype.split(":", 1)[0]
        separator = "|" if "|" in genotype else "/" if "/" in genotype else None
        if separator is None:
            continue
        fields = genotype.split(separator)
        if len(fields) != 2 or "." in fields:
            continue
        try:
            first, second = int(fields[0]), int(fields[1])
        except ValueError:
            continue
        if separator == "/" and first != second:
            continue
        alleles[2 * sample_index] = first
        alleles[2 * sample_index + 1] = second
    return alleles


def _r2(values_ref: np.ndarray, values_alt: np.ndarray) -> tuple[float, float]:
    def calculate(first: np.ndarray, second: np.ndarray) -> float:
        combined = np.r_[first, second]
        total = float(np.square(combined - combined.mean()).sum())
        if not np.isfinite(total) or total <= 0:
            return 0.0
        grand = float(combined.mean())
        between = float(
            len(first) * (float(first.mean()) - grand) ** 2
            + len(second) * (float(second.mean()) - grand) ** 2
        )
        return between / total

    return calculate(values_ref, values_alt), calculate(
        np.log1p(values_ref), np.log1p(values_alt)
    )


def _variant_type(ref: str, alt: str) -> str:
    if alt.startswith("<") or "[" in alt or "]" in alt:
        return "SV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"


def _nearest_index(positions: np.ndarray, position: int) -> int:
    right = int(np.searchsorted(positions, position))
    if right == 0:
        return 0
    if right == len(positions):
        return len(positions) - 1
    left = right - 1
    return left if position - positions[left] <= positions[right] - position else right


def _score_region_variants(
    *,
    bcftools: str | Path,
    bcf_path: str | Path,
    sample_list_path: str | Path,
    contig: str,
    start_0based: int,
    end_0based_exclusive: int,
    population: str,
    chromosome: int,
    region_id: str,
    tmrca: np.ndarray,
    tmrca_positions: np.ndarray,
    pairs: np.ndarray,
    n_samples: int,
    minimum_genotype_pairs: int,
) -> pd.DataFrame:
    command = [
        str(bcftools),
        "query",
        "--samples-file",
        str(sample_list_path),
        "--regions",
        f"{contig}:{start_0based + 1}-{end_0based_exclusive}",
        "--format",
        "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n",
        str(bcf_path),
    ]
    rows = []
    pair_i, pair_j = pairs[:, 0], pairs[:, 1]
    # A large warning stream must not fill a stderr pipe while stdout is streamed.
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as error_file:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=error_file,
            text=True,
        )
        assert process.stdout is not None
        for line_number, line in enumerate(process.stdout, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7 + n_samples:
                process.kill()
                raise ValueError(
                    f"bcftools genotype row {line_number} has {len(fields) - 7} samples; "
                    f"expected {n_samples}"
                )
            chrom, pos1_text, identifier, ref, alt_text, qual, filter_value = fields[:7]
            if not alt_text or alt_text == ".":
                continue
            pos1 = int(pos1_text)
            pos0 = pos1 - 1
            target_index = _nearest_index(tmrca_positions, pos0)
            target = np.asarray(tmrca[target_index], dtype=float)
            alleles = _parse_genotypes(fields[7:])
            called_haplotypes = alleles[alleles >= 0]
            if not len(called_haplotypes):
                continue
            first = alleles[pair_i]
            second = alleles[pair_j]
            finite = np.isfinite(target) & (target >= 0)
            for alt_index, alt in enumerate(alt_text.split(","), 1):
                ref_selector = finite & (first == 0) & (second == 0)
                alt_selector = finite & (first == alt_index) & (second == alt_index)
                n_ref = int(ref_selector.sum())
                n_alt = int(alt_selector.sum())
                if min(n_ref, n_alt) < minimum_genotype_pairs:
                    continue
                values_ref = target[ref_selector]
                values_alt = target[alt_selector]
                r2, log_r2 = _r2(values_ref, values_alt)
                rows.append({
                    "population": population,
                    "chromosome": int(chromosome),
                    "region_id": region_id,
                    "chrom": chrom,
                    "position_0based": pos0,
                    "position_1based": pos1,
                    "id": identifier,
                    "ref": ref,
                    "alt": alt,
                    "alt_index": int(alt_index),
                    "variant_type": _variant_type(ref, alt),
                    "qual": qual,
                    "filter": filter_value,
                    "alternate_allele_frequency": float(
                        np.count_nonzero(called_haplotypes == alt_index)
                        / len(called_haplotypes)
                    ),
                    "tmrca_position_0based": int(tmrca_positions[target_index]),
                    "n_hom_ref_pairs": n_ref,
                    "n_hom_alt_pairs": n_alt,
                    "mean_hom_ref_tmrca_generations": float(values_ref.mean()),
                    "mean_hom_alt_tmrca_generations": float(values_alt.mean()),
                    "hom_alt_minus_hom_ref_generations": float(
                        values_alt.mean() - values_ref.mean()
                    ),
                    "r2_tmrca": float(r2),
                    "r2_log1p_tmrca": float(log_r2),
                })
        return_code = process.wait()
        error_file.seek(0)
        stderr = error_file.read()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command, stderr=stderr)
    return pd.DataFrame(rows, columns=VARIANT_COLUMNS)


def _confidence_profile(
    values: np.ndarray,
    positions: np.ndarray,
    selector: np.ndarray,
    *,
    population: str,
    chromosome: int,
    region_id: str,
    representative: pd.Series,
    pair_class: str,
) -> pd.DataFrame:
    rows = []
    for index, position in enumerate(positions):
        observed = np.asarray(values[index, selector], dtype=float)
        observed = observed[np.isfinite(observed) & (observed > 0)]
        n_pairs = int(len(observed))
        if not n_pairs:
            continue
        mean = float(observed.mean())
        median = float(np.median(observed))
        if n_pairs > 1:
            standard_error = float(observed.std(ddof=1) / np.sqrt(n_pairs))
            log_half_width = float(
                t.ppf(0.975, df=n_pairs - 1) * standard_error / mean
            )
            lower = mean * np.exp(-log_half_width)
            upper = mean * np.exp(log_half_width)
        else:
            standard_error = lower = upper = np.nan
        rows.append({
            "population": population,
            "chromosome": int(chromosome),
            "region_id": region_id,
            "representative_position_1based": int(representative["position_1based"]),
            "representative_ref": representative["ref"],
            "representative_alt": representative["alt"],
            "position_0based": int(position),
            "position_1based": int(position) + 1,
            "pair_class": pair_class,
            "n_pairs": n_pairs,
            "mean_tmrca_generations": mean,
            "median_tmrca_generations": median,
            "standard_error_generations": standard_error,
            "ci95_lower_generations": lower,
            "ci95_upper_generations": upper,
        })
    return pd.DataFrame(rows, columns=PROFILE_COLUMNS)


def _plot_pair_distribution(
    summary: pd.DataFrame,
    distribution: pd.DataFrame,
    *,
    called_column: str,
    population: str,
    chromosome: int,
    region_id: str,
    peak: int,
    threshold_years: float,
    minimum_fraction: float,
    output_stem: Path,
) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True, constrained_layout=True)
    x = summary["position_0based"].to_numpy(dtype=float) / 1e6
    axes[0].plot(x, summary[called_column], color="#5b21b6", lw=1.4)
    axes[0].axhline(
        minimum_fraction,
        color="0.35",
        lw=0.9,
        ls=":",
        label=f"{minimum_fraction:.1%} screen",
    )
    axes[0].set_ylabel(
        f"Fraction posterior-mean TMRCA < {threshold_years:,.0f} years"
    )
    axes[0].legend(loc="best")
    qx = distribution["position_0based"].to_numpy(dtype=float) / 1e6
    axes[1].fill_between(
        qx,
        distribution["tmrca_p05_generations"],
        distribution["tmrca_p95_generations"],
        color="#bfdbfe",
        alpha=0.45,
        label="5th-95th percentile",
    )
    axes[1].fill_between(
        qx,
        distribution["tmrca_p25_generations"],
        distribution["tmrca_p75_generations"],
        color="#60a5fa",
        alpha=0.45,
        label="25th-75th percentile",
    )
    axes[1].plot(
        qx,
        distribution["tmrca_median_generations"],
        color="#1d4ed8",
        lw=1.4,
        label="pair median",
    )
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Posterior-mean pair TMRCA (generations)")
    axes[1].set_xlabel(f"Chromosome {chromosome} position (Mb)")
    axes[1].legend(loc="best")
    for axis in axes:
        axis.axvline(peak / 1e6, color="#e69f00", lw=1.0, ls="--")
        axis.grid(alpha=0.16, linewidth=0.6)
    figure.suptitle(f"{population} {region_id}: decoded-pair TMRCA around scan peak")
    figure.savefig(Path(f"{output_stem}.png"), dpi=190, bbox_inches="tight")
    figure.savefig(Path(f"{output_stem}.pdf"), bbox_inches="tight")
    plt.close(figure)


def _plot_representative_variant(
    scores: pd.DataFrame,
    profile: pd.DataFrame,
    representative: pd.Series,
    *,
    population: str,
    chromosome: int,
    region_id: str,
    peak: int,
    output_stem: Path,
) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=False, constrained_layout=True)
    axes[0].scatter(
        scores["position_0based"] / 1e6,
        scores["r2_tmrca"],
        s=9,
        alpha=0.55,
        color="#6b7280",
        rasterized=True,
    )
    axes[0].scatter(
        [float(representative["position_0based"]) / 1e6],
        [float(representative["r2_tmrca"])],
        s=55,
        color="#d97706",
        label="representative variant",
        zorder=4,
    )
    axes[0].set_ylabel("Variance explained in pair TMRCA (R²)")
    axes[0].set_xlabel(f"Chromosome {chromosome} position (Mb)")
    axes[0].legend(loc="best")
    styles = {
        "hom_ref": ("hom ref", "#595959", "#bdbdbd"),
        "hom_alt": ("hom alt", "#0072b2", "#56b4e9"),
    }
    for pair_class, (label, line_color, fill_color) in styles.items():
        group = profile.loc[profile["pair_class"].eq(pair_class)].sort_values(
            "position_0based"
        )
        x = group["position_0based"].to_numpy(dtype=float) / 1e6
        axes[1].fill_between(
            x,
            group["ci95_lower_generations"],
            group["ci95_upper_generations"],
            color=fill_color,
            alpha=0.28,
            linewidth=0,
        )
        axes[1].plot(
            x,
            group["mean_tmrca_generations"],
            color=line_color,
            lw=1.4,
            label=f"{label} mean",
        )
    axes[1].axvline(
        float(representative["position_0based"]) / 1e6,
        color="#d97706",
        lw=1.0,
        ls="--",
        label="representative variant",
    )
    axes[1].axvline(peak / 1e6, color="#7c3aed", lw=1.0, ls=":", label="scan peak")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Mean posterior pair TMRCA (generations)")
    axes[1].set_xlabel(f"Chromosome {chromosome} position (Mb)")
    axes[1].legend(loc="best")
    for axis in axes:
        axis.grid(alpha=0.16, linewidth=0.6)
    figure.suptitle(
        f"{population} {region_id}: {representative['chrom']}:{int(representative['position_1based'])} "
        f"{representative['ref']}>{representative['alt']}\n"
        f"hom-alt vs hom-ref among the exact decoded pair set"
    )
    figure.savefig(Path(f"{output_stem}.png"), dpi=190, bbox_inches="tight")
    figure.savefig(Path(f"{output_stem}.pdf"), bbox_inches="tight")
    plt.close(figure)


def _artifact_records(output_dir: Path, excluded: set[Path]) -> list[dict]:
    records = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file() or path in excluded:
            continue
        records.append({
            "path": str(path.relative_to(output_dir)),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        })
    return records


def analyze_candidate_regions(
    *,
    summary_path: str | Path,
    regions_path: str | Path,
    raw_path: str | Path,
    pairs_manifest_path: str | Path,
    bcf_path: str | Path,
    sample_list_path: str | Path,
    bcftools: str | Path,
    contig: str,
    population: str,
    chromosome: int,
    output_dir: str | Path,
    mutation_rate: float = 1.29e-8,
    threshold_years: float = 4500,
    profile_half_width: int = 500_000,
    variant_half_width: int = 100_000,
    minimum_genotype_pairs: int = 20,
    minimum_fraction: float = 0.02,
    merge_gap: int = 20_000,
) -> dict:
    population = population.upper()
    if population not in POPULATIONS:
        raise ValueError(f"unsupported population: {population}")
    if min(profile_half_width, variant_half_width, minimum_genotype_pairs) < 1:
        raise ValueError("profile/variant widths and minimum pair count must be positive")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    scores_dir = output_dir / "variant_scores"
    figures_dir.mkdir(parents=True, exist_ok=True)
    scores_dir.mkdir(parents=True, exist_ok=True)
    regions = pd.read_csv(regions_path, sep="\t")
    if regions.empty:
        raise ValueError("candidate analysis was invoked with no regions")
    if set(REGION_COLUMNS).difference(regions.columns):
        raise ValueError("candidate-region table has an unexpected schema")
    summary, _ = validate_summary(summary_path, threshold_years=threshold_years)
    called_column = called_fraction_column(threshold_years)
    sample_names = [
        value.strip()
        for value in Path(sample_list_path).read_text(encoding="utf-8").splitlines()
        if value.strip()
    ]
    raw_sample_names = _sample_names_from_raw_meta(Path(raw_path))
    if raw_sample_names != sample_names:
        raise ValueError("raw posterior sample order does not match the BCF sample list")
    pairs = read_pair_manifest(pairs_manifest_path)
    if int(pairs.max()) >= 2 * len(sample_names):
        raise ValueError("pair manifest references a haplotype outside the sample list")

    ordered_array = output_dir / "candidate_tmrca.f32.zst"
    ordered_meta = output_dir / "candidate_tmrca.meta.json"
    working_array = output_dir / f".candidate_tmrca.{os.getpid()}.f32"
    tmrca, tmrca_meta = convert_raw_to_ordered_tmrca(
        raw_path,
        pairs_manifest_path=pairs_manifest_path,
        mutation_rate=mutation_rate,
        output_path=ordered_array,
        metadata_path=ordered_meta,
        working_path=working_array,
    )
    positions = np.asarray(tmrca_meta["positions_0based"], dtype=np.int64)
    representatives = []
    profiles = []
    distributions = []
    try:
        for region in regions.itertuples(index=False):
            peak = int(region.peak_position_0based)
            profile_start = max(0, peak - profile_half_width)
            profile_end = peak + profile_half_width
            position_selector = np.flatnonzero(
                (positions >= profile_start) & (positions <= profile_end)
            )
            if not len(position_selector):
                raise ValueError(f"no pair TMRCA positions overlap {region.region_id}")
            region_values = np.asarray(tmrca[position_selector], dtype=float)
            region_positions = positions[position_selector]
            quantiles = np.nanpercentile(region_values, [5, 25, 50, 75, 95], axis=1)
            distribution = pd.DataFrame({
                "population": population,
                "chromosome": int(chromosome),
                "region_id": region.region_id,
                "position_0based": region_positions,
                "position_1based": region_positions + 1,
                "n_pairs": np.isfinite(region_values).sum(axis=1),
                "tmrca_p05_generations": quantiles[0],
                "tmrca_p25_generations": quantiles[1],
                "tmrca_median_generations": quantiles[2],
                "tmrca_p75_generations": quantiles[3],
                "tmrca_p95_generations": quantiles[4],
            }, columns=DISTRIBUTION_COLUMNS)
            distributions.append(distribution)
            summary_region = summary.loc[
                summary["position_0based"].between(
                    profile_start, profile_end, inclusive="both"
                )
            ]
            _plot_pair_distribution(
                summary_region,
                distribution,
                called_column=called_column,
                population=population,
                chromosome=chromosome,
                region_id=region.region_id,
                peak=peak,
                threshold_years=threshold_years,
                minimum_fraction=minimum_fraction,
                output_stem=figures_dir / f"{region.region_id}.pair_tmrca_distribution",
            )

            variant_start = max(0, peak - variant_half_width)
            variant_end = peak + variant_half_width + 1
            scores = _score_region_variants(
                bcftools=bcftools,
                bcf_path=bcf_path,
                sample_list_path=sample_list_path,
                contig=contig,
                start_0based=variant_start,
                end_0based_exclusive=variant_end,
                population=population,
                chromosome=chromosome,
                region_id=region.region_id,
                tmrca=tmrca,
                tmrca_positions=positions,
                pairs=pairs,
                n_samples=len(sample_names),
                minimum_genotype_pairs=minimum_genotype_pairs,
            )
            score_path = scores_dir / f"{region.region_id}.variant_scores.tsv.gz"
            if scores.empty:
                _atomic_frame(score_path, scores, gzip=True)
                continue
            scores = scores.sort_values(
                ["r2_tmrca", "n_hom_alt_pairs", "position_0based", "alt_index"],
                ascending=[False, False, True, True],
            ).reset_index(drop=True)
            _atomic_frame(score_path, scores, gzip=True)
            representative = scores.iloc[0].copy()
            representatives.append(representative)

            # Re-query only the representative record to reconstruct phased alleles.
            query = [
                str(bcftools),
                "query",
                "--samples-file",
                str(sample_list_path),
                "--regions",
                f"{contig}:{int(representative['position_1based'])}-{int(representative['position_1based'])}",
                "--format",
                "%POS\t%REF\t%ALT[\t%GT]\n",
                str(bcf_path),
            ]
            completed = subprocess.run(query, check=True, capture_output=True, text=True)
            matching_alleles = None
            for line in completed.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) != 3 + len(sample_names):
                    continue
                if int(fields[0]) != int(representative["position_1based"]):
                    continue
                if fields[1] != representative["ref"]:
                    continue
                alternate = str(representative["alt"])
                alternatives = fields[2].split(",")
                if int(representative["alt_index"]) > len(alternatives):
                    continue
                if alternatives[int(representative["alt_index"]) - 1] != alternate:
                    continue
                matching_alleles = _parse_genotypes(fields[3:])
                break
            if matching_alleles is None:
                raise ValueError("could not re-query the representative variant")
            pair_i, pair_j = pairs[:, 0], pairs[:, 1]
            first, second = matching_alleles[pair_i], matching_alleles[pair_j]
            alt_index = int(representative["alt_index"])
            ref_selector = (first == 0) & (second == 0)
            alt_selector = (first == alt_index) & (second == alt_index)
            profile = pd.concat([
                _confidence_profile(
                    region_values,
                    region_positions,
                    ref_selector,
                    population=population,
                    chromosome=chromosome,
                    region_id=region.region_id,
                    representative=representative,
                    pair_class="hom_ref",
                ),
                _confidence_profile(
                    region_values,
                    region_positions,
                    alt_selector,
                    population=population,
                    chromosome=chromosome,
                    region_id=region.region_id,
                    representative=representative,
                    pair_class="hom_alt",
                ),
            ], ignore_index=True)
            profiles.append(profile)
            _plot_representative_variant(
                scores,
                profile,
                representative,
                population=population,
                chromosome=chromosome,
                region_id=region.region_id,
                peak=peak,
                output_stem=figures_dir / f"{region.region_id}.hom_alt_vs_hom_ref",
            )
    finally:
        mapping = tmrca._mmap
        del tmrca
        mapping.close()
        working_array.unlink(missing_ok=True)

    representative_frame = pd.DataFrame(representatives, columns=VARIANT_COLUMNS)
    profile_frame = (
        pd.concat(profiles, ignore_index=True)
        if profiles
        else pd.DataFrame(columns=PROFILE_COLUMNS)
    )
    distribution_frame = pd.concat(distributions, ignore_index=True)
    _atomic_frame(output_dir / "representative_variants.tsv", representative_frame)
    _atomic_frame(output_dir / "pair_tmrca_profiles.tsv", profile_frame)
    _atomic_frame(output_dir / "pair_tmrca_distributions.tsv", distribution_frame)
    manifest_path = output_dir / "candidate_analysis.json"
    result = {
        "schema_version": 1,
        "population": population,
        "chromosome": int(chromosome),
        "screen": {
            "called_fraction_column": called_column,
            "strictly_greater_than": float(minimum_fraction),
            "merge_gap_bp": int(merge_gap),
            "profile_half_width_bp": int(profile_half_width),
            "variant_half_width_bp": int(variant_half_width),
            "minimum_hom_genotype_pairs": int(minimum_genotype_pairs),
        },
        "pair_contract": {
            "n_pairs": int(len(pairs)),
            "pairs_manifest": str(Path(pairs_manifest_path).resolve()),
            "pairs_manifest_sha256": sha256_file(pairs_manifest_path),
            "tmrca_array_order": "column k is data row k of pairs_manifest",
        },
        "n_regions": int(len(regions)),
        "n_regions_with_representative_variant": int(len(representative_frame)),
        "regions_sha256": sha256_file(regions_path),
        "artifacts": [],
    }
    result["artifacts"] = _artifact_records(output_dir, {manifest_path})
    _atomic_text(manifest_path, json.dumps(result, indent=2) + "\n")
    return result


def write_empty_candidate_analysis(
    *,
    regions_path: str | Path,
    population: str,
    chromosome: int,
    pairs_manifest_path: str | Path,
    output_dir: str | Path,
) -> dict:
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    _atomic_frame(output_dir / "representative_variants.tsv", pd.DataFrame(columns=VARIANT_COLUMNS))
    _atomic_frame(output_dir / "pair_tmrca_profiles.tsv", pd.DataFrame(columns=PROFILE_COLUMNS))
    _atomic_frame(
        output_dir / "pair_tmrca_distributions.tsv",
        pd.DataFrame(columns=DISTRIBUTION_COLUMNS),
    )
    manifest_path = output_dir / "candidate_analysis.json"
    result = {
        "schema_version": 1,
        "population": population.upper(),
        "chromosome": int(chromosome),
        "n_regions": 0,
        "n_regions_with_representative_variant": 0,
        "regions_sha256": sha256_file(regions_path),
        "pair_contract": {
            "pairs_manifest_sha256": sha256_file(pairs_manifest_path),
            "tmrca_array_order": "not emitted because no region exceeded the screen",
        },
        "artifacts": [],
    }
    result["artifacts"] = _artifact_records(output_dir, {manifest_path})
    _atomic_text(manifest_path, json.dumps(result, indent=2) + "\n")
    return result
