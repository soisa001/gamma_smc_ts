"""run6: a spatial, genotype-stratified TMRCA scan.

run5 found that the within-individual TMRCA at a swept introgressed locus is
trimodal, and that the three genotype classes pull in opposite directions:
hom-carrier pairs are far *younger* than neutral, heterozygous pairs are far
*older* and can never be recent, and hom-non-carrier pairs are neutral-like.
Pooling them into a single ``P(TMRCA < x)`` throws that structure away, and below
a carrier frequency of 2/3 it actively inverts the signal.

This module asks the natural follow-up: **at each 10 kb stride, does the
carrier class look much younger than a neutrally evolving region?**  It scans
the whole contig, stratifying by the genotype at the focal site, and compares
each class against a null built from the neutral replicates' *unstratified*
statistic -- which is the right reference because a neutral region has no allele
and therefore no classes.

Two things make the comparison fair:

* **Matched pair counts.** A class holding 30 pairs has a noisier statistic than
  one holding 100. Since the statistic is a proportion out of a known total, a
  matched-size null is drawn exactly by hypergeometric sampling from each
  neutral replicate's counts rather than approximately by resampling.
* **A pooled marginal null.** Under a constant demography the neutral statistic
  is stationary along the genome, so pooling every neutral position and replicate
  estimates the marginal null far more finely than 100 replicates alone. Note
  this is the *marginal* null: positions within a replicate are linked, so it
  does not license genome-wide multiple-testing claims.
"""

from __future__ import annotations

import glob
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

CLASSES = ("hom_carrier", "het", "hom_noncarrier", "overall")
_CLASS_CODE = {"hom_carrier": 2, "het": 1, "hom_noncarrier": 0}
DEFAULT_THRESHOLDS_YEARS: tuple[float, ...] = (
    1_000.0,
    4_500.0,
    10_000.0,
    20_000.0,
    30_000.0,
    50_000.0,
    100_000.0,
)
#: A class smaller than this is too noisy to score at a single position.
MIN_CLASS_PAIRS = 10


@dataclass(frozen=True)
class Replicate:
    positions: np.ndarray
    tmrca: np.ndarray  # (n_positions, n_pairs), generations
    genotypes: np.ndarray | None  # (n_pairs,) in {0, 1, 2}
    census_af: float
    replicate_index: int


def load_replicates(arm_dir: str | Path, mode: str) -> list[Replicate]:
    """Load the stored raw TMRCA matrices and focal genotypes for one mode."""
    records: list[Replicate] = []
    pattern = str(Path(arm_dir) / mode / "replicates" / "*")
    for path in sorted(glob.glob(pattern)):
        directory = Path(path)
        matrix = directory / "tmrca_matrix.npz"
        endpoint = directory / "endpoint.json"
        if not (matrix.is_file() and endpoint.is_file()):
            continue
        record = json.loads(endpoint.read_text(encoding="utf-8"))
        if record.get("status") != "completed":
            continue
        with np.load(matrix) as data:
            positions = np.asarray(data["positions"], dtype=float)
            tmrca = np.asarray(data["tmrca_generations"], dtype=np.float32)
        genotype_path = directory / "focal_genotypes.npy"
        genotypes = np.load(genotype_path) if genotype_path.is_file() else None
        census = record.get("final_allele_frequency", {}).get("census_af")
        records.append(
            Replicate(
                positions=positions,
                tmrca=tmrca,
                genotypes=genotypes,
                census_af=float(census) if census is not None else float("nan"),
                replicate_index=int(record["replicate_index"]),
            )
        )
    return records


def _below_counts(tmrca: np.ndarray, generations: float) -> np.ndarray:
    """Per-position count of pairs coalescing more recently than the cutoff."""
    return (tmrca < generations).sum(axis=1)


def build_matched_null(
    neutral: Sequence[Replicate],
    generations: float,
    n_pairs: int,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Pooled neutral statistic, drawn at a matched pair count.

    The statistic is a proportion out of a known total, so a class of ``n_pairs``
    is emulated exactly by hypergeometric sampling from each neutral position's
    counts -- no approximation and no need to revisit the raw matrices.
    """
    draws = []
    for replicate in neutral:
        total = replicate.tmrca.shape[1]
        if n_pairs > total:
            continue
        counts = _below_counts(replicate.tmrca, generations)
        sampled = rng.hypergeometric(counts, total - counts, n_pairs)
        draws.append(sampled / float(n_pairs))
    if not draws:
        return np.empty(0, dtype=float)
    return np.sort(np.concatenate(draws))


def _percentile(values: np.ndarray, sorted_null: np.ndarray) -> np.ndarray:
    """Fraction of the null below each value, ties split."""
    if sorted_null.size == 0:
        return np.full(values.shape, np.nan)
    left = np.searchsorted(sorted_null, values, side="left")
    right = np.searchsorted(sorted_null, values, side="right")
    return (left + 0.5 * (right - left)) / sorted_null.size


def scan_arm(
    arm_dir: str | Path,
    *,
    generation_time: float,
    thresholds_years: Sequence[float] = DEFAULT_THRESHOLDS_YEARS,
    classes: Iterable[str] = CLASSES,
    seed: int = 20260823,
) -> pd.DataFrame:
    """Return the per-position, per-class, per-threshold scan for one arm."""
    arm_dir = Path(arm_dir)
    selected = load_replicates(arm_dir, "selected")
    neutral = load_replicates(arm_dir, "neutral")
    if not selected or not neutral:
        raise RuntimeError(f"missing replicates under {arm_dir}")

    rng = np.random.default_rng(seed)
    positions = selected[0].positions
    n_positions = positions.size
    rows: list[dict[str, Any]] = []

    for years in thresholds_years:
        generations = float(years) / generation_time
        null_cache: dict[int, np.ndarray] = {}

        for name in classes:
            statistic = np.full((len(selected), n_positions), np.nan)
            percentile = np.full((len(selected), n_positions), np.nan)

            for row, replicate in enumerate(selected):
                if name == "overall":
                    mask = np.ones(replicate.tmrca.shape[1], dtype=bool)
                elif replicate.genotypes is None:
                    continue
                else:
                    mask = replicate.genotypes == _CLASS_CODE[name]
                size = int(mask.sum())
                if size < MIN_CLASS_PAIRS:
                    continue
                values = (replicate.tmrca[:, mask] < generations).mean(axis=1)
                statistic[row] = values
                if size not in null_cache:
                    null_cache[size] = build_matched_null(
                        neutral, generations, size, rng=rng
                    )
                percentile[row] = _percentile(values, null_cache[size])

            usable = np.isfinite(statistic).any(axis=1)
            if not usable.any():
                continue
            rows.append(
                pd.DataFrame(
                    {
                        "position_0based": positions,
                        "threshold_years": float(years),
                        "genotype_class": name,
                        "n_replicates": int(usable.sum()),
                        "mean_statistic": np.nanmean(statistic[usable], axis=0),
                        "median_statistic": np.nanmedian(statistic[usable], axis=0),
                        "auc_vs_neutral": np.nanmean(percentile[usable], axis=0),
                    }
                )
            )

    # The neutral arm's own unstratified profile, for reference.
    for years in thresholds_years:
        generations = float(years) / generation_time
        values = np.vstack(
            [(r.tmrca < generations).mean(axis=1) for r in neutral]
        )
        rows.append(
            pd.DataFrame(
                {
                    "position_0based": positions,
                    "threshold_years": float(years),
                    "genotype_class": "neutral_reference",
                    "n_replicates": len(neutral),
                    "mean_statistic": values.mean(axis=0),
                    "median_statistic": np.median(values, axis=0),
                    "auc_vs_neutral": np.full(n_positions, 0.5),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def detection_width(
    scan: pd.DataFrame,
    *,
    focal_position: float,
    genotype_class: str = "hom_carrier",
    auc_floor: float = 0.75,
) -> pd.DataFrame:
    """Half-width of the region whose AUC stays above ``auc_floor``.

    Measured as the distance from the focal base to the first position, moving
    outwards in each direction, where the AUC drops below the floor.
    """
    rows = []
    subset = scan[scan["genotype_class"] == genotype_class]
    for years, group in subset.groupby("threshold_years"):
        group = group.sort_values("position_0based")
        pos = group["position_0based"].to_numpy()
        auc = group["auc_vs_neutral"].to_numpy()
        focal_index = int(np.argmin(np.abs(pos - focal_position)))
        if not np.isfinite(auc[focal_index]):
            continue

        left = focal_index
        while left > 0 and np.isfinite(auc[left - 1]) and auc[left - 1] >= auc_floor:
            left -= 1
        right = focal_index
        while (
            right < len(pos) - 1
            and np.isfinite(auc[right + 1])
            and auc[right + 1] >= auc_floor
        ):
            right += 1
        rows.append(
            {
                "threshold_years": float(years),
                "genotype_class": genotype_class,
                "auc_floor": auc_floor,
                "auc_at_focal": float(auc[focal_index]),
                "left_bp": float(focal_position - pos[left]),
                "right_bp": float(pos[right] - focal_position),
                "width_bp": float(pos[right] - pos[left]),
            }
        )
    return pd.DataFrame(rows)


__all__ = [
    "CLASSES",
    "DEFAULT_THRESHOLDS_YEARS",
    "MIN_CLASS_PAIRS",
    "Replicate",
    "build_matched_null",
    "detection_width",
    "load_replicates",
    "scan_arm",
]
