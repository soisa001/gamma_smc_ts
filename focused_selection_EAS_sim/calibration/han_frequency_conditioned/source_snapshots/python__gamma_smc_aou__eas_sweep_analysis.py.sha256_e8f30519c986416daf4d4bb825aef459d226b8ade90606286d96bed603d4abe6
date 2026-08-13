"""Analysis primitives for the EAS selected-sweep simulation studies.

The functions in this module deliberately describe *internal genotype
contrasts*.  A normalized CDF area compares within-individual TMRCA
distributions for hom-alt and hom-ref simulated samples; it is not a ROC AUC,
and no neutral-null p-value is calculated here.
"""

from __future__ import annotations

import os
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import PercentFormatter


DEFAULT_GENERATION_TIME_YEARS = 25.0
DEFAULT_THRESHOLDS_YEARS = (1_000, 4_500, 10_000, 20_000, 30_000, 40_000, 50_000)
GENOTYPE_CLASSES = ("overall", "hom_ref", "heterozygous", "hom_alt")

_GENOTYPE_CLASS_BY_COUNT = {0: "hom_ref", 1: "heterozygous", 2: "hom_alt"}
_CLASS_COLORS = {
    "overall": "#000000",
    "hom_ref": "#595959",
    "heterozygous": "#e69f00",
    "hom_alt": "#0072b2",
}
_CLASS_LABELS = {
    "overall": "Overall",
    "hom_ref": "Ref/ref",
    "heterozygous": "Ref/alt",
    "hom_alt": "Alt/alt",
}
_SOURCE_LINESTYLES = ("-", "--", ":", "-.")
_PLOT_FONT = {"title": 20, "axis": 18, "tick": 15, "legend": 13}


def _sampled_diploids(ts) -> list[tuple[int, tuple[int, int]]]:
    """Return sampled individuals and nodes in tskit's VCF individual order."""
    sample_nodes = {int(node) for node in ts.samples()}
    records: list[tuple[int, tuple[int, int]]] = []
    assigned: list[int] = []
    for individual in ts.individuals():
        nodes = tuple(
            int(node) for node in individual.nodes if int(node) in sample_nodes
        )
        if not nodes:
            continue
        if len(nodes) != 2:
            raise ValueError(
                f"individual {individual.id} has {len(nodes)} sampled nodes; expected 2"
            )
        records.append((int(individual.id), (nodes[0], nodes[1])))
        assigned.extend(nodes)
    if not records:
        raise ValueError("tree sequence has no sampled diploid individuals")
    if len(assigned) != len(set(assigned)):
        raise ValueError("a sample node is assigned to more than one diploid")
    if set(assigned) != sample_nodes:
        missing = sorted(sample_nodes.difference(assigned))
        raise ValueError(
            "every sample node must belong to a diploid individual; "
            f"unassigned sample nodes: {missing[:5]}"
        )
    return records


def _site_at_coordinate(ts, focal_coordinate: float):
    matches = [
        site
        for site in ts.sites()
        if np.isclose(
            float(site.position),
            float(focal_coordinate),
            rtol=0.0,
            atol=1e-9,
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one focal site at {focal_coordinate:g}; found {len(matches)}"
        )
    return matches[0]


def _focal_genotype_counts_from_site(
    ts,
    focal_coordinate: float,
    *,
    selected_allele_index: int = 1,
) -> np.ndarray:
    records = _sampled_diploids(ts)
    site = _site_at_coordinate(ts, focal_coordinate)
    ordered_nodes = [node for _, nodes in records for node in nodes]
    variant = tskit.Variant(ts, samples=ordered_nodes)
    variant.decode(site.id)

    allele_index = int(selected_allele_index)
    if allele_index < 0 or allele_index >= len(variant.alleles):
        raise ValueError(
            f"selected allele index {allele_index} is absent from focal alleles "
            f"{variant.alleles}"
        )
    genotypes = np.asarray(variant.genotypes, dtype=np.int64)
    if np.any(genotypes < 0):
        raise ValueError("missing focal genotypes cannot be assigned to pair classes")
    counts = (genotypes.reshape(len(records), 2) == allele_index).sum(axis=1)
    return counts.astype(np.int8, copy=False)


def focal_genotype_counts(
    ts,
    focal_coordinate: float,
    *,
    selected_allele_index: int = 1,
) -> np.ndarray:
    """Count focal selected alleles for sampled diploids in VCF/Gamma order.

    The returned array has one value (0, 1, or 2) per sampled diploid.  The
    ordering is the individual ordering used by ``TreeSequence.write_vcf`` and
    therefore by Gamma-SMC's consecutive haplotype numbering.
    """
    return _focal_genotype_counts_from_site(
        ts,
        focal_coordinate,
        selected_allele_index=selected_allele_index,
    )


def _coerce_focal_counts(
    ts,
    records: Sequence[tuple[int, tuple[int, int]]],
    counts,
) -> np.ndarray:
    if isinstance(counts, Mapping):
        missing = [
            individual_id for individual_id, _ in records if individual_id not in counts
        ]
        if missing:
            raise ValueError(
                "focal genotype counts do not cover sampled individual IDs "
                f"{missing[:5]}"
            )
        array = np.asarray([counts[individual_id] for individual_id, _ in records])
    else:
        array = np.asarray(counts)
        if array.ndim != 1:
            raise ValueError("focal genotype counts must be one-dimensional")
        if len(array) == ts.num_individuals and len(array) != len(records):
            array = np.asarray([array[individual_id] for individual_id, _ in records])
        elif len(array) != len(records):
            raise ValueError(
                "focal genotype counts must have one value per sampled diploid "
                f"({len(records)} values required; found {len(array)})"
            )
    numeric = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise ValueError("focal genotype counts must be finite integers")
    integers = numeric.astype(np.int8)
    if not np.isin(integers, [0, 1, 2]).all():
        raise ValueError("focal genotype counts must contain only 0, 1, or 2")
    return integers


def build_diploid_pair_table(
    ts,
    focal_coordinate: float,
    focal_genotype_counts: Sequence[int] | Mapping[int, int] | None = None,
    *,
    selected_allele_index: int = 1,
) -> pd.DataFrame:
    """Map sampled diploids to tree nodes and deterministic Gamma haplotypes.

    If ``focal_genotype_counts`` is omitted, counts are read from the unique
    tree-sequence site at ``focal_coordinate``.  Explicit arrays are interpreted
    in sampled VCF diploid order (or in individual-ID order when the tree
    sequence also contains nonsampled individuals).  A mapping is keyed by
    tree-sequence individual ID.
    """
    coordinate = float(focal_coordinate)
    if not 0 <= coordinate < float(ts.sequence_length):
        raise ValueError("focal coordinate lies outside the tree sequence")
    records = _sampled_diploids(ts)
    if focal_genotype_counts is None:
        counts = _focal_genotype_counts_from_site(
            ts,
            coordinate,
            selected_allele_index=selected_allele_index,
        )
    else:
        counts = _coerce_focal_counts(ts, records, focal_genotype_counts)

    rows = []
    for vcf_index, ((individual_id, nodes), count) in enumerate(
        zip(records, counts, strict=True)
    ):
        rows.append(
            {
                "vcf_diploid_index": vcf_index,
                "tree_sequence_individual_id": individual_id,
                "sample_node_0": nodes[0],
                "sample_node_1": nodes[1],
                "gamma_smc_haplotype_0": 2 * vcf_index,
                "gamma_smc_haplotype_1": 2 * vcf_index + 1,
                "focal_position_0based": coordinate,
                "focal_selected_allele_index": int(selected_allele_index),
                "focal_selected_allele_count": int(count),
                "genotype_class": _GENOTYPE_CLASS_BY_COUNT[int(count)],
            }
        )
    return pd.DataFrame(rows)


def _validate_pair_table(pair_table: pd.DataFrame) -> None:
    required = {
        "sample_node_0",
        "sample_node_1",
        "gamma_smc_haplotype_0",
        "gamma_smc_haplotype_1",
        "genotype_class",
    }
    missing = sorted(required.difference(pair_table.columns))
    if missing:
        raise ValueError(f"pair table is missing columns: {', '.join(missing)}")
    observed = set(pair_table["genotype_class"].astype(str))
    unexpected = sorted(observed.difference(GENOTYPE_CLASSES[1:]))
    if unexpected:
        raise ValueError(f"unexpected genotype classes: {', '.join(unexpected)}")
    pairs = pair_table[["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]].to_numpy(
        dtype=np.int64
    )
    if len(pairs) and np.any(pairs[:, 0] >= pairs[:, 1]):
        raise ValueError("Gamma-SMC pairs must have increasing haplotype indices")
    if len({tuple(pair) for pair in pairs}) != len(pairs):
        raise ValueError("Gamma-SMC pair table contains duplicate pairs")


def gamma_pair_tables(
    pair_table: pd.DataFrame,
    *,
    require_nonempty: bool = True,
) -> dict[str, pd.DataFrame]:
    """Return two-column, header-free-ready Gamma pair tables for each class."""
    _validate_pair_table(pair_table)
    result: dict[str, pd.DataFrame] = {}
    for genotype_class in GENOTYPE_CLASSES:
        selected = (
            pair_table
            if genotype_class == "overall"
            else pair_table[pair_table["genotype_class"] == genotype_class]
        )
        if require_nonempty and selected.empty:
            raise ValueError(
                f"no {genotype_class} diploid pairs; reject this simulation draw"
            )
        result[genotype_class] = (
            selected[["gamma_smc_haplotype_0", "gamma_smc_haplotype_1"]]
            .rename(
                columns={
                    "gamma_smc_haplotype_0": "haplotype_0",
                    "gamma_smc_haplotype_1": "haplotype_1",
                }
            )
            .reset_index(drop=True)
        )
    return result


def write_gamma_pair_files(
    pair_table: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "",
    require_nonempty: bool = True,
) -> dict[str, Path]:
    """Atomically write explicit Gamma-SMC pair files for all genotype classes.

    Each non-comment row contains exactly two zero-based haplotype indices and
    can be supplied directly to Gamma-SMC's ``--pairs_file`` option.
    """
    tables = gamma_pair_tables(pair_table, require_nonempty=require_nonempty)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for genotype_class, table in tables.items():
        path = destination / f"{prefix}{genotype_class}.pairs.tsv"
        temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
        try:
            with temporary.open("w", encoding="utf-8", newline="\n") as handle:
                handle.write("# Gamma-SMC explicit haplotype pairs\n")
                handle.write(f"# genotype_class\t{genotype_class}\n")
                handle.write(f"# n_pairs\t{len(table)}\n")
                table.to_csv(handle, sep="\t", header=False, index=False)
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        paths[genotype_class] = path
    return paths


def _validated_thresholds(thresholds_years: Sequence[float]) -> np.ndarray:
    thresholds = np.asarray(thresholds_years, dtype=float)
    if thresholds.ndim != 1 or not len(thresholds):
        raise ValueError("thresholds_years must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(thresholds)) or np.any(thresholds <= 0):
        raise ValueError("TMRCA thresholds must be finite and positive")
    if len(np.unique(thresholds)) != len(thresholds):
        raise ValueError("TMRCA thresholds must be unique")
    return np.sort(thresholds)


def _validated_positions(ts, positions: Sequence[float]) -> np.ndarray:
    values = np.asarray(positions, dtype=float)
    if values.ndim != 1 or not len(values):
        raise ValueError("positions must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("positions must be finite")
    if np.any(values < 0) or np.any(values >= float(ts.sequence_length)):
        raise ValueError("profile position lies outside the tree sequence")
    return np.unique(values)


def tree_truth_profiles(
    ts,
    positions: Sequence[float],
    pair_table: pd.DataFrame,
    *,
    thresholds_years: Sequence[float] = DEFAULT_THRESHOLDS_YEARS,
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
    require_nonempty: bool = True,
) -> pd.DataFrame:
    """Compute long-form within-diploid tree-truth TMRCA profiles.

    ``P(TMRCA < x)`` uses a strict inequality after converting years to
    generations with the declared generation time.  Mean and median TMRCA are
    repeated on threshold rows to match the decoded long-form schema.
    """
    _validate_pair_table(pair_table)
    generation_time = float(generation_time_years)
    if not np.isfinite(generation_time) or generation_time <= 0:
        raise ValueError("generation_time_years must be finite and positive")
    thresholds = _validated_thresholds(thresholds_years)
    profile_positions = _validated_positions(ts, positions)

    node_pairs: dict[str, np.ndarray] = {}
    for genotype_class in GENOTYPE_CLASSES:
        selected = (
            pair_table
            if genotype_class == "overall"
            else pair_table[pair_table["genotype_class"] == genotype_class]
        )
        if require_nonempty and selected.empty:
            raise ValueError(
                f"no {genotype_class} diploid pairs; reject this simulation draw"
            )
        node_pairs[genotype_class] = selected[
            ["sample_node_0", "sample_node_1"]
        ].to_numpy(dtype=np.int64)

    rows = []
    for position in profile_positions:
        tree = ts.at(float(position))
        for genotype_class in GENOTYPE_CLASSES:
            pairs = node_pairs[genotype_class]
            if not len(pairs):
                continue
            times = np.fromiter(
                (tree.tmrca(int(left), int(right)) for left, right in pairs),
                dtype=float,
                count=len(pairs),
            )
            if not np.all(np.isfinite(times)):
                raise ValueError(f"nonfinite pairwise TMRCA at position {position:g}")
            mean_generations = float(times.mean())
            median_generations = float(np.median(times))
            for threshold_years in thresholds:
                threshold_generations = threshold_years / generation_time
                probability = float(np.mean(times < threshold_generations))
                rows.append(
                    {
                        "source": "tree_truth",
                        "position_0based": float(position),
                        "position_1based": float(position) + 1,
                        "genotype_class": genotype_class,
                        "n_pairs": int(len(times)),
                        "generation_time_years": generation_time,
                        "threshold_years": float(threshold_years),
                        "threshold_generations": float(threshold_generations),
                        "p_tmrca_lt_threshold": probability,
                        "mean_p_tmrca_lt_threshold": probability,
                        "frac_recent": probability,
                        "mean_tmrca_generations": mean_generations,
                        "median_tmrca_generations": median_generations,
                        "mean_tmrca_years": mean_generations * generation_time,
                        "median_tmrca_years": median_generations * generation_time,
                    }
                )
    return pd.DataFrame(rows)


_MEAN_PROBABILITY_COLUMN = re.compile(r"^mean_p_lt_(.+)$")
_RECENT_FRACTION_COLUMN = re.compile(r"^frac_recent_(.+)$")


def _threshold_columns(
    columns: Sequence[str],
    pattern: re.Pattern[str],
) -> dict[float, str]:
    result: dict[float, str] = {}
    for column in columns:
        match = pattern.match(str(column))
        if match is None:
            continue
        try:
            threshold = float(match.group(1))
        except ValueError as error:
            raise ValueError(
                f"cannot parse threshold from Gamma-SMC column {column}"
            ) from error
        if not np.isfinite(threshold) or threshold <= 0:
            raise ValueError(f"invalid threshold in Gamma-SMC column {column}")
        if threshold in result:
            raise ValueError(f"duplicate Gamma-SMC threshold {threshold:g} years")
        result[threshold] = str(column)
    return result


def gamma_summary_to_long(
    summary: pd.DataFrame,
    *,
    genotype_class: str | None = None,
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
    source: str = "gamma_smc",
    id_columns: Sequence[str] = (),
) -> pd.DataFrame:
    """Convert a multi-threshold Gamma-SMC summary to the truth-long schema.

    ``mean_p_lt_X`` is retained as the soft posterior CDF estimate, while
    ``frac_recent_X`` remains a separate posterior-call fraction.  These two
    quantities are never silently substituted for one another.
    """
    if summary.empty:
        raise ValueError("Gamma-SMC summary is empty")
    required = {"position_0based", "mean_tmrca_generations"}
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError(f"Gamma-SMC summary is missing columns: {', '.join(missing)}")
    if genotype_class is None:
        if "genotype_class" not in summary:
            raise ValueError(
                "genotype_class must be provided or present in the summary"
            )
    elif genotype_class not in GENOTYPE_CLASSES:
        raise ValueError(f"unknown genotype class: {genotype_class}")
    if genotype_class is not None and "genotype_class" in summary:
        observed = set(summary["genotype_class"].astype(str))
        if observed != {genotype_class}:
            raise ValueError(
                "summary genotype_class conflicts with the requested class"
            )

    generation_time = float(generation_time_years)
    if not np.isfinite(generation_time) or generation_time <= 0:
        raise ValueError("generation_time_years must be finite and positive")
    mean_columns = _threshold_columns(summary.columns, _MEAN_PROBABILITY_COLUMN)
    fraction_columns = _threshold_columns(summary.columns, _RECENT_FRACTION_COLUMN)
    thresholds = sorted(set(mean_columns).union(fraction_columns))
    if not thresholds:
        raise ValueError(
            "Gamma-SMC summary has no mean_p_lt_X or frac_recent_X threshold columns"
        )
    absent_ids = sorted(set(id_columns).difference(summary.columns))
    if absent_ids:
        raise ValueError(
            f"Gamma-SMC summary is missing ID columns: {', '.join(absent_ids)}"
        )

    rows = []
    for _, input_row in summary.iterrows():
        row_class = (
            genotype_class
            if genotype_class is not None
            else str(input_row["genotype_class"])
        )
        if row_class not in GENOTYPE_CLASSES:
            raise ValueError(f"unknown genotype class: {row_class}")
        position = float(input_row["position_0based"])
        mean_generations = float(input_row["mean_tmrca_generations"])
        median_generations = (
            float(input_row["median_tmrca_generations"])
            if "median_tmrca_generations" in summary
            else np.nan
        )
        for threshold in thresholds:
            probability = (
                float(input_row[mean_columns[threshold]])
                if threshold in mean_columns
                else np.nan
            )
            recent_fraction = (
                float(input_row[fraction_columns[threshold]])
                if threshold in fraction_columns
                else np.nan
            )
            output_row = {column: input_row[column] for column in id_columns}
            output_row.update(
                {
                    "source": source,
                    "position_0based": position,
                    "position_1based": (
                        float(input_row["position_1based"])
                        if "position_1based" in summary
                        else position + 1
                    ),
                    "genotype_class": row_class,
                    "n_pairs": (
                        int(input_row["n_pairs"])
                        if "n_pairs" in summary and pd.notna(input_row["n_pairs"])
                        else pd.NA
                    ),
                    "generation_time_years": generation_time,
                    "threshold_years": float(threshold),
                    "threshold_generations": float(threshold / generation_time),
                    "p_tmrca_lt_threshold": probability,
                    "mean_p_tmrca_lt_threshold": probability,
                    "frac_recent": recent_fraction,
                    "mean_tmrca_generations": mean_generations,
                    "median_tmrca_generations": median_generations,
                    "mean_tmrca_years": mean_generations * generation_time,
                    "median_tmrca_years": median_generations * generation_time,
                }
            )
            rows.append(output_row)
    return pd.DataFrame(rows)


def gamma_summaries_to_long(
    summaries: Mapping[str, pd.DataFrame],
    *,
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
    source: str = "gamma_smc",
    id_columns: Sequence[str] = (),
) -> pd.DataFrame:
    """Convert one Gamma-SMC summary per genotype class and concatenate them."""
    unknown = sorted(set(summaries).difference(GENOTYPE_CLASSES))
    if unknown:
        raise ValueError(f"unknown genotype classes: {', '.join(unknown)}")
    if not summaries:
        raise ValueError("no Gamma-SMC summaries were provided")
    return pd.concat(
        [
            gamma_summary_to_long(
                frame,
                genotype_class=genotype_class,
                generation_time_years=generation_time_years,
                source=source,
                id_columns=id_columns,
            )
            for genotype_class, frame in summaries.items()
        ],
        ignore_index=True,
    )


def normalized_cdf_auc(
    profile: pd.DataFrame,
    *,
    value_column: str = "p_tmrca_lt_threshold",
    group_columns: Sequence[str] = (
        "source",
        "position_0based",
        "genotype_class",
    ),
) -> pd.DataFrame:
    """Integrate each TMRCA CDF from ``(0, 0)`` and divide by max time.

    This is a normalized area under a *within-genotype TMRCA CDF*.  It is not a
    classifier/ROC AUC.  Larger values indicate more probability mass at recent
    coalescence times over the supplied threshold range.
    """
    required = set(group_columns).union({"threshold_years", value_column})
    missing = sorted(required.difference(profile.columns))
    if missing:
        raise ValueError(f"profile is missing columns: {', '.join(missing)}")
    if not group_columns:
        raise ValueError("at least one AUC grouping column is required")

    rows = []
    grouper = list(group_columns)
    for keys, group in profile.groupby(grouper, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        ordered = group.sort_values("threshold_years")
        x = ordered["threshold_years"].to_numpy(dtype=float)
        y = ordered[value_column].to_numpy(dtype=float)
        if len(np.unique(x)) != len(x):
            raise ValueError("CDF groups must have one row per threshold")
        if not np.all(np.isfinite(x)) or np.any(x <= 0):
            raise ValueError("CDF thresholds must be finite and positive")
        if not np.all(np.isfinite(y)):
            raise ValueError(f"CDF column {value_column} contains missing values")
        if np.any(y < -1e-9) or np.any(y > 1 + 1e-9):
            raise ValueError("CDF probabilities must lie between zero and one")
        if np.any(np.diff(y) < -1e-8):
            raise ValueError("CDF probabilities must be nondecreasing with threshold")
        y = np.clip(np.maximum.accumulate(y), 0.0, 1.0)
        x_with_origin = np.r_[0.0, x]
        y_with_origin = np.r_[0.0, y]
        widths = np.diff(x_with_origin)
        area = float(np.sum(widths * (y_with_origin[:-1] + y_with_origin[1:]) / 2))
        maximum = float(x[-1])
        output = dict(zip(group_columns, keys, strict=True))
        output.update(
            {
                "cdf_auc_kind": "within_genotype_tmrca_cdf",
                "cdf_auc_includes_origin_0_0": True,
                "cdf_auc_max_threshold_years": maximum,
                "cdf_area_probability_years": area,
                "normalized_cdf_auc": area / maximum,
            }
        )
        rows.append(output)
    return pd.DataFrame(rows)


def _one_tmrca_value_per_class(
    group: pd.DataFrame,
    value_column: str,
) -> dict[str, float]:
    values: dict[str, float] = {}
    for genotype_class, class_group in group.groupby("genotype_class", sort=False):
        unique = class_group[value_column].dropna().to_numpy(dtype=float)
        if not len(unique):
            continue
        if not np.allclose(unique, unique[0], rtol=1e-12, atol=1e-12):
            raise ValueError(
                f"{value_column} changes across threshold rows for {genotype_class}"
            )
        values[str(genotype_class)] = float(unique[0])
    return values


def internal_genotype_contrasts(
    profile: pd.DataFrame,
    *,
    value_column: str = "p_tmrca_lt_threshold",
    group_columns: Sequence[str] = ("source", "position_0based"),
) -> pd.DataFrame:
    """Compute alt-minus-ref internal CDF-area and mean-TMRCA contrasts.

    Under a successful recent sweep, the anticipated signs are positive for
    normalized CDF area and negative for mean TMRCA.  These descriptive
    contrasts have no neutral null and are neither ROC statistics nor p-values.
    """
    required = set(group_columns).union(
        {
            "genotype_class",
            "mean_tmrca_generations",
            "generation_time_years",
        }
    )
    missing = sorted(required.difference(profile.columns))
    if missing:
        raise ValueError(f"profile is missing columns: {', '.join(missing)}")
    auc = normalized_cdf_auc(
        profile,
        value_column=value_column,
        group_columns=tuple(group_columns) + ("genotype_class",),
    )
    auc_lookup = {
        tuple(row[column] for column in group_columns)
        + (row["genotype_class"],): float(row["normalized_cdf_auc"])
        for _, row in auc.iterrows()
    }

    rows = []
    grouper = list(group_columns)
    for keys, group in profile.groupby(grouper, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        tmrca = _one_tmrca_value_per_class(group, "mean_tmrca_generations")
        missing_classes = {"hom_ref", "hom_alt"}.difference(tmrca)
        if missing_classes:
            raise ValueError(
                "internal genotype contrast requires hom_ref and hom_alt; missing "
                + ", ".join(sorted(missing_classes))
            )
        ref_auc = auc_lookup.get(tuple(keys) + ("hom_ref",))
        alt_auc = auc_lookup.get(tuple(keys) + ("hom_alt",))
        if ref_auc is None or alt_auc is None:
            raise ValueError("CDF AUC is missing hom_ref or hom_alt")
        generation_times = group["generation_time_years"].dropna().unique()
        if len(generation_times) != 1:
            raise ValueError("each contrast must use one declared generation time")
        generation_time = float(generation_times[0])
        mean_difference = tmrca["hom_alt"] - tmrca["hom_ref"]
        output = dict(zip(group_columns, keys, strict=True))
        output.update(
            {
                "contrast_scope": "internal_genotype_contrast_no_neutral_null",
                "alt_minus_ref_normalized_cdf_auc": alt_auc - ref_auc,
                "alt_minus_ref_mean_tmrca_generations": mean_difference,
                "alt_minus_ref_mean_tmrca_years": mean_difference * generation_time,
                "expected_sweep_sign_normalized_cdf_auc": "positive",
                "expected_sweep_sign_mean_tmrca": "negative",
            }
        )
        rows.append(output)
    return pd.DataFrame(rows)


def _output_stem(path: str | Path) -> Path:
    stem = Path(path)
    if stem.suffix.lower() in {".png", ".pdf"}:
        stem = stem.with_suffix("")
    stem.parent.mkdir(parents=True, exist_ok=True)
    return stem


def _save_png_pdf(fig, output_stem: str | Path) -> dict[str, Path]:
    stem = _output_stem(output_stem)
    paths = {"png": stem.with_suffix(".png"), "pdf": stem.with_suffix(".pdf")}
    fig.savefig(paths["png"], dpi=220, bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)
    return paths


def _ordered_classes(values: Sequence[str]) -> list[str]:
    observed = set(map(str, values))
    return [name for name in GENOTYPE_CLASSES if name in observed]


def plot_center_cdf(
    profile: pd.DataFrame,
    focal_coordinate: float,
    output_stem: str | Path,
    *,
    value_column: str = "p_tmrca_lt_threshold",
    title: str | None = None,
) -> dict[str, Path]:
    """Plot genotype-class TMRCA CDFs at the focal coordinate."""
    required = {
        "source",
        "position_0based",
        "genotype_class",
        "threshold_years",
        value_column,
    }
    missing = sorted(required.difference(profile.columns))
    if missing:
        raise ValueError(f"profile is missing columns: {', '.join(missing)}")
    tolerance = max(1e-9, abs(float(focal_coordinate)) * 1e-12)
    center = profile[
        np.isclose(
            profile["position_0based"].to_numpy(dtype=float),
            float(focal_coordinate),
            rtol=0.0,
            atol=tolerance,
        )
    ]
    if center.empty:
        raise ValueError(f"profile has no row at focal coordinate {focal_coordinate:g}")

    fig, axis = plt.subplots(figsize=(11, 8.5), constrained_layout=True)
    sources = list(dict.fromkeys(center["source"].astype(str)))
    classes = _ordered_classes(center["genotype_class"].astype(str))
    for source_index, source in enumerate(sources):
        for genotype_class in classes:
            group = center[
                (center["source"].astype(str) == source)
                & (center["genotype_class"].astype(str) == genotype_class)
            ].sort_values("threshold_years")
            if group.empty:
                continue
            x = np.r_[0.0, group["threshold_years"].to_numpy(dtype=float) / 1_000]
            y = np.r_[0.0, group[value_column].to_numpy(dtype=float)]
            axis.plot(
                x,
                y,
                color=_CLASS_COLORS[genotype_class],
                linestyle=_SOURCE_LINESTYLES[source_index % len(_SOURCE_LINESTYLES)],
                linewidth=2.8,
                marker="o",
                markersize=5,
                label=f"{source}: {_CLASS_LABELS[genotype_class]}",
            )
    axis.set_xlim(left=0)
    axis.set_ylim(0, 1.02)
    axis.set_xlabel("TMRCA threshold (kya)", fontsize=_PLOT_FONT["axis"])
    axis.set_ylabel("P(TMRCA < threshold)", fontsize=_PLOT_FONT["axis"])
    axis.set_title(
        title or f"Focal-site TMRCA CDF at {float(focal_coordinate) / 1e6:g} Mb",
        fontsize=_PLOT_FONT["title"],
    )
    axis.tick_params(axis="both", labelsize=_PLOT_FONT["tick"])
    axis.grid(alpha=0.2)
    axis.legend(fontsize=_PLOT_FONT["legend"], ncol=2, frameon=False)
    return _save_png_pdf(fig, output_stem)


def _centers_to_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if len(values) == 1:
        half_width = max(abs(float(values[0])) * 0.05, 0.5)
        return np.asarray([values[0] - half_width, values[0] + half_width])
    midpoints = (values[:-1] + values[1:]) / 2
    return np.r_[
        values[0] - (midpoints[0] - values[0]),
        midpoints,
        values[-1] + (values[-1] - midpoints[-1]),
    ]


def plot_spatial_probability_heatmaps(
    profile: pd.DataFrame,
    output_stem: str | Path,
    *,
    source: str | None = None,
    value_column: str = "p_tmrca_lt_threshold",
    title: str | None = None,
) -> dict[str, Path]:
    """Plot spatial position-by-threshold heatmaps for all genotype classes."""
    required = {
        "source",
        "position_0based",
        "genotype_class",
        "threshold_years",
        value_column,
    }
    missing = sorted(required.difference(profile.columns))
    if missing:
        raise ValueError(f"profile is missing columns: {', '.join(missing)}")
    sources = list(dict.fromkeys(profile["source"].astype(str)))
    if source is None:
        if len(sources) != 1:
            raise ValueError("select one source for a spatial heatmap figure")
        source = sources[0]
    data = profile[profile["source"].astype(str) == str(source)]
    if data.empty:
        raise ValueError(f"profile has no rows for source {source}")
    classes = _ordered_classes(data["genotype_class"].astype(str))
    if not classes:
        raise ValueError("profile has no recognized genotype classes")

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(11, 8.5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    image = None
    for axis, genotype_class in zip(axes.flat, classes, strict=False):
        group = data[data["genotype_class"].astype(str) == genotype_class]
        matrix = (
            group.pivot(
                index="threshold_years",
                columns="position_0based",
                values=value_column,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        x = matrix.columns.to_numpy(dtype=float) / 1e6
        y = matrix.index.to_numpy(dtype=float) / 1_000
        image = axis.pcolormesh(
            _centers_to_edges(x),
            _centers_to_edges(y),
            matrix.to_numpy(dtype=float),
            shading="flat",
            cmap="viridis",
            vmin=0,
            vmax=1,
        )
        axis.set_title(_CLASS_LABELS[genotype_class], fontsize=_PLOT_FONT["title"])
        axis.tick_params(axis="both", labelsize=_PLOT_FONT["tick"])
    for axis in axes.flat[len(classes) :]:
        axis.set_visible(False)
    for axis in axes[-1, :]:
        if axis.get_visible():
            axis.set_xlabel("Position (Mb)", fontsize=_PLOT_FONT["axis"])
    for axis in axes[:, 0]:
        if axis.get_visible():
            axis.set_ylabel("TMRCA threshold (kya)", fontsize=_PLOT_FONT["axis"])
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes, shrink=0.86)
        colorbar.set_label("P(TMRCA < threshold)", fontsize=_PLOT_FONT["axis"])
        colorbar.ax.tick_params(labelsize=_PLOT_FONT["tick"])
    fig.suptitle(
        title or f"Spatial genotype-class TMRCA CDF: {source}",
        fontsize=_PLOT_FONT["title"],
    )
    return _save_png_pdf(fig, output_stem)


def plot_class_tmrca_profiles(
    profile: pd.DataFrame,
    output_stem: str | Path,
    *,
    focal_coordinate: float | None = None,
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
    title: str | None = None,
) -> dict[str, Path]:
    """Plot spatial mean (and available median) TMRCA by genotype class."""
    required = {
        "source",
        "position_0based",
        "genotype_class",
        "mean_tmrca_generations",
    }
    missing = sorted(required.difference(profile.columns))
    if missing:
        raise ValueError(f"profile is missing columns: {', '.join(missing)}")
    generation_time = float(generation_time_years)
    if not np.isfinite(generation_time) or generation_time <= 0:
        raise ValueError("generation_time_years must be finite and positive")
    sources = list(dict.fromkeys(profile["source"].astype(str)))
    fig, axes = plt.subplots(
        1,
        len(sources),
        figsize=(11, 8.5),
        squeeze=False,
        sharey=True,
        constrained_layout=True,
    )
    legend_entries: dict[str, object] = {}
    for axis, source in zip(axes.flat, sources, strict=True):
        source_data = profile[profile["source"].astype(str) == source]
        classes = _ordered_classes(source_data["genotype_class"].astype(str))
        for genotype_class in classes:
            group = source_data[
                source_data["genotype_class"].astype(str) == genotype_class
            ].sort_values("position_0based")
            columns = ["position_0based", "mean_tmrca_generations"]
            if "median_tmrca_generations" in group:
                columns.append("median_tmrca_generations")
            unique = group[columns].drop_duplicates()
            if unique["position_0based"].duplicated().any():
                raise ValueError("TMRCA summaries change across threshold rows")
            x = unique["position_0based"].to_numpy(dtype=float) / 1e6
            mean_kya = (
                unique["mean_tmrca_generations"].to_numpy(dtype=float)
                * generation_time
                / 1_000
            )
            axis.plot(
                x,
                mean_kya,
                color=_CLASS_COLORS[genotype_class],
                linewidth=2.5,
                label=f"{_CLASS_LABELS[genotype_class]} mean",
            )
            if "median_tmrca_generations" in unique:
                median = unique["median_tmrca_generations"].to_numpy(dtype=float)
                if np.isfinite(median).any():
                    axis.plot(
                        x,
                        median * generation_time / 1_000,
                        color=_CLASS_COLORS[genotype_class],
                        linewidth=1.8,
                        linestyle=":",
                        label=f"{_CLASS_LABELS[genotype_class]} median",
                    )
        if focal_coordinate is not None:
            axis.axvline(
                float(focal_coordinate) / 1e6,
                color="#cc79a7",
                linewidth=1.8,
                linestyle="--",
                label="Focal site",
            )
        axis.set_yscale("log")
        axis.set_xlabel("Position (Mb)", fontsize=_PLOT_FONT["axis"])
        axis.set_title(source, fontsize=_PLOT_FONT["title"])
        axis.tick_params(axis="both", labelsize=_PLOT_FONT["tick"])
        axis.grid(alpha=0.2)
        handles, labels = axis.get_legend_handles_labels()
        for handle, label in zip(handles, labels, strict=True):
            legend_entries.setdefault(label, handle)
    axes[0, 0].set_ylabel("Within-diploid TMRCA (kya)", fontsize=_PLOT_FONT["axis"])
    fig.suptitle(title or "TMRCA profiles by focal genotype class", fontsize=22)
    fig.legend(
        list(legend_entries.values()),
        list(legend_entries),
        loc="outside lower center",
        ncol=3,
        fontsize=11,
        frameon=False,
    )
    return _save_png_pdf(fig, output_stem)


def _safe_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "source"


def write_publication_plots(
    profile: pd.DataFrame,
    focal_coordinate: float,
    output_dir: str | Path,
    *,
    prefix: str = "eas_sweep",
    value_column: str = "p_tmrca_lt_threshold",
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
) -> dict[str, dict[str, Path]]:
    """Write the core publication figures, each in PNG and PDF format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "center_cdf": plot_center_cdf(
            profile,
            focal_coordinate,
            output_dir / f"{prefix}_center_cdf",
            value_column=value_column,
        ),
        "tmrca_profiles": plot_class_tmrca_profiles(
            profile,
            output_dir / f"{prefix}_class_tmrca_profiles",
            focal_coordinate=focal_coordinate,
            generation_time_years=generation_time_years,
        ),
    }
    for source in dict.fromkeys(profile["source"].astype(str)):
        key = f"spatial_heatmaps_{_safe_filename(source)}"
        outputs[key] = plot_spatial_probability_heatmaps(
            profile,
            output_dir / f"{prefix}_{key}",
            source=source,
            value_column=value_column,
        )
    return outputs


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns).difference(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def _resolve_scenario_column(frame: pd.DataFrame, requested: str) -> str:
    if requested in frame.columns:
        return requested
    if requested == "scenario" and "study_type" in frame.columns:
        return "study_type"
    raise ValueError(
        f"study table has no {requested!r} column"
        + (" (or 'study_type' fallback)" if requested == "scenario" else "")
    )


def _selection_label(value) -> str:
    try:
        return f"s={float(value):g}"
    except (TypeError, ValueError):
        return f"s={value}"


def _curve_order(value) -> tuple[int, str]:
    text = str(value).lower()
    if "lower" in text or text in {"lo", "low"}:
        return 0, text
    if "median" in text or text in {"med", "mid"}:
        return 1, text
    if "upper" in text or text in {"hi", "high"}:
        return 2, text
    return 3, text


def _specification_layout(
    frame: pd.DataFrame,
    *,
    specification_column: str,
    selection_column: str,
    demography_curve_column: str,
    target_af_column: str,
    label_column: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """Order and label one-row-per-specification study metadata."""
    required = [
        specification_column,
        selection_column,
        demography_curve_column,
        target_af_column,
    ]
    if label_column is not None:
        required.append(label_column)
    _require_columns(frame, required, "study specification table")
    if frame[specification_column].isna().any():
        raise ValueError("specification IDs cannot be missing")
    if frame[selection_column].isna().any():
        raise ValueError("selection coefficients cannot be missing")
    if label_column is not None and frame[label_column].isna().any():
        raise ValueError("plot labels cannot be missing")
    if frame[specification_column].duplicated().any():
        duplicates = frame.loc[
            frame[specification_column].duplicated(keep=False), specification_column
        ].astype(str)
        raise ValueError(
            "study table must have one row per specification; duplicate IDs: "
            + ", ".join(sorted(set(duplicates)))
        )

    has_curve = frame[demography_curve_column].notna()
    has_target = frame[target_af_column].notna()
    ordered = frame.copy()
    if has_curve.any():
        if not has_curve.all():
            raise ValueError(
                "one scenario cannot mix rows with and without demography curves"
            )
        ordered["_dimension_order"] = ordered[demography_curve_column].map(_curve_order)
        dimension = "demography curve"
    elif has_target.any():
        if not has_target.all():
            raise ValueError(
                "target_allele_frequency must be present for every introgression specification"
            )
        target = pd.to_numeric(ordered[target_af_column], errors="coerce")
        if target.isna().any() or not target.between(0, 1, inclusive="both").all():
            raise ValueError("target allele frequencies must lie between zero and one")
        ordered["_dimension_order"] = target
        dimension = "target allele frequency"
    else:
        ordered["_dimension_order"] = 0
        dimension = "specification"

    selection_numeric = pd.to_numeric(ordered[selection_column], errors="coerce")
    ordered["_selection_order"] = np.where(
        selection_numeric.notna(), selection_numeric, np.inf
    )
    ordered["_specification_order"] = ordered[specification_column].astype(str)
    ordered = ordered.sort_values(
        ["_dimension_order", "_selection_order", "_specification_order"],
        kind="stable",
    ).reset_index(drop=True)

    if label_column is not None:
        ordered["_plot_label"] = ordered[label_column].astype(str)
    else:
        labels = []
        for row in ordered.to_dict(orient="records"):
            if dimension == "demography curve":
                dimension_label = f"curve={row[demography_curve_column]}"
            elif dimension == "target allele frequency":
                dimension_label = f"target AF={float(row[target_af_column]):.0%}"
            else:
                dimension_label = "specification"
            labels.append(
                f"{dimension_label}; {_selection_label(row[selection_column])}\n"
                f"{row[specification_column]}"
            )
        ordered["_plot_label"] = labels
    return ordered, dimension


def plot_acceptance_effort(
    acceptance_rows: pd.DataFrame,
    output_stem: str | Path,
    *,
    scenario_label: str | None = None,
    specification_column: str = "specification_id",
    selection_column: str = "selection_coefficient",
    demography_curve_column: str = "demography_curve",
    target_af_column: str = "target_allele_frequency",
    achieved_af_column: str = "achieved_allele_frequency",
    attempts_column: str = "attempts_to_accept",
    label_column: str | None = None,
) -> dict[str, Path]:
    """Plot achieved allele frequency and seed-ordered acceptance effort.

    ``acceptance_rows`` must contain exactly one accepted row per specification.
    Attempts are one-based counts through the accepted draw, not zero-based
    attempt IDs.
    """
    required = [achieved_af_column, attempts_column]
    _require_columns(acceptance_rows, required, "acceptance table")
    if acceptance_rows.empty:
        raise ValueError("acceptance table is empty")
    ordered, dimension = _specification_layout(
        acceptance_rows,
        specification_column=specification_column,
        selection_column=selection_column,
        demography_curve_column=demography_curve_column,
        target_af_column=target_af_column,
        label_column=label_column,
    )
    if label_column is not None:
        compact_labels = ordered["_plot_label"].astype(str).tolist()
        tick_rotation = 35
        tick_alignment = "right"
    else:
        compact_labels = []
        for row in ordered.to_dict(orient="records"):
            if dimension == "demography curve":
                compact_labels.append(str(row[demography_curve_column]))
            elif dimension == "target allele frequency":
                compact_labels.append(f"{float(row[target_af_column]):.0%}")
            else:
                compact_labels.append(str(row[specification_column]))
        tick_rotation = 0 if dimension != "specification" else 35
        tick_alignment = "center" if tick_rotation == 0 else "right"
    achieved = pd.to_numeric(ordered[achieved_af_column], errors="coerce")
    if achieved.isna().any() or not achieved.between(0, 1, inclusive="both").all():
        raise ValueError("achieved allele frequencies must lie between zero and one")
    attempts = pd.to_numeric(ordered[attempts_column], errors="coerce")
    if (
        attempts.isna().any()
        or (attempts < 1).any()
        or not np.equal(attempts, np.floor(attempts)).all()
    ):
        raise ValueError("attempts_to_accept must contain positive integer counts")

    selection_values = list(dict.fromkeys(ordered[selection_column]))
    colors = {
        value: plt.get_cmap("tab10")(index % 10)
        for index, value in enumerate(selection_values)
    }
    x = np.arange(len(ordered), dtype=float)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11, 8.5),
        sharex=True,
        constrained_layout=True,
    )
    seen_selection = set()
    for index, row in ordered.iterrows():
        selection = row[selection_column]
        label = (
            _selection_label(selection)
            if selection not in seen_selection
            else "_nolegend_"
        )
        seen_selection.add(selection)
        axes[0].scatter(
            x[index],
            achieved.iloc[index],
            s=95,
            color=colors[selection],
            edgecolor="black",
            linewidth=0.6,
            label=label,
            zorder=3,
        )
        axes[1].vlines(
            x[index],
            1,
            attempts.iloc[index],
            color=colors[selection],
            linewidth=2.6,
            alpha=0.85,
        )
        axes[1].scatter(
            x[index],
            attempts.iloc[index],
            s=80,
            color=colors[selection],
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
        )
    target = pd.to_numeric(ordered[target_af_column], errors="coerce")
    if target.notna().any():
        axes[0].scatter(
            x[target.notna()],
            target[target.notna()],
            marker="D",
            s=70,
            facecolor="white",
            edgecolor="black",
            linewidth=1.3,
            label="Target AF",
            zorder=4,
        )

    axes[0].set_ylim(-0.02, 1.02)
    axes[0].yaxis.set_major_formatter(PercentFormatter(1.0))
    axes[0].set_ylabel("Achieved allele frequency", fontsize=_PLOT_FONT["axis"])
    axes[0].set_title(
        "Accepted draw allele frequency",
        fontsize=_PLOT_FONT["title"],
    )
    axes[0].legend(fontsize=_PLOT_FONT["legend"], frameon=False, ncol=2)
    axes[1].set_yscale("log")
    axes[1].set_ylim(bottom=0.9)
    axes[1].set_ylabel(
        "Attempts through accepted draw\n(log scale)",
        fontsize=_PLOT_FONT["axis"],
    )
    axes[1].set_title("Seed-ordered acceptance effort", fontsize=_PLOT_FONT["title"])
    axes[1].set_xticks(
        x,
        compact_labels,
        rotation=tick_rotation,
        ha=tick_alignment,
    )
    axes[1].set_xlabel(
        f"Accepted specification ({dimension}; color denotes s)",
        fontsize=_PLOT_FONT["axis"],
    )
    for axis in axes:
        axis.tick_params(axis="both", labelsize=_PLOT_FONT["tick"])
        axis.grid(axis="y", alpha=0.2)
    axes[1].tick_params(axis="x", labelsize=12)
    fig.suptitle(
        f"{scenario_label + ': ' if scenario_label else ''}"
        "Achieved AF and simulation acceptance effort",
        fontsize=22,
    )
    return _save_png_pdf(fig, output_stem)


def _center_specification_metadata(
    center_summaries: pd.DataFrame,
    *,
    specification_column: str,
    selection_column: str,
    demography_curve_column: str,
    target_af_column: str,
    label_column: str | None,
) -> pd.DataFrame:
    columns = [
        specification_column,
        selection_column,
        demography_curve_column,
        target_af_column,
    ]
    if label_column is not None:
        columns.append(label_column)
    metadata = center_summaries[columns].drop_duplicates()
    counts = metadata.groupby(specification_column, dropna=False).size()
    inconsistent = counts[counts != 1]
    if len(inconsistent):
        raise ValueError(
            "center summary metadata changes across sources for specification IDs: "
            + ", ".join(map(str, inconsistent.index))
        )
    return metadata


def _source_order(values: Sequence[str]) -> list[str]:
    observed = list(dict.fromkeys(map(str, values)))
    preferred = [source for source in ("tree_truth", "gamma_smc") if source in observed]
    return preferred + [source for source in observed if source not in preferred]


def plot_center_internal_contrast_heatmaps(
    center_summaries: pd.DataFrame,
    output_stem: str | Path,
    *,
    scenario_label: str | None = None,
    specification_column: str = "specification_id",
    selection_column: str = "selection_coefficient",
    demography_curve_column: str = "demography_curve",
    target_af_column: str = "target_allele_frequency",
    source_column: str = "source",
    cdf_contrast_column: str = "alt_minus_ref_normalized_cdf_auc",
    mean_tmrca_contrast_column: str = "alt_minus_ref_mean_tmrca_generations",
    label_column: str | None = None,
) -> dict[str, Path]:
    """Plot center alt-minus-ref contrasts across simulation specifications.

    Cells are descriptive, within-simulation genotype contrasts.  They have no
    neutral null, are not p-values, and the normalized CDF area is not ROC AUC.
    """
    required = [
        specification_column,
        selection_column,
        demography_curve_column,
        target_af_column,
        source_column,
        cdf_contrast_column,
        mean_tmrca_contrast_column,
    ]
    if label_column is not None:
        required.append(label_column)
    _require_columns(center_summaries, required, "center contrast table")
    if center_summaries.empty:
        raise ValueError("center contrast table is empty")
    if center_summaries[source_column].isna().any():
        raise ValueError("center contrast sources cannot be missing")
    if center_summaries[[specification_column, source_column]].duplicated().any():
        duplicate = center_summaries.loc[
            center_summaries[[specification_column, source_column]].duplicated(
                keep=False
            ),
            [specification_column, source_column],
        ]
        formatted = [f"{spec}/{source}" for spec, source in duplicate.to_numpy()]
        raise ValueError(
            "center contrast table must have one row per specification/source; "
            "duplicates: " + ", ".join(sorted(set(formatted)))
        )
    metadata = _center_specification_metadata(
        center_summaries,
        specification_column=specification_column,
        selection_column=selection_column,
        demography_curve_column=demography_curve_column,
        target_af_column=target_af_column,
        label_column=label_column,
    )
    ordered, dimension = _specification_layout(
        metadata,
        specification_column=specification_column,
        selection_column=selection_column,
        demography_curve_column=demography_curve_column,
        target_af_column=target_af_column,
        label_column=label_column,
    )
    sources = _source_order(center_summaries[source_column].astype(str))
    spec_ids = ordered[specification_column].tolist()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11, 8.5),
        sharey=True,
        constrained_layout=True,
    )
    definitions = (
        (
            cdf_contrast_column,
            "Normalized TMRCA-CDF area\nalt/alt - ref/ref (expected > 0)",
            ".3f",
        ),
        (
            mean_tmrca_contrast_column,
            "Mean TMRCA (generations)\nalt/alt - ref/ref (expected < 0)",
            ",.0f",
        ),
    )
    for axis, (metric, metric_title, number_format) in zip(
        axes, definitions, strict=True
    ):
        numeric = pd.to_numeric(center_summaries[metric], errors="coerce")
        if numeric.isna().any() or not np.isfinite(numeric).all():
            raise ValueError(f"center contrast metric {metric} must be finite")
        if metric == cdf_contrast_column and (numeric.abs() > 1 + 1e-9).any():
            raise ValueError("alt-minus-ref normalized CDF AUC must lie in [-1, 1]")
        plotting = center_summaries.assign(_metric=numeric)
        matrix = plotting.pivot(
            index=specification_column,
            columns=source_column,
            values="_metric",
        ).reindex(index=spec_ids, columns=sources)
        values = matrix.to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        magnitude = max(float(np.max(np.abs(finite))), np.finfo(float).eps)
        color_map = plt.get_cmap("RdBu_r").with_extremes(bad="#e6e6e6")
        image = axis.imshow(
            np.ma.masked_invalid(values),
            aspect="auto",
            cmap=color_map,
            norm=TwoSlopeNorm(vmin=-magnitude, vcenter=0.0, vmax=magnitude),
        )
        for row_index in range(values.shape[0]):
            for column_index in range(values.shape[1]):
                value = values[row_index, column_index]
                if not np.isfinite(value):
                    text = "pending"
                    color = "#595959"
                else:
                    text = format(value, f"+{number_format}")
                    color = "white" if abs(value) > 0.55 * magnitude else "black"
                axis.text(
                    column_index,
                    row_index,
                    text,
                    ha="center",
                    va="center",
                    fontsize=9 if len(spec_ids) <= 14 else 7,
                    color=color,
                )
        axis.set_xticks(range(len(sources)), sources, rotation=25, ha="right")
        axis.set_yticks(range(len(spec_ids)), ordered["_plot_label"])
        axis.set_title(metric_title, fontsize=16)
        axis.tick_params(axis="x", labelsize=12)
        axis.tick_params(axis="y", labelsize=10)
        colorbar = fig.colorbar(image, ax=axis, shrink=0.72)
        colorbar.ax.tick_params(labelsize=11)
    axes[0].set_ylabel(
        f"Simulation specification ({dimension})", fontsize=_PLOT_FONT["axis"]
    )
    fig.suptitle(
        f"{scenario_label + ': ' if scenario_label else ''}"
        "Center internal genotype contrasts\n"
        "Descriptive only: no neutral null or p-values; CDF area is not ROC AUC",
        fontsize=19,
    )
    return _save_png_pdf(fig, output_stem)


def write_cross_spec_study_plots(
    center_summaries: pd.DataFrame,
    acceptance_rows: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "eas_sweep",
    scenario_column: str = "scenario",
    specification_column: str = "specification_id",
    selection_column: str = "selection_coefficient",
    demography_curve_column: str = "demography_curve",
    target_af_column: str = "target_allele_frequency",
    achieved_af_column: str = "achieved_allele_frequency",
    attempts_column: str = "attempts_to_accept",
    source_column: str = "source",
    cdf_contrast_column: str = "alt_minus_ref_normalized_cdf_auc",
    mean_tmrca_contrast_column: str = "alt_minus_ref_mean_tmrca_generations",
    label_column: str | None = None,
    require_matching_specifications: bool = True,
) -> dict[str, dict[str, dict[str, Path]]]:
    """Write cross-spec acceptance and center-contrast figures per scenario.

    ``scenario`` is the canonical scenario column.  For compatibility with a
    compact table using ``study_type``, the default automatically falls back to
    that name.  Every figure is written as publication-sized PNG and PDF.
    """
    if center_summaries.empty or acceptance_rows.empty:
        raise ValueError("center and acceptance study tables must both be nonempty")
    center_scenario_column = _resolve_scenario_column(center_summaries, scenario_column)
    acceptance_scenario_column = _resolve_scenario_column(
        acceptance_rows, scenario_column
    )
    _require_columns(center_summaries, [specification_column], "center contrast table")
    _require_columns(acceptance_rows, [specification_column], "acceptance table")
    if center_summaries[center_scenario_column].isna().any():
        raise ValueError("center contrast scenarios cannot be missing")
    if acceptance_rows[acceptance_scenario_column].isna().any():
        raise ValueError("acceptance scenarios cannot be missing")
    center_scenarios = set(center_summaries[center_scenario_column].astype(str))
    acceptance_scenarios = set(acceptance_rows[acceptance_scenario_column].astype(str))
    if center_scenarios != acceptance_scenarios:
        raise ValueError(
            "center and acceptance tables cover different scenarios: "
            f"center={sorted(center_scenarios)}, "
            f"acceptance={sorted(acceptance_scenarios)}"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, dict[str, dict[str, Path]]] = {}
    safe_names = [_safe_filename(scenario) for scenario in center_scenarios]
    if len(set(safe_names)) != len(safe_names):
        raise ValueError("scenario names collide after conversion to output filenames")
    for scenario in sorted(center_scenarios):
        center = center_summaries[
            center_summaries[center_scenario_column].astype(str) == scenario
        ]
        acceptance = acceptance_rows[
            acceptance_rows[acceptance_scenario_column].astype(str) == scenario
        ]
        if require_matching_specifications:
            center_specs = set(center[specification_column].astype(str))
            acceptance_specs = set(acceptance[specification_column].astype(str))
            if center_specs != acceptance_specs:
                raise ValueError(
                    f"scenario {scenario} has different center and acceptance "
                    f"specifications: center={sorted(center_specs)}, "
                    f"acceptance={sorted(acceptance_specs)}"
                )
        safe_scenario = _safe_filename(scenario)
        outputs[scenario] = {
            "acceptance_effort": plot_acceptance_effort(
                acceptance,
                output_dir / f"{prefix}_{safe_scenario}_acceptance_effort",
                scenario_label=scenario,
                specification_column=specification_column,
                selection_column=selection_column,
                demography_curve_column=demography_curve_column,
                target_af_column=target_af_column,
                achieved_af_column=achieved_af_column,
                attempts_column=attempts_column,
                label_column=label_column,
            ),
            "center_internal_contrasts": plot_center_internal_contrast_heatmaps(
                center,
                output_dir / f"{prefix}_{safe_scenario}_center_internal_contrasts",
                scenario_label=scenario,
                specification_column=specification_column,
                selection_column=selection_column,
                demography_curve_column=demography_curve_column,
                target_af_column=target_af_column,
                source_column=source_column,
                cdf_contrast_column=cdf_contrast_column,
                mean_tmrca_contrast_column=mean_tmrca_contrast_column,
                label_column=label_column,
            ),
        }
    return outputs


def _step_value_at_age(
    times: np.ndarray,
    values: np.ndarray,
    age_generations: float,
) -> float:
    index = int(np.searchsorted(times, age_generations, side="right") - 1)
    index = int(np.clip(index, 0, len(values) - 1))
    return float(values[index])


def plot_eas_ne_origin_design(
    trajectories: pd.DataFrame,
    origins: pd.DataFrame,
    output_stem: str | Path,
    *,
    curve_column: str = "demography_curve",
    time_column: str = "time_generations",
    ne_column: str = "effective_population_size",
    selection_column: str = "selection_coefficient",
    origin_age_column: str = "origin_age_generations",
    origin_ne_column: str = "ne_at_origin",
    expected_curves: Sequence[str] = ("q025", "median", "q975"),
    expected_selection_coefficients: Sequence[float] = (0.01, 0.005, 0.001),
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
) -> dict[str, Path]:
    """Plot pointwise EAS Ne quantiles and all nine de novo origin ages.

    Both inputs are tidy tables. ``trajectories`` has one row per
    curve/time/Ne point; ``origins`` has one row per curve/selection
    combination. If ``ne_at_origin`` is absent, its value is recovered from the
    corresponding piecewise-constant trajectory.
    """
    _require_columns(
        trajectories,
        [curve_column, time_column, ne_column],
        "EAS trajectory table",
    )
    _require_columns(
        origins,
        [curve_column, selection_column, origin_age_column],
        "EAS origin table",
    )
    if trajectories.empty or origins.empty:
        raise ValueError("EAS trajectory and origin tables must both be nonempty")
    generation_time = float(generation_time_years)
    if not np.isfinite(generation_time) or generation_time <= 0:
        raise ValueError("generation_time_years must be finite and positive")

    expected_curve_names = tuple(map(str, expected_curves))
    if len(expected_curve_names) != 3 or len(set(expected_curve_names)) != 3:
        raise ValueError("expected_curves must contain three unique curve names")
    observed_curves = set(trajectories[curve_column].astype(str))
    if observed_curves != set(expected_curve_names):
        raise ValueError(
            "trajectory table must contain exactly the requested pointwise curves: "
            f"expected={list(expected_curve_names)}, observed={sorted(observed_curves)}"
        )
    origin_curves = set(origins[curve_column].astype(str))
    if origin_curves != set(expected_curve_names):
        raise ValueError(
            "origin table must cover every pointwise curve exactly: "
            f"expected={list(expected_curve_names)}, observed={sorted(origin_curves)}"
        )
    if trajectories[[curve_column, time_column]].duplicated().any():
        raise ValueError("trajectory table has duplicate curve/time rows")

    clean_trajectories = trajectories.copy()
    clean_trajectories[time_column] = pd.to_numeric(
        clean_trajectories[time_column], errors="coerce"
    )
    clean_trajectories[ne_column] = pd.to_numeric(
        clean_trajectories[ne_column], errors="coerce"
    )
    if (
        clean_trajectories[[time_column, ne_column]].isna().any().any()
        or not np.isfinite(clean_trajectories[[time_column, ne_column]]).all().all()
        or (clean_trajectories[time_column] < 0).any()
        or (clean_trajectories[ne_column] <= 0).any()
    ):
        raise ValueError("trajectory times must be nonnegative and Ne must be positive")

    curve_data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    positive_time_grid = None
    for curve in expected_curve_names:
        group = clean_trajectories[
            clean_trajectories[curve_column].astype(str) == curve
        ].sort_values(time_column)
        times = group[time_column].to_numpy(dtype=float)
        sizes = group[ne_column].to_numpy(dtype=float)
        if np.any(np.diff(times) <= 0):
            raise ValueError(f"trajectory times must increase for {curve}")
        positive = times > 0
        if np.count_nonzero(positive) < 2:
            raise ValueError(
                f"trajectory {curve} needs at least two positive time points"
            )
        if positive_time_grid is None:
            positive_time_grid = times[positive]
        elif not np.array_equal(times[positive], positive_time_grid):
            raise ValueError("pointwise EAS curves must share one positive time grid")
        curve_data[curve] = (times, sizes)

    expected_s = np.asarray(expected_selection_coefficients, dtype=float)
    if (
        expected_s.shape != (3,)
        or not np.isfinite(expected_s).all()
        or (expected_s <= 0).any()
        or len(np.unique(expected_s)) != 3
    ):
        raise ValueError(
            "expected_selection_coefficients must contain three unique positive values"
        )
    clean_origins = origins.copy()
    clean_origins[selection_column] = pd.to_numeric(
        clean_origins[selection_column], errors="coerce"
    )
    clean_origins[origin_age_column] = pd.to_numeric(
        clean_origins[origin_age_column], errors="coerce"
    )
    if (
        clean_origins[[selection_column, origin_age_column]].isna().any().any()
        or not np.isfinite(clean_origins[[selection_column, origin_age_column]])
        .all()
        .all()
        or (clean_origins[origin_age_column] <= 0).any()
    ):
        raise ValueError("origin ages and selection coefficients must be positive")
    observed_s = np.sort(clean_origins[selection_column].unique())
    if observed_s.shape != (3,) or not np.allclose(
        observed_s, np.sort(expected_s), rtol=1e-12, atol=0.0
    ):
        raise ValueError(
            "origin table must contain exactly the requested selection coefficients"
        )
    if len(clean_origins) != 9:
        raise ValueError("origin table must contain nine curve-by-selection rows")
    for curve in expected_curve_names:
        for coefficient in expected_s:
            matches = clean_origins[
                (clean_origins[curve_column].astype(str) == curve)
                & np.isclose(
                    clean_origins[selection_column].to_numpy(dtype=float),
                    coefficient,
                    rtol=1e-12,
                    atol=0.0,
                )
            ]
            if len(matches) != 1:
                raise ValueError(
                    "origin table must have one row for every curve-by-selection cell"
                )

    if origin_ne_column in clean_origins:
        clean_origins[origin_ne_column] = pd.to_numeric(
            clean_origins[origin_ne_column], errors="coerce"
        )
        if (
            clean_origins[origin_ne_column].isna().any()
            or not np.isfinite(clean_origins[origin_ne_column]).all()
            or (clean_origins[origin_ne_column] <= 0).any()
        ):
            raise ValueError("origin Ne values must be finite and positive")
    else:
        origin_sizes = []
        for row in clean_origins.to_dict(orient="records"):
            curve = str(row[curve_column])
            times, sizes = curve_data[curve]
            origin_sizes.append(
                _step_value_at_age(times, sizes, float(row[origin_age_column]))
            )
        clean_origins[origin_ne_column] = origin_sizes

    curve_colors = {
        expected_curve_names[0]: "#56b4e9",
        expected_curve_names[1]: "#222222",
        expected_curve_names[2]: "#d55e00",
    }
    curve_labels = {
        expected_curve_names[0]: "q025 (pointwise 2.5%)",
        expected_curve_names[1]: "median (pointwise 50%)",
        expected_curve_names[2]: "q975 (pointwise 97.5%)",
    }
    marker_by_s = {
        float(coefficient): marker
        for coefficient, marker in zip(expected_s, ("o", "s", "^"), strict=True)
    }

    fig, axis = plt.subplots(figsize=(11, 8.5), constrained_layout=True)
    curve_handles = []
    for curve in expected_curve_names:
        times, sizes = curve_data[curve]
        positive = times > 0
        (line,) = axis.step(
            times[positive],
            sizes[positive],
            where="post",
            color=curve_colors[curve],
            linewidth=2.7,
            label=curve_labels[curve],
        )
        curve_handles.append(line)
    for row in clean_origins.to_dict(orient="records"):
        coefficient = float(row[selection_column])
        matched_s = float(expected_s[np.argmin(np.abs(expected_s - coefficient))])
        curve = str(row[curve_column])
        axis.scatter(
            float(row[origin_age_column]),
            float(row[origin_ne_column]),
            marker=marker_by_s[matched_s],
            s=125,
            facecolor=curve_colors[curve],
            edgecolor="black",
            linewidth=0.9,
            zorder=4,
        )
    selection_handles = [
        Line2D(
            [0],
            [0],
            marker=marker_by_s[float(coefficient)],
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=10,
            label=_selection_label(coefficient),
        )
        for coefficient in expected_s
    ]
    curve_legend = axis.legend(
        handles=curve_handles,
        title="EAS Ne trajectory",
        loc="upper left",
        fontsize=_PLOT_FONT["legend"],
        title_fontsize=_PLOT_FONT["legend"],
        frameon=False,
    )
    axis.add_artist(curve_legend)
    axis.legend(
        handles=selection_handles,
        title="De novo origin marker",
        loc="lower right",
        fontsize=_PLOT_FONT["legend"],
        title_fontsize=_PLOT_FONT["legend"],
        frameon=False,
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(
        f"Generations ago (log scale; {generation_time:g} years/generation)",
        fontsize=_PLOT_FONT["axis"],
    )
    axis.set_ylabel("Diploid effective population size, Ne (log scale)", fontsize=18)
    axis.set_title(
        "EAS pointwise Ne sensitivity curves and nine de novo origin ages",
        fontsize=_PLOT_FONT["title"],
    )
    axis.tick_params(axis="both", which="both", labelsize=_PLOT_FONT["tick"])
    axis.grid(which="both", alpha=0.18)
    axis.text(
        0.98,
        0.97,
        "Design provenance\n"
        "3 pointwise Ne paths x 3 selection coefficients\n"
        "Time zero omitted from the logarithmic x-axis",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=11,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85},
    )
    return _save_png_pdf(fig, output_stem)


def _timeline_label(
    name: str,
    generations_ago: float,
    generation_time_years: float,
) -> str:
    kya = generations_ago * generation_time_years / 1_000
    return f"{name}\n{generations_ago:,.0f} g ({kya:,.1f} kya)"


def plot_ancient_eurasia_introgression_timeline(
    output_stem: str | Path,
    *,
    human_neanderthal_split_generations: float = 27_840,
    mutation_origin_generations: float = 2_400,
    introgression_pulse_generations: float = 2_272,
    realized_introgression_pulse_generations: float | None = None,
    han_split_generations: float = 2_016,
    realized_han_split_generations: float | None = None,
    introgression_fraction: float = 0.0296,
    generation_time_years: float = DEFAULT_GENERATION_TIME_YEARS,
) -> dict[str, Path]:
    """Draw the AncientEurasia archaic-variant and selection-onset timeline."""
    split = float(human_neanderthal_split_generations)
    mutation = float(mutation_origin_generations)
    pulse = float(introgression_pulse_generations)
    realized_pulse = (
        pulse
        if realized_introgression_pulse_generations is None
        else float(realized_introgression_pulse_generations)
    )
    han_split = float(han_split_generations)
    realized_han_split = (
        han_split
        if realized_han_split_generations is None
        else float(realized_han_split_generations)
    )
    fraction = float(introgression_fraction)
    generation_time = float(generation_time_years)
    values = np.asarray(
        [
            split,
            mutation,
            pulse,
            realized_pulse,
            han_split,
            realized_han_split,
            generation_time,
            fraction,
        ]
    )
    if not np.isfinite(values).all():
        raise ValueError("timeline parameters must be finite")
    if not split > mutation > pulse > han_split > 0:
        raise ValueError(
            "timeline must satisfy human-Neanderthal split > mutation origin > "
            "introgression pulse > Han split > present"
        )
    if not pulse > realized_han_split > 0:
        raise ValueError(
            "realized Han split must occur after the introgression pulse and "
            "before the present"
        )
    if not mutation > realized_pulse > realized_han_split:
        raise ValueError(
            "realized introgression pulse must occur after the mutation origin "
            "and before the realized Han split"
        )
    if generation_time <= 0:
        raise ValueError("generation_time_years must be positive")
    if not 0 < fraction < 1:
        raise ValueError("introgression_fraction must lie between zero and one")

    zoom_max = max(3_000.0, mutation * 1.12)
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(11, 8.5),
        gridspec_kw={"height_ratios": [1.0, 2.1, 0.30]},
    )
    fig.subplots_adjust(left=0.17, right=0.98, top=0.89, bottom=0.04, hspace=0.58)
    full, recent, note = axes
    human_color = "#0072b2"
    archaic_color = "#d55e00"
    selected_color = "#cc79a7"

    full.plot([split, 0], [1, 1], color=human_color, linewidth=3)
    full.plot([split, 0], [2, 2], color=archaic_color, linewidth=3)
    full.vlines(split, 1, 2, color="#333333", linewidth=2)
    full.scatter([split], [1.5], marker="D", s=75, color="#333333", zorder=3)
    full.axvspan(zoom_max, 0, color="#f0f0f0", alpha=0.9)
    full.annotate(
        _timeline_label("Human-Neanderthal split", split, generation_time),
        xy=(split, 1.5),
        xytext=(split * 0.70, 1.55),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        fontsize=12,
        va="center",
    )
    full.text(
        zoom_max / 2,
        2.35,
        "Recent-event interval enlarged below",
        ha="center",
        fontsize=12,
    )
    full.set_xlim(split * 1.03, 0)
    full.set_ylim(0.65, 2.55)
    full.set_yticks([1, 2], ["Modern-human lineage", "Neanderthal lineage"])
    full.set_xlabel("Generations ago", fontsize=15)
    full.tick_params(axis="both", labelsize=12)
    full.set_title("Full demographic timescale", fontsize=17)
    full.grid(axis="x", alpha=0.18)

    neanderthal_y = 2.05
    ancestral_y = 1.18
    loschbour_y = 0.78
    han_y = 0.35
    recent.plot(
        [zoom_max, realized_pulse],
        [neanderthal_y, neanderthal_y],
        color=archaic_color,
        linewidth=3,
    )
    recent.plot(
        [realized_pulse, 0],
        [neanderthal_y, neanderthal_y],
        color=archaic_color,
        linewidth=1.5,
        linestyle=":",
        alpha=0.7,
    )
    recent.plot(
        [zoom_max, realized_han_split],
        [ancestral_y, ancestral_y],
        color=human_color,
        linewidth=3,
    )
    recent.vlines(
        realized_han_split,
        han_y,
        ancestral_y,
        color=human_color,
        linewidth=2.2,
    )
    recent.plot(
        [realized_han_split, 0],
        [han_y, han_y],
        color=human_color,
        linewidth=3,
        label="Han",
    )
    recent.plot(
        [realized_han_split, 0],
        [loschbour_y, loschbour_y],
        color="#56b4e9",
        linewidth=3,
        label="Loschbour",
    )

    recent.scatter(
        [mutation],
        [neanderthal_y],
        marker="*",
        s=260,
        color="#f0e442",
        edgecolor="black",
        linewidth=0.9,
        zorder=5,
    )
    recent.annotate(
        _timeline_label("Archaic-specific mutation origin", mutation, generation_time),
        xy=(mutation, neanderthal_y),
        xytext=(zoom_max * 0.98, 2.53),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
        fontsize=11,
        ha="left",
    )
    recent.annotate(
        "",
        xy=(realized_pulse, ancestral_y + 0.03),
        xytext=(realized_pulse, neanderthal_y - 0.03),
        arrowprops={"arrowstyle": "-|>", "color": selected_color, "lw": 3},
    )
    if np.isclose(realized_pulse, pulse):
        pulse_label = _timeline_label(
            f"{fraction:.2%} introgression pulse\nselection onset",
            pulse,
            generation_time,
        )
    else:
        requested_pulse_kya = pulse * generation_time / 1_000
        realized_pulse_kya = realized_pulse * generation_time / 1_000
        pulse_label = (
            f"{fraction:.2%} introgression pulse\nselection onset\n"
            f"catalog: {pulse:,.0f} g ({requested_pulse_kya:,.1f} kya)\n"
            f"Q-scaled: {realized_pulse:,.0f} g ({realized_pulse_kya:,.1f} kya)"
        )
    recent.annotate(
        pulse_label,
        xy=(realized_pulse, 1.60),
        xytext=(1_760, 1.62),
        arrowprops={"arrowstyle": "->", "color": selected_color},
        fontsize=11,
        ha="center",
        va="center",
    )
    recent.plot(
        [realized_pulse, realized_han_split],
        [ancestral_y, ancestral_y],
        color=selected_color,
        linewidth=5,
        solid_capstyle="butt",
    )
    recent.vlines(
        realized_han_split,
        han_y,
        ancestral_y,
        color=selected_color,
        linewidth=4,
    )
    recent.plot(
        [realized_han_split, 0],
        [han_y, han_y],
        color=selected_color,
        linewidth=5,
        solid_capstyle="butt",
    )
    recent.scatter(
        [realized_han_split], [ancestral_y], marker="D", s=70, color=human_color
    )
    if np.isclose(realized_han_split, han_split):
        han_split_label = _timeline_label("Han split", han_split, generation_time)
    else:
        requested_kya = han_split * generation_time / 1_000
        realized_kya = realized_han_split * generation_time / 1_000
        han_split_label = (
            "Han split\n"
            f"catalog: {han_split:,.0f} g ({requested_kya:,.1f} kya)\n"
            f"Q-scaled: {realized_han_split:,.0f} g ({realized_kya:,.1f} kya)"
        )
    recent.annotate(
        han_split_label,
        xy=(realized_han_split, ancestral_y),
        xytext=(1_700, 0.02),
        arrowprops={"arrowstyle": "->", "color": human_color},
        fontsize=11,
        ha="center",
    )
    recent.axvline(0, color="#333333", linewidth=1.4)
    recent.text(
        25,
        2.42,
        "Present\n0 g (0 kya)",
        fontsize=11,
        ha="right",
    )
    recent.text(80, han_y + 0.08, "Han", color=selected_color, fontsize=12, ha="right")
    recent.text(
        80, loschbour_y + 0.08, "Loschbour", color="#0072b2", fontsize=12, ha="right"
    )
    recent.set_xlim(zoom_max, 0)
    recent.set_ylim(-0.15, 2.80)
    recent.set_yticks(
        [han_y, loschbour_y, ancestral_y, neanderthal_y],
        ["Han", "Loschbour", "Shared Eurasian ancestor", "Neanderthal"],
    )
    recent.set_xlabel(
        f"Generations ago ({generation_time:g} years/generation)", fontsize=15
    )
    recent.tick_params(axis="both", labelsize=11)
    recent.set_title(
        "Recent-event zoom (chronological and generation-scaled)", fontsize=17
    )
    recent.grid(axis="x", alpha=0.18)

    fig.suptitle(
        "AncientEurasia introgression design: archaic mutation to selected Han allele",
        fontsize=20,
    )
    note.set_axis_off()
    note.text(
        0.5,
        0.5,
        "No ILS by construction: the focal mutation arises on the Neanderthal branch "
        "after the human-Neanderthal split and reaches the shared Loschbour/Han "
        "ancestral lineage only through the 2.96% pulse.",
        transform=note.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "#fff4cc",
            "edgecolor": "#555555",
        },
    )
    return _save_png_pdf(fig, output_stem)


_SIMULATION_STATUS_STYLES = {
    "complete": ("COMPLETE", "#b8e3d3", "Complete"),
    "cached": ("CACHED", "#b9d7ea", "Cached (complete artifact)"),
    "partial": ("PARTIAL", "#c9d7f0", "Partially complete"),
    "timed_out": ("TIMED OUT", "#f9d98c", "Timed out"),
    "exhausted": ("EXHAUSTED", "#e8b66f", "Search budget exhausted"),
    "failed": ("FAILED", "#f3b5a5", "Failed"),
    "not_requested": ("NOT REQ.", "#dddddd", "Not requested"),
}


def _format_elapsed_seconds(value: float) -> str:
    if value < 60:
        return f"t={value:.0f}s"
    if value < 3_600:
        return f"t={value / 60:.1f}m"
    return f"t={value / 3_600:.1f}h"


def _match_design_value(
    value: float,
    expected: np.ndarray,
    *,
    label: str,
) -> float:
    matches = np.flatnonzero(np.isclose(value, expected, rtol=1e-12, atol=1e-12))
    if len(matches) != 1:
        raise ValueError(f"unexpected {label}: {value:g}")
    return float(expected[matches[0]])


def plot_simulation_status_grid(
    simulation_status: pd.DataFrame,
    output_stem: str | Path,
    *,
    scenario_label: str | None = None,
    specification_column: str = "specification_id",
    selection_column: str = "selection_coefficient",
    demography_curve_column: str = "demography_curve",
    target_af_column: str = "target_allele_frequency",
    status_column: str = "status",
    elapsed_column: str = "elapsed_seconds",
    completed_draws_column: str = "completed_external_draws",
    sample_af_column: str = "last_completed_draw_sample_af",
    expected_selection_coefficients: Sequence[float] = (0.01, 0.005, 0.001),
    expected_demography_curves: Sequence[str] = ("q025", "median", "q975"),
    expected_introgression_target_frequencies: Sequence[float] = tuple(
        value / 10 for value in range(1, 10)
    ),
) -> dict[str, Path]:
    """Plot one scenario's complete conditioning-status design matrix.

    Rows are selection coefficients. Columns are the three pointwise EAS Ne
    curves for a no-introgression scenario or the nine target allele
    frequencies for an introgression scenario. Every expected design cell must
    occur exactly once; the function never silently drops or aggregates cells.
    """
    required = [
        specification_column,
        selection_column,
        demography_curve_column,
        target_af_column,
        status_column,
        elapsed_column,
        completed_draws_column,
        sample_af_column,
    ]
    _require_columns(simulation_status, required, "simulation status table")
    if simulation_status.empty:
        raise ValueError("simulation status table is empty")
    if simulation_status[specification_column].isna().any():
        raise ValueError("simulation status specification IDs cannot be missing")
    if simulation_status[specification_column].duplicated().any():
        raise ValueError("simulation status specification IDs must be unique")

    expected_s = np.asarray(expected_selection_coefficients, dtype=float)
    if (
        expected_s.shape != (3,)
        or not np.isfinite(expected_s).all()
        or (expected_s <= 0).any()
        or len(np.unique(expected_s)) != 3
    ):
        raise ValueError(
            "expected_selection_coefficients must contain three unique positive values"
        )
    curves = tuple(map(str, expected_demography_curves))
    if len(curves) != 3 or len(set(curves)) != 3:
        raise ValueError("expected_demography_curves must contain three unique labels")
    target_frequencies = np.asarray(
        expected_introgression_target_frequencies, dtype=float
    )
    if (
        target_frequencies.shape != (9,)
        or not np.isfinite(target_frequencies).all()
        or not ((target_frequencies > 0) & (target_frequencies < 1)).all()
        or len(np.unique(target_frequencies)) != 9
    ):
        raise ValueError(
            "expected introgression targets must contain nine unique frequencies in (0, 1)"
        )

    frame = simulation_status.copy()
    frame[selection_column] = pd.to_numeric(frame[selection_column], errors="coerce")
    if (
        frame[selection_column].isna().any()
        or not np.isfinite(frame[selection_column]).all()
    ):
        raise ValueError("selection coefficients must be finite")
    frame["_selection_key"] = [
        _match_design_value(value, expected_s, label="selection coefficient")
        for value in frame[selection_column].to_numpy(dtype=float)
    ]

    curve_present = frame[demography_curve_column].notna() & frame[
        demography_curve_column
    ].astype(str).str.strip().ne("")
    if curve_present.any():
        if not curve_present.all():
            raise ValueError(
                "one simulation-status scenario cannot mix curve and non-curve rows"
            )
        mode = "no_introgression"
        frame["_column_key"] = frame[demography_curve_column].astype(str)
        observed_curves = set(frame["_column_key"])
        if observed_curves != set(curves):
            raise ValueError(
                "no-introgression status grid must contain q025, median, and q975 "
                f"curves; observed={sorted(observed_curves)}"
            )
        column_keys: list[str | float] = list(curves)
        column_labels = [
            "q025\n(2.5%)",
            "median\n(50%)",
            "q975\n(97.5%)",
        ]
        x_label = "Pointwise EAS Ne sensitivity curve"
        default_scenario_label = "No-introgression de novo sweep"
    else:
        mode = "introgression"
        target = pd.to_numeric(frame[target_af_column], errors="coerce")
        if target.isna().any() or not target.between(0, 1, inclusive="neither").all():
            raise ValueError(
                "introgression target allele frequencies must lie strictly between zero and one"
            )
        frame["_column_key"] = [
            _match_design_value(
                value, target_frequencies, label="target allele frequency"
            )
            for value in target.to_numpy(dtype=float)
        ]
        column_keys = list(map(float, target_frequencies))
        column_labels = [f"{value:.0%}" for value in target_frequencies]
        x_label = "Target final allele-frequency threshold"
        default_scenario_label = "AncientEurasia introgression sweep"

    normalized_status = frame[status_column].astype(str).str.strip().str.lower()
    unknown_statuses = sorted(
        set(normalized_status).difference(_SIMULATION_STATUS_STYLES)
    )
    if unknown_statuses:
        raise ValueError(
            "unsupported simulation statuses: " + ", ".join(unknown_statuses)
        )
    frame["_status"] = normalized_status

    for column, label in (
        (elapsed_column, "elapsed seconds"),
        (completed_draws_column, "completed external draws"),
        (sample_af_column, "last completed draw sample AF"),
    ):
        raw = frame[column]
        present = raw.notna() & raw.astype(str).str.strip().ne("")
        converted = pd.to_numeric(raw, errors="coerce")
        if converted[present].isna().any():
            raise ValueError(f"{label} must be numeric when present")
        frame[column] = converted
        finite = frame[column].dropna()
        if not np.isfinite(finite).all():
            raise ValueError(f"{label} must be finite when present")
    if (frame[elapsed_column].dropna() < 0).any():
        raise ValueError("elapsed seconds cannot be negative")
    draws = frame[completed_draws_column].dropna()
    if (draws < 0).any() or not np.equal(draws, np.floor(draws)).all():
        raise ValueError("completed external draws must be nonnegative integers")
    sample_af = frame[sample_af_column].dropna()
    if not sample_af.between(0, 1, inclusive="both").all():
        raise ValueError("last completed draw sample AF must lie between zero and one")

    expected_cells = {
        (float(selection), column_key)
        for selection in expected_s
        for column_key in column_keys
    }
    observed_cell_list = list(
        zip(frame["_selection_key"], frame["_column_key"], strict=True)
    )
    observed_cells = set(observed_cell_list)
    if len(observed_cell_list) != len(observed_cells):
        raise ValueError("simulation status table has duplicate design cells")
    if observed_cells != expected_cells:
        missing = sorted(expected_cells.difference(observed_cells), key=str)
        extra = sorted(observed_cells.difference(expected_cells), key=str)
        raise ValueError(
            "simulation status design grid is incomplete or unexpected: "
            f"missing={missing}, extra={extra}"
        )
    row_lookup = {
        (float(row["_selection_key"]), row["_column_key"]): row
        for row in frame.to_dict(orient="records")
    }

    n_columns = len(column_keys)
    fig, axis = plt.subplots(figsize=(11, 8.5))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.76, bottom=0.18)
    status_font = 12 if mode == "no_introgression" else 8.5
    detail_font = 10 if mode == "no_introgression" else 7.5
    for row_index, selection in enumerate(expected_s):
        for column_index, column_key in enumerate(column_keys):
            row = row_lookup[(float(selection), column_key)]
            status = str(row["_status"])
            display_status, face_color, _ = _SIMULATION_STATUS_STYLES[status]
            axis.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1,
                    1,
                    facecolor=face_color,
                    edgecolor="white",
                    linewidth=3,
                )
            )
            details = []
            if pd.notna(row[elapsed_column]):
                details.append(_format_elapsed_seconds(float(row[elapsed_column])))
            if pd.notna(row[completed_draws_column]):
                details.append(f"draws={int(row[completed_draws_column]):,}")
            if pd.notna(row[sample_af_column]):
                details.append(f"last AF={float(row[sample_af_column]):.1%}")
            axis.text(
                column_index,
                row_index - 0.18,
                display_status,
                ha="center",
                va="center",
                fontsize=status_font,
                fontweight="bold",
            )
            if details:
                axis.text(
                    column_index,
                    row_index + 0.20,
                    "\n".join(details),
                    ha="center",
                    va="center",
                    fontsize=detail_font,
                    linespacing=1.1,
                )
    axis.add_patch(
        Rectangle(
            (-0.5, -0.5),
            n_columns,
            len(expected_s),
            fill=False,
            edgecolor="#333333",
            linewidth=1.8,
        )
    )
    axis.set_xlim(-0.5, n_columns - 0.5)
    axis.set_ylim(len(expected_s) - 0.5, -0.5)
    axis.set_xticks(range(n_columns), column_labels)
    axis.set_yticks(
        range(len(expected_s)), [_selection_label(value) for value in expected_s]
    )
    axis.set_xlabel(x_label, fontsize=_PLOT_FONT["axis"], labelpad=12)
    axis.set_ylabel("Selection coefficient", fontsize=_PLOT_FONT["axis"], labelpad=10)
    axis.tick_params(axis="x", labelsize=14 if mode == "no_introgression" else 12)
    axis.tick_params(axis="y", labelsize=_PLOT_FONT["tick"])
    for spine in axis.spines.values():
        spine.set_visible(False)

    title = scenario_label or default_scenario_label
    fig.suptitle(
        f"{title}: simulation status and conditioning diagnostics",
        fontsize=21,
        y=0.97,
    )
    fig.text(
        0.5,
        0.915,
        f"Every expected design cell is shown ({len(expected_cells)} total)",
        ha="center",
        fontsize=14,
    )
    legend_handles = [
        Patch(facecolor=color, edgecolor="#555555", label=legend_label)
        for _, color, legend_label in _SIMULATION_STATUS_STYLES.values()
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.88),
        ncol=5,
        fontsize=12,
        frameon=False,
    )
    fig.text(
        0.5,
        0.045,
        "Cell annotations: recorded final invocation time (not cumulative across retries); "
        "completed external draws; "
        "sample AF from the last completed draw (when available).",
        ha="center",
        fontsize=11,
    )
    return _save_png_pdf(fig, output_stem)
