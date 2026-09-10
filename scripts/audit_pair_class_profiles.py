#!/usr/bin/env python3
"""Independently check class partitions and recompute profile averages."""

import argparse
import csv
from collections import defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.out / "manifest.json").read_text())
    cutoffs = manifest["config"]["tmrca_cutoffs_years"]
    metrics = ["mean_tmrca_years", *(f"frac_recent_{t}" for t in cutoffs)]
    sums, n = defaultdict(float), defaultdict(int)
    row_count, partition_count = 0, 0
    for key in manifest["sample_ids"]:
        folder = args.out / "profiles" / key
        receipt = json.loads((folder / "complete.json").read_text())
        profile = folder / "profile.csv.gz"
        assert sha(profile) == receipt["profile_sha256"]
        expected_pairs = dict(receipt["pair_counts"])
        expected_pairs["all pairs"] = sum(expected_pairs.values())
        assert expected_pairs["all pairs"] == manifest["config"]["haplotype_pairs"]
        by_position = defaultdict(dict)
        with gzip.open(profile, "rt", newline="") as stream:
            for row in csv.DictReader(stream):
                row_count += 1
                assert int(row["replicate"]) == receipt["replicate"]
                assert int(row["onset_years"]) == receipt["onset_years"]
                label = row["pair_class"]
                size = int(row["n_pairs"])
                assert size == expected_pairs[label]
                identity = row["source"], row["position_0based"]
                assert label not in by_position[identity]
                by_position[identity][label] = row
                counts = [int(row[f"count_recent_{t}"]) for t in cutoffs]
                assert all(0 <= c <= size for c in counts)
                assert counts == sorted(counts)
                for cutoff, count in zip(cutoffs, counts):
                    value = row[f"frac_recent_{cutoff}"]
                    if size:
                        assert abs(float(value) - count / size) < 1e-12
                    else:
                        assert value == ""
                for metric in metrics:
                    value = row[metric]
                    if not size:
                        assert value == ""
                        continue
                    value = float(value)
                    assert math.isfinite(value) and value >= 0
                    group = (
                        row["onset_years"],
                        row["source"],
                        label,
                        row["position_0based"],
                        metric,
                    )
                    sums[group] += value
                    n[group] += 1
        for identity, classes in by_position.items():
            assert set(classes) == {"alt/alt", "ref/ref", "alt/ref", "all pairs"}
            total = classes["all pairs"]
            subsets = [classes[label] for label in ("alt/alt", "ref/ref", "alt/ref")]
            assert sum(int(r["n_pairs"]) for r in subsets) == int(total["n_pairs"])
            for t in cutoffs:
                assert sum(int(r[f"count_recent_{t}"]) for r in subsets) == int(
                    total[f"count_recent_{t}"]
                )
            age_sum = sum(
                int(r["n_pairs"]) * float(r["mean_tmrca_years"])
                for r in subsets
                if int(r["n_pairs"])
            )
            assert math.isclose(
                age_sum / int(total["n_pairs"]),
                float(total["mean_tmrca_years"]),
                rel_tol=1e-12,
                abs_tol=1e-8,
            )
            partition_count += 1
    summaries_checked = 0
    with (args.out / "analysis/mean_profiles.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            for metric in metrics:
                group = (
                    row["onset_years"],
                    row["source"],
                    row["pair_class"],
                    row["position_0based"],
                    metric,
                )
                assert n[group] == int(row[metric + "_n_regions"])
                if n[group]:
                    assert math.isclose(
                        sums[group] / n[group],
                        float(row[metric]),
                        rel_tol=1e-8,
                        abs_tol=1e-8,
                    ), group
                else:
                    assert row[metric] == ""
                summaries_checked += 1
    provenance = json.loads((args.out / "analysis/provenance.json").read_text())
    for name, spec in provenance["outputs"].items():
        p = args.out / "analysis" / name
        assert p.stat().st_size == spec["bytes"] and sha(p) == spec["sha256"]
    npos = manifest["config"]["scored_length_bp"] // manifest["config"]["stride_bp"]
    assert row_count == len(manifest["sample_ids"]) * 8 * npos
    assert partition_count == len(manifest["sample_ids"]) * 2 * npos
    result = dict(
        passed=True,
        regions=len(manifest["sample_ids"]),
        profile_rows_checked=row_count,
        pair_class_partitions_checked=partition_count,
        summary_means_recomputed=summaries_checked,
        analysis_provenance_sha256=sha(args.out / "analysis/provenance.json"),
        audit_source_sha256=sha(Path(__file__)),
    )
    (args.out / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
