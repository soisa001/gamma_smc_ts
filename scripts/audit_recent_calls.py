#!/usr/bin/env python3
"""Independently audit complete regional predictions and baseline preservation."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--baseline-analysis", type=Path, required=True)
    args = parser.parse_args()
    out = args.analysis
    provenance = json.loads((out / "provenance.json").read_text())
    for name, spec in provenance["outputs"].items():
        path = out / name
        assert path.stat().st_size == spec["bytes"], name
        assert file_hash(path) == spec["sha256"], name
    keys = [
        "onset_years",
        "source",
        "method",
        "cutoff_years",
        "policy",
        "target",
        "mode",
    ]
    count, called, identities, groups = Counter(), Counter(), set(), defaultdict(set)

    def key(row):
        return tuple(row[k] for k in keys)

    with gzip.open(out / "held_out_regions.csv.gz", "rt", newline="") as stream:
        for row in csv.DictReader(stream):
            group = key(row)
            identity = group + (row["task_id"],)
            assert identity not in identities, identity
            identities.add(identity)
            groups[group].add(row["task_id"])
            count[group] += 1
            called[group] += int(row["called"])
            assert int(row["called"]) in (0, 1)
            assert 0 < float(row["p"]) <= 1
            if row["policy"] == "neutral_target":
                assert int(row["called"]) == int(
                    float(row["p"]) <= float(row["target"])
                )
            assert row["task_id"].startswith(
                "neutral/"
                if row["mode"] == "neutral"
                else "onset" + row["onset_years"] + "/"
            )
    with (out / "region_calls.csv").open(newline="") as stream:
        summaries = list(csv.DictReader(stream))
    assert len(summaries) == len(count) == 960
    for row in summaries:
        k = key(row)
        expected = 1000 if row["mode"] == "neutral" else 100
        assert count[k] == int(row["n"]) == len(groups[k]) == expected
        assert called[k] == int(row["called"])
        assert abs(called[k] / count[k] - float(row["region_call_rate"])) < 1e-12
    baseline_keys = ["source", "method", "cutoff_years", "policy", "target", "mode"]
    current = {}
    for row in summaries:
        if row["onset_years"] == "50000" and row["source"] in ("mean", "truth"):
            item = dict(row, source="decoded" if row["source"] == "mean" else "truth")
            current[tuple(item[k] for k in baseline_keys)] = item
    matched = 0
    with (args.baseline_analysis / "region_calls.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["haplotypes"] != "400":
                continue
            other = current[tuple(row[k] for k in baseline_keys)]
            assert (
                abs(float(row["region_call_rate"]) - float(other["region_call_rate"]))
                < 1e-12
            )
            assert int(row["n"]) == int(other["n"])
            matched += 1
    assert matched == 192
    result = dict(
        passed=True,
        unique_held_out_predictions=len(identities),
        summaries_recomputed=len(summaries),
        existing_baseline_rows_reproduced=matched,
        output_hashes_checked=len(provenance["outputs"]),
        audit_source_sha256=file_hash(Path(__file__)),
    )
    (out.parent / "analysis_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
