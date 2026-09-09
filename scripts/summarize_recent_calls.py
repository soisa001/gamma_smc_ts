#!/usr/bin/env python3
"""Render audited regional summaries as reproducible Markdown tables.

This reads completed analysis only. It never selects or fits a detector, reads
genotypes, reruns simulation, or emits a figure.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics


SOURCES = ("mean", "median", "prob80", "prob90", "truth")
CUTOFFS = (5000, 10000, 20000, 30000, 40000, 50000)


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one(rows, **criteria):
    matches = [r for r in rows if all(r[k] == str(v) for k, v in criteria.items())]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for {criteria}; found {len(matches)}")
    return matches[0]


def cell(rows, **criteria):
    selected = one(rows, mode="selected", **criteria)
    neutral = one(rows, mode="neutral", **criteria)
    assert int(selected["n"]) == 100 and int(neutral["n"]) == 1000
    return (
        f"{float(selected['region_call_rate']):.0%} / "
        f"{float(neutral['region_call_rate']):.1%}"
    )


def table(header, rows):
    return "\n".join(
        [
            "| " + " | ".join(header) + " |",
            "|" + "|".join("---" for _ in header) + "|",
            *("| " + " | ".join(map(str, row)) + " |" for row in rows),
        ]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.analysis
    audit = json.loads((root.parent / "analysis_audit.json").read_text())
    assert audit["passed"] and audit["summaries_recomputed"] == 960
    provenance_path = root / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    for name, expected in provenance["outputs"].items():
        path = root / name
        assert path.stat().st_size == expected["bytes"]
        assert digest(path) == expected["sha256"], name
    calls = read_rows(root / "region_calls.csv")
    grid = read_rows(root / "paired_grid.csv")
    meta = read_rows(root / "replicates.csv")
    assert len(calls) == 960 and len(grid) == 10800 and len(meta) == 1200
    sections = [
        "# Audited tables: posterior calls and selection onset",
        "Every rate cell is **selected regions called / neutral regions called** "
        "(100 selected and 1,000 shared neutral regions). A call anywhere in the "
        "isolated 10 Mb region counts. All selected replicates are retained.",
        "## Allele frequencies",
    ]
    rows = []
    for onset in (50000, 10000):
        values = [
            float(r["sample_af"])
            for r in meta
            if r["mode"] == "selected" and int(r["onset_years"]) == onset
        ]
        assert len(values) == 100
        rows.append(
            (
                onset // 1000,
                len(values),
                f"{statistics.mean(values):.6f}",
                f"{statistics.median(values):.5f}",
                f"{min(values):.4f}",
                f"{max(values):.4f}",
                sum(v == 1 for v in values),
                sum(v < 0.1 for v in values),
            )
        )
    sections.append(
        table(
            ["Onset (kya)", "n", "Mean", "Median", "Min", "Max", "Fixed", "AF < 0.1"],
            rows,
        )
    )
    for target in (0.05, 0.1):
        sections.append(f"## Regional maximum at nominal alpha {target:g}")
        for onset in (50000, 10000):
            rows = []
            for cutoff in CUTOFFS:
                rows.append(
                    [cutoff // 1000]
                    + [
                        cell(
                            calls,
                            onset_years=onset,
                            source=source,
                            method="maximum_frac_recent",
                            cutoff_years=cutoff,
                            policy="neutral_target",
                            target=target,
                        )
                        for source in SOURCES
                    ]
                )
            sections.extend(
                [
                    f"Selection onset **{onset // 1000} kya**.",
                    table(["T (kya)", *SOURCES], rows),
                ]
            )
    sections.extend(
        [
            "## Thresholds trained for 70% selected sensitivity",
            "These fixed T=10 and T=50 kya examples use the raw regional maximum. "
            "Each threshold is trained on 80 selected regions, then evaluated on "
            "held-out regions. The table reports realized held-out rates, which "
            "need not equal the training target. All other cutoffs and the "
            "standardized score remain available in region_calls.csv.",
        ]
    )
    rows = []
    for source in SOURCES:
        for cutoff in (10000, 50000):
            rows.append(
                [source, cutoff // 1000]
                + [
                    cell(
                        calls,
                        onset_years=onset,
                        source=source,
                        method="maximum_frac_recent",
                        cutoff_years=cutoff,
                        policy="selected_target",
                        target=0.7,
                    )
                    for onset in (50000, 10000)
                ]
            )
    sections.append(table(["Rule", "T (kya)", "Onset 50 kya", "Onset 10 kya"], rows))
    sections.extend(
        [
            "## Stride p <= 0.001 with consecutive evidence",
            "These examples use 10 kb reporting strides, with five or ten "
            "significant strides required and no gaps. The pointwise p threshold "
            "does not imply a regional 0.001 error rate. The full exploratory "
            "grid is in paired_grid.csv; no winner from that grid has received "
            "separate confirmation.",
        ]
    )
    for cutoff in (10000, 50000):
        rows = []
        for source in SOURCES:
            for minimum in (5, 10):
                rows.append(
                    [source, minimum]
                    + [
                        cell(
                            grid,
                            onset_years=onset,
                            source=source,
                            cutoff_years=cutoff,
                            alpha=0.001,
                            stride_bp=10000,
                            max_gap_strides=0,
                            min_significant_strides=minimum,
                        )
                        for onset in (50000, 10000)
                    ]
                )
        sections.extend(
            [
                f"TMRCA cutoff **{cutoff // 1000} kya**.",
                table(
                    ["Rule", "Required strides", "Onset 50 kya", "Onset 10 kya"], rows
                ),
            ]
        )
    sections.extend(
        [
            "## Provenance",
            f"Analysis source: `{root.resolve()}`. Input/output inventory SHA-256: "
            f"`{digest(provenance_path)}`. Table-rendering script SHA-256: "
            f"`{digest(Path(__file__))}`. The independent audit recomputed all "
            "960 regional summaries from 528,000 unique held-out predictions "
            "and reproduced the 192 existing mean/truth baseline rows.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text("\n\n".join(sections) + "\n", encoding="utf-8", newline="\n")
    temporary.replace(args.output)
    print(
        json.dumps(
            {"output": str(args.output.resolve()), "sha256": digest(args.output)}
        )
    )


if __name__ == "__main__":
    main()
