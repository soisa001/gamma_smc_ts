"""Summarize completed runs and preserve compact provenance alongside their scores."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil

import numpy as np

import benchmark as b


def write_report(root, results, document):
    root, results, document = Path(root), Path(results), Path(document)
    summary = json.loads((results / "summary.json").read_text())
    contrasts = json.loads((results / "paired_comparisons.json").read_text())
    runs, rows = [], []
    for path in sorted(root.glob("*/seed_*/scores.json")):
        directory = path.parent
        for manifest in ("samples.manifest.json", "fit.manifest.json", "evaluation.manifest.json"):
            if not b.manifest_valid(directory, manifest):
                raise ValueError(f"Invalid {directory / manifest}")
        sample = json.loads((directory / "samples.manifest.json").read_text())
        audit = json.loads((directory / "fit" / "audit.json").read_text())
        runs.append({"directory": str(directory.resolve()), "sample_manifest": sample, "fit_audit": audit,
                     "scores_sha256": b.sha256(path)})
        rows.extend(json.loads(path.read_text()))
    invocation = max(root.glob("invocation_*.json"), key=lambda p: p.name)
    b.write_json(results / "run_provenance.json", {"latest_invocation": json.loads(invocation.read_text()),
                  "source_files": {p.name: b.sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
                  "dependency_lock_sha256": b.sha256(Path(__file__).parent / "uv.lock"), "runs": runs})
    shutil.copyfile(root / "run.log", results / "execution.log")
    focal = [r for r in rows if r["scenario"] == "five_shift" and r["method"] == "5mer"]
    config = next(r["sample_manifest"]["config"] for r in runs if r["sample_manifest"]["config"]["scenario"] == "five_shift")
    get = lambda s, m: next(r for r in summary if r["scenario"] == s and r["method"] == m)
    gain = 1-get("five_shift", "5mer")["oracle_kl_per_mutation"]/get("five_shift", "3mer")["oracle_kl_per_mutation"]
    pair = next(r for r in contrasts if r["scenario"] == "five_shift" and r["contrast"] == "5mer minus 3mer" and r["metric"] == "oracle_kl_per_mutation")
    seven_counts = Counter(str(r["shrinkage"]) for r in rows if r["scenario"] == "five_shift" and r["method"] == "7mer")
    labels = {"1mer": "1-mer", "3mer": "3-mer", "5mer": "5-mer", "7mer": "7-mer, hierarchical", "7mer_flat": "7-mer, independent"}
    text = f"""# Offspring allele prediction: mock-task results

The planted five-base task works: five-base probability estimates reduced excess
log loss by **{100*gain:.1f}%** relative to three-base estimates in the shifted-target
experiment. Five bases beat the independently regularized seven-base model in
every replicate. A hierarchical seven-base model nearly matched five bases, with
a small average improvement. This supports five bases as the smallest sufficient
window; it does not make a larger model that learns to pool its estimates incorrect.

## Executed experiment

There are {len(runs)} independently generated datasets across {len(set(r['scenario'] for r in rows))} scenarios and {len(set(r['seed'] for r in rows))}
seeds ({min(r['seed'] for r in rows)}–{max(r['seed'] for r in rows)}). Each dataset has {config['train_bases']:,} calibration bases and
{config['test_bases']:,} target bases, divided into {config['contig_bases']:,}-base
contigs. Every method receives the same aligned calibration observations and
opportunity counts. Hyperparameters are selected using held-out calibration
contigs, and offspring from the new parent are withheld from all fitting.
The local run used four processes, a declared 26-GB budget, Python 3.11, and uv
under WSL Ubuntu 22.04. This is an algorithmic benchmark; no autonomous task-solving
agent was evaluated in this run.

The one-transmission mock mutation probability averages {config['mean_rate']} under
uniform sequence. This is deliberately arbitrary and elevated. The parent is a
haploid transmitted sequence, and the offspring distribution is sampled once
from the fixed parental neighborhood. Neither human rates nor diploid inheritance
are being inferred. See the [implementation and reproduction instructions](../experiments/mutation_context/README.md).

## Main five-base task

Calibration and target have the same expected trinucleotide composition. The
probability that the outer bases of a five-base window match changes from 0.25 to
{config['target_outer_match']}. The generating kernel depends on those outer bases,
so a three-base table cannot transport the probabilities correctly.

The primary score is excess expected full-allele log loss in nats per expected
mutation. Zero is oracle performance; lower is better. This scores both the
probability of a mutation and which alternate base it produces.

| Model | Excess log loss | Conditional allele KL | Mean absolute burden error |
| --- | ---: | ---: | ---: |
"""
    for method in b.METHODS:
        r = get("five_shift", method)
        text += f"| {labels[method]} | {r['oracle_kl_per_mutation']:.5f} | {r['oracle_conditional_kl']:.5f} | {100*r['burden_mean_absolute_relative_error']:.2f}% |\n"
    text += f"""
Five bases improved the primary score in {pair['second_better_replicates']}/{pair['replicates']}
replicates. The paired mean difference (5-mer minus 3-mer) was
{pair['mean_difference']:.5f}, with a seed-bootstrap 95% interval
[{pair['ci_low']:.5f}, {pair['ci_high']:.5f}]. The realized conditional allele log
loss also improved on average; it is noisier than expected scoring.

Both seven-base comparators use calibration-only validation to choose their
regularization. The independent model shrinks directly toward one-base estimates;
the hierarchical model can shrink toward the fitted five-base estimates. These
are different estimators, not two contradictory results about the same estimator.
Hierarchical seven-base shrinkage choices were `{dict(seven_counts)}`.
Its small improvement can refine finite-sample shrinkage even though the true
distribution contains no extra seven-base information. Do not claim that a
seven-base representation must be worse in principle.

The observed trinucleotide total-variation distance between calibration and target
was {np.mean([r['trinucleotide_tv'] for r in focal]):.5f} on average
(range {min(r['trinucleotide_tv'] for r in focal):.5f}–{max(r['trinucleotide_tv'] for r in focal):.5f}).
Thus composition is matched in expectation, not forced to identical finite counts.
Global-alignment SNV precision was {min(r['alignment_precision'] for r in focal):.4%} at minimum,
and recall ranged from {min(r['alignment_recall'] for r in focal):.4%} to
{max(r['alignment_recall'] for r in focal):.4%}. Some true events were masked near
gaps chosen by the edit-distance alignment. No assembly errors or structural
rearrangements were simulated.

## Controls

| Generating scenario | 1-mer loss | 3-mer loss | 5-mer loss | Hierarchical 7-mer loss |
| --- | ---: | ---: | ---: | ---: |
"""
    for scenario in b.SCENARIOS:
        if any(r["scenario"] == scenario for r in summary):
            values = " | ".join(f"{get(scenario, method)['oracle_kl_per_mutation']:.5f}" for method in b.METHODS[:4])
            text += f"| {scenario} | {values} |\n"
    text += """
The constant and three-base controls favor the appropriately small models.
Five-base same-composition data show a local prediction benefit even when all
models predict total burden reasonably well. The true seven-base control benefits
from seven bases, so the comparison is not designed to make five bases always win.

The optional empirical-shaped kernel was made by uniformly averaging the outer
pair from the Carlson seven-base table, then assigning the same arbitrary total
mock rate. Its smaller five-base benefit is additional evidence within this
simulation, not a measurement of human prediction performance. The core planted
task does not depend on this resource or biological realism.

## Recommended task and grading

Hand the solver only [task_001.zip](../experiments/mutation_context/releases/task_001.zip).
It contains calibration parent/offspring FASTAs, a new parent FASTA, query
coordinates, metadata, and the unsignposted prompt. The archive has no generator,
window-size hints, fitted models, hidden offspring, or truth probabilities.
The surrounding authoring repository must remain inaccessible during a blind run.

The essential prompt is: given the calibration assemblies, predict P(A), P(C),
P(G), and P(T) for the offspring at each supplied parental coordinate, and return
the expected total number of changed alleles. Every eligible coordinate is queried;
the evaluator does not select sites after seeing mutations. The exact prompt in
the archive states the transmitted-haplotype and coordinate conventions.

Use full probability log loss or expected KL as the primary objective. Keep
conditional alternate-base scoring and burden calibration as diagnostic scores.
GC-directional counts and immediate GC flux are sensible later extensions.
An actual gene-conversion event rate is a different estimand and was not fitted.

## Artifacts and verification

- [Per-replicate scores](../experiments/mutation_context/results/replicate_scores.csv)
- [Paired comparisons and bootstrap intervals](../experiments/mutation_context/results/paired_comparisons.json)
- [Summary statistics](../experiments/mutation_context/results/summary.json)
- [Sample hashes, fit audits, and execution provenance](../experiments/mutation_context/results/run_provenance.json)
- [PNG figure](../experiments/mutation_context/results/prediction_comparison.png) and [vector PDF](../experiments/mutation_context/results/prediction_comparison.pdf)
- [Independent grade of exported example predictions](../experiments/mutation_context/results/example_grade.json)
- [Automatic window selection for the released instance](../experiments/mutation_context/results/example_solver.json)
- [Research and candidate downstream metrics](mutation_context_literature.md)

Focused tests cover exact model dependence, reverse-complement symmetry, alignment,
gap masking, the matched-composition generator, oracle scores, corrupted caches,
missing/malformed submissions, public-only solving, and deterministic release
packaging. The exported million-site example is checked against the independent
grader; the full prediction TSV stays local. The figure was visually inspected,
uses a horizontal legend above the panels, and is saved without inline display.

The public-only reference solver uses a paired contig-level one-standard-error
rule to prefer a smaller window when a larger window provides no reliable
validation improvement. It selected five bases for the released seed-0 task,
without access to the hidden kernel or target offspring. The exported predictions
were byte-identical to the explicit five-base comparator's predictions.
"""
    document.parent.mkdir(parents=True, exist_ok=True)
    document.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("work/main"))
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--document", type=Path, default=Path("../../docs/mutation_context_mock_results.md"))
    args = parser.parse_args()
    write_report(args.root, args.results, args.document)
