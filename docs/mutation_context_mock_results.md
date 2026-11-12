# Offspring allele prediction: mock-task results

The planted five-base task works: five-base probability estimates reduced excess
log loss by **79.7%** relative to three-base estimates in the shifted-target
experiment. Five bases beat the independently regularized seven-base model in
every replicate. A hierarchical seven-base model nearly matched five bases, with
a small average improvement. This supports five bases as the smallest sufficient
window; it does not make a larger model that learns to pool its estimates incorrect.

## Executed experiment

There are 72 independently generated datasets across 6 scenarios and 12
seeds (0–11). Each dataset has 4,000,000 calibration bases and
1,000,000 target bases, divided into 50,000-base
contigs. Every method receives the same aligned calibration observations and
opportunity counts. Hyperparameters are selected using held-out calibration
contigs, and offspring from the new parent are withheld from all fitting.
The local run used four processes, a declared 26-GB budget, Python 3.11, and uv
under WSL Ubuntu 22.04. This is an algorithmic benchmark; no autonomous task-solving
agent was evaluated in this run.

The one-transmission mock mutation probability averages 0.005 under
uniform sequence. This is deliberately arbitrary and elevated. The parent is a
haploid transmitted sequence, and the offspring distribution is sampled once
from the fixed parental neighborhood. Neither human rates nor diploid inheritance
are being inferred. See the [implementation and reproduction instructions](../experiments/mutation_context/README.md).

## Main five-base task

Calibration and target have the same expected trinucleotide composition. The
probability that the outer bases of a five-base window match changes from 0.25 to
0.8. The generating kernel depends on those outer bases,
so a three-base table cannot transport the probabilities correctly.

The primary score is excess expected full-allele log loss in nats per expected
mutation. Zero is oracle performance; lower is better. This scores both the
probability of a mutation and which alternate base it produces.

| Model | Excess log loss | Conditional allele KL | Mean absolute burden error |
| --- | ---: | ---: | ---: |
| 1-mer | 0.25009 | 0.04702 | 28.36% |
| 3-mer | 0.15432 | 0.02028 | 28.36% |
| 5-mer | 0.03136 | 0.01492 | 6.55% |
| 7-mer, hierarchical | 0.03099 | 0.01505 | 6.40% |
| 7-mer, independent | 0.15914 | 0.05421 | 19.25% |

Five bases improved the primary score in 12/12
replicates. The paired mean difference (5-mer minus 3-mer) was
-0.12296, with a seed-bootstrap 95% interval
[-0.12807, -0.11675]. The realized conditional allele log
loss also improved on average; it is noisier than expected scoring.

Both seven-base comparators use calibration-only validation to choose their
regularization. The independent model shrinks directly toward one-base estimates;
the hierarchical model can shrink toward the fitted five-base estimates. These
are different estimators, not two contradictory results about the same estimator.
Hierarchical seven-base shrinkage choices were `{'100000.0': 2, '1000000000000.0': 8, '10000.0': 2}`.
Its small improvement can refine finite-sample shrinkage even though the true
distribution contains no extra seven-base information. Do not claim that a
seven-base representation must be worse in principle.

The observed trinucleotide total-variation distance between calibration and target
was 0.00606 on average
(range 0.00473–0.00742).
Thus composition is matched in expectation, not forced to identical finite counts.
Global-alignment SNV precision was 100.0000% at minimum,
and recall ranged from 99.3997% to
99.6505%. Some true events were masked near
gaps chosen by the edit-distance alignment. No assembly errors or structural
rearrangements were simulated.

## Controls

| Generating scenario | 1-mer loss | 3-mer loss | 5-mer loss | Hierarchical 7-mer loss |
| --- | ---: | ---: | ---: | ---: |
| constant | 0.00020 | 0.00020 | 0.00029 | 0.00053 |
| three | 0.16548 | 0.00262 | 0.00269 | 0.00269 |
| five_iid | 0.23239 | 0.09866 | 0.03476 | 0.03464 |
| five_shift | 0.25009 | 0.15432 | 0.03136 | 0.03099 |
| seven_shift | 0.33929 | 0.25696 | 0.16083 | 0.12673 |
| empirical_five | 0.46013 | 0.02532 | 0.01670 | 0.01670 |

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
