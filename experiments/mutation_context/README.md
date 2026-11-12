# Predict offspring alleles from assembled genomes

This mock task asks for P(A), P(C), P(G), and P(T) at supplied coordinates in a new
parent assembly. The solver learns from calibration parent–offspring assemblies.
The public prompt never names sequence context or the intended model size.

The main synthetic process has an exact five-base dependency. A three-base model
misses part of the signal. An independently fitted seven-base table has more
parameters and can lose through sampling variance; a hierarchical seven-base
model can pool to five bases and match its performance. The task should reward
good probabilities, not require a particular implementation or forbid a model
that learns an equivalent solution.

See [the research note](../../docs/mutation_context_literature.md) for background
and alternative downstream metrics. Biological realism is not a requirement of
the core mock task. The optional published-table experiment is an empirically
shaped synthetic control, not an estimate of human germline rates.

## Public task and evaluator boundary

Each generated instance has:

```text
public/
  TASK.md
  metadata.json
  calibration_parent.fasta
  calibration_offspring.fasta
  new_parent.fasta
  queries.tsv
private/
  new_offspring.fasta
  kernel.npy
```

Only `public/` belongs in the solver's environment. Do not expose this repository,
the generator, solver, research note, scenario names, private directory, fitted
models, or results in a blind task. The public release archive contains only those
six public files. The evaluator and reference solver here are authoring tools.

Coordinates are one-based in the parent assembly. All eligible sites are queried;
query selection does not depend on which sites mutated. The target child is
entirely withheld, including its flanks. The parent sequence supplies the known
flanks. This is a haploid transmitted-allele task, so Mendelian segregation,
diploid phasing, and recombination are outside the observation model.

## Generator and comparison

Each seed generates 80 calibration contigs of 50 kb (4 Mb) and 20 target contigs
(1 Mb). Seeds are 0–11. The simulator independently samples one offspring allele
per eligible site from a categorical distribution based on the fixed parent.
It is an exact one-transmission mock model, with no recurrent evolution or
changing neighborhoods during that transmission. The mean mutation probability
under uniform sequence is 0.005. This deliberately elevated synthetic rate is
not a human mutation-rate claim. The first/last three bases of each contig are
excluded identically for all models.

The synthetic channels begin with transition weight 2 and transversion weight 1.
The 3-base rule adds a sixfold CpG transition multiplier. The 5-base rule adds a
fourfold multiplier to A/T→G/C changes when the outer two bases of the centered
five-base window match. A separately normalized 7-base positive control adds a
twelvefold multiplier to C→A/G→T when the outer bases of the seven-base window
match. All rules are reverse-complement symmetric and normalized over uniform
7-mer opportunities to the supplied mean rate.

Calibration sequences have uniform bases. In shifted targets, positions four
bases apart match with probability 0.8 instead of 0.25. Four independent stationary
Markov chains, interleaved by position modulo four, implement this construction.
Every contiguous four-base window is still uniform in expectation. Thus shorter
composition cannot explain the planted fifth-base effect. Finite empirical
trinucleotide distributions are not exactly identical; their total variation
distance is recorded for each replicate.

The six scenarios are constant, 3-base, 5-base with unchanged composition, 5-base
with shifted composition, a true 7-base control, and an optional empirical-shaped
5-base kernel. The latter uniformly marginalizes the outermost pair of bases from
the Carlson table, then renormalizes to the same arbitrary mock mean. It is not
the original paper's human opportunity-weighted 5-mer estimate.

The reference estimator globally aligns each calibration pair with Edlib NW,
extracts mismatches, masks sites near inferred gaps, and counts opportunities in
the parent. Both variant counts and opportunities are required. Reverse
complements are pooled. Contig indices divisible by five form a validation set;
all other contigs fit the candidate models. Smoothing is chosen only using these
calibration data. Final models use all calibration contigs. Target children and
simulator probabilities never enter fitting or tuning.

The `1mer`, `3mer`, `5mer`, and `7mer` models use hierarchical shrinkage toward
the next smaller odd window. The `7mer_flat` comparator instead shrinks directly
toward the one-base probabilities, so it must learn its longer-range effects from
seven-base cells. Both seven-base models choose their shrinkage from the same
validation grid (10, 100, 1,000, 10,000, 100,000, or 10^12 effective opportunities).
Thus the standalone seven-base comparison is regularized rather than an
intentionally unsmoothed straw man. The hierarchical model is the stronger
comparator for whether extra flanks add information beyond five bases.

These are error-free collinear assemblies. This pilot does not validate general
assembly comparison with rearrangements, uncertain homology, or assembly error.
Edlib may choose a gapped optimum even when the generator made only substitutions;
the resulting mask can remove real nearby mutations. Per-run precision and recall
against generator events are recorded as evaluator diagnostics.

## Scores

The primary score is expected full four-allele excess log loss:

`sum_i KL(p_true_i || p_estimated_i) / sum_i P_true(mutation at i)`.

Its unit is nats per expected mutation; lower is better and zero is the oracle.
The normalization improves readability without changing model rankings within a
task. Expected scoring uses evaluator-only truth and reduces noise from rare
realized events. We also compute observed full log loss on the withheld child,
conditional alternate-base log loss, conditional oracle KL, Brier excess, and
relative error in total expected mutation burden. The full log score evaluates
both mutation incidence and alternate-allele identity. Conditional scoring alone
does not test the mutation rate. Raw base accuracy is unsuitable because an
unchanged-base guess is correct at almost every coordinate.

GC-changing prediction is summarized by log error in predicted A/T→G/C versus
G/C→A/T expected counts. This is a mutation-spectrum diagnostic, not an inferred
biological gene-conversion event rate. Five-base context can improve it, but the
primary task remains offspring probability prediction.

Paired comparisons use the same 12 seeds, with 10,000 seed-level bootstrap samples
and fixed bootstrap seed 20261112. These intervals describe simulation-repeat
variability, not uncertainty in every individual predicted probability.

## Run in separate phases

From this directory, use Python 3.11 and `uv`. Dependencies are pinned in `uv.lock`.
The alignment dependency has Linux wheels; the local Windows run used WSL Ubuntu
22.04 rather than requiring a Windows C++ build toolchain.

```bash
uv sync --locked --group dev
uv run --locked pytest -q
uv run --locked benchmark.py prepare --workers 4 --replicates 12
uv run --locked benchmark.py fit --workers 4 --replicates 12
uv run --locked benchmark.py evaluate --workers 4 --replicates 12
uv run --locked benchmark.py report --output results
```

The defaults run the five wholly synthetic scenarios. To reproduce the additional
empirical-shaped control, pass this same option to prepare, fit, and evaluate:

```bash
--scenarios constant three five_iid five_shift seven_shift empirical_five
```

Defaults are four worker processes and a logged 26-GB memory budget; numerical
libraries use one thread per worker to avoid oversubscription. The memory option
records the allocation, not an OS-enforced limit. Up to 24 local workers are
accepted. The implementation processes short contigs and uses far less than the
default memory allowance. Logs contain parameters, seeds, paths, environment,
source hash, and Git revision. Sample manifests and output checksums control
reuse; a Git revision change alone does not trigger regeneration.

`prepare` verifies or regenerates an instance if its declared inputs or file
checksums differ. `fit` and `evaluate` reject corrupted upstream files. Semantic
changes to the generation or fitting algorithm must explicitly bump the relevant
schema/version; cache keys deliberately do not depend on the raw Git commit.

## Solve and grade one instance

The reference solver requires only public files:

```bash
uv run --locked solve.py work/main/five_shift/seed_000/public --output work/predictions.tsv
```

Its automatic choice uses calibration validation log loss and a paired
contig-level one-standard-error rule: select the smallest window whose disadvantage
relative to the best candidate is within one standard error. This favors simpler
models when extra flanks provide no reliable improvement. To export a comparator, use
`--method 1mer`, `3mer`, `5mer`, `7mer`, or `7mer_flat`. The task prompt does not
tell the task solver to use this reference implementation.

The evaluator runs with a separate private directory:

```bash
uv run --locked grade.py work/main/five_shift/seed_000 work/predictions.tsv --output work/grade.json
```

The grader checks every coordinate and the probability schema, rejects missing
or reordered queries and malformed probabilities, and reports expected and
realized scores. Exact zero probabilities have infinite oracle loss here because
every allele has positive true probability; they produce an explicit failed-score
message rather than being silently clipped.

Figures are saved as PNG and PDF with editable vector text, using the Agg backend.
They are not displayed inline. Tests cover alignment and gap masking, exact kernel
dependencies, the sequence-composition control, probability scoring, cache
corruption, and public-only solving followed by independent grading.
