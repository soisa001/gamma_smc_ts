# Blind discovery of an arbitrary eight-base mutation spectrum

## Fixed experiment

This trial extends the arbitrary six-base benchmark to a forward-oriented window
at offsets **−4 through +3**, including the mutated base at offset zero. There are
65,536 contexts and 196,608 applicable substitution channels. Every channel starts
with an independent lognormal draw (log SD 1.25, log mean −1.25²/2), followed by one
common normalization to mean mutation probability 0.005 under uniform DNA. The
unchanged-base probability completes each row. No equality, CpG, transition, or
reverse-complement rules are imposed. All eight positions matter.

The experiment retains seed 0, four million calibration positions, one million
target positions, 50,000-base circular scaffolds, and independent uniform parental
DNA. The parental sequences, offspring-sampling random uniforms, query order, task
prose, and public metadata match case 003. Only the generating rate table and its
encoding change. The largest context mutation probability is 0.203033; no rates
were clipped or resampled. This is an arbitrary synthetic spectrum, not a model
of human absolute mutation rates.

Every target position is queried exactly once, including both ends of every
scaffold, in independently shuffled order. Circular sequence semantics are stated
in the public prompt without mentioning context, window lengths, or asymmetry.
The six-file public release contains no simulator, private truth, evaluator, or
reference solution. In particular, `kernel_encoding.json` remains private.

The public archive is
[task_004.zip](../experiments/mutation_context/releases/task_004.zip), SHA-256
`64eb0f9f03da8da71f0e00866dd6645e5f9d739a93c354a038eac974bc4434da`.
The [protocol](../experiments/mutation_context/blind_runs/case_004/protocol.json)
was written before launching the fresh solver. There is one blind attempt, no
generator retuning, and no feedback of hidden scores before final submission.

## Information available in calibration

The calibration has 19,901 observed mutations. Each context has 61.04 site
opportunities and 0.304 mutations on average. There are only 0.101 events per
substitution channel; **179,956 of 196,608 channels (91.53%) have no observed
events**. All contexts nevertheless occur among the calibration opportunities.

For comparison, case 003 had 4.88 mutations per six-base context and 45.26% of
channels without events. The two extra bases multiply the number of contexts by
16 while retaining the same sample size. Roughly 64 million calibration bases
would match the earlier average opportunities per context. This is a support
calculation, not a claim that such a run has been performed or would guarantee
exact rate recovery.

Window identification and rate estimation are therefore distinct outcomes. A
shared pattern in predictive validation can identify the window even when most
individual replacement rates have never produced an observed event.

## What the blind agent discovered

The agent independently identified **the exact −4..+3 eight-base window**. It
started from uniform and source-base models, explored local composition, spatial
patterns and circular offset relationships, then compared contiguous windows of
lengths 3–10. The strongest transferable context signal localized to the correct
eight-base window. Removing any one position substantially reduced validation
performance. This meets the preregistered full-discovery criterion.

Detailed discovery used 40 scaffolds, model selection another 20, and the final
reserved assessment the remaining 20, assigned with seed 0. The initial global
count audit included all scaffolds, but detailed feature development excluded the
reserved assessment partition. The agent selected a logistic-normal prior for
per-context mutation rates and a symmetric Dirichlet prior for replacement-base
probabilities. Both use fitted shared hyperparameters to regularize the sparse
individual contexts. These are useful approximations; they are not the exact
distribution used to generate the arbitrary channel rates.

After locking the model, it fitted 60 scaffolds and evaluated the reserved million
sites. Mean four-allele negative log likelihood fell from 0.036793077 for a uniform
model to 0.035903488, a gain of **889.589 nats**. Of this, 466.622 nats came from
mutation probability and 422.967 from conditional replacement probabilities. A
5,000-resample scaffold bootstrap gave a descriptive 95% interval of
792.936–991.619 nats for the gain. This interval is conditional on the trained
model and does not correct for all exploratory comparisons. The same model
specification was then refitted to all four million calibration positions.

The agent's [report](../experiments/mutation_context/blind_runs/case_004/agent_submission/REPORT.md)
and [chronological log](../experiments/mutation_context/blind_runs/case_004/agent_submission/decision_log.md)
preserve its independent hypotheses and decisions.

## Author-side probability comparisons

Reference tables use an ungapped alignment justified by the public SNV-only
homology guarantee, and circular context extraction at every position. They fit
on calibration scaffolds with index modulo five nonzero and choose smoothing on
the remaining scaffolds. The smoothing grid is 21 logarithmically spaced values
from 1 to 100,000, plus 10¹². Each table can borrow from the best already-fitted
contained smaller window. No strand pooling is imposed. All model choices precede
access to withheld target outcomes and rates.

| Reference estimator | Excess log loss |
| --- | ---: |
| Parent-base only | 0.778248 |
| Centered 5-mer, −2..+2 | 0.771628 |
| Centered 7-mer, −3..+3 | 0.720219 |
| Opposite 8-mer, −3..+4 | 0.720219 |
| Correct 8-mer, −4..+3 | 0.642483 |
| Hierarchical 9-mer, −4..+4 | 0.642463 |
| **Blind agent, discovered −4..+3** | **0.570586** |

Loss is expected KL from the true four-allele probabilities divided by the true
expected number of mutations, in nats per expected mutation; zero is oracle
performance. The correctly oriented eight-base reference improves this loss by
10.8% relative to the centered seven-base reference. The opposite orientation
shrinks almost entirely to its seven-base parent. A larger table can pool toward
the correct eight-base table, so the task does not guarantee that every nine-base
estimator must perform worse. Reference scores were never sent to the blind agent.

The [reference output](../experiments/mutation_context/blind_runs/case_004/reference_scores.json)
also records binary mutation loss, conditional replacement loss, burden, and
support counts. Aggregate burden alone is insufficient: the parent-base model
predicts 4,976.45 changes versus a true expectation of 5,005.44 (−0.58%), despite
its poor site-specific probabilities. The realized target contains 5,037 changes.

## Frozen blind submission and hidden grading

The final million-row submission was frozen and hashed before private grading.
All input hashes, query coordinates and order, row counts, and probability checks
passed. No scientific follow-up hints or hidden feedback were supplied to the
agent, and no predictions were revised after scoring.

The agent achieved excess KL **0.570586 nats per expected mutation**, a 20.8%
reduction relative to the centered seven-base reference and a 26.1% reduction
relative to the centered five-base reference. It also improved 11.2% over the
author's correctly oriented eight-base table. That last comparison includes
different regularization choices and calibration splits; it is not solely a
window-length effect. Conditional replacement-allele KL was 0.303903, with the
remaining 0.266682 attributable to mutation/no-mutation probabilities.

The predicted burden was **4,974.7287**, compared with the true expectation
**5,005.4441**, an error of **−0.6136%**. Its burden estimate is slightly less
accurate than the parent-base reference's, despite substantially better sitewise
probabilities. This reinforces why the full probabilistic score is the primary
metric. Realized target mean allele NLL was 0.036233348.

The residual excess loss is much larger than the prior arbitrary six-base agent's
0.171698. The same four-million-base sample supports identifying the longer window,
but does not support accurate recovery of every arbitrary channel rate. This is
one run per kernel, not a replicated effect-size or success-rate estimate.

The [hidden grade](../experiments/mutation_context/blind_runs/case_004/grade.json),
[submission manifest](../experiments/mutation_context/blind_runs/case_004/submission_manifest.json),
and unchanged agent code/reports are committed. The full 85 MB prediction TSV
remains at `D:/agent_evaluations/case_004/predictions.tsv`; a short preview is
committed. Its SHA-256 is
`a69acbae4504c51c47111368c63c79e5bf05ead53cf7215e89173e600eacb9f1`.

## Blinding and reproducibility

The solver started with `fork_turns="none"` in `D:/agent_evaluations/case_004`,
using the same generic prompt as prior blind runs. Only the six public files and
standard runtime/package facilities were allowed. This is a fresh conversation
with instruction-based restrictions on a shared filesystem, not an OS sandbox.
No OS-level access trace was collected, and a single run does not estimate an
agent success rate.

The new private encoding allows the grader to handle the eight-base kernel while
retaining the legacy seven-base default for earlier tasks. Fourteen focused tests
passed, including legacy cases, deterministic kernel generation, dependence on
every planted position, shuffled complete queries, independent modular encoding
at circular boundaries, oracle zero-loss grading, reference-score agreement on
the target mutation expectation, and corruption detection. No plots were made or
displayed. See the [README](../experiments/mutation_context/README.md) for commands.
The frozen estimator was replayed from the public archive with Python 3.13.15
and the agent's pinned dependencies. It regenerated a byte-identical prediction
TSV; the [replay check](../experiments/mutation_context/blind_runs/case_004/replay_check.json)
also verified archived source/report hashes and all public input hashes. This
reproduces the fitted method, not the agent's stochastic reasoning process.
