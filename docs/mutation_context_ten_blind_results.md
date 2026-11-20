# Ten-base dependence under sparse calibration

**The 4M blind agent did not discover the dependency; the independent 64M blind
agent recovered the exact −5..+4 window.** Additional calibration made the signal
much clearer in both the blind comparison and a fixed-estimator learning curve.
This is one attempt per size on one fixed kernel, not a sharp sample-size threshold
or an estimate of discovery success rates.

## Fixed design and paired experiments

The planted mutation spectrum is an arbitrary forward-oriented ten-base lookup
table at offsets **−5 through +4**, including the target base at zero. It has
1,048,576 contexts and 3,145,728 applicable substitution channels. Each applicable
channel starts with an independent lognormal draw with log SD 1.25 and log mean
−1.25²/2; one common normalization sets the mean mutation probability to 0.005
under uniform DNA. The unchanged-base probability completes each row. No equality,
CpG, transition, or reverse-complement constraints are imposed. The largest context
mutation probability is 0.405060; no rates were clipped or resampled.

The primary trial retains **4 million calibration sites** and **1 million target
sites**, with seed 0 and 50,000-base circular scaffolds. Parental sequences, public
wording, metadata, and shuffled query coordinates match the earlier arbitrary
eight-base experiment. Every target position is queried once, including circular
boundaries. The public files contain no context-length or orientation clues.

A separate fresh agent received **64 million calibration sites** under the same
kernel. This second trial was launched before the first agent completed, after
author-side comparisons suggested a much stronger signal with additional data.
It was not a revision of the first agent's submission. The larger calibration set
contains the exact original four-million-site prefix; the target parent, target
offspring, rate table, and query order are identical. The
[input audit](../experiments/mutation_context/blind_runs/case_005/nested_input_audit.json)
verifies byte-identical shared files and both calibration FASTA prefixes.

Each solver started with no conversation inheritance and received only the same
generic task instructions and six public files in a neutral directory. Neither
solver was given scientific hints, model suggestions, reference results, or hidden
scores. They were instructed not to access other directories or communicate with
other agents. This is instruction-based isolation on a shared filesystem, not an
OS sandbox; no OS access trace was collected. Both inherited the current parent
model without an override. Their model/backend should not be assumed identical
to runs conducted on earlier dates.

The protocols distinguish identifying the exact asymmetric span from merely
detecting some longer-context signal:
[4M protocol](../experiments/mutation_context/blind_runs/case_005/protocol.json) and
[64M protocol](../experiments/mutation_context/blind_runs/case_006/protocol.json).

## How sparse are the data?

| Calibration bases | Observed changes | Sites per 10-mer, mean | Mutations per 10-mer, mean | Channels with zero events |
| ---: | ---: | ---: | ---: | ---: |
| 4,000,000 | 19,832 | 3.815 | 0.0189 | 99.3786% |
| 16,000,000 | 79,846 | 15.259 | 0.0761 | 97.5960% |
| 64,000,000 | 320,288 | 61.035 | 0.3055 | 91.4792% |

At 4M, 2.21% of possible contexts never occur even as unchanged opportunities;
most contexts occur, but almost all individual replacement channels have no
events. All contexts occur by 16M in this realized sample. Consequently,
discovering the window and estimating its individual rates are separate goals.

The 64M sample restores the average opportunities per context of the earlier
eight-base task with 4M calibration bases. Matching the earlier six-base task's
opportunities per context would take approximately **1.024 billion bases**. This
last figure is arithmetic scaling, not an experiment performed here or a guarantee
of exact rate recovery.

## Four-million-base blind result: no window discovery

The primary blind agent selected a homogeneous process: mutation probability
0.004958, probability 0.995042 for the parental allele, and 0.0016526667 for each
alternative. It therefore **did not discover the planted ten-base dependency**.

Its initial marginal tests looked uninformative. It evaluated rate models through
seven-base contexts, including forward and reverse-complement representations,
and full-spectrum models through reverse-complement seven-base contexts. Five-fold
whole-scaffold and contiguous-region validation gave small gains for strongly
shrunk seven-base models. The best rate model gained 6.61 nats in scaffold
validation, with a paired standard error of 3.77 nats; the agent judged this weak
after screening many alternatives and chose a conservative constant model.

That decision is documented in the unchanged
[agent report](../experiments/mutation_context/blind_runs/case_005/agent_submission/REPORT.md)
and [decision log](../experiments/mutation_context/blind_runs/case_005/agent_submission/decision_log.md).
Crucially, it did **not** search the true ten-base span, and reverse-complement
pooling is not a symmetry of this generator. This run cannot establish an
information-theoretic limit. Sparsity weakened the short-context evidence that
might motivate further search, while the agent's search scope and pooling choices
also contributed. The correct-window author reference still extracts a small
predictive benefit from these data.

The frozen submission's hidden excess KL was **0.785108 nats per expected
mutation**. It predicted **4,958.0 changes**, versus a true expectation of
**4,996.2174**, an error of **−0.7649%**. The actual offspring contained 4,977
changes. Accurate aggregate burden therefore did not indicate recovery of the
context-specific spectrum. Full output integrity and all six input hashes passed.
The [hidden grade](../experiments/mutation_context/blind_runs/case_005/grade.json)
was computed only after submission was final, with no feedback or revision.

## Sixty-four-million-base blind result: exact window recovered

The second agent initially found flat marginal substitution patterns, but its
nine-base tables gave substantial validation gains. It explored longer and
asymmetric windows, found reverse-complement pooling harmful, and identified the
two best shifted nine-base windows as −5..+3 and −4..+4. Combining them into
**−5..+4** was substantially better, while extra flanking positions degraded its
independently estimated tables. It explicitly selected this oriented ten-base
window, satisfying the full-discovery criterion.

Its final estimator uses empirical Bayes regularization of alternate-allele
probabilities: a shared discrete prior on 100 probability-grid points, fitted
through grouped binomial likelihoods, followed by posterior mean predictions for
each context and alternate allele. This is an approximate composite model for
multinomial outcomes; the three alternate fits are not literally independent
events. All emitted probabilities were positive and normalized. The EM procedure
hit its 3,000-iteration cap, and the agent did not claim exact convergence.

The selected method gained 10,202.26 nats over pooled uniform substitutions on
12.8M calibration validation sites, with scaffold-clustered standard error 184.36
nats. The agent inspected its initial assessment split and then continued model
development, so these are **adaptive validation results**, not an untouched
confirmatory test. Its
[report](../experiments/mutation_context/blind_runs/case_006/agent_submission/REPORT.md)
and [decision log](../experiments/mutation_context/blind_runs/case_006/agent_submission/decision_log.md)
explicitly acknowledge this. The hidden target remained independent throughout.

After freezing the submission, its hidden excess KL was **0.566496 nats per
expected mutation**, versus 0.785108 for the smaller-data blind submission, a
**27.8% reduction**. Conditional replacement-allele KL was 0.301009. It predicted
**5,005.3708 changes**, versus the shared true expectation of **4,996.2174**,
an error of **+0.1832%**. Realized target mean allele NLL was 0.035725581, versus
0.036825097 for the smaller-data model. All output and input-integrity checks
passed; see the [hidden grade](../experiments/mutation_context/blind_runs/case_006/grade.json).

The blind-agent comparison includes different search paths and estimators as well
as different sample sizes. It supports the practical value of more data, but does
not isolate a purely causal effect of sample size on agent reasoning. The fixed
reference comparisons below isolate the data increase for one estimation procedure.

## Fixed-estimator learning curve

Author-side comparisons used nested prefixes of the same 64M calibration dataset.
Whole scaffolds with index modulo five zero formed validation; the others formed
training. A grid of 21 strengths from 1 to 100,000, plus 10¹², chose smoothing using
calibration data only. Each table could borrow from the best already-fitted
contained smaller window. Models were refitted on each full prefix before hidden
scoring. The target genome and true probabilities were held fixed across sizes.

| Calibration bases | Parent-base model | Centered 9-mer, −4..+4 | Correct 10-mer, −5..+4 | Opposite 10-mer, −4..+5 | Hierarchical 11-mer |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4,000,000 | 0.785498 | 0.778384 | 0.768926 | 0.778398 | 0.768676 |
| 16,000,000 | 0.785159 | 0.764776 | 0.735819 | 0.764807 | 0.735819 |
| 64,000,000 | 0.785079 | 0.723991 | 0.642819 | 0.723991 | 0.642785 |

Scores are expected KL from the true four-base probabilities, divided by expected
mutations: nats per expected mutation, lower is better, zero is oracle performance.
The correct ten-base reference reduces loss relative to the parent-base model by
2.1%, 6.3%, and 18.1% at these three calibration sizes. At 4M its improvement over
the centered nine-base reference is only about 1.2%.

The 11-mer can shrink toward its correct ten-base parent and therefore need not
perform worse. Its small validation gains at 4M and 64M do not establish that the
extra position is biologically relevant. The best raw validation score selected
the 11-mer at those sizes, and the ten-base model at 16M. These author comparisons
include the true window among candidates and are not themselves blind discovery
trials. They demonstrate increased estimability with more data, not a success
probability for arbitrary agents.

For the same correctly oriented reference estimator, going from 4M to 64M reduces
excess KL by **16.4%**. Thus the data increase improves probability estimation even
without changing the chosen window or the fitting procedure. The 16M point was
evaluated only with reference models; no blind agent was run at that size.

For the next benchmark version, 64M calibration bases are supported by an observed
blind discovery success and a stronger predictive signal. Sixteen million is a
reasonable intermediate size to examine with replicated blind trials. Neither
64M nor successful window identification guarantees accurate recovery of all
3.1 million arbitrary channel rates; 91.5% still have zero observed events at 64M.

Full per-size scores, support counts, smoothing choices, timings, source hashes,
and cache manifests are in the
[learning-curve directory](../experiments/mutation_context/blind_runs/case_005/learning_curve).

## Verification and reproduction

Five focused tests passed: the new ten-base tests and adjacent eight-base
regressions. They cover kernel validity and dependence on every position,
determinism, complete shuffled queries, independent circular encoding,
zero-loss oracle grading, independent count aggregation for each reference
window, target-score agreement, corruption detection, nested data prefixes, and
worker-count-independent generation. No plots were generated or displayed.

The public archives contain only the six solver inputs:
[4M task](../experiments/mutation_context/releases/task_005.zip) and
[64M task](../experiments/mutation_context/releases/task_006.zip). All author-side
code, protocols, reports, private kernels and saved solutions must remain outside
a fresh blind solver's environment. Reproduction commands are in the experiment
[README](../experiments/mutation_context/README.md).

Both submissions were replayed with Python 3.13.15 and their own pinned packages;
their prediction TSVs were byte-identical to the frozen originals. See the
[4M replay check](../experiments/mutation_context/blind_runs/case_005/replay_check.json)
and [64M replay check](../experiments/mutation_context/blind_runs/case_006/replay_check.json).
Replay regenerates the selected estimator, not the agent's reasoning. The larger
agent's two large, regenerable count caches remain local rather than in Git;
their hashes and reconstruction scripts are recorded in
[omitted artifacts](../experiments/mutation_context/blind_runs/case_006/omitted_artifacts.json).
The required eleven-base count cache was rebuilt from public calibration for
replay. Full million-row predictions remain in `D:/agent_evaluations/case_005`
and `D:/agent_evaluations/case_006`; their hashes and short previews are committed.
