# Blind discovery of an arbitrary asymmetric six-base spectrum

The fresh agent discovered the exact asymmetric window, **−3 through +2**, even
with arbitrary context-specific rates and no query-boundary exclusion. It tested
all 28 contiguous windows of lengths 1–7 that contain the focal base and selected
the correct six-base window. Its final model is a regularized 4,096-context table,
not a compact equality or CpG rule.

This is one fixed-seed task and one fresh agent run. No generator retuning,
score-based retry, or scientific hints were provided during the run.

## Harder generator and boundary removal

There are 4,096 ordered contexts with the focal base fourth, covering offsets
`[-3, -2, -1, 0, +1, +2]`. Each of the 12,288 possible substitution channels gets
an independent lognormal draw, with log standard deviation 1.25 and log mean
−1.25²/2, before one common normalization to mean mutation probability 0.005 under
uniform sequence. The unchanged probability is one minus the three substitution
probabilities. There is no planted CpG effect, equality rule, transition bias, or
reverse-complement symmetry. Tests verify that every required offset can affect
the spectrum and +3 cannot.

The calibration remains 4 Mb, with 1 Mb of target sequence and 50-kb scaffolds.
Both sets are independent uniform DNA; no motif-composition shift was planted.
All numerical randomness uses seed 0 with separate SeedSequence streams for the
kernel, each contig, and query order. Generation defaults to four workers and
records a 26-GB resource budget.

All scaffolds are circular, and the public prompt explicitly explains coordinate
wraparound. **Every target base is queried**, including position 1 and the final
base of every scaffold. The one million rows are independently shuffled. Thus
there is no excluded boundary margin from which to guess a context radius.
The public prompt does not mention context, window length, or asymmetry.

The [public task archive](../experiments/mutation_context/releases/task_003.zip),
[generator](../experiments/mutation_context/arbitrary_task.py), and
[generation provenance](../experiments/mutation_context/blind_runs/case_003/generation_provenance.json)
are saved separately. Only the archive's six files were supplied to the agent.

## What the agent independently discovered

Its early analysis found that highly mutable seven-base words often shared their
first six letters while varying in the final letter. It treated this as a
hypothesis and systematically compared asymmetric windows. The −3..+2 window won
on development scaffolds; extending it on either side worsened its tested table
models. It then compared several regularization approaches rather than searching
for a low-dimensional motif rule.

The selected estimator uses an empirical-Bayes logit-normal prior on the total
mutation probability and an estimated symmetric Dirichlet prior on the three
replacement-allele probabilities. The agent did not recover or claim the exact
generating prior: its model factorization differs from the independent lognormal
channel draws. The recovered window and improved probabilities are the relevant
successes.

It used a seed-0 split of 48 training, 16 development, and 16 reserved audit
scaffolds. After freezing the model family, it evaluated 800,000 audit sites,
obtaining allele log loss 0.03438183 versus 0.03688424 for its constant baseline.
It then refitted on all calibration data for submission. Some initial descriptive
summaries used all calibration data; its report acknowledges that the audit was
not completely untouched by exploration. The private target grader remained
unseen throughout model construction.

- [Chronological decision log](../experiments/mutation_context/blind_runs/case_003/agent_submission/decision_log.md).
- [Agent report and model comparisons](../experiments/mutation_context/blind_runs/case_003/agent_submission/REPORT.md).
- [All tested window results](../experiments/mutation_context/blind_runs/case_003/agent_submission/window_validation.json).

## Frozen submission and hidden scores

All 1,000,000 submitted rows passed independent coordinate, ordering, probability,
and input-integrity checks. Predictions were frozen and hashed before accessing
private target outcomes or the true rate table. There was no post-score revision.

| Estimator | Excess full-allele log loss |
| --- | ---: |
| Parent-base-only reference | 0.753061 |
| Centered 5-mer reference, −2..+2 | 0.522728 |
| Oppositely offset 6-mer reference, −2..+3 | 0.522728 |
| Correctly offset 6-mer reference, −3..+2 | 0.240108 |
| Hierarchical 7-mer reference | 0.237611 |
| **Blind agent, discovered −3..+2** | **0.171698** |

Loss is expected KL from the true four-allele probabilities, divided by expected
mutations, in nats per expected mutation. Lower is better; zero is oracle
performance. The agent reduced this loss by about 67% versus the centered 5-mer
reference and 28% versus the correctly offset reference table. The author-side
tables use a separate fixed calibration validation split and a coarser smoothing
grid, so the extra improvement reflects estimator choice as well as correct
window identification. No reference scores were shown to the agent.

The agent predicted **5,000.441** mutations against a true expectation of
**4,996.795**, a **+0.07296%** error. The realized offspring happened to contain
4,936 changes; that random count is not the ground truth for expected burden.
Conditional alternate-allele KL was 0.103501.

Aggregate burden is deliberately not a sufficient success metric here. With
matched training and target composition, even the parent-base-only reference has
only about 0.0403% burden error while its sitewise probabilities are much worse.
The full probability score detects the missed spectrum.

## Why rate estimation is harder even though the window was found

The calibration contains 19,987 mutations: 4.88 per six-base context and 1.63 per
substitution channel on average. **5,561 of 12,288 channels (45.26%) have no observed
events.** There is ample evidence for the correct window but limited information
for many individual rates. Regularization is essential, and exact table recovery
is not supported by these data.

The remaining excess loss, 0.1717, is much larger than the roughly 0.00026 obtained
in the earlier simple-rule six-base trial. This is a descriptive comparison across
different generators, not a controlled effect estimate: topology and target
composition also changed. The practical result is that the agent still locates
the asymmetric window, while arbitrary-spectrum estimation remains substantially
imperfect. The prior task's compact-rule shortcut is gone.

## Blinding and reproducibility

The agent started with `fork_turns="none"` in a new neutral directory. It received
the same generic [agent prompt](../experiments/mutation_context/blind_runs/agent_prompt_template.txt)
as earlier trials, substituting `case_003`, and no follow-up hints. The
[preregistered protocol](../experiments/mutation_context/blind_runs/case_003/protocol.json)
distinguished full asymmetric-window discovery from merely finding that larger
centered windows help. The full-discovery criterion was met.

As before, the runtime shares the host filesystem: access restriction was
instruction-based, not an OS sandbox. The agent reported compliance; its saved
code uses its own task directory and standard packages. No OS-level access trace
was collected. This one successful run does not estimate an agent success rate.

The [grade](../experiments/mutation_context/blind_runs/case_003/grade.json),
[submission hashes](../experiments/mutation_context/blind_runs/case_003/submission_manifest.json),
[reference scores](../experiments/mutation_context/blind_runs/case_003/reference_scores.json),
and [support diagnostics](../experiments/mutation_context/blind_runs/case_003/support_diagnostics.json)
are committed with the unchanged agent source and reports. The full prediction
TSV remains at `D:/agent_evaluations/case_003/predictions.tsv`.

The circular grader passed oracle tests at all positions of small test scaffolds,
including wraparound sites, and regression checks for earlier linear tasks. Twelve
focused tests passed initially; the five affected generator/baseline tests passed
again after the circular baseline was finalized. The agent independently checked
19,898 sampled and boundary context encodings. No plots were generated or displayed.
See the [README](../experiments/mutation_context/README.md#arbitrary-six-base-trial)
for generation and solver-replay commands.
The frozen solver was replayed with its pinned dependencies and produced a
byte-identical prediction TSV; hashes are recorded in the
[replay check](../experiments/mutation_context/blind_runs/case_003/replay_check.json).
