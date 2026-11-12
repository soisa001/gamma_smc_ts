# Blind agent discovery of five-base and asymmetric six-base processes

Both fresh agents discovered the intended dependency from the public task data.
The first recovered the centered five-base interaction; the second recovered the
asymmetric six-base span, including both planted flank interactions. Neither
received the generator, reference solver, original conversation, window-size
labels, or hidden-score feedback before submitting final predictions.

This is two descriptive trials, one agent on one seed for each process. It is not
an estimated success rate across agents, models, or datasets.

## What the agents discovered

The first agent began with base-specific and CpG rates. It found unexplained
heterogeneity in five-base contexts, inspected high-residual motifs, and identified
`parent[i-2] == parent[i+2]` as the relevant interaction. It correctly localized
the increase to A/T→C/G substitutions and detected the target's shift from roughly
25% to 80% matching outer bases. Its final model has five fitted rate parameters;
it represents a five-base dependency without estimating an unrestricted five-base
table. Its own five-fold contig validation favored this compact model.

The second agent initially found little improvement from centered five-base
tables but a real improvement from seven-base tables. It scanned flank
interactions and discovered `parent[i-3] == parent[i+2]`, followed by the smaller
`parent[i-2] == parent[i-1]` effect on A↔T. It rediscovered both pairs separately
within each training fold. It also recognized that reversing the coordinates
before evaluating these interactions would be wrong. Its nine-rate model uses
the exact asymmetric span `[-3, -2, -1, 0, +1, +2]`; it does not require +3.

These discoveries appear in the agents' chronological logs and final reports,
written before hidden grading:

- [Five-base decision log](../experiments/mutation_context/blind_runs/case_001/agent_submission/decision_log.md) and [report](../experiments/mutation_context/blind_runs/case_001/agent_submission/REPORT.md).
- [Six-base decision log](../experiments/mutation_context/blind_runs/case_002/agent_submission/decision_log.md) and [report](../experiments/mutation_context/blind_runs/case_002/agent_submission/REPORT.md).

## Hidden evaluation

Each submission contains all 999,880 requested sites in order with valid
four-allele probabilities. Input hashes were verified, source/report files were
copied without modification, and prediction hashes were frozen before the
independent grader accessed private outcomes or simulator probabilities. No
submission was revised in response to a hidden score.

| Trial | True expected changes | Agent prediction | Burden error | Excess log loss |
| --- | ---: | ---: | ---: | ---: |
| Centered 5-base truth | 6,941.157 | 6,897.168 | −0.634% | 0.00073138 |
| Asymmetric 6-base truth | 6,914.795 | 6,974.538 | +0.864% | 0.00026326 |

Excess log loss is the sum of KL divergences between true and submitted allele
probabilities divided by the true expected mutation count, in nats per expected
mutation. Zero is oracle performance. Burden errors compare with the expectation,
not the random realized count. Conditional alternate-allele KL was 0.00054997 and
0.00019203, respectively.

On the identical seed-0 five-base dataset, the previous 3-mer table scored
0.153013 and the regularized 5-mer table scored 0.049967. The blind agent's
0.000731 score is much better because it shares information through the inferred
low-dimensional rule instead of estimating many sparse rates.

On the six-base dataset, author-side directional tables scored 0.152949 for a
centered 5-mer, 0.124004 for the correctly offset 6-mer, 0.152981 for the oppositely
offset 6-mer, and 0.119636 for a hierarchical 7-mer. These baselines were not shown
to the agents. They tune smoothing on calibration data and assume no
reverse-complement symmetry. The blind agent's 0.000263 score again reflects its
compact inferred rule. A larger representation that learns equivalent
probabilities is a valid solution; success does not require a particular model name.

- [Five-base grade](../experiments/mutation_context/blind_runs/case_001/grade.json) and [same-dataset references](../experiments/mutation_context/blind_runs/case_001/reference_scores.json).
- [Six-base grade](../experiments/mutation_context/blind_runs/case_002/grade.json) and [directional references](../experiments/mutation_context/blind_runs/case_002/reference_scores.json).

## Blinding protocol and limits

Both agents used `fork_turns="none"`, with no inherited conversation turns and no
model override. Each received a neutral directory containing only six public
files. The prompt supplied the task path, output requirements, numerical seed and
resource defaults, and runtime instructions. It prohibited other projects,
repository history, chats/agents, hidden outcomes, scientific web sources,
generator code, and reference solutions. No follow-up hints or process answers
were provided. The [recorded prompt](../experiments/mutation_context/blind_runs/agent_prompt_template.txt)
is identical apart from the neutral case identifier. Public `TASK.md`, metadata,
and query coordinates are byte-identical between cases; only the FASTAs differ.

The collaboration runtime shares the host filesystem. Isolation was therefore a
fresh conversation plus explicit file-access restrictions, not an OS-enforced
sandbox. Both agents reported compliance; their saved solution code uses the
supplied directory and standard libraries. No independent OS-level file-access
trace was collected. Protocols and assessment criteria were saved before final
submissions: [case 001](../experiments/mutation_context/blind_runs/case_001/protocol.json)
and [case 002](../experiments/mutation_context/blind_runs/case_002/protocol.json).

Both agents used exploratory calibration analysis to invent features and then
reused calibration folds for comparisons. Their reports acknowledge that these
cross-validation scores are not untouched tests of the entire discovery process.
The private target evaluation remains independent and occurred only after final
submission. Given the SNV-only, homologous, error-free assemblies, both agents used
direct ungapped comparisons. This tests inference rather than general alignment
robustness.

## Asymmetric generator and task-design implications

The new task uses 4 Mb of calibration sequence, 1 Mb of target sequence, seed 0,
50-kb contigs, and an arbitrary mean mutation probability of 0.005 under uniform
sequence. It samples one offspring allele from a categorical distribution based
on the fixed parent. Rates include CpG transitions, a fourfold A/T→C/G multiplier
when −3 and +2 match, and a twofold A↔T multiplier when −2 and −1 match. Tests
verify that all six positions from −3 through +2 can matter and +3 cannot.
This is deliberately a directional mock process, without ordinary
reverse-complement symmetry.

Five stationary chains interleaved by position modulo five produce uniform
consecutive five-base words in expectation. Lag-five matching changes from 0.25
in calibration to 0.8 in targets. This changes the six-base interaction while
preserving shorter-window composition in expectation. The public prompt does not
disclose the construction. [task_002.zip](../experiments/mutation_context/releases/task_002.zip)
contains only public inputs; [asymmetric_task.py](../experiments/mutation_context/asymmetric_task.py)
is the author-only generator.

The intended discovery occurred even with an unfamiliar asymmetric offset.
However, low-dimensional equality rules are easier to infer than arbitrary tables
of thousands of unrelated context rates. Increasing interaction complexity is a
more meaningful difficulty adjustment than changing window length alone.

One residual cue should be removed in a stricter next version: the first agent
explicitly noticed that queries start at position 4, suggesting local sequence
context. This did not reveal the five-base rule, but it is an indirect hint. A
query panel or callable mask chosen independently of the generating window would
make the test cleaner. These trial archives remain unchanged for auditability.

Frozen agent source, reports, diagnostics, and small public-calibration count
arrays live under `experiments/mutation_context/blind_runs/case_*/agent_submission`.
Large predictions and per-site caches remain local under
`D:/agent_evaluations/case_001` and `D:/agent_evaluations/case_002`; prediction hashes
are committed. No plots were needed. Three focused generator/baseline tests passed,
and both full submissions passed independent grading and integrity checks.
See the [README](../experiments/mutation_context/README.md#blind-agent-trials) for
reproduction. Replaying saved solver code is not a new blind agent trial.
The first solver was replayed from its frozen source and public archive with
pinned dependencies; its prediction TSV was byte-identical to the submitted file.
The [replay check](../experiments/mutation_context/blind_runs/replay_check.json)
records the hashes.
