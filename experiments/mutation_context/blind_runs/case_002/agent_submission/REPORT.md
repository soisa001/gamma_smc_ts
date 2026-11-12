# Offspring allele prediction report

## Deliverables and main result

`predictions.tsv` contains **999,880 rows**, in exactly the supplied query order, with columns `contig`, `position`, `p_A`, `p_C`, `p_G`, `p_T`. The expected number of changed alleles is **6,974.537599**, or **0.697537%** of queried sites. `summary.json` records the numeric expectation, fitted parameters, verification results, input/output checksums, runtime configuration, and uncertainty estimates.

These are marginal probabilities for one transmitted haploid offspring. The parent allele remains the most probable outcome at every site; predicting every site deterministically as unchanged would discard the calibrated uncertainty and give incorrect expected mutation burden. The expectation is not a claim that precisely this many changes will occur.

## Scope, data integrity, and assumptions

Only the six initially supplied files in this case directory were used. No external references, project directories, repository history, hidden outcomes, or other agents were consulted. The original six files were not modified; their SHA-256 hashes were checked again after prediction generation. All sequences contained only A/C/G/T. Calibration parent and offspring identifiers and lengths matched exactly: 80 homologous pairs, each 50,000 bases, for 4 million aligned sites. The new parent contains 20 contigs of 50,000 bases. All queries are valid one-based coordinates, positions 4–49,997 on each new contig, totaling 999,880 sites.

The task specifies error-free assemblies, one transmission, single-nucleotide changes only, and a shared process between calibration and prediction. I therefore used direct aligned parent/offspring substitutions as labels, without alignment ambiguity, segregation modeling, or assembly-error correction. Calibration fitting and comparable validation use the 3,999,520 sites with complete three-base flanks, excluding only 480 edge sites. There are 20,189 changes at these eligible sites.

## Hypotheses and alternatives tested

1. **Base-only substitution frequencies.** This is a useful baseline but cannot explain the marked elevation of C-to-T when the right parental neighbor is G, and its G-to-A counterpart when the left neighbor is C.
2. **CpG-conditioned rates.** This substantially improves held-out log likelihood. Complementary substitutions had similar estimates, motivating pooling while testing strand-specific alternatives.
3. **Full trinucleotide, pentanucleotide, and heptanucleotide tables.** Both oriented and fully reverse-complement-pooled versions were tested. Tables were shrunk toward lower-dimensional models using prior strengths 100, 1,000, 10,000, 100,000, and 1,000,000. Sparse contexts received positive probabilities instead of zero estimates. Neither trinucleotide nor pentanucleotide detail provided convincing improvement. The best heptanucleotide table improved over CpG alone by 56.69 total NLL, with paired contig standard error 6.37. This led to investigating long-range interactions instead of dismissing sparse extreme counts.
4. **Compact flank relations.** Residual aggregation over individual flank positions and pairs/triples identified a strong relation between positions -3 and +2. Their equality increases A/T-to-C/G substitutions approximately fourfold. Both hot and cold contexts were observed sufficiently often to estimate rates directly. Equality of particular nucleotides, full pair identities, separate strand rates, and adding equality effects to all C/G outcomes were compared. Their additional parameters did not produce persuasive improvements over a compact pooled model. Equality identity improved only 1.35 NLL over the compact model, with standard error 5.50.
5. **Substitution-specific residuals.** After the first compact model, equality of positions -2 and -1 remained strongly associated specifically with A-to-T/T-to-A. These rates approximately double. Adding this rule improved held-out NLL by 105.76, with paired contig standard error 14.61. The largest residual statistic across 90 remaining equality/substitution comparisons was 5.98, insufficient after accounting for the number of comparisons. I stopped adding context rules.
6. **Contig variation.** Out-of-fold expected counts were compared with observed counts on each contig. The approximate dispersion statistic was 98.91 on 79 degrees of freedom (nominal p=0.064). This does not supply convincing evidence for a new-contig latent rate adjustment. No contig identifier was used as a predictive feature.

## Validation design and evidence

A NumPy generator seeded with 0 permuted the 80 calibration contigs into five groups of 16. Every validation fold held out entire contigs; fitting counts and all probability estimates for that fold excluded those contigs. This avoids neighboring-context leakage caused by randomly splitting adjacent sites. Every model was scored on the same 3,999,520 eligible sites. The primary proper scoring rule was multiclass negative log likelihood (NLL; smaller is better), with empirical calibration and multiclass Brier score as supporting checks. Every candidate includes the unchanged parental allele.

| Model | Total held-out NLL | NLL per site |
|---|---:|---:|
| Base only, complement pooled | 145,355.326 | 0.036343193 |
| CpG only, complement pooled | 142,637.781 | 0.035663725 |
| Oriented heptanucleotide, shrinkage 10,000 | 142,581.093 | 0.035649551 |
| CpG plus oriented long-range equality, compact model | 140,740.415 | 0.035189326 |
| Final three-rule model | 140,634.659 | 0.035162884 |

The final multiclass Brier sum was 40153.596242, equivalent to a mean of 0.010039604. Its out-of-fold expected count was 20,194.388, versus 20,189 observed. These totals are a calibration check, not a separate tuning objective.

Crucially, this was exploratory model discovery. The initial residual inspection used the whole calibration set; consequently the quoted cross-validation is **not a pristine, untouched confirmation set**, and small differences between many candidates could be optimistic. As a reproducibility check, the (-3,+2) equality was reselected from all 15 flank pairs using only each training fold, with statistics 2,868–2,936 and the next candidate at most 45.4. After fitting that rule, each training fold independently selected the (-2,-1) A/T-interchange relation as the strongest substitution-specific residual (statistics 176–206). The large, repeated improvements support these effects; they do not remove all adaptive-selection uncertainty. The withheld new offspring was never available and no claim about its realized accuracy is possible.

Calibration by parental base class and trigger state is shown below. Groups 0–3 are A/T centers, ordered by the (-3,+2) equality flag and then the (-2,-1) equality flag. Groups 4–7 are C/G centers, ordered by the CpG flag and then the (-2,-1) equality flag. The latter equality does not affect final C/G probabilities and is retained here as a diagnostic split.

| Group | Sites | Observed changes | Out-of-fold expected changes |
|---|---:|---:|---:|
| 0 | 1,122,074 | 3,531 | 3524.84 |
| 1 | 374,355 | 1,454 | 1471.59 |
| 2 | 375,632 | 3,808 | 3804.03 |
| 3 | 125,424 | 1,373 | 1369.22 |
| 4 | 1,124,227 | 3,381 | 3429.74 |
| 5 | 375,530 | 1,155 | 1145.66 |
| 6 | 376,354 | 4,097 | 4083.10 |
| 7 | 125,924 | 1,390 | 1366.20 |

## Final probability model

For a queried parent base at position `i`, compute all conditions from the **parent** sequence:

- `E1 = (parent[i-3] == parent[i+2])` in the supplied assembly orientation.
- `E2 = (parent[i-2] == parent[i-1])` in the supplied assembly orientation.
- `CpG = (parent[i] == C and parent[i+1] == G) or (parent[i] == G and parent[i-1] == C)`.

For A/T centers, the transition A→G or T→C and the transversion A→C or T→G each have separate cold/hot rates according to E1. The remaining A→T or T→A event has separate cold/hot rates according to E2. For C/G centers, C→T or G→A has separate rates according to CpG; the other two substitutions share one common rate. Complementary outcomes are pooled to stabilize estimation. **The E1 and E2 coordinate offsets are never reverse-complemented.** This distinction is supported by the calibration data and is required to preserve the oriented long-range signal.

There are nine fitted rate parameters. For each event class, with `m` observed events over `N` applicable site-event exposures, estimate `p = (m + 0.5)/(N + 1)`. This is a weak half-count smoothing of a binomial event rate. For the common C/G transversion rate, there are two possible events per eligible site, so the exposure is twice the site count. The fit uses a marginal event-rate estimator rather than claiming an exact jointly constrained multinomial maximum-likelihood solution. Each queried site's three substitution probabilities are assigned these rates; its unchanged probability is `1 - sum(substitution probabilities)`. All resulting vectors were verified to be valid.

| Parameter | Events | Calibration exposure | Estimated probability | Query exposure |
|---|---:|---:|---:|---:|
| `at_other_cold` | 1,186 | 1,497,706 | 0.000792211 | 373,443 |
| `at_other_hot` | 790 | 499,779 | 0.001581696 | 125,307 |
| `at_transition_cold` | 2,376 | 1,496,429 | 0.001588113 | 99,347 |
| `at_transition_hot` | 3,067 | 501,056 | 0.006122058 | 399,403 |
| `at_transversion_cold` | 1,138 | 1,496,429 | 0.000760811 | 99,347 |
| `at_transversion_hot` | 1,609 | 501,056 | 0.003212209 | 399,403 |
| `cg_transversions` | 3,004 | 4,004,070 | 0.000750361 | 1,002,260 |
| `cg_transition_cold` | 2,324 | 1,499,757 | 0.001549917 | 374,152 |
| `cg_transition_hot` | 4,695 | 502,278 | 0.009348390 | 126,978 |

The expected total is `sum_i (1 - p_i[parent_i])`. It was also independently reconstructed as the sum of each parameter's rate times its query exposure, and matched the file-derived expectation.

## Why the new-parent expectation differs from the calibration average

The new parent has a substantial, directly observed sequence-composition shift in the learned long-range feature. Among A/T-centered sites, E1 is true at **399,403 of 498,750 query sites (80.08%)**, compared with **25.08%** in calibration. The final expectation therefore rises to about 0.6975% even though the calibration mutation fraction is about 0.5048%. Scaling the calibration-average mutation rate by the query count would miss this predictable change in burden. The hot and cold conditional states are both well represented in calibration; this is reweighting observed states, not extrapolating rates beyond their observed feature support. The new sequence was used for prediction and this diagnostic only, not for discovering the rate rules or choosing their parameters.

## Uncertainty and practical limits

The approximate standard error of the fitted expected query count is 61.03 changes, using marginal rate-posterior variances and a delta calculation. Under conditionally independent site outcomes, the standard deviation of the realized count around its fitted expectation is 83.14; combining those contributions approximately gives 103.14. These are scale estimates, not rigorously calibrated coverage guarantees: they omit cross-rate covariance, exploratory model-selection uncertainty, and any remaining dependence between sites. Site independence is unnecessary for the expected-count calculation itself.

The compact rules are empirical predictive findings, not assertions about the hidden generating code. More subtle effects, spatial dependence, and parameter drift cannot be ruled out. Higher-order tables were noisy, and the new offspring is stochastic. Predictions should be assessed with probabilistic scores or calibration, rather than expecting exact identification of the changed positions.

## Verification and reproducibility

The entire saved TSV was independently re-read alongside `queries.tsv`. All 999,880 coordinates matched in order, all probabilities were finite and between zero and one, and the maximum absolute row-sum error was 4.44e-16. The expected count reconstructed from the saved probabilities matched the computed expectation. A second implementation using explicit seven-character strings checked all 16,384 possible contexts against the numeric prediction table. All six original input hashes and the output checksum matched recorded provenance.

Computations used Python 3.13.15, NumPy 2.5.4, seed 0, a four-thread runtime cap, and a Linux process address-space ceiling of 26 GiB. The dataset is small enough for vectorized in-memory counting; no large model or plot was required. No plots were displayed or produced. `requirements.txt` records the installed package versions. All code and artifacts are confined to the case directory, apart from the explicitly authorized dedicated Python environment.

Scripts are split into repeatable phases. Run from WSL using the dedicated environment; outputs are regenerated, and the six source inputs remain unchanged. No repository clone or push is used because this is an isolated blind case with no authorized repository destination.

```bash
%%bash
set -euo pipefail
cd /mnt/d/agent_evaluations/case_002
test -x /home/mew/.cache/case_002_env/bin/python || uv venv /home/mew/.cache/case_002_env
uv pip install --python /home/mew/.cache/case_002_env/bin/python -r requirements.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
PYTHON=/home/mew/.cache/case_002_env/bin/python
# Phase 1: inspect inputs and build sufficient statistics.
"$PYTHON" analyze.py > exploration.log
# Phase 2: broad context comparison and compact-rule validation.
"$PYTHON" models.py > model_validation.log
"$PYTHON" residuals.py > residuals.log
"$PYTHON" rule_models.py > rule_validation.log
"$PYTHON" final_validation.py > final_validation.log
# Phase 3: produce and independently verify deliverables.
"$PYTHON" produce.py > production.log
"$PYTHON" verify_model.py > model_verification.log
"$PYTHON" write_report.py
```

Detailed evidence is retained in `input_diagnostics.json`, `validation.json`, `rule_validation.json`, `final_validation.json`, the corresponding logs, and the chronological `decision_log.md`. The two primary machine-readable outputs are `predictions.tsv` and `summary.json`.
