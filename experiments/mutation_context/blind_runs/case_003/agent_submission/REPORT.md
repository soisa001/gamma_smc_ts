# Offspring-allele prediction report

## Deliverables and main result

`predictions.tsv` contains the requested 1,000,000 rows, with columns `contig`, `position`, `p_A`, `p_C`, `p_G`, and `p_T`, in exactly the query order. All probabilities include retention of the parent allele. The expected total number of changed alleles is **5,000.441186784352**, or an average change probability of **0.005000441186784352**. The corresponding expectation recomputed from the serialized TSV is 5,000.441186801279; the insignificant difference comes from decimal serialization. Both values are recorded in `summary.json`.

The selected predictor uses a circular six-base parent-sequence context, from offsets -3 through +2 relative to the queried base. It estimates both the context-specific probability of change and the distribution of the three possible replacement alleles, with empirical-Bayes shrinkage. The frozen model performed well on reserved calibration scaffolds; it does not identify which individual sites will change in the withheld offspring.

## Scope, inputs, and integrity

Only the six supplied files in this case directory were consulted for task information and data. No external scientific sources, other project directories, repository history, other chats, other agents, hidden outcomes, or evaluation machinery were accessed. The general-purpose numerical dependencies were installed into the explicitly allowed dedicated environment. All generated code and case outputs are in this case directory. None of the six supplied files was edited. Input SHA-256 hashes and sizes are recorded in `input_manifest.json` and were rechecked after finalization. Because that manifest was created during finalization, it verifies that phase explicitly; earlier scripts also opened supplied inputs only for reading.

There are 80 calibration scaffolds of 50,000 bases each, for 4,000,000 paired bases, and 20 target scaffolds of the same length. Matching calibration contig identifiers and lengths were checked, and every calibration base was verified as A/C/G/T. There are 19,987 calibration changes, a raw average rate of 0.00499675. The million query rows are all unique valid sites and cover the million target bases. Target and calibration sequences share no identical whole-scaffold sequences. Contig identifiers were treated as labels, not transferable biological features.

The supplied metadata describe error-free haploid assemblies, one transmission, circular topology, and one-based coordinates. Accordingly, observations are direct parent-to-offspring allele outcomes. Context extraction always wraps separately within each scaffold, including both ends. No diploid segregation or assembly-error model was added.

## Hypotheses and chronological investigation

The chronological reasoning is recorded separately in `decision_log.md`. The main alternatives were:

1. **Constant or parent-base-specific mutation and substitution rates.** The overall spectrum is fairly balanced, but these predictors omit a substantial context effect.
2. **Scaffold-specific or coordinate-phase rates.** Calibration scaffold change counts are consistent with sampling variation under a common average rate (exploratory chi-square p=0.685). Coordinate phases with periods 2, 3, 4, 5, 10, and 100 were checked. Isolated nominal phase p-values were not persuasive after considering the number of exploratory tests. No transferable position rule was supported.
3. **Centered contexts.** Contexts of length 1, 3, 5, 7, and 9 were compared with smoothed mutation counts. The 5- and 7-base windows improved development mutation log loss. Hot 7-base contexts often shared their first six letters while tolerating all four final letters. This suggested testing an asymmetric rule.
4. **All contiguous windows containing the focal base, lengths 1 through 7.** All 28 windows were compared. For each, mutation-count pseudocounts of 10, 30, 100, 300, 1,000, 3,000, 10,000, 30,000, and 100,000 were crossed with destination pseudocounts of 0.3, 3, 10, 30, 100, 1,000, and 1e9. The clear best window was -3..+2. A length-seven extension on either side performed worse, consistent with dilution into sparse subcontexts. No reverse-complement symmetry was imposed: directional context was directly supported by the observed comparisons; symmetry was not itself formally tested.
5. **Different shrinkage distributions for the selected six-base window.** Fixed global shrinkage, an estimated beta-binomial rate prior, an estimated logit-normal rate prior, and a flexible nonparametric mixture rate prior were compared. Replacement-allele probabilities were modeled with a fixed symmetric Dirichlet prior, an estimated symmetric Dirichlet prior, or a flexible mixture of symmetric Dirichlet priors. The smooth logit-normal rate plus estimated Dirichlet destination model was chosen.

The context model is an empirical predictive description. The data do not prove that the generating mechanism literally samples parameters from the chosen prior families.

## Validation design and development evidence

A seed-0 random permutation of the 80 scaffolds fixed 48 training scaffolds, 16 development scaffolds, and 16 audit scaffolds. Exact membership is in `split.json`. All sites from a scaffold remain in one partition, avoiding neighboring sequence overlap across partitions. The development set was used for model-family and hyperparameter selection. After selection, the model family was frozen, refitted on the 64 training-plus-development scaffolds, and evaluated once on the 16 audit scaffolds. Finally it was refitted on all calibration data for the requested predictions.

Initial descriptive summaries and the first display of context-enrichment examples used all calibration data. Therefore the audit is protected from model-selection scoring and parameter fitting, but is not a completely untouched dataset in the strictest exploratory sense. This limitation is recorded explicitly rather than presenting the audit as externally blinded. New-parent offspring outcomes were never available.

All log losses below are natural-log four-allele losses per site; smaller is better.

| Development model | Log loss |
| --- | ---: |
| Parent-base-only, smoothed | 0.0378710040 |
| Centered five-base context, best tested smoothing | 0.0367828492 |
| Six-base context -3..+2, fixed rate/destination smoothing | 0.0352038271 |
| Seven-base context -4..+2, best tested smoothing | 0.0362926665 |
| Seven-base context -3..+3, best tested smoothing | 0.0362809705 |
| Six-base, beta-binomial rate plus estimated Dirichlet destination | 0.0351987843 |
| Six-base, logit-normal rate plus estimated Dirichlet destination | 0.0351630315 |
| Six-base, nonparametric rate plus estimated Dirichlet destination | 0.0351624187 |

The nonparametric advantage over the logit-normal model is approximately 0.0000006 per site, too small to justify its additional flexibility here. The fitted flexible destination mixture was slightly worse than a single estimated Dirichlet. The selected training-only prior had logit-rate mean -5.66687, logit-rate standard deviation 0.85305, and per-alternative Dirichlet concentration 0.54216. Full comparison records are in `context_validation.json`, `window_validation.json`, and `bayes_validation.json`, with readable exploration logs.

## Reserved audit results

The audit comprises 800,000 sites across 16 scaffolds, with 3,989 observed changes. The selected model predicted **3,994.1391** changes.

| Audit model | Four-allele log loss | Four-allele Brier score |
| --- | ---: | ---: |
| Constant change rate, uniform alternative allele | 0.0368842368 | 0.00993934998 |
| Six-base model with fixed smoothing | 0.0343811477 | 0.00990064089 |
| Frozen selected empirical-Bayes model | 0.0343818301 | 0.00990011488 |

The selected model reduces log loss by approximately 6.8% relative to the constant baseline. A seed-0 bootstrap resampling the 16 audit scaffolds, using 20,000 replicates, gives a 95% interval of [0.00237082, 0.00263620] for baseline-minus-selected log loss. This interval describes audit sampling variation under scaffold resampling, not all model-selection or generative uncertainty.

The fixed-smoothed six-base model has a marginally smaller audit log loss than the selected model, by only 0.00000068, while the selected model has a smaller Brier score. No selection change was made after audit inspection. The reliable result is the benefit of the six-base motif; evidence for a particular smooth prior is much weaker.

Audit calibration by ascending predicted-rate decile is retained in `audit_validation.json`. The lowest decile predicted 131.2 changes and observed 147, and the highest predicted 1,249.3 and observed 1,171. Thus aggregate calibration is strong, but the distribution of predicted rates is not perfectly calibrated: there is some apparent overprediction in the highest bins and underprediction in several lower bins. I did not fit a correction on audit outcomes. The final all-calibration fit has more information per motif than the audit fit.

## Final statistical model

Let c be the ordered parent sequence at offsets [-3,-2,-1,0,+1,+2]. There are 4^6 = 4,096 contexts. Let N_c be the calibration exposure count, M_c the number of changed offspring alleles, and A_cj the count of changes to alternative base j. The central parent base is encoded as the fourth of the six letters.

For the mutation probability r_c, the prior is:

    logit(r_c) ~ Normal(mu, sigma^2)
    M_c | r_c, N_c ~ Binomial(N_c, r_c)

The hyperparameters are fitted by maximizing the marginal calibration likelihood across contexts. Integrals are approximated on a fixed 500-point grid of logit values from -15 to 0, using normalized Gaussian prior weights and stable log-sum-exp arithmetic. The site-specific rate prediction is the posterior mean E[r_c | N_c, M_c]. The grid covers probabilities from about 3e-7 to 0.5, well outside the observed relevant rate range; it is a numerical approximation, not evidence for sharp biological bounds.

Conditional on mutation, the three possible replacement probabilities theta_c follow a symmetric Dirichlet prior with per-allele concentration a. The hyperparameter a is fitted by maximizing the Dirichlet-multinomial marginal likelihood. Its posterior mean is:

    E[theta_cj | A_c] = (A_cj + a) / (M_c + 3*a)

The prior treats rate and destination proportions as independent. Under the binomial/multinomial factorization this yields:

    P(offspring = parent | c) = 1 - E[r_c | data]
    P(offspring = j | c) = E[r_c | data] * E[theta_cj | data]

The final hyperparameters fitted to all 4,000,000 calibration bases are:

- Logit-rate mean mu: **-5.6398513724**.
- Logit-rate standard deviation sigma: **0.8290602418**.
- Implied grid prior mean mutation probability: **0.00496172298**.
- Per-alternative Dirichlet concentration a: **0.6191702523**.

The fitted prior mean is not the target expectation. Target motif frequencies and posterior context-specific rates determine the requested total. For query i, with parent allele b_i, its change probability is 1-p_i(b_i). Summing this quantity over query rows gives 5,000.441186784352. Linearity of expectation makes this sum valid without assuming independence among offspring mutations.

## Output verification and reproducibility

`fitted_model.npz` stores context exposures, alternative counts, the complete 4,096-by-4 probability table, and prior integration grid/weights. `summary.json` records the final model, parameters, environment versions, checks, expectation, and prediction SHA-256.

All one million output rows were read back and checked for exact query-order and coordinate identity, expected column shape, finite nonnegative probabilities, and row sums equal to one within floating-point serialization precision. Maximum absolute row-sum error is **5.21e-13**. A separate verification script compared the vectorized context extractor with direct scalar modular-indexing calculations at **19,898** sampled and boundary positions, including every target scaffold's first three and last two positions. The prediction checksum was verified independently. Details are in `verification.json`.

The prediction SHA-256 is:

    e40d92df593618303edfa14e2e40ba5b7f1e8cd1134a09a1bb23bda54fbfdd7f

Scripts are idempotent, recomputing deterministic outputs rather than appending prediction rows. Random sampling is seeded with 0; NumPy/BLAS/OpenMP thread settings are fixed at four. Data structures are far below the 26 GB memory budget. No plots were needed or displayed. The only case logs with append semantics are the intentionally chronological decision notes.

The primary scripts are `analysis.py` (input handling and initial statistics), `context_explore.py` (centered contexts and partition), `window_explore.py` (asymmetric-window selection), `model.py` (prior-family comparison and reusable model functions), `finalize.py` (frozen audit and full-data predictions), and `verify.py` (separate correctness checks).

From an environment where the supplied case directory is already available, final outputs can be regenerated with these phased commands. No repository is present or required for this blind case, so no repository URL, push, or clone operation is invented.

```bash
%%bash
set -euo pipefail
cd /mnt/d/agent_evaluations/case_003
uv venv --allow-existing /home/mew/.cache/case_003_env
uv pip install --python /home/mew/.cache/case_003_env/bin/python -r requirements.txt
# Phase 1: fit the frozen selected model, audit, and produce outputs.
/home/mew/.cache/case_003_env/bin/python finalize.py > final_run.log
# Phase 2: independent format and circular-context verification.
/home/mew/.cache/case_003_env/bin/python verify.py
```

The existing `split.json` is the fixed, seed-0 partition. For a full rerun of exploration, run `analysis.py`, `context_explore.py`, `window_explore.py`, and `model.py` in that order before `finalize.py`. Exploration logs and tabular candidate results are already preserved.

## Remaining uncertainty

The withheld offspring contains stochastic outcomes; even accurately calibrated probabilities cannot predict its exact mutation count or sites. Roughly five observed changes per six-base context on average leave substantial uncertainty about rare contexts and replacement preferences. Shrinkage addresses this statistically but does not remove it. Empirical-Bayes hyperparameters are fitted as point estimates, so the deliverable does not integrate uncertainty over prior-family selection or hyperparameters.

The model assumes the same directional motif process applies to target scaffolds, as stated by the task. It does not model unobserved lineage effects, long-distance interactions, mutation clustering, or dependencies among changes. The available comparisons support a six-base signal but cannot exclude additional weak effects. No very broad genomic-feature learner or reverse-complement pooling model was tested, and no claims are made that all possible mechanisms were ruled out.

Finally, the audit is internal calibration evidence with the exploratory limitation stated above. Only scoring against the withheld target offspring can establish performance on this exact requested transmission.
