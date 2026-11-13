# Offspring allele predictions for case 004

The delivered `predictions.tsv` contains 1,000,000 rows, in the exact supplied query order, with columns `contig`, `position`, `p_A`, `p_C`, `p_G`, and `p_T`. The expected number of changed alleles is **4,974.728737409099**, a mean change probability of **0.004974728737409099**. The model retains a high probability of the parent's allele at every site. `summary.json` contains the numerical result and model parameters.

## Inputs and scope

Only the six initially supplied files in this task directory were consulted. No external scientific sources, other project directories, repository history, agents, withheld outcomes, or evaluator code were used. General-purpose numerical packages were installed with `uv` into the designated environment. All code and outputs are in this task directory; the six supplied files remain byte-for-byte unchanged, as verified by SHA-256.

Calibration consists of 80 homologous parent/offspring scaffold pairs, each 50,000 bases long: 4,000,000 transmission opportunities and 19,901 changes. The empirical overall rate is 0.00497525. Sequences contain only A/C/G/T, and homologous lengths agree. The new parent has 20 scaffolds of 50,000 bases. The queries cover all 1,000,000 new-parent sites exactly once, in shuffled order. Assemblies are interpreted as haploid, error-free, circular sequences, as specified. Contig names in the new parent are used for coordinate lookup, not as links to calibration outcomes.

## Hypotheses and evidence

I began with a uniform mutation process and considered source-base effects, short local context, longer context, compositional windows, coordinate periodicity, scaffold variation, and circular offset templating. The detailed development chronology is recorded in `decision_log.md`.

A seeded scaffold permutation (`numpy.random.default_rng(0)`) allocated 40 scaffolds to discovery, 20 to model selection, and 20 to final evaluation. The exact indices are in `splits.json`. The initial basic audit reported aggregate counts from all calibration scaffolds; subsequent detailed feature discovery and model selection excluded the final evaluation set. Entire scaffolds were held out to prevent neighboring-site or repeated-molecule leakage.

The discovery set's change rate was 0.004951. Source-base and centered three-base context tests did not indicate rate heterogeneity (nominal p-values 0.458 and 0.572). Scaffold and coordinate-bin tests were also unremarkable. Circular cross-correlations scanned all nonzero offsets for matches between changed offspring alleles and neighboring parent bases, complementary matches, and mutation-pair clustering. Extreme standardized values of about four are unsurprising across roughly 50,000 offsets per scan and were not used as rules. GC-window and periodicity diagnostics were treated as exploratory, with attention to multiple comparisons.

Centered seven-base contexts showed excess dispersion. Its sparse cells make the nominal chi-square p-value unsuitable as definitive evidence; independent predictive validation was the deciding test. Smoothed seven-base tables gained 99.76 nats of mutation log likelihood on the separate selection set. Centered nine-base tables gained 104.28 nats.

Scanning contiguous windows of lengths 3–10 containing the target localized the strongest signal to **eight parent bases at offsets -4, -3, -2, -1, 0, +1, +2, +3**. A two-component rate model gained 333.67 validation nats. Removing any one position from this eight-base identity reduced gains to about 89–107 nats. Adding flanking positions reduced evidence sharing across occurrences. This pattern supports a context-specific lookup process, with all eight positions relevant, rather than a simple CpG rule, broad GC effect, or nearby-base copying.

Mutation destinations also share information across identical eight-base contexts. Among 1,689 pairs of discovery mutations with the same context, 1,165 (68.98%) share the same alternative allele, versus one third under uniform independent alternatives. A context-specific destination model therefore has strong justification even though the aggregate substitution matrix is close to uniform.

The data support these predictive dependencies; they do not identify the exact hidden generating algorithm. No particular pseudorandom hash or inaccessible generating seed was assumed or reconstructed.

## Alternatives compared during selection

The following rate gains are binary mutation log-likelihood improvements over a fitted uniform rate on the 1,000,000-site selection set. All corresponding priors were fitted using the 2,000,000-site discovery set.

| Eight-base rate prior | Validation improvement, nats |
|---|---:|
| Beta prior, fitted by marginal likelihood | 331.099 |
| Two fixed-rate mixture components | 333.667 |
| Flexible discrete mixture | 350.493 |
| Logistic-normal prior | 350.035 |

The logistic-normal prior was selected because its two parameters capture essentially all the flexible mixture's predictive benefit. The 0.458-nat difference across one million validation sites does not justify a substantially more flexible mixing distribution.

For the conditional alternative-allele distribution, gains over uniform alternatives were:

| Destination prior | Validation improvement, nats |
|---|---:|
| Symmetric Dirichlet, fitted concentration | 356.721 |
| One preferred allele with fitted preference strength | 353.272 |
| Symmetric discrete simplex mixture | 355.551 |

The symmetric Dirichlet was selected. Its discovery-fit total concentration was 0.977995. A broad pseudocount grid independently favored a value near one. These gains add to the binary-rate gains because the likelihood factors into the occurrence of a change and its destination conditional on a change.

Some exploratory two-point mixtures reached their iteration cap, and unsupported components could numerically collapse; they were not used in the final model. The helper was subsequently guarded against division by zero for reproducibility. The final logistic-normal and Dirichlet optimization paths are separate from that exploratory EM implementation.

## Final probabilistic model

For each context c, let n_c be its calibration occurrence count, m_c its changed-allele count, and d_cb its count of changes to alternative base b. The context determines the central parental allele a_c.

1. Its change rate p_c has a logistic-normal prior: logit(p_c) follows Normal(mu, sigma squared).
2. Conditional on a change, the three alternative-base probabilities q_c have a symmetric Dirichlet prior with each parameter tau/3.
3. The two prior quantities are independent conditional on the fitted hyperparameters. Contexts share these hyperparameters but have separate latent rates and destination distributions.

The rate posterior is proportional to p_c raised to m_c times (1-p_c) raised to (n_c-m_c), multiplied by its logistic-normal prior density. Its mean r_c is computed by 80-point Gaussian quadrature. The prediction for the parental allele is 1-r_c. For an alternative allele b, it is:

`r_c * (d_cb + tau/3) / (m_c + tau)`.

Thus a rare or unobserved context is automatically pulled toward the population prior and never receives an unjustified zero mutation probability. A context with no observed changes receives uniform alternatives conditional on a future change. Repeat observations update both rate and destination probabilities. The alternative probabilities and unchanged probability sum to one by construction.

Hyperparameters maximize the calibration marginal likelihood, integrating out context-specific latent values. After the specification was locked, it was fitted to the combined 60 discovery/selection scaffolds for final evaluation. The same specification was then refitted to all 80 scaffolds for delivery. Final fitted parameters are:

- Logistic-normal logit mean: **-5.69395284739004**.
- Logistic-normal logit standard deviation: **0.8946945036380947**.
- Symmetric Dirichlet total concentration: **1.048019815816202**.

`model.json` stores the hyperparameters, quadrature representation, offsets, and defaults. `model_context_table.npz` stores all 65,536 context predictions, counts, and conditional posterior rate variances.

Circular indexing is applied separately within each scaffold. An offset before its first base uses its final bases; an offset after its final base uses its first bases. Coordinates in the input and output remain one-based. Contexts never cross between scaffolds.

## Locked final evaluation

The final evaluation set contains 20 complete scaffolds, 1,000,000 sites, and 4,972 observed changes. No model choices were changed after this evaluation.

| Model trained on the other 60 scaffolds | Mean allele negative log likelihood | Mean multiclass Brier score | Expected changes |
|---|---:|---:|---:|
| Uniform change rate and alternatives | 0.036793077 | 0.009911039 | 4,976.333 |
| Source-base substitution matrix | 0.036794850 | 0.009911069 | 4,978.611 |
| Eight-base rate; uniform alternatives | 0.036326454 | 0.009898753 | 4,973.801 |
| Selected rate and destination model | **0.035903488** | **0.009893091** | **4,973.801** |

Lower scores are better. The selected model improves allele log likelihood by **889.589 nats** relative to the uniform baseline: 466.622 nats from rate prediction and 422.967 from destination prediction. A scaffold bootstrap with 5,000 resamples and seed 0 gives a 95% percentile interval of **792.936–991.619 nats** for the total gain. This interval is descriptive and conditional on the trained model, rather than a correction for every earlier exploratory comparison. Mutation-discrimination AUC is **0.584708**. Decile calibration and scaffold-level scores are available in `validation.json`.

The modest AUC is expected for a stochastic process with low rates and sparse context observations. A high unchanged-allele accuracy alone would be uninformative. Proper probability scores and expected versus observed mutation totals provide a more useful check.

## Output verification and reproducibility

`solve.py` writes predictions through a temporary file and then atomically replaces `predictions.tsv`. It re-reads the delivered TSV, checks all column names and exact query order, verifies finite nonnegative probabilities and row sums, and recalculates the expected change count from the rounded on-disk probabilities. The maximum row-sum error is **5.500044863993026e-13**. All six input hashes are checked again.

`verify.py` separately checks 4,000 circular contexts against direct modular indexing at scaffold boundaries and seeded interior coordinates, verifies that the encoded central allele is correct, checks all output probabilities, and confirms exhaustive unique query coverage. Numerical provenance, input/output hashes, source-code hashes, package versions, and all final model parameters are recorded in `provenance.json`, `model.json`, `summary.json`, and `run.log`. Detailed exploration outputs accompany their scripts. The fixed solve pipeline was rerun unchanged to record final source hashes and confirm deterministic reproduction; no models were adapted after viewing final evaluation results.

Commands below run from the permitted task directory in WSL. They do not access any other project, require Git, alter the supplied inputs, access controlled data, or publish results. The environment is the dedicated path authorized for this task. The scripts set numerical thread counts to four before importing numerical libraries. Their working arrays are comfortably below 26 GB, with the largest dense model tables having only 65,536 contexts and four alleles; the exploratory ten-base tables have 1,048,576 entries. Every stochastic operation uses seed 0.

```bash
%%bash
set -euo pipefail
cd /mnt/d/agent_evaluations/case_004
uv venv --allow-existing /home/mew/.cache/case_004_env --python 3.13
uv pip install --python /home/mew/.cache/case_004_env/bin/python -r requirements.txt
PY=/home/mew/.cache/case_004_env/bin/python

# Phase 1: input audit and independent discovery/model selection.
"$PY" explore.py > exploration.log
"$PY" diagnose.py >> exploration.log
"$PY" context_models.py >> exploration.log
"$PY" window_scan.py >> exploration.log
"$PY" refine_model.py >> exploration.log
"$PY" destination_models.py >> exploration.log
"$PY" check_structure.py >> exploration.log

# Phase 2: fixed-model held-out evaluation, all-data fit, and delivery.
"$PY" solve.py > run.log 2>&1
"$PY" verify.py > verification.log
```

The environment and output paths are reused, and outputs are deterministically replaced. There are no plots, notebook displays of sequence data, or inline figures.

## Remaining uncertainty

The expected total of 4,974.73 is an expectation, not a claim that the next offspring has exactly that many changes. An individual site's future allele remains stochastic, and the new offspring's outcomes were unavailable throughout.

Conditional on the selected model and fitted hyperparameters, uncertainty in the latent context rates gives a posterior standard deviation of **18.35** for the expected total. Including new transmission randomness gives a posterior predictive count standard deviation of **72.50**. These calculations account for shared uncertain rates across repeated occurrences of a context; they do not assume those posterior uncertainties are independent per occurrence. They exclude hyperparameter and model uncertainty and therefore are not advertised as complete frequentist confidence intervals.

The posterior sharing assumptions could miss subtle interactions, dependence among realized mutations, or a different distribution of context-specific rates. Circular offset and spatial diagnostics did not reveal a compelling alternative, and the fixed held-out evaluation supports generalization, but they cannot prove absence of every possible structure. The exact eight-base lookup is well supported; treating its latent rates as logistic-normal and its conditional alternatives as Dirichlet remains a statistical approximation. The final fitted priors were selected for predictive performance and parsimony, without claiming they are the hidden generator's exact distributions.
