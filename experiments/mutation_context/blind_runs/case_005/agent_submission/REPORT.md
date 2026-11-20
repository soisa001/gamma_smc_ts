# Offspring allele prediction report

## Result

The selected model assigns probability **0.995042** to the parent's allele and **0.001652666666667** to each of the other three alleles at every queried site. The expected number of changed alleles across **1,000,000** query rows is **4,958.000000**. `predictions.tsv` contains the requested six columns in the original query order. `summary.json` contains the estimate, model parameters, uncertainty estimates, hashes and audit results.

The estimate is learned from the calibration data; it is not rounded to a conjectured synthetic generating parameter. The withheld offspring was never accessed. This is a probability forecast: individual realized mutations cannot be reliably identified under the selected model.

## Scope and input integrity

Only the six supplied files in this case directory were used. No other project files, repository history, hidden outcomes, evaluator code, chats, agents, scientific references or external datasets were consulted. General-purpose Python dependencies were installed with uv in the permitted dedicated environment. All generated scripts and result files are in this case directory; the environment itself is at `/home/mew/.cache/case_005_env`.

The calibration has 80 homologous parent/offspring scaffold pairs, each 50,000 bases long, totaling 4,000,000 sites. The new parent has 20 scaffolds of 50,000 bases, totaling 1,000,000 sites. Homologous names and lengths were checked. All sequences contain only A, C, G and T. Assemblies were compared directly at matching coordinates, as warranted by the task's error-free, single-nucleotide-change assumptions. The new assembly is treated as a separate parent: repeated scaffold names do not imply matching sequence or offspring outcomes across calibration and prediction assemblies.

There are 19,832 calibration changes. Queries have 1,000,000 rows, 1,000,000 unique coordinates, and 0 duplicate rows. Every coordinate is valid. All new-parent sites are queried exactly once. SHA-256 checks before and after output generation verify that all six supplied inputs remained unchanged; their hashes appear in `exploration.json` and `summary.json`.

## Hypotheses and empirical evidence

The leading hypothesis is a common mutation probability with a uniform replacement among the three alternative nucleotides. Alternatives considered were source-base-specific risk, GC versus AT risk, local motifs and reverse-complement motifs, sequence composition, short homopolymers, coordinate-dependent rates, periodic coordinate effects, scaffold heterogeneity, mutation clustering and nonuniform/context-dependent replacement spectra.

The calibration transition counts (rows parent, columns offspring, in A/C/G/T order) are:

| Parent | A | C | G | T |
|---|---:|---:|---:|---:|
| A | 996512 | 1721 | 1687 | 1698 |
| C | 1703 | 995895 | 1667 | 1636 |
| G | 1678 | 1653 | 994050 | 1580 |
| T | 1601 | 1627 | 1581 | 993711 |

The overall change rate is 19832/4000000 = 0.00495800. Source-base rates show a weak nominal association (chi-square 8.778, 3 df, p=0.0324), but this is one of many examined hypotheses. The combined test of all 12 directed substitution rates against the common-rate uniform-replacement model gives chi-square 14.244 on 11 df, p=0.220. Symmetric replacement classes A↔C/G↔T, A↔G/C↔T, and A↔T/C↔G have counts 6,585, 6,628 and 6,619 (p=0.925 against equality). GC-versus-AT risk has p=0.981.

Scaffold change counts range from 207 to 287, averaging 247.9. Their observed variance is 237.84; a common-rate binomial model gives heterogeneity chi-square 76.17 on 79 df, p=0.569. Ten equal coordinate bins have p=0.444. Adjacent identical-base counts have p=0.817. Coarse GC composition bins in a 1,001-base circular window have p=0.657. These tests support pooling but do not prove absence of every possible effect.

Circular offsets -30 through +30 were checked separately, excluding zero. Their source-stratified mutation-risk association statistics show no compelling signal after multiple comparisons. Among mutations, destination-base association with each offset was also tested conditional on the source base; the smallest p-value was 0.0255 at +15, which does not survive correction over 60 offsets. Positional moduli 2 through 40 were screened. Mutation pair counts at lags 1, 2, 3, 5, 10, 50, 100 and 1,000 were 96, 94, 100, 77, 110, 91, 85 and 80, against approximately 98.33 expected at each under independence. These exploratory summaries do not suggest a substantial burst process. Pair counts are diagnostics only and do not use unavailable offspring information as prediction features.

## Validation design

Every calibration site received an out-of-fold prediction in each of two five-fold validation schemes:

1. **Whole scaffolds:** seed-0 permutation assigns 16 entire scaffolds to each fold; each fit uses the other 64. This tests transfer to unseen scaffolds.
2. **Contiguous regions:** each 50,000-base scaffold is divided into five consecutive 10,000-base regions; each fold holds out the corresponding region of every scaffold. This checks sensitivity to positional structure.

No held-out offspring allele contributes to training counts or priors. Context uses parent sequence only; circular windows wrap within the same scaffold and never across scaffolds. Parent contexts can span a held-out boundary because the whole parent assembly is available at prediction time. Initial descriptive screening used all calibration pairs, and the two CV schemes reuse the same outcomes. These are internal model-comparison checks, not independent confirmatory experiments or performance estimates on the withheld new offspring.

The primary metric is four-allele negative log likelihood (NLL), including the parent allele. Brier score is also recorded. Accuracy would be dominated by unchanged sites and was not used for model selection. The pooled model's scaffold-CV log loss is 0.036703648643 per site and Brier score is 0.009883226415 (sum over four classes).

Candidate categorical risk models estimate the mutation probability in each category and divide it equally among alternatives. A prior of strength alpha is centered on the training fold's pooled rate: (mutations + alpha × pooled rate)/(exposure + alpha). Full-spectrum models estimate the four outcomes coded by parent XOR offspring, with a four-category Dirichlet prior centered on the same homogeneous model. Reverse-complement context categories are canonicalized; XOR replacement classes are invariant under simultaneous complementation. Global source-spectrum models use a 0.5 pseudocount per cell. Global XOR-spectrum models use a 0.5 pseudocount per alternative; transition-spectrum models pool the two transversion classes. Source-rate-plus-XOR-spectrum models combine independently estimated source-specific risks and a shared replacement spectrum. Probabilities are clipped only where necessary for numerical safety in validation.

Risk candidates cover source/GC classes, immediate dinucleotides, full and reverse-complement 3-, 5- and 7-mers, 10/50/100 positional bins, selected screened positional moduli, GC windows of 11/51/201/1,001 bases, and zero/one/two identical immediate neighbors. Full-spectrum candidates cover immediate dinucleotides, 3- and 5-mers, and reverse-complement 7-mers. Alpha is an exposure pseudocount, not a count of mutations. All alternatives and smoothing settings evaluated are listed below; none were silently discarded.

## Complete model comparison

Values are total held-out NLL **minus pooled NLL** across all four million sites. Positive values are worse. The last column is an approximate standard error of the scaffold-split total NLL difference, computed as sqrt(80) × the sample standard deviation of paired per-scaffold NLL differences. It describes heterogeneity of the validation difference and is not a formal independent-CV uncertainty estimate because training folds overlap.

| Candidate | Alpha | Scaffold ΔNLL | Region ΔNLL | Paired scaffold SE |
|---|---:|---:|---:|---:|
| pooled | 0 | 0.000 | 0.000 | 0.000 |
| source | 0 | -1.982 | -2.881 | 3.076 |
| source_gc | 0 | 0.606 | 0.546 | 0.343 |
| context_3 | 0 | 29.064 | 29.303 | 8.069 |
| context_rc_3 | 0 | 17.670 | 20.480 | 5.859 |
| context_5 | 1000 | 404.257 | 411.133 | 27.426 |
| context_5 | 10000 | 47.444 | 50.639 | 8.360 |
| context_5 | 100000 | 2.181 | 2.624 | 1.065 |
| context_rc_5 | 1000 | 266.360 | 263.265 | 21.128 |
| context_rc_5 | 10000 | 57.679 | 57.035 | 9.190 |
| context_rc_5 | 100000 | 2.689 | 2.666 | 1.397 |
| context_7 | 1000 | 200.089 | 216.343 | 22.165 |
| context_7 | 10000 | -3.998 | -1.389 | 2.701 |
| context_7 | 100000 | -0.768 | -0.492 | 0.277 |
| context_rc_7 | 1000 | 296.550 | 315.902 | 27.576 |
| context_rc_7 | 10000 | -6.613 | -3.697 | 3.771 |
| context_rc_7 | 100000 | -1.379 | -1.071 | 0.394 |
| dinuc_-1 | 0 | 7.205 | 1.588 | 5.094 |
| dinuc_1 | 0 | -0.090 | 3.286 | 4.698 |
| position_bins_10 | 0 | 10.552 | 0.000 | 3.671 |
| position_bins_50 | 0 | 44.131 | 0.000 | 7.986 |
| position_bins_100 | 0 | 73.937 | 0.000 | 11.497 |
| position_mod_2 | 0 | 0.649 | 0.892 | 1.208 |
| position_mod_3 | 0 | 1.440 | 1.792 | 0.621 |
| position_mod_7 | 0 | -0.595 | 0.651 | 3.565 |
| position_mod_13 | 0 | 8.780 | -1.916 | 4.651 |
| position_mod_39 | 0 | 14.936 | 12.659 | 8.871 |
| gc_window_11 | 0 | 7.323 | 9.383 | 2.601 |
| gc_window_51 | 0 | 2.627 | 2.378 | 3.448 |
| gc_window_201 | 0 | 2.124 | 2.741 | 1.183 |
| gc_window_1001 | 0 | 0.710 | 0.397 | 0.705 |
| homopolymer_neighbor | 0 | 1.400 | 2.571 | 0.760 |
| source_spectrum | — | 5.650 | 8.258 | 3.853 |
| xor_spectrum | — | 4.039 | 1.225 | 1.005 |
| transition_spectrum | — | 1.173 | 0.622 | 0.495 |
| source_rate_xor_spectrum | — | 2.057 | -1.656 | 3.183 |
| dinuc_-1_full_spectrum | 0 | 26.850 | 22.141 | 8.138 |
| dinuc_-1_full_spectrum | 1000 | 26.539 | 21.857 | 8.097 |
| dinuc_-1_full_spectrum | 10000 | 23.953 | 19.495 | 7.747 |
| dinuc_1_full_spectrum | 0 | 26.082 | 25.619 | 8.535 |
| dinuc_1_full_spectrum | 1000 | 25.808 | 25.346 | 8.492 |
| dinuc_1_full_spectrum | 10000 | 23.518 | 23.063 | 8.121 |
| context_3_full_spectrum | 0 | 114.941 | 93.365 | 13.350 |
| context_3_full_spectrum | 1000 | 110.213 | 89.068 | 13.078 |
| context_3_full_spectrum | 10000 | 78.025 | 60.100 | 11.051 |
| context_rc_3_full_spectrum | 0 | 58.411 | 44.720 | 9.816 |
| context_rc_3_full_spectrum | 1000 | 57.183 | 43.641 | 9.716 |
| context_rc_3_full_spectrum | 10000 | 47.639 | 35.302 | 8.901 |
| context_5_full_spectrum | 10000 | 124.770 | 117.192 | 12.964 |
| context_5_full_spectrum | 100000 | 4.093 | 3.141 | 1.652 |
| context_rc_5_full_spectrum | 10000 | 145.915 | 142.507 | 14.661 |
| context_rc_5_full_spectrum | 100000 | 3.562 | 3.194 | 2.236 |
| context_rc_7_full_spectrum | 10000 | -5.748 | -5.019 | 6.358 |
| context_rc_7_full_spectrum | 100000 | -2.661 | -2.562 | 0.677 |

The small source-rate gains (about 2–3 total nats) are smaller than their approximately 3.1-nat paired standard errors. Strongly shrunk reverse-complement 7-mer risk models gain 6.61 nats in scaffold CV and 3.70 nats in region CV; their respective paired standard errors are 3.77 and 3.95 nats, with approximate 95% intervals including zero. The full-spectrum version also has a small gain. These gains are tiny after screening many candidate models and smoothing strengths; there is no stable, identifiable context rule. Less-shrunk context models markedly overfit, and direct destination-offset tests provide no corroboration. A conservative common-rate model is therefore selected. The choice explicitly favors parsimony and stability over the absolute minimum internal CV score. The best context model is not claimed to have been disproved.

## Final fitted model and expected count

Fit the chosen model on all calibration sites:

```
q = number of observed changes / number of calibration sites
  = 19,832 / 4,000,000 = 0.004958
P(offspring = parent) = 1 - q = 0.995042
P(offspring = b) = q / 3 for each b different from the parent
E[number changed over query rows] = sum_i q = 4958.000000
```

Although the estimated hazard is constant, the placement of the high probability among p_A/p_C/p_G/p_T depends on the observed new-parent allele at each queried coordinate. Probabilities are output with 17 significant digits. The output expected count is the sum of non-parent probability over the original query rows; with this query set every physical site is counted once.

The 95% exact binomial interval for q is [0.0048894013, 0.0050273169]. Multiplying by the query count gives [4889.401, 5027.317] for uncertainty in the expected count under the homogeneous model. A normal predictive interval for the realized number of changes, adding approximate rate-estimation variance to binomial outcome variance, is [4804.1, 5111.9]. This is distinct from the expected number and assumes independent changes and a stable process. It is not a guarantee about the withheld outcome.

## Output audit and reproducibility

`finalize.py` independently re-read every exported row alongside `queries.tsv`, checked exact coordinate and order agreement, checked the six column names, checked finite nonnegative probabilities bounded by one, checked row sums, recomputed the expected number of changes from serialized probabilities, and rechecked all original input hashes. All 1,000,000 rows passed. Maximum absolute row-sum error was 0. The expected count recomputed from the file was 4958.000000000018. The final prediction SHA-256 is `4bfe3ca7ed8f550d7ab8bd059df33aff3f9586ff2b75d3531f59f5687bd6e7e7`.

The scripts are idempotent and are divided into phases. `explore.py` reads and audits the inputs and writes descriptive summaries; `validate.py` evaluates risk and global spectrum models; `validate_spectra.py` evaluates context-dependent replacement probabilities; `finalize.py` fits the selected model, writes predictions atomically, audits them, and writes this report and summary. NumPy RNG seed is 0; only fold assignment requires randomness. BLAS/OpenMP thread limits are four. No large model or training operation requires anywhere near the 26 GB memory budget. No figures or notebook data displays were produced. LightGBM was installed as an available general-purpose dependency but was not used to fit the final or comparison models.

Run inside WSL/Linux (the existing environment was provisioned with uv):

```bash
%%bash
set -euo pipefail
cd /mnt/d/agent_evaluations/case_005
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
uv venv --allow-existing /home/mew/.cache/case_005_env
uv pip install --python /home/mew/.cache/case_005_env/bin/python -r requirements.txt
PY=/home/mew/.cache/case_005_env/bin/python
"$PY" explore.py > exploration_stdout.txt
"$PY" validate.py > validation_stdout.txt
"$PY" validate_spectra.py > validation_spectra_stdout.txt
"$PY" finalize.py > finalization_stdout.txt
```

No clone, pull or push was performed: no repository is supplied within the authorized blind-evaluation inputs. `requirements.txt` pins the installed package versions. Complete numerical results are retained in `exploration.json`, `diagnostics.json`, `validation.json`, `validation_spectra.json` and `model_comparison.json`; `decision_log.md` records the chronological decisions.

## Remaining uncertainty

Only one transmission is observed for each parent, so individual-site propensities cannot be estimated without a shared pattern. The observed data cannot rule out subtle long-range effects, rare motifs, latent state, or dependencies that do not substantially alter marginal probabilities. Multiple exploratory comparisons can produce small apparent gains by chance. Conversely, pooling could discard a genuine weak source or context effect; the validation evidence does not justify confidently resolving effects of that size. Binomial uncertainty intervals are conditional on the selected model and omit this structural uncertainty. The task states that the same process applies to the new parent, but its offspring outcomes remain unavailable; no claim is made about realized test-set performance.
