from analyze import *
def main():
    s=json.loads((ROOT/'summary.json').read_text());v=s['validation'];basic=json.loads((ROOT/'validation.json').read_text());pars=json.loads((ROOT/'final_model_parameters.json').read_text())
    table='\n'.join(f"| {b['lower']:.3f}–{b['upper']:.3f} | {b['sites']:,} | {b['expected_changes']:.1f} | {b['observed_changes']:,} |" for b in v['calibration_bins'])
    report=f'''# Offspring allele prediction: case 006

## Delivered result

`predictions.tsv` contains {s['number_of_queries']:,} query rows in the exact supplied order, with columns `contig`, `position`, `p_A`, `p_C`, `p_G`, and `p_T`. The expected number of changed alleles is **{s['expected_number_of_changed_alleles']:.6f}**. This is the sum of one minus the predicted parent-allele probability across the queried rows; it is not a rounded mutation-count prediction. There are {s['number_of_unique_queried_sites']:,} distinct queried coordinates. `summary.json` contains the numerical result, validation summaries, and input/output hashes.

The final estimator uses the oriented parent sequence from five bases before the queried base through four bases after it, with circular wrapping. A flexible empirical Bayes prior shares statistical strength across ten-base contexts while preserving context-specific evidence for each alternate allele. No withheld new-offspring outcome was available or used.

## Authorized inputs and integrity

Only the six supplied case files were inspected for task information: TASK.md, metadata.json, calibration_parent.fasta, calibration_offspring.fasta, new_parent.fasta, and queries.tsv. No other project directories, repository history, external scientific references, evaluator material, or other agents were consulted. Supplied files were read without modification. General-purpose numerical dependencies were installed using uv in the permitted dedicated environment.

The calibration consists of 1,280 matching parent/offspring scaffolds, each 50,000 bases long, totaling 64,000,000 haploid sites. Matching names and lengths and the A/C/G/T alphabet were checked. There are 320,288 observed differences, a pooled change rate of 0.0050045. The new assembly contains twenty 50,000-base scaffolds. All query contigs and one-based coordinates were checked against that assembly. Input SHA-256 hashes are recorded in summary.json.

The substitution count matrix, with rows as parent A/C/G/T and columns as offspring A/C/G/T, is:

```
15924008    26808    26575    27173
   26484 15916557    26754    26648
   26533    26619 15913577    26624
   26637    26627    26806 15925570
```

## Hypotheses and chronological model development

1. **Stationary uniform substitutions.** The near-equal off-diagonal counts initially favored a pooled rate and equal alternate probabilities. Parent-specific mutation rates, transition/transversion bias, and the complete substitution matrix provided negligible held-out improvement.
2. **Marginal sequence, position, and scaffold effects.** Association checks at offsets −30 through +30, scaffold GC composition, relative scaffold position, and small coordinate periods showed little evidence of useful marginal predictors. The smallest offset association p-value was 0.00480, but 0.288 after correction for sixty offsets. Scaffold rate heterogeneity gave chi-square 1147.65 with 1279 degrees of freedom (p=0.9963); GC-rate correlation was 0.00958 (p=0.732). The 100-bin relative-position check gave p=0.470. These tests are exploratory, not proof of complete independence.
3. **High-order contexts.** Whole-context models contradicted the simple pooled hypothesis. Smoothed nine-base contexts improved prediction substantially despite flat marginal effects. Seven-base contexts improved modestly and three-base contexts did not help. Thus averaging over missing context bases hides the useful signal.
4. **Directionality and window extent.** Reverse-complement pooling worsened held-out predictions. The ten-base window −5..+4 outperformed −4..+5 and symmetric eleven-base windows. A scan of nine-base windows shifted from −6 through +6 peaked for windows −5..+3 and −4..+4, the two nine-base subwindows of the selected ten-base window. This is consistent with an oriented ten-base lookup effect. It does not identify the hidden implementation or prove that every possible longer-range dependence is absent.
5. **Sparse-context shrinkage.** Ordinary Dirichlet shrinkage, beta-distributed change rates with Dirichlet alternate proportions, and symmetric mixtures of background plus one to five hotspot classes were tested. A nonparametric binomial-mixture prior for each alternate probability performed best on both validation portions. This model was selected and refitted using all 64 million calibration sites.

`decision_log.md` records these decisions as they occurred, including the initial incorrect impression and the high-order-context finding.

## Validation design and evidence

A seed-0 permutation split whole scaffolds into 768 fitting scaffolds (38.4 million sites), 256 tuning scaffolds (12.8 million), and 256 assessment scaffolds (12.8 million). The split is recorded per scaffold in exploration.json. Whole-scaffold separation avoids using outcomes from the evaluation scaffold when building the context lookup table. Repeated contexts across separate scaffolds are intentional: learning those shared probabilities is the prediction task.

Both evaluation portions were inspected during iterative development. Consequently, these are **adaptive validation results, not an untouched final test**. No label from either validation portion was included in the original fitting counts for the comparisons below. Model architecture and model-family decisions did use the observed validation scores. The truly withheld new offspring remains unseen.

All log likelihoods below are natural-log sums; larger is better. Each evaluation portion has 12.8 million sites.

| Model fitted on the 60% portion | Tuning log likelihood | Second validation log likelihood |
|---|---:|---:|
| Pooled uniform substitutions | −474291.501 | −474003.866 |
| Full parent-to-offspring substitution matrix | −474292.641 | −474002.508 |
| Centered 7-mer, Dirichlet strength 10,000 | −473994.426 | −473661.300 |
| Centered 9-mer, Dirichlet strength 1,000 | −470752.890 | −470394.139 |
| Centered 11-mer, Dirichlet strength 300 | −471467.755 | −471098.983 |
| Ten-base −5..+4, beta/Dirichlet model | −465025.388 | −464297.060 |
| Ten-base −5..+4, five-hot-class mixture | −464530.997 | −463862.729 |
| **Ten-base −5..+4, selected binomial-mixture prior** | **−464501.205** | **−463801.610** |

The selected model gains {v['assessment_log_likelihood_gain']:.3f} nats over the pooled model on the second portion. The standard error of this gain, clustering by its 256 scaffolds, is {v['gain_standard_error_scaffold_clustered']:.3f} nats. Its multiclass Brier score is {v['assessment_multiclass_brier']:.10f}. Its expected change count on that portion is {v['assessment_expected_changes']:.2f}, versus {v['assessment_observed_changes']:,} observed changes. Brier score is strongly dominated by unchanged alleles, so held-out log likelihood and calibration are more informative for this comparison.

Calibration by predicted change-probability range on the second validation portion:

| Predicted change range | Sites | Expected changes | Observed changes |
|---|---:|---:|---:|
{table}

Pair counts of changes at distances 1–100 were close to the independent-process expectation of about 1,603 pairs per distance; the first five counts were 1598, 1560, 1640, 1696, and 1613. This is supportive of a mostly local context mechanism but does not prove site independence.

## Exact final estimator

Encode A/C/G/T as 0/1/2/3 and extract the length-ten parent word at offsets −5 through +4. Each word has four accumulated outcome counts. Express outcomes relative to the parent base: zero for unchanged and one/two/three for cyclic alternate alleles. Pool the three alternate count distributions when learning the prior, reflecting the approximately symmetric global substitution spectrum.

For each context and alternate allele let `n` be the calibration number of occurrences of that context and `x` the count of that alternate allele. Model the alternate probability `p` with a discrete prior on 100 grid points: zero plus 99 logarithmically spaced values from 0.00001 to 0.3. Estimate the shared prior weights by expectation maximization using the marginal binomial likelihood, grouping identical `(n, x)` observations for efficiency. There are at most 3,000 iterations and a marginal-log-likelihood stopping tolerance of 1e-7. The final full-data fit ran {pars['iterations']} iterations; hitting the iteration cap is not claimed to be exact convergence. All fitted grid points and weights are stored in final_model_parameters.json and final_model.npz.

For a requested alternate allele, compute the posterior mean

```
E[p | n, x] = sum_g grid[g] * weight[g] * grid[g]**x * (1-grid[g])**(n-x)
              / sum_g weight[g] * grid[g]**x * (1-grid[g])**(n-x).
```

Calculations use log-sum-exp for numerical stability, with an infinitesimal numerical substitute for log zero. Compute the three alternate posterior means separately and set the unchanged probability to one minus their sum. The independent binomial fits constitute a composite model: competing alternate outcomes are not literally independent multinomial events. The resulting four probabilities were verified positive and normalized for all learned contexts and for every output row. This approximation was selected because its held-out predictions improved over the tested joint categorical alternatives.

The final lookup is fitted to all calibration data, not just the fitting subset. Query positions are converted from one-based to zero-based only for indexing; output identifiers and positions are copied unchanged. Circular extraction uses modular shifts. For every query, the relative probabilities are mapped back to absolute A/C/G/T columns.

## Output verification

- Serialized predictions were reopened and checked against the original query table for exact order and identifier/position equality.
- All probabilities are finite and nonnegative; maximum row-sum error after text serialization is {s['maximum_row_sum_error_after_serialization']:.3g}.
- The minimum emitted probability is {s['min_probability']:.8g}, and the maximum is {s['max_probability']:.8g}.
- Expected changes were recomputed from the reread prediction file and actual new-parent alleles.
- Input and output SHA-256 hashes and software versions are preserved in summary.json; dependencies are pinned in requirements.lock.txt.
- No plots were generated or displayed, and no sequence data or predictions were sent outside this case directory.

## Remaining uncertainty

Individual offspring outcomes remain stochastic. An expected count is not a guarantee of the realized count. The model estimates context-conditioned marginals and does not provide a validated joint distribution across queried sites. Parameter uncertainty and adaptive selection uncertainty are not fully propagated into the output probabilities. In particular, extreme contexts are rare, and a probability grid capped at 0.3 is an approximation. The learned prior is not uniquely identifiable at very small probabilities from the available context counts; posterior predictive probabilities are the practical target. New contexts are treated through the learned prior rather than assigned zero risk. The approximate binomial factorization may miss dependencies among competing alternates, and the initial independence diagnostics do not rule out hidden long-range structure. Nevertheless, the large, replicated validation gain and calibration checks support the selected model over all tested simpler alternatives.

## Reproduction and retained files

All code and artifacts are in this case directory. Random sampling is seeded with 0. Numerical thread environment variables are set to four; no parallel worker pool was created. Arrays remain well below the 26 GB memory allowance. The supplied inputs are unchanged. The workflow is separated into idempotent phases, each overwriting only derived files:

```bash
%%bash
cd /mnt/d/agent_evaluations/case_006
uv venv --allow-existing /home/mew/.cache/case_006_env
uv pip install --python /home/mew/.cache/case_006_env/bin/python -r requirements.lock.txt
/home/mew/.cache/case_006_env/bin/python analyze.py
/home/mew/.cache/case_006_env/bin/python validate.py
/home/mew/.cache/case_006_env/bin/python context_probe.py
/home/mew/.cache/case_006_env/bin/python extended.py
/home/mew/.cache/case_006_env/bin/python validate11.py
/home/mew/.cache/case_006_env/bin/python mixture.py
/home/mew/.cache/case_006_env/bin/python asymmetric.py
/home/mew/.cache/case_006_env/bin/python shift_scan.py
/home/mew/.cache/case_006_env/bin/python beta_model.py
/home/mew/.cache/case_006_env/bin/python nonparam.py
/home/mew/.cache/case_006_env/bin/python finalize.py
/home/mew/.cache/case_006_env/bin/python write_report.py
```

The main production path after input inspection is `analyze.py`, `extended.py`, and `finalize.py`; the other scripts document the investigated alternatives and regenerate supporting evidence. For repeated production, the saved final_model.npz avoids refitting, though finalize.py deliberately recomputes from calibration for an auditable full rerun. No Git repository was accessed or created because this task was explicitly restricted to the supplied case directory.
'''
    (ROOT/'REPORT.md').write_text(report)
    with open(ROOT/'decision_log.md','a') as f:f.write(f'\n## 6. Completed output verification\nWrote {s["number_of_queries"]:,} probabilities in query order. Reopened predictions.tsv and verified coordinates, finite nonnegative probabilities, and normalization (maximum row-sum error {s["maximum_row_sum_error_after_serialization"]:.3g}). The expected changed-allele count from the serialized output is {s["expected_number_of_changed_alleles"]:.9f}. Saved full model, software pins, input/output hashes, detailed validation metrics, summary.json, and REPORT.md. No supplied input file was modified.\n')
    print('REPORT.md complete')
if __name__=='__main__':main()
