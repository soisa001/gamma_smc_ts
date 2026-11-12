import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'summary.json').read_text());a=json.loads((ROOT/'advanced_results.json').read_text());cv=json.loads((ROOT/'cv_results.json').read_text());d=json.loads((ROOT/'diagnostics.json').read_text())
rows=[]
for name,st in [('base',0),('symmetric_base',0),('cpg_unpooled',0),('cpg',0),('context3',1000000),('context5',2000),('context7',10000)]:
    z=next(x for x in cv if x['model']==name and x['strength']==st)
    rows.append(f"| {name}, CpG shrinkage {st:g} | {z['mean_log_loss']:.10f} |")
for name,st in [('tied5',0),('tied6',0),('motif4',0),('motif6',0),('motif8',0),('equal1',0),('equal3',0),('context5',50000),('context7',1000000)]:
    z=next(x for x in a['cv'] if x['model']==name and x['strength']==st)
    rows.append(f"| {name}, motif shrinkage {st:g} | {z['mean_log_loss']:.10f} |")
p=s['xor_outcome_probabilities'];names=['A/T, unequal bases at -2 and +2','A/T, equal bases at -2 and +2','C/G, outside CpG','C/G, in CpG']
probrows=['| '+names[i]+' | '+' | '.join(f'{v:.12f}' for v in p[i])+f" | {s['query_group_exposures'][i]:,} |" for i in range(4)]
diagrows=[f"| {k} | {v['pearson_stat']:.2f} | {v['df']} | {v['dispersion']:.4f} | {v['p_approx']:.4g} |" for k,v in s['residual_diagnostics'].items()]
report=f'''# Offspring allele probability model

The output predicts **{s['expected_total_changed_alleles']:.6f} changed alleles** across **{s['query_count']:,} queried sites**. `predictions.tsv` contains the requested six columns in exactly the query order. The final model uses a CpG transition effect and an A/T motif defined by equality of the bases two positions to either side, with reverse-complement symmetry and shared background substitution probabilities. It estimates allele probabilities; it does not pretend to identify which individual stochastic events will happen.

## Inputs, scope, and safeguards

Only the six supplied files were used for task information and data. No repository history, other case directories, hidden outcomes, other agents, or external scientific/reference sources were inspected. General-purpose Python dependencies were installed using uv. All supplied inputs were kept unchanged and their SHA-256 digests verified after prediction generation; see `provenance.json`. Code and outputs stay in this case directory, apart from the explicitly authorized dedicated Python environment.

Calibration contains 80 homologous contig pairs, each 50,000 bases long. All sequences contain only A/C/G/T and homologous lengths agree. The assemblies differ at 20,138 of 4,000,000 positions (0.0050345). The new parent has 20 contigs of 50,000 bases. Queries cover 999,880 positions with complete seven-base windows. For consistent comparison of candidate models, fitting and validation used 3,999,520 calibration positions, excluding three bases at each contig end. All 20,138 observed changes remain in this subset. Matching scaffold identifiers across the calibration and new assembly were not treated as sample-specific predictors.

## Hypotheses and empirical evidence

The initial model was a parent-base-specific categorical distribution, including the unchanged allele. Directional mutation counts suggested reverse-complement symmetry. C followed by G, and G preceded by C, showed change rates around 0.011, compared with approximately 0.0032 for C/G outside that context. This motivated the CpG model.

Single positions farther from the center showed no strong marginal signal. However, five-base contexts had much larger residual heterogeneity than expected from CpG alone: the pooled 5-mer Pearson dispersion was 3.851, versus 1.170 for 3-mers. The largest residual contexts shared equal bases at offsets -2 and +2. For A/T centers, this equality changed the observed rate from 0.0031794 to 0.0099575. In A-oriented coordinates it increased A-to-C and A-to-G, while A-to-T stayed near its background probability. Equality at offsets +/-1 or +/-3 did not produce a comparable signal. This is a statistical rule inferred from calibration; it is not an assertion about a specific biochemical mechanism.

I tested parent-base models; strand-pooled and unpooled CpG models; 3-, 5-, and 7-mer categorical models with shrinkage; a four-group motif model; a six-group version allowing the equality effect in C/G; strand-specific motif models; extra equality indicators at +/-1 and +/-3; and two models that share baseline substitution rates across groups. The five-parameter model ties the background transversions; the six-parameter alternative allows them to differ.

## Validation design and results

Five-fold cross-validation assigns whole contigs to folds using `numpy.random.default_rng(0)`. Each fold contains 16 contigs. Every candidate is fitted without that fold's offspring outcomes. The proper scoring rule used for selection is mean negative log probability of the actual offspring allele, in natural-log units. Lower is better. Accuracy would be misleading because predicting the parent allele everywhere already succeeds at about 99.5% of sites.

| Candidate | Held-out mean log loss |
|---|---:|
{chr(10).join(rows)}

The five-parameter model's cross-validated Brier score is 0.0100146781 (sum over the four allele outcomes, then mean over sites). Its held-out expected count is 20,138.34 versus 20,138 observed. This pooled agreement is reassuring but, because each fold estimates similar class frequencies, is weaker evidence than the comparative log scores. Fold-by-group expected and observed counts are provided in `advanced_results.json`.

The six-parameter model has the best mean log score, improving on the five-parameter model by 0.0000004502 per site. The paired standard error across 80 contigs is 0.0000005729, so the simpler model is within one paired standard error. A heavily regularized 5-mer model is also statistically close. I selected five parameters to reduce variance and avoid unsupported context-specific variation. A likelihood-ratio comparison of the full four-group model against five parameters gives p=0.131. Allowing unequal background transversions gives nominal p=0.0198, so the equality assumption has some uncertainty; this marginal distinction has negligible influence on the total expected count. No probabilities were rounded to presumed generating constants.

**New-parent context shift:** equal flanks occur at 398,545 of 498,027 A/T query sites (80.025%), compared with 500,124 of 1,998,519 calibration A/T sites (25.025%). Thus applying the calibration average change rate directly would be inappropriate. The fitted motif model explicitly accounts for this increased exposure. As a further check, I reweighted held-out losses to the new parent's four-group frequencies. The five-parameter loss is 0.0460042584, versus 0.0460038260 for six parameters; the difference 4.32e-7 is smaller than its paired-contig standard error 7.27e-7. Splitting groups by the exact base at offsets -2, -1, +1, -3, or +3 worsens weighted held-out loss, as does the unconstrained four-group model. See `shift_diagnostics.json`. This check supports retaining the same model under the observed shift without using any withheld offspring outcome.

**Selection limitation:** exploratory motif discovery used the entire calibration set, and the same folds were used to compare several candidate models. Therefore the scores are cross-validated parameter-estimation results, not a completely untouched confirmatory test of the discovery procedure. The wide margin over base/CpG models and simplicity of the rule support it, but a new independent calibration batch would be stronger validation. Nominal tests are descriptive and not adjusted for all exploratory comparisons.

## Final probability model

Encode A=0, C=1, G=2, T=3. Define an outcome class by XOR of the parent and offspring codes: 0 is unchanged; 1 is a non-complement transversion (A/C or G/T); 2 is a transition (A/G or C/T); 3 is a complement transversion (A/T or C/G). This representation automatically preserves reverse-complement symmetry. CpG means parent C immediately followed by G, or parent G immediately preceded by C.

The final parameters maximize the multinomial likelihood of the calibration group/outcome counts, using analytic gradients. A Newton fallback with an exact Hessian handles rare floating-point line-search failures during bootstrap; it requires a small likelihood gradient. Each unchanged probability is one minus the three mutation probabilities. The five free probabilities are background transversion, background transition, motif non-complement transversion, motif transition, and CpG transition.

| Parent context | Unchanged | XOR 1 | XOR 2 | XOR 3 | Query sites |
|---|---:|---:|---:|---:|---:|
{chr(10).join(probrows)}

To emit `p_A`, `p_C`, `p_G`, `p_T`, map each candidate allele through XOR with the actual parent base. The two background rows share parameters. All site predictions use parent sequence alone; neighboring offspring outcomes are never used.

After the motif is accounted for, residual diagnostics are:

| Feature partition | Pearson statistic | Approximate df | Dispersion | Approximate p |
|---|---:|---:|---:|---:|
{chr(10).join(diagrows)}

These are descriptive diagnostics. Sparse seven-base classes have low expected mutation counts, so asymptotic chi-square p-values there are approximate. Richer-context cross-validation is the stronger check. Contig and position comparisons do not establish perfect homogeneity, and no finite calibration can rule out weak unmodeled effects. Neighbor-change frequencies at lags 1,2,3,4,5,10,20,100 were broadly near the marginal rate, with no compelling large clustering signal. Marginal allele predictions and their expected total do not require independent offspring changes.

The residual five-base dispersion remains modestly elevated (1.133, nominal p=0.000203), which is evidence against claiming an exact description of every context effect. The tested richer models did not provide a reliable out-of-sample gain, so I retain the compact rule and report this as remaining model uncertainty. The bootstrap interval below does not cover such potential misspecification.

## Expected count and uncertainty

For each query i with parent allele r_i, the expected change probability is `1 - p_i(r_i)`. Their sum is **{s['expected_total_changed_alleles']:.9f}**. The expected unchanged count is {s['expected_unchanged_alleles']:.9f}.

A 500-replicate whole-contig bootstrap, seed 0, refits the fixed five-parameter model. The estimated standard error of the expected count is **{s['parameter_uncertainty_bootstrap_se']:.3f}**, with percentile 95% interval **[{s['parameter_uncertainty_bootstrap_95_percent_interval'][0]:.3f}, {s['parameter_uncertainty_bootstrap_95_percent_interval'][1]:.3f}]**. This describes uncertainty in the estimated expectation conditional on model choice; it does not include all model-discovery uncertainty.

Separately, if offspring site outcomes are conditionally independent, the realized number of changes would have standard deviation **{s['conditional_independent_site_count_sd']:.3f}** around its expectation with fixed model parameters. This is process variation, not the bootstrap uncertainty. A combined predictive interval is deliberately not asserted because unknown dependence and model error are not fully identified by one offspring per calibration parent.

## Correctness and reproducibility

Validation re-read the serialized TSV and checked the exact column names, all {s['query_count']:,} coordinate pairs in original query order, finite nonnegative probabilities, and row sums. The largest row-sum deviation after decimal serialization is {s['max_roundtrip_row_sum_error']:.3g}. Recomputing the expected count from the written probabilities gives {s['expected_changes_from_serialized_predictions']:.9f}, agreeing with the in-memory result. `summary.json` stores these checks, fitted probabilities, exposure counts, software versions, input/output provenance references, and residual diagnostics. All sampling uses seed 0. Numeric libraries and bootstrap worker count are set to four; data fit comfortably within 26 GB. No plots were produced or displayed.

The idempotent workflow is split into exploration, validation, and prediction phases. Run from this directory in WSL/Linux; dependencies are version-pinned and installed by uv:

```bash
%%bash
cd /mnt/d/agent_evaluations/case_001
bash run.sh explore
bash run.sh validate
bash run.sh predict
/home/mew/.cache/case_001_env/bin/python shift_check.py > shift_output.txt
/home/mew/.cache/case_001_env/bin/python write_report.py
```

The standalone `%%bash` line is a notebook cell magic and should be omitted when pasting into a terminal. For prediction-only reruns after the supplied exploration provenance exists, `bash run.sh predict` is sufficient. No repository clone/push step is provided because this blind task supplied no authorized repository and explicitly restricted access to this directory.

Primary deliverables: `predictions.tsv`, `summary.json`, this `REPORT.md`, and chronological `decision_log.md`. Supporting scripts are `explore.py`, `model.py`, `diagnostics.py`, `motif.py`, `advanced.py`, `finalize.py`, `write_report.py`, and `run.sh`. Machine-readable comparison results and logs are retained for review.
'''
(ROOT/'REPORT.md').write_text(report)
print('Wrote REPORT.md')
