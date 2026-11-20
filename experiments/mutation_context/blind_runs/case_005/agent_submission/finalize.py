"""Fit the selected pooled model, export probabilities, and audit every row."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='4'
import csv, hashlib, itertools, json, math, platform
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy, sklearn, lightgbm
from explore import fasta, encoded, ROOT

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

explore=json.loads((ROOT/'exploration.json').read_text())
diagnostics=json.loads((ROOT/'diagnostics.json').read_text())
validation=json.loads((ROOT/'validation.json').read_text())
spectra=json.loads((ROOT/'validation_spectra.json').read_text())
for name, metadata in explore['inputs'].items():
    assert sha(ROOT/name)==metadata['sha256'], f'Input changed: {name}'
p=encoded(fasta('calibration_parent.fasta')); o=encoded(fasta('calibration_offspring.fasta')); new=encoded(fasta('new_parent.fasta'))
assert p.keys()==o.keys()
n=sum(len(v) for v in p.values()); nchanged=sum(int(np.count_nonzero(v!=o[k])) for k,v in p.items())
rate=nchanged/n
assert n==explore['n'] and nchanged==explore['changed']
assert all(v.max()<4 for v in new.values())
prob=np.full((4,4),rate/3); np.fill_diagonal(prob,1-rate)
seen={k:np.zeros(len(v),bool) for k,v in new.items()}
source_counts=np.zeros(4,np.int64); count=0; duplicates=0
tmp=ROOT/'predictions.tsv.tmp'
with (ROOT/'queries.tsv').open(newline='') as fi, tmp.open('w',newline='') as fo:
    reader=csv.DictReader(fi,delimiter='\t'); assert reader.fieldnames==['contig','position']
    writer=csv.writer(fo,delimiter='\t',lineterminator='\n'); writer.writerow(['contig','position','p_A','p_C','p_G','p_T'])
    for row in reader:
        key=row['contig']; pos=int(row['position']); assert key in new and 1<=pos<=len(new[key])
        base=int(new[key][pos-1]); source_counts[base]+=1
        duplicates+=int(seen[key][pos-1]); seen[key][pos-1]=True
        writer.writerow([key,row['position']]+[format(v,'.17g') for v in prob[base]])
        count+=1
os.replace(tmp,ROOT/'predictions.tsv')

# Re-read the serialized output independently, paired with the original query file.
max_error=0.; audited=0; change_by_base=[0.,0.,0.,0.]; pmin=1.; pmax=0.
with (ROOT/'queries.tsv').open(newline='') as fq, (ROOT/'predictions.tsv').open(newline='') as fp:
    queries=csv.DictReader(fq,delimiter='\t'); predictions=csv.DictReader(fp,delimiter='\t')
    assert predictions.fieldnames==['contig','position','p_A','p_C','p_G','p_T']
    for q,r in itertools.zip_longest(queries,predictions):
        assert q is not None and r is not None
        assert (q['contig'],q['position'])==(r['contig'],r['position'])
        ps=[float(r['p_'+b]) for b in 'ACGT']
        assert all(math.isfinite(v) and 0<=v<=1 for v in ps)
        err=abs(math.fsum(ps)-1); max_error=max(max_error,err); assert err<1e-12
        b=int(new[q['contig']][int(q['position'])-1])
        assert np.allclose(ps,prob[b],rtol=0,atol=1e-15)
        change_by_base[b]=1-ps[b]; pmin=min(pmin,min(ps)); pmax=max(pmax,max(ps)); audited+=1
assert count==audited
expected=count*rate
expected_from_export=math.fsum(int(source_counts[b])*change_by_base[b] for b in range(4))
assert abs(expected_from_export-expected)<1e-7
for name, metadata in explore['inputs'].items(): assert sha(ROOT/name)==metadata['sha256']

allmodels=validation['models']+spectra['models']
comparison=[]
for r in allmodels:
    base=next(b for b in validation['models'] if b['scheme']==r['scheme'] and b['model']=='pooled')
    delta=np.array(r['group_nll'])-np.array(base['group_nll'])
    se=float(delta.std(ddof=1)*math.sqrt(len(delta)))
    comparison.append({k:r[k] for k in ['scheme','model','prior_strength','log_loss','nll_delta_vs_pooled','brier']}|{'paired_scaffold_se_total_nll':se})
(ROOT/'model_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
rate_ci=diagnostics['rate_ci95']
pred_sd=math.sqrt(count*rate*(1-rate)*(1+count/n))
summary={
    'model':'Homogeneous Bernoulli mutation rate; uniform among the three non-parent alleles',
    'seed':0,'threads':4,'memory_limit_gb':26,
    'calibration_sites':n,'calibration_changed_alleles':nchanged,'estimated_mutation_probability':rate,
    'probability_parent_allele':1-rate,'probability_each_nonparent_allele':rate/3,
    'query_count':count,'unique_query_count':sum(int(v.sum()) for v in seen.values()),'duplicate_query_rows':duplicates,
    'expected_changed_alleles':expected,
    'expected_changed_alleles_from_serialized_predictions':expected_from_export,
    'mutation_probability_95_percent_binomial_ci':rate_ci,
    'expected_changed_alleles_95_percent_parameter_ci':[count*v for v in rate_ci],
    'approximate_realized_count_95_percent_predictive_interval':[expected-1.96*pred_sd,expected+1.96*pred_sd],
    'uncertainty_assumptions':'Homogeneous independent changes and unchanged process; predictive interval adds approximate rate estimation variance. These intervals do not cover model misspecification.',
    'parent_base_counts_queried':dict(zip('ACGT',source_counts.tolist())),
    'audit':{'rows_checked':audited,'query_order_preserved':True,'all_probabilities_finite_nonnegative':True,'maximum_absolute_row_sum_error':max_error,'minimum_probability':pmin,'maximum_probability':pmax,'six_original_inputs_unchanged':True},
    'predictions_file':'predictions.tsv','predictions_sha256':sha(ROOT/'predictions.tsv'),
    'inputs':explore['inputs'],
    'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'scikit_learn':sklearn.__version__,'lightgbm':lightgbm.__version__},
}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

rows=[]
keys=list(dict.fromkeys((r['model'],r['prior_strength']) for r in allmodels))
for model,a in keys:
    r1=next(r for r in comparison if r['model']==model and r['prior_strength']==a and r['scheme']=='scaffold')
    r2=next(r for r in comparison if r['model']==model and r['prior_strength']==a and r['scheme']=='contiguous_region')
    rows.append(f"| {model} | {a if a is not None else '—'} | {r1['nll_delta_vs_pooled']:.3f} | {r2['nll_delta_vs_pooled']:.3f} | {r1['paired_scaffold_se_total_nll']:.3f} |")
base=next(r for r in allmodels if r['model']=='pooled' and r['scheme']=='scaffold')
report=f'''# Offspring allele prediction report

## Result

The selected model assigns probability **{1-rate:.6f}** to the parent's allele and **{rate/3:.15f}** to each of the other three alleles at every queried site. The expected number of changed alleles across **{count:,}** query rows is **{expected:,.6f}**. `predictions.tsv` contains the requested six columns in the original query order. `summary.json` contains the estimate, model parameters, uncertainty estimates, hashes and audit results.

The estimate is learned from the calibration data; it is not rounded to a conjectured synthetic generating parameter. The withheld offspring was never accessed. This is a probability forecast: individual realized mutations cannot be reliably identified under the selected model.

## Scope and input integrity

Only the six supplied files in this case directory were used. No other project files, repository history, hidden outcomes, evaluator code, chats, agents, scientific references or external datasets were consulted. General-purpose Python dependencies were installed with uv in the permitted dedicated environment. All generated scripts and result files are in this case directory; the environment itself is at `/home/mew/.cache/case_005_env`.

The calibration has 80 homologous parent/offspring scaffold pairs, each 50,000 bases long, totaling {n:,} sites. The new parent has 20 scaffolds of 50,000 bases, totaling 1,000,000 sites. Homologous names and lengths were checked. All sequences contain only A, C, G and T. Assemblies were compared directly at matching coordinates, as warranted by the task's error-free, single-nucleotide-change assumptions. The new assembly is treated as a separate parent: repeated scaffold names do not imply matching sequence or offspring outcomes across calibration and prediction assemblies.

There are {nchanged:,} calibration changes. Queries have {count:,} rows, {summary['unique_query_count']:,} unique coordinates, and {duplicates} duplicate rows. Every coordinate is valid. All new-parent sites are queried exactly once. SHA-256 checks before and after output generation verify that all six supplied inputs remained unchanged; their hashes appear in `exploration.json` and `summary.json`.

## Hypotheses and empirical evidence

The leading hypothesis is a common mutation probability with a uniform replacement among the three alternative nucleotides. Alternatives considered were source-base-specific risk, GC versus AT risk, local motifs and reverse-complement motifs, sequence composition, short homopolymers, coordinate-dependent rates, periodic coordinate effects, scaffold heterogeneity, mutation clustering and nonuniform/context-dependent replacement spectra.

The calibration transition counts (rows parent, columns offspring, in A/C/G/T order) are:

| Parent | A | C | G | T |
|---|---:|---:|---:|---:|
'''+'\n'.join('| '+b+' | '+' | '.join(str(v) for v in row)+' |' for b,row in zip('ACGT',explore['spectrum']))+f'''

The overall change rate is {nchanged}/{n} = {rate:.8f}. Source-base rates show a weak nominal association (chi-square 8.778, 3 df, p=0.0324), but this is one of many examined hypotheses. The combined test of all 12 directed substitution rates against the common-rate uniform-replacement model gives chi-square 14.244 on 11 df, p=0.220. Symmetric replacement classes A↔C/G↔T, A↔G/C↔T, and A↔T/C↔G have counts 6,585, 6,628 and 6,619 (p=0.925 against equality). GC-versus-AT risk has p=0.981.

Scaffold change counts range from 207 to 287, averaging 247.9. Their observed variance is 237.84; a common-rate binomial model gives heterogeneity chi-square 76.17 on 79 df, p=0.569. Ten equal coordinate bins have p=0.444. Adjacent identical-base counts have p=0.817. Coarse GC composition bins in a 1,001-base circular window have p=0.657. These tests support pooling but do not prove absence of every possible effect.

Circular offsets -30 through +30 were checked separately, excluding zero. Their source-stratified mutation-risk association statistics show no compelling signal after multiple comparisons. Among mutations, destination-base association with each offset was also tested conditional on the source base; the smallest p-value was 0.0255 at +15, which does not survive correction over 60 offsets. Positional moduli 2 through 40 were screened. Mutation pair counts at lags 1, 2, 3, 5, 10, 50, 100 and 1,000 were 96, 94, 100, 77, 110, 91, 85 and 80, against approximately 98.33 expected at each under independence. These exploratory summaries do not suggest a substantial burst process. Pair counts are diagnostics only and do not use unavailable offspring information as prediction features.

## Validation design

Every calibration site received an out-of-fold prediction in each of two five-fold validation schemes:

1. **Whole scaffolds:** seed-0 permutation assigns 16 entire scaffolds to each fold; each fit uses the other 64. This tests transfer to unseen scaffolds.
2. **Contiguous regions:** each 50,000-base scaffold is divided into five consecutive 10,000-base regions; each fold holds out the corresponding region of every scaffold. This checks sensitivity to positional structure.

No held-out offspring allele contributes to training counts or priors. Context uses parent sequence only; circular windows wrap within the same scaffold and never across scaffolds. Parent contexts can span a held-out boundary because the whole parent assembly is available at prediction time. Initial descriptive screening used all calibration pairs, and the two CV schemes reuse the same outcomes. These are internal model-comparison checks, not independent confirmatory experiments or performance estimates on the withheld new offspring.

The primary metric is four-allele negative log likelihood (NLL), including the parent allele. Brier score is also recorded. Accuracy would be dominated by unchanged sites and was not used for model selection. The pooled model's scaffold-CV log loss is {base['log_loss']:.12f} per site and Brier score is {base['brier']:.12f} (sum over four classes).

Candidate categorical risk models estimate the mutation probability in each category and divide it equally among alternatives. A prior of strength alpha is centered on the training fold's pooled rate: (mutations + alpha × pooled rate)/(exposure + alpha). Full-spectrum models estimate the four outcomes coded by parent XOR offspring, with a four-category Dirichlet prior centered on the same homogeneous model. Reverse-complement context categories are canonicalized; XOR replacement classes are invariant under simultaneous complementation. Global source-spectrum models use a 0.5 pseudocount per cell. Global XOR-spectrum models use a 0.5 pseudocount per alternative; transition-spectrum models pool the two transversion classes. Source-rate-plus-XOR-spectrum models combine independently estimated source-specific risks and a shared replacement spectrum. Probabilities are clipped only where necessary for numerical safety in validation.

Risk candidates cover source/GC classes, immediate dinucleotides, full and reverse-complement 3-, 5- and 7-mers, 10/50/100 positional bins, selected screened positional moduli, GC windows of 11/51/201/1,001 bases, and zero/one/two identical immediate neighbors. Full-spectrum candidates cover immediate dinucleotides, 3- and 5-mers, and reverse-complement 7-mers. Alpha is an exposure pseudocount, not a count of mutations. All alternatives and smoothing settings evaluated are listed below; none were silently discarded.

## Complete model comparison

Values are total held-out NLL **minus pooled NLL** across all four million sites. Positive values are worse. The last column is an approximate standard error of the scaffold-split total NLL difference, computed as sqrt(80) × the sample standard deviation of paired per-scaffold NLL differences. It describes heterogeneity of the validation difference and is not a formal independent-CV uncertainty estimate because training folds overlap.

| Candidate | Alpha | Scaffold ΔNLL | Region ΔNLL | Paired scaffold SE |
|---|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

The small source-rate gains (about 2–3 total nats) are smaller than their approximately 3.1-nat paired standard errors. Strongly shrunk reverse-complement 7-mer risk models gain 6.61 nats in scaffold CV and 3.70 nats in region CV; their respective paired standard errors are 3.77 and 3.95 nats, with approximate 95% intervals including zero. The full-spectrum version also has a small gain. These gains are tiny after screening many candidate models and smoothing strengths; there is no stable, identifiable context rule. Less-shrunk context models markedly overfit, and direct destination-offset tests provide no corroboration. A conservative common-rate model is therefore selected. The choice explicitly favors parsimony and stability over the absolute minimum internal CV score. The best context model is not claimed to have been disproved.

## Final fitted model and expected count

Fit the chosen model on all calibration sites:

```
q = number of observed changes / number of calibration sites
  = 19,832 / 4,000,000 = 0.004958
P(offspring = parent) = 1 - q = 0.995042
P(offspring = b) = q / 3 for each b different from the parent
E[number changed over query rows] = sum_i q = {expected:.6f}
```

Although the estimated hazard is constant, the placement of the high probability among p_A/p_C/p_G/p_T depends on the observed new-parent allele at each queried coordinate. Probabilities are output with 17 significant digits. The output expected count is the sum of non-parent probability over the original query rows; with this query set every physical site is counted once.

The 95% exact binomial interval for q is [{rate_ci[0]:.10f}, {rate_ci[1]:.10f}]. Multiplying by the query count gives [{count*rate_ci[0]:.3f}, {count*rate_ci[1]:.3f}] for uncertainty in the expected count under the homogeneous model. A normal predictive interval for the realized number of changes, adding approximate rate-estimation variance to binomial outcome variance, is [{expected-1.96*pred_sd:.1f}, {expected+1.96*pred_sd:.1f}]. This is distinct from the expected number and assumes independent changes and a stable process. It is not a guarantee about the withheld outcome.

## Output audit and reproducibility

`finalize.py` independently re-read every exported row alongside `queries.tsv`, checked exact coordinate and order agreement, checked the six column names, checked finite nonnegative probabilities bounded by one, checked row sums, recomputed the expected number of changes from serialized probabilities, and rechecked all original input hashes. All {audited:,} rows passed. Maximum absolute row-sum error was {max_error:.3g}. The expected count recomputed from the file was {expected_from_export:.12f}. The final prediction SHA-256 is `{summary['predictions_sha256']}`.

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
'''
(ROOT/'REPORT.md').write_text(report)
print(json.dumps({k:summary[k] for k in ['query_count','expected_changed_alleles','audit','predictions_sha256']},indent=2))
