# Posterior calling rules and delayed selection: completed results

Changing the hard-call definition does not resolve the low-power problem for
s=0.005. With selection beginning at 50 kya and a nominal 5% regional error
target, the highest observed decoded sensitivity across the tested regional
maximum scores and time cutoffs is 21%; true TMRCA reaches 45%. With selection
beginning at 10 kya, both decoded and true-TMRCA regional maxima have little
discrimination: their highest observed sensitivities are at most 8% at the same
nominal error level. These are descriptive extrema across the tested methods,
not independently confirmed winning detectors.

All 100 new later-onset selected simulations and all 1,200 posterior replays
completed without failures. The earlier fresh 400-haplotype study supplies
1,000 neutrals and 100 selected replicates with onset at 50 kya. No old handoff
simulation data were added. The independent audit passed, recomputing all 960
regional summaries from 528,000 held-out predictions and reproducing all 192
previous mean/truth baseline rows exactly.

## Design and interpretation

Both selected arms use EAS, 400 haplotypes, 10,000 fixed unique pairs,
s=0.005, dominance 0.5, and a 2% introgression pulse at 50 kya. Mutation rate
1.25e-8, recombination rate 1e-8, generation time 25 years, decoder parameters,
and the constant-Ne flow field remain unchanged. The new arm starts selection
at 10 kya. It is an illustrative timing experiment, not a fitted EPAS1 model.

Every statistic is `frac_recent_T`, a fraction of pairs passing a hard criterion:
posterior mean below T, posterior median at most T, posterior mass below T at
least 0.8 (`prob80`), or posterior mass below T at least 0.9 (`prob90`). True
pairwise TMRCA supplies a separate reference. None uses `mean_p_lt_T`. Cutoffs
are 5, 10, 20, 30, 40 and 50 kya.

A qualifying call anywhere in the isolated 10 Mb region counts. The analysis
retains all selected replicates and permits fixation. Selected alleles must be
observed; neutrals align to the nearest segregating site without an AF floor.
These are sensitivity estimates conditional on allele observation. The two
onset cohorts have reproducible but distinct genealogies; the four posterior
rules are paired on identical input sequences and haplotype pairs.

Five folds split whole regions. In each fold, 400 neutrals fit normalization,
another 400 calibrate the regional score, and 200 are held out. Selected regions
split 80 for threshold training and 20 held out. The same folds are used across
rules. The shared neutral panel is not an additional independent panel for
each onset. Each rule and cutoff receives separate null calibration.

The neutral fractions below are measured regional false-positive probabilities
under this simulation model. They are not estimates of the fraction of real
genomic discoveries that would be false; that also depends on the prevalence
of selection. The matched simulation calibration addresses fixed-Ne decoder
misspecification for these null call rates. It does not restore information
lost in decoding or establish genome-wide error control.

## Regional maxima at a low neutral-call rate

This table shows the highest observed selected sensitivity across the six
cutoffs for each rule using the **raw regional maximum**, with nominal regional
alpha=0.05. Ties are resolved by the lower neutral rate and then the lower T.
The choice of T in this summary is exploratory. Rate cells are **selected /
neutral**, with denominators 100 / 1,000.

| Rule | Onset 50 kya: T | Onset 50 kya: rates | Onset 10 kya: T | Onset 10 kya: rates |
|---|---:|---:|---:|---:|
| Mean | 20 kya | 19% / 5.5% | 30 kya | 8% / 5.5% |
| Median | 10 kya | 18% / 5.8% | 30 kya | 7% / 5.2% |
| P(TMRCA < T) >= 0.8 | 10 kya | 18% / 4.6% | 50 kya | 7% / 5.3% |
| P(TMRCA < T) >= 0.9 | 50 kya | 19% / 5.4% | 50 kya | 6% / 5.4% |
| True TMRCA | 50 kya | 45% / 5.5% | 10 kya | 8% / 4.3% |

Position standardization does not change the conclusion. For onset at 50 kya,
its highest observed decoded sensitivities are 18% for mean, 18% for median,
20% for prob80 and 21% for prob90. The latter uses T=50 kya and calls 5.6% of
neutrals. For onset at 10 kya, standardized decoded maxima call at most 7% of
selected regions at the nominal 5% level. With only 100 selected replicates,
differences of a few percentage points are not persuasive evidence of a general
improvement, especially after comparing several cutoffs.

At nominal regional alpha=0.10, raw maxima reach at most 29% decoded sensitivity
for onset at 50 kya, compared with 47% for true TMRCA. The later-onset raw maxima
remain near background, reaching at most 13% decoded sensitivity and 12% for
true TMRCA. The [complete tables](results/recent_calls/tables.md) show all six
cutoffs, including unfavorable results.

## What happens near the 70% sensitivity target

Thresholds trained to call at least 70% of training selected regions produce
the following held-out examples for selection onset at 50 kya:

| Rule and regional score | T | Selected called | Neutral called |
|---|---:|---:|---:|
| Mean, raw maximum | 10 kya | 68% | 54.0% |
| Median, raw maximum | 5 kya | 70% | 50.8% |
| Prob80, raw maximum | 10 kya | 69% | 48.8% |
| Prob80, standardized maximum | 10 kya | 71% | 46.1% |
| True TMRCA, raw maximum | 50 kya | 68% | 25.9% |

The prob80 standardized result has the lowest observed neutral rate among the
decoded trained-threshold comparisons, so it is an optimistic exploratory
summary. It still calls almost half of neutral regions. For onset at 10 kya,
the corresponding most favorable decoded example is 69% selected / 62.0%
neutral. The experiment does not establish a useful 70%-power detector at a
low neutral-region call rate.

## Stronger consecutive evidence

At T=50 kya, stride p<=0.001, 10 kb reporting strides, and five consecutive
significant strides with no gaps, onset-50k results are:

| Rule | Selected called | Neutral called |
|---|---:|---:|
| Mean | 17% | 3.4% |
| Median | 7% | 1.8% |
| Prob80 | 22% | 5.7% |
| Prob90 | 22% | 7.6% |
| True TMRCA | 34% | 3.7% |

Requiring ten consecutive strides lowers neutral calls further, but selected
calls fall to 3%, 1%, 8%, 10% and 15%, respectively. At T=10 kya, the mean
five-stride rule gives 25% selected / 7.7% neutral; ten strides give 17% / 3.0%.
Thus a larger selected rate under a posterior-mass rule can accompany a larger
neutral rate. The full grid contains 540 rules per source and onset. None,
including the truth rules, reaches both at least 70% selected and at most 5%
neutral calls. This is an exploratory grid, not independent confirmation of
the best entry. Its p=0.001 pointwise tail is resolution-limited by 1,000 neutrals.

## Why stricter posterior calls are not enough

At T=50 kya, changing mean to prob90 reduces the neutral SD across replicates,
averaged over positions, from 0.08219 to 0.03273. However, the onset-50k focal
selected-minus-neutral mean also falls from 0.23995 to 0.11402. Stricter calls
reduce signal as well as variation. Matched regional calibration exposes the
remaining overlap between selected and neutral regions.

Different rules also shift the cutoff where the decoded signal appears.
At young T, a decoded statistic can have more detection power than the
true-TMRCA statistic at that same T while being less accurate about ages.
True TMRCA with one specified summary and cutoff is a reference, not an
information-theoretic upper bound on every alternative detector.

The later-onset cohort is substantially less swept. Mean selected AF falls
from 0.783925 for onset at 50 kya to 0.212400 for onset at 10 kya; medians are
0.82125 and 0.17250. In the later cohort, 31/100 have AF below 0.1 and 54/100
below 0.2; none in either cohort is sample-fixed. This comparison changes both
time under selection and the resulting allele-frequency distribution. The
weak true-TMRCA regional results show that decoding is not the only limitation
of this statistic in the later-onset conditions.

## Flow-field investigation

A fixed demographic history Ne(t) can support a two-dimensional flow table,
but replacing the current table alone would be inconsistent. The current
constant-Ne prior also determines initialization, the algebra combining
forward and reverse filters, and parts of clipping and posterior summaries.
The source construction projects the infinitesimal SMC-prime density change
onto the gamma family's two tangent directions by weighted least squares;
it does not simply match mean and variance.

The [flow-field research note](variable_ne_flow_field.md) derives the
variable-history transition kernel and proposes a two-parameter family
q_N(t;a,b) proportional to pi_N(t) t^a exp(-bt), where pi_N is the demographic
coalescent prior. It retains exact prior initialization, simple mutation
updates, and prior-corrected filter combination within the family, and reduces
to ordinary gamma when Ne is constant. Recombination still requires projection.

The algebra/quadrature checks passed with absolute errors below 2e-15 at the
tested points. An independent source/table check also found a recombination
factor-of-two convention discrepancy. The note documents an entropy-clipping
heuristic that can exclude valid exact posteriors. Neither observation has been
shown to cause the measured power loss. A demographic decoder has not yet been
implemented or benchmarked; the current experiment keeps the original binary
and field unchanged. A useful next step would be an accurate small SMC-prime
reference, followed by a validated demographic family and field construction.
A scalar Ne(t) would still approximate the structured archaic-pulse model.

## Files, verification, and rerun

- [Complete comparison tables](results/recent_calls/tables.md).
- [Regional summaries](results/recent_calls/analysis/region_calls.csv),
  [stride/run grid](results/recent_calls/analysis/paired_grid.csv), and
  [diagnostics](results/recent_calls/analysis/diagnostics.csv).
- [Nominal 5% figure, vector PDF](results/recent_calls/analysis/power_neutral_0.05.pdf)
  and [nominal 10% figure, vector PDF](results/recent_calls/analysis/power_neutral_0.1.pdf).
  PNG versions are saved alongside them. Figures are not displayed in chat.
- [Independent audit](results/recent_calls/analysis_audit.json),
  [analysis provenance](results/recent_calls/analysis/provenance.json), and
  [verified artifact inventory](results/recent_calls/artifact_manifest.json).

Sixteen focused tests passed, including native hard-call comparisons with
direct gamma CDF calculations, raw-posterior replay, padded pair lanes,
truncated-stream rejection, delayed-onset events, and held-out evaluation.
All 1,200 native and replayed mean summaries reproduce their original integer
counts. The full comparison passes posterior-mass nesting checks. The actual
generated SLiM event arrays were checked for a 10-kya fitness start and
unchanged survival conditioning after the 50-kya pulse.

Figure checks run on the headless canvas: all text lies within the canvas,
legends lie outside panels, PDFs contain no raster images, and PDF fonttype is
42. Figures were not visually displayed. Raw posteriors, detailed held-out
predictions, and simulation data remain under D:/phase2simselection/sim. The
repository contains compact summaries and provenance, with byte-for-byte copy
verification. Work used 20 production workers and stayed below the resource
allowances.

On the configured WSL/Linux environment, using the existing staged simulation
inputs and native dependencies described in [the experiment specification](recent_call_experiment.md):

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
bash scripts/launch_recent_calls.sh run
```

The workflow reuses valid simulations and posterior receipts based on scientific
inputs, pair/position manifests and output hashes. The `analyse` phase reruns
only analysis and its audit. A fresh clone does not itself contain the large
simulation inputs. `scripts/summarize_recent_calls.py` recreates the tables
from a completed, audited analysis without rerunning simulation or decoding.
