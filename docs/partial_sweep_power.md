# Partial-sweep power and neutral false calls

The analysis target is approximately 70% of s=0.005 selected regions called
anywhere while reducing the fraction of wholly neutral, isolated 10 Mb regions
called by the identical rule. Correct localization is not required for this
primary region-level endpoint. Of the 100 fresh
s=0.005 simulations, 99 have sample allele frequency below one. Fixed samples
are retained in evaluation and reported separately; they are excluded when
training a detector specifically for partial sweeps.

## Error definitions

The main error measure is the fraction of neutral regions with at least one
called peak. A second measure is the number of false peaks per neutral region.
Every call in a wholly neutral region is false, so its false-discovery proportion
is one when any peak is called and zero otherwise. The mean of those proportions
therefore equals the regional any-call rate. This is also the family-wise error
rate for an all-neutral 10 Mb region.

A prespecified focal-stride false-positive rate is reported separately. It asks
a different question from scanning the region and selecting its strongest peak.
The scan's peak-finding and evidence filters must be applied to the neutral test
regions as well as selected regions.

For a mixture of selected and neutral genomic regions, FDR also depends on sweep
prevalence and the definition of one discovery. Neutral-region error rates alone
do not establish genome-wide FDR. A null-calibrated position is not a null-calibrated
maximum over a region, and a reported cluster does not imply every base in it is
affected by selection.

Within selected regions, peak summits within ±100 kb of the selected allele are
used only for a separate localization diagnostic. At most one called peak matches the one
selected locus. This is an explicit tolerance rather than a causal footprint.
Interval-overlap results and wider radii are retained as secondary diagnostics.
Calls outside this tolerance can belong to a decaying sweep tail and should not
automatically be interpreted as biologically neutral.

The selected-allele carrier label is defined at the focal variant. It does not
assign a binary causal label to every distal genomic position: recombination
can change the ancestry and tree as position changes. If the focal allele is
fixed, every sampled haplotype is a carrier, so carrier status alone cannot
define the footprint or identify false peaks.

## Paired truth and decoded grid

`paired_peak_grid.py` reads matched saved `frac_recent_T` profiles from the same
simulations and pair manifest, using either true pair TMRCA or Gamma-SMC posterior
mean TMRCA. Each source is calibrated against its own matched neutral statistic;
truth thresholds are never applied directly to decoded values.

The fixed grid is:

- TMRCA cutoffs: 5, 10, 20, 30, 40, and 50 kya.
- Inclusive p ≤ 0.05, 0.01, or 0.001, including all ties in the upper tail.
- Reporting strides: 10, 20, and 50 kb, by subsampling the saved 10 kb profiles.
  This does not change or rerun the underlying decoder.
- Minimum evidence: 1, 3, 5, 10, or 20 significant strides.
- Allowed interruption: zero or one nonsignificant reporting stride.

Selected ranks use `(1 + neutral exceedances) / 1001`; neutral test ranks leave
that replicate out, giving `(1 + other-neutral exceedances) / 1000`.
The inclusive comparison is explicit because p=0.001 is the lowest possible
neutral leave-one-out rank. It avoids interpreting the unresolvable strict
p<0.001 neutral comparison as zero observed error. Selected ranks do not hit
exactly 0.001. No p-value is zero and linked positions are not counted as
independent null simulations. The p=0.001 estimate has limited tail resolution.

The grid is exploratory, evaluated on 100 selected s=0.005 replicates and 1,000
neutrals. Its best settings are not independent confirmation of an optimized
rule. The output includes paired stride-call agreement; disagreement with the
true-TMRCA statistic measures decoder fidelity, not causal-selection FDR.

Mean profiles are also plotted against distance from the allele. Distance-bin
averages are calculated within each simulation first, then averaged over
replicates. Bands show approximate Monte Carlo uncertainty for the difference
between selected and neutral means, with one independent simulation as the unit.

## Internally cross-validated detector comparison

`partial_sweep_cv.py` uses the fixed specification in
`configs/eas_partial_sweep_analysis.json`. Five deterministic outer folds split
whole simulations, separately within each selection arm and the neutral group.
Each test fold is excluded from null normalization, supervised template fitting,
threshold selection, and method choice. All randomness uses seed 20380101.

Candidate scores all derive from decoded `frac_recent_T`:

- Individual cutoffs, with several smoothing widths.
- Local contrasts between the center and a broader surrounding average.
- A regularized linear combination of the six time cutoffs.
- A regularized combination of time cutoffs and spatial contrasts.
- Cluster mass above a training-neutral threshold.

Candidates scan the whole region using the same rule at every coordinate. The
selected coordinate supplies training labels and a template example, but no
prediction is favored because it is near the center. Peak summits have a fixed
200 kb minimum separation. Both all-peak calling and at-most-one-peak calling
are considered, with the latter explicitly using the isolated-region scenario.

Training thresholds must call at least 70% or, in a separate policy, 80%
of nonfixed training s=0.005 regions, anywhere in each region. The 80% policy tests whether a
training sensitivity margin preserves about 70% on held-out simulations. Method
and threshold selection minimize the training neutral-region any-call rate,
then mean neutral peak count, then favor higher power. Spatial FDR is only a
final tie-breaker and an evaluation diagnostic.

An adaptive result selects the candidate independently within each training
fold and reports only that candidate's predictions on the held-out fold. Scores
for individual candidates are also saved, but choosing a winner from those
pooled held-out scores introduces another model-selection step. Infeasible
candidate/fold combinations are not silently treated as complete evaluations.

These are internal cross-validation estimates from previously summarized
simulations, not a new independent confirmatory dataset. The s=0.005 sample
size is only 100; a power estimate near 70% has substantial Monte Carlo error.
None of the learned rules carries an automatic formal FDR-control guarantee.

## Integrity, outputs and invocation

Both workflows check all input profile SHA-256s, receipt sizes, coordinates,
cutoff monotonicity, and integer pair-count fractions against the fresh study
manifest. The cross-validation workflow has idempotent per-fold caches keyed by
scientific inputs, analysis specification and analysis source, with payload
checksums. No simulations are regenerated.

The paired grid writes to `analysis/paired_peak_grid`; the primary region-level
cross-validation writes to `analysis/partial_sweep_region_cv`. An earlier
localization-constrained run is retained separately in `analysis/partial_sweep_cv`
and is not the primary region-call result. Per-replicate results, parameter choices, fold
assignments, model weights, summaries, figures and provenance are retained on D.
Figures use large text, PNG and vector PDF with editable text.

Run in Linux/WSL with the existing pinned uv environment. Raw simulations must
already be available locally; they are not distributed through Git.

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
uv --no-config run --no-project --python .venv/bin/python \
  python -m gamma_smc_aou.paired_peak_grid \
  --study /mnt/d/phase2simselection/sim/eas_q02 --workers 20
uv --no-config run --no-project --python .venv/bin/python \
  python -m gamma_smc_aou.partial_sweep_cv \
  --study /mnt/d/phase2simselection/sim/eas_q02 \
  --spec configs/eas_partial_sweep_analysis.json
```
