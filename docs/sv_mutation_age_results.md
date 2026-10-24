# SV mutation-age pilot results, 2026-10-24

The pilot successfully generated and decoded neutral SV-surrogate data under
the project's introgression demography. The important result is a limitation:
modern pairwise TMRCA estimates constrain mutation ages, but this initial
estimator comparison does not demonstrate a useful incremental improvement
over an allele-frequency baseline for precise age inference. Uncalibrated
products of pairwise probabilities give substantially overconfident intervals.

See [the protocol and runnable command](sv_mutation_age_pilot.md) and
[the motivating method sketch](sv_origin_tmrca_method.md). Numerical summaries
and figures are stored in `docs/results/sv_mutation_age_20261024/`; raw synthetic
trees, VCFs and posterior arrays remain in the local simulation directory.

## Dataset and separation of training from evaluation

- 160 independent 500-kb regions, totaling 80 Mb of independent simulated
  sequence, with 400 modern EAS haplotypes per region.
- The EAS PHLASH size trajectory, 2% pulse at 50 kya, donor Ne=3,600, split
  at 500 kya, 25-year generations, and heterogeneous recombination.
- 277 observed candidate events: 60 old imported, 57 post-introgression, and
  160 ancestral-retention examples. Missing history classes were retained as
  missing rather than filled by rejection sampling whole regions.
- 80 training regions containing 137 events; 40 calibration regions containing
  68 events; 40 test regions containing 72 events.
- The test set has 16 old imported, 16 post-introgression and 40 ancestral
  events. The first two groups constitute the central, 32-event comparison.
- Each region was decoded with the focal base masked and again with 2,001 bp
  masked around each event. The principal table below uses the focal-base mask.

Candidate SV events and observed SNPs were independent mutation overlays.
Focal mutation ages and trees were withheld from inference. Archaic samples
were used for ascertainment/truth only. This is a neutral, correctly phased,
single-origin marker experiment, not a structural read/genotyping benchmark.

## Exact-tree bounds expose an asymmetry

Among all 60 old imported events, 59 have a modern carrier-clade MRCA younger
than 50 kya. In the held-out old imported group, 15/16 have age bounds spanning
introgression, even when the modern tree is known exactly. Their median ratio
of upper to lower bound is 32.8. Thus a young carrier crown cannot exclude an
old mutation, and precise age estimation requires additional information or
strong model assumptions about placement on the long ancestral branch.

New mutations are more identifiable in this simulation: 14/16 held-out events
have an exact carrier/noncarrier upper bound younger than 50 kya. The same
14/16 have sampled introgressed noncarriers. None of the 16 held-out old
imported events has introgressed noncarriers. The median exact-tree bound
width for new mutations is 2.65-fold.

These figures are conditional on the pilot's donor history and event
ascertainment. They are not expected to hold unchanged under different donor
diversity, multiple pulses, selection or sampling schemes.

## Performance on the 32 imported-versus-new test events

Fold error is max(estimated age / true age, true age / estimated age), so one
is exact. Intervals target 90% coverage, but uncalibrated score intervals do
not attain it. A confident timing call requires the whole interval to lie
before or after the 50-kya event.

| Method | Median fold error | Point estimate on correct side of 50 kya | Interval coverage | Confident timing calls |
| --- | ---: | ---: | ---: | ---: |
| Carrier mean TMRCA | 2.34 | 65.6% | Not supplied | Not supplied |
| Raw carrier/carrier + carrier/noncarrier CDF score | 2.17 | 71.9% | 37.5% | 26/32; 4 wrong |
| CDF score plus oracle-ancestry pair targeting | 2.56 | 87.5% | 28.1% | 30/32; 3 wrong |
| Trained attenuation + calibrated CDF interval | 4.63 | 50.0% | 93.8% | 0/32 |
| TMRCA-feature kNN + calibrated interval | 3.14 | 78.1% | 96.9% | 1/32; 0 wrong |
| Frequency-only kNN + calibrated interval | 2.86 | 71.9% | 96.9% | 5/32; 0 wrong |

The calibrated TMRCA-feature kNN intervals have median upper/lower ratio
169-fold in this subset. Calibration achieves coverage largely by retaining
substantial uncertainty, rather than resolving mutation ages. Zero observed
errors among one or five calls is not evidence of a negligible error rate.

The 32 focal events occupy 16 independent test regions. Region-bootstrap 95%
intervals for median fold error are 1.46-3.37 for carrier mean, 1.68-3.53 for
the raw CDF score, 2.69-4.05 for TMRCA-feature kNN and 2.55-4.53 for
frequency-only kNN. These are pilot-scale uncertainty estimates, not evidence
for a decisive ranking among small point-estimate differences.

The oracle-ancestry comparison points toward a better sampling design:
carrier/noncarrier comparisons within introgressed ancestry improve timing
classification. However, its poor age coverage means its raw intervals should
not be used for scientific claims. Exact pair duplicates were removed before
combining pair sets; correlation among distinct pairs remains.

Masking 2,001 bp around the SV did not change the main conclusion. In the focal
32-event subset, raw CDF median fold error changed from 2.17 to 2.49 and
calibrated TMRCA-kNN from 3.14 to 3.25. The latter's coverage changed from
96.9% to 93.8%. This small sensitivity analysis does not establish robustness
to mapping, phasing or genotype errors.

## Why the overall score looks more favorable

Across all 72 test events, including much older ancestral-retention variants,
median fold errors are 5.24 for the carrier mean, 2.50 for the raw CDF score,
1.75 for TMRCA-feature kNN and 1.78 for frequency-only kNN. The nearly identical
last two values are why improvement over carrier mean alone is insufficient
evidence that Gamma-SMC adds useful age information.

Across the same full set, raw score 90% intervals cover only 18.1% of mutation
ages. Calibrated kNN coverage is 95.8%, with a 66.4-fold median interval width.
Report the full and restricted comparisons together; the older ancestral
examples substantially change the apparent difficulty of the task.

## Interpretation and next method changes

1. Preserve explicit sampled-tree age bounds and an unresolved category.
   Selection is not required for old imported mutations to have young carrier
   crowns: this neutral pilot already shows the phenomenon.
2. Prioritize informative carrier/noncarrier comparisons within local archaic
   ancestry, or sequence-based nearest noncarriers, rather than predominantly
   random cross-ancestry pairs. Evaluate inferred ancestry rather than assuming
   truth in the next benchmark.
3. For age intervals, move toward a coherent local-tree posterior with mutation
   placement on its branches. The current independent pair-product construction
   loses calibration; merely weakening its weights does not restore precision.
4. Add reliable archaic SV observations or an age-aware modern/archaic genealogy
   to constrain old imported mutations. The modern-only pilot deliberately
   excludes that information; contemporary Gamma-SMC emissions should not be
   applied to ancient samples without an appropriate sampling-time model.
5. Expand demographic, selection, phasing, genotyping and recurrence scenarios
   before interpreting a trained model as transferable to real SVs. The kNN
   model is conditional on this simulation and ascertainment distribution.

## Verification

All 160 simulation filesets passed checksum and semantic audits: tree/SNP
sample consistency, mutation on the recorded branch, age within carrier bounds,
carrier/ancestry labels, and pair membership. Exact decoder output coordinates,
pair order, finite positive Gamma parameters and complete raw posterior streams
were checked. Nine focused tests passed, including mutation-age ambiguity,
mask boundaries, cache corruption, pair-label correctness, duplicate-pair
removal, region split consistency, and absence of held-out-age leakage into
training. The full repository test suite was not run.

The output tables include 500-replicate region-bootstrap intervals. These
resample independent regions and preserve dependence between co-located SVs;
they should be used when assessing the small point-estimate differences.

All three figures were visually inspected: labels are readable, legends sit
above the panels, and no legend overlaps the data. Both PNG and editable
vector PDF are saved; no notebook display is required.

| Figure | PNG | PDF |
| --- | --- | --- |
| Exact modern-tree age bounds | [PNG](results/sv_mutation_age_20261024/oracle_age_bounds.png) | [PDF](results/sv_mutation_age_20261024/oracle_age_bounds.pdf) |
| True versus estimated ages | [PNG](results/sv_mutation_age_20261024/age_estimates.png) | [PDF](results/sv_mutation_age_20261024/age_estimates.pdf) |
| Interval coverage and call fraction | [PNG](results/sv_mutation_age_20261024/interval_calibration.png) | [PDF](results/sv_mutation_age_20261024/interval_calibration.pdf) |
