# EAS time-resolved selection pilot: measured performance

Simulation and scoring finished in 21.52 minutes including setup, within the
requested 30-minute limit. All 239 analyzed runs completed: 40 development
nulls, 99 independent calibration nulls, 40 independent neutral tests, 40
selected tests at s=0.01, and 20 additional strong-selection controls at s=0.05.
The strong controls share onset-ancestry seeds with the first 20 selected tests;
they are a paired diagnostic, not 20 additional independent demographic histories.
Two earlier setup simulations were excluded from evaluation.

The main finding is that posterior mass improves the interval-overlap criterion
for a strong sweep, but does not yet give reliable selection dates. The weaker
s=0.01 episode was not temporally recovered by any method.

## What ran

Unscaled EAS PHLASH median demography, no introgression, selection active from
30 to 20 kya, dominance h=0.5, and standing variation near 5% frequency at onset.
There was no final-frequency conditioning. Each region was 1 Mb with 100 sampled
diploids, 1,000 fixed sampled haplotype pairs, 10-kb scoring stride, and cutoffs
1 through 50 kya at one-kyr steps. These smaller region/panel sizes fit the
interactive budget; they are not the planned 10-Mb/200-diploid/10,000-pair study.

The mean final frequency was 0.290 in the s=0.01 tests (range 0.016-0.684) and
0.999 in the s=0.05 controls (range 0.987-1.000). The latter were added as an
exploratory positive control after weak recovery appeared in the initial case.
No primary thresholds or metrics were changed in response to the results.

Each method was calibrated separately at nominal region-level alpha=0.05 using
the maximum over the entire position-by-time search. Thus time bins, pairs, and
positions were not treated as independent observations.

## Primary one-kyr interval scan

"Near-site detection" means any significant cell within 100 kb of the selected
coordinate, at any time. "Joint overlap TPR" additionally requires that cell to
overlap the true [20,30) kya episode. Center overlap and temporal F1 use only
the exact selected-coordinate profile. No-call temporal F1 is zero.

| Strength | Method | Any call in region | Near-site detection | Joint overlap TPR | Center overlap TPR | Mean temporal F1 |
|---|---|---:|---:|---:|---:|---:|
| s=0.01, n=40 | True TMRCA | 5/40 (12.5%) | 2/40 (5%) | 0/40 | 0/40 | 0.000 |
| s=0.01, n=40 | Posterior median | 5/40 (12.5%) | 0/40 | 0/40 | 0/40 | 0.000 |
| s=0.01, n=40 | Posterior mass | 2/40 (5%) | 1/40 (2.5%) | 0/40 | 0/40 | 0.000 |
| s=0.05, n=20 | True TMRCA | 20/20 (100%) | 19/20 (95%) | 0/20 | 0/20 | 0.000 |
| s=0.05, n=20 | Posterior median | 20/20 (100%) | 19/20 (95%) | 0/20 | 0/20 | 0.000 |
| s=0.05, n=20 | Posterior mass | 19/20 (95%) | 18/20 (90%) | 15/20 (75%) | 10/20 (50%) | 0.078 |

The same held-out neutral tests gave FPR 0/40 for truth, 4/40 (10%) for median,
and 0/40 for posterior mass. Wilson 95% intervals are 0-8.8% for 0/40 and
4.0-23.1% for 4/40; these are small-sample estimates, not proof of zero false
positives or a definitive violation of nominal 5% calibration.

For strong selection, posterior mass's joint-overlap TPR of 75% has a Wilson
95% interval of 53.1-88.8%. Its paired improvement over median is 75 percentage
points (replicate-bootstrap interval 55-90 points), conditional on this fitted
calibration. This does not include uncertainty from regenerating calibration
cohorts. There were no temporal successes in 40 weaker selected tests; the
corresponding Wilson upper bound is 8.8%.

## Outside-interval calls and wider windows

For the user's temporal false-positive criterion, at the exact selected
coordinate in the strong controls, at least one wholly outside-episode interval
was called in 17/20 truth, 9/20 median, and 7/20 mass runs for the one-kyr scan.
This is separate from neutral FPR. An interval that partially overlaps truth
does not count as wholly outside, which is why duration-based F1 is essential.

The prespecified secondary scan combined window widths 1,2,5,10 kyr, all stepping
by one kyr, and recalibrated the entire 182-window-by-position search:

| Strong-control method | Joint overlap TPR | Center overlap TPR | Mean temporal F1 | Event F1 | Any wholly outside interval at center |
|---|---:|---:|---:|---:|---:|
| True TMRCA | 16/20 (80%) | 14/20 (70%) | 0.290 | 0.348 | 4/20 (20%) |
| Posterior median | 0/20 | 0/20 | 0.000 | 0.000 | 11/20 (55%) |
| Posterior mass | 18/20 (90%) | 16/20 (80%) | 0.263 | 0.480 | 1/20 (5%) |

Mass's 90% joint-overlap TPR has a Wilson 95% interval of 69.9-97.2%. Its average
called support at the center spans 20.25 kyr for a true 10-kyr episode, while
average temporal recall is only 45.5%. The high any-overlap TPR therefore
overstates precise temporal recovery. Event F1 uses 20 selected and 40 neutral
test regions, one-to-one episode matching, and all unmatched spatial/time
components as false events. It is not directly interchangeable with temporal F1.

## Interpretation

In the strong controls, the mean true-tree excess at the selected coordinate
peaks at 33.5 kya, before selection started. The decoded median excess peaks at
10.5 kya, while posterior-mass excess is broad and peaks at 2.5 kya. These are
coalescence-profile modes, not estimates of selection onset. Standing-variation
genealogies and the decoder's temporal behavior both matter here.

Posterior mass is worth retaining as a scoring arm: it preserves fractional
evidence and improves overlap in this strong-selection case. However, the low
temporal F1 and shifted profiles mean we should not label the called support a
recovered selection episode yet. The next timing experiment should vary onset,
cessation, and starting frequency using true genealogies first, then quantify
how the decoder changes the mapping. Additional null cohorts are needed for
precise FPR estimates and for genome-wide claims.

## Artifacts and verification

- [Primary s=0.01 results](s0p01_one_kyr/README.md)
- [Primary s=0.05 results](s0p05_one_kyr/README.md)
- [Multiscale s=0.01 results](s0p01_multiscale/README.md)
- [Multiscale s=0.05 results](s0p05_multiscale/README.md)
- [Execution design and commands](../../time_resolved_pilot_execution.md)
- [Timing](timing.json) and [integrity checks](integrity_checks.json)

Each result directory contains per-replicate metrics, calibration parameters,
Wilson intervals, paired comparisons, center coalescence summaries, cumulative
cutoff baselines, PNG and editable vector PDF plots, and a checksum manifest.
Raw simulated trees and posteriors remain in the local simulation directories.
All score/tree/report hashes and raw posterior hashes were verified. Reconstructed
median fractions agree exactly with native summaries at 10 and 50 kya in every
run. All trajectories have the expected 1,201 generations/time entries. Six
focused tests passed. The plots were inspected locally; PDFs contain embedded
TrueType fonts and no raster image objects. Plots are not displayed inline.
