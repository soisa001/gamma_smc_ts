# EAS selection beginning 10 kya

Selection coefficients 0.005, 0.01, 0.02, 0.03; h=0.5; 25 years/generation; unscaled WF forward simulation. Selected focal allele must be strictly segregating in the whole population at present.
De novo: introduce one copy at 10 kya. Introgression: a 1% Neanderthal pulse at 50 kya, with allele fixed in donors, human–Neanderthal split 500 kya and Neanderthal Ne=3,600. The introgressed allele is neutral until 10 kya. Pre-onset losses/fixations are rejected because they cannot pass the final segregation condition.
Neutral regions are unconditioned, with focal label at the nearest present-day sample-polymorphic variant. Nulls match the origin-specific demographic simulation and forward duration. EAS PHLASH median history; 1 Mb; 100 diploids; 1,000 unfiltered pairs; 10 kb stride; cutoffs 1–50 kya.
Independent development and calibration nulls; conservative empirical max calibration over all positions and time windows at alpha=0.05. Primary intervals are one kyr. Secondary widths 1/2/5/10 kyr and cumulative 10/50-kya comparators are separately calibrated.
Spatial TPR allows any time within 100 kb of the known center. Joint TPR additionally requires overlap with [1,10) kya. Temporal F1 is computed at the exact selected coordinate against those nine cells. Neutral FPR is any call anywhere in a held-out neutral region. Center mistiming means a center call with no temporal overlap; this differs from neutral FPR.

Completed proposals: 3505; analyzed selected regions: 201.

| Origin | s | Accepted | Proposals in prefix | Lost | Fixed |
|---|---:|---:|---:|---:|---:|
| de_novo | 0.005 | 11 | 396 | 385 | 0 |
| de_novo | 0.01 | 12 | 396 | 384 | 0 |
| de_novo | 0.02 | 15 | 396 | 381 | 0 |
| de_novo | 0.03 | 16 | 396 | 380 | 0 |
| introgression | 0.005 | 30 | 396 | 366 | 0 |
| introgression | 0.01 | 37 | 396 | 359 | 0 |
| introgression | 0.02 | 38 | 396 | 358 | 0 |
| introgression | 0.03 | 42 | 395 | 353 | 0 |

| Origin | s | Mass spatial TPR | Mass joint TPR | Neutral FPR | Temporal F1 |
|---|---:|---:|---:|---:|---:|
| de_novo | 0.005 | 0/11 | 0/11 | 2/50 | 0.000 |
| de_novo | 0.01 | 0/12 | 0/12 | 2/50 | 0.000 |
| de_novo | 0.02 | 0/15 | 0/15 | 2/50 | 0.000 |
| de_novo | 0.03 | 4/16 | 4/16 | 2/50 | 0.144 |
| introgression | 0.005 | 0/30 | 0/30 | 2/50 | 0.000 |
| introgression | 0.01 | 1/37 | 1/37 | 2/50 | 0.000 |
| introgression | 0.02 | 25/38 | 25/38 | 2/50 | 0.406 |
| introgression | 0.03 | 39/42 | 39/42 | 2/50 | 0.656 |

## One-kyr estimator comparison

| Origin | s | True-tree joint TPR | Median joint TPR | Mass joint TPR | Mass region F1 | Mass event F1 |
|---|---:|---:|---:|---:|---:|---:|
| de_novo | 0.005 | 0/11 | 0/11 | 0/11 | 0.000 | 0.000 |
| de_novo | 0.01 | 0/12 | 0/12 | 0/12 | 0.000 | 0.000 |
| de_novo | 0.02 | 0/15 | 0/15 | 0/15 | 0.000 | 0.000 |
| de_novo | 0.03 | 7/16 | 4/16 | 4/16 | 0.364 | 0.308 |
| introgression | 0.005 | 0/30 | 0/30 | 0/30 | 0.000 | 0.000 |
| introgression | 0.01 | 0/37 | 0/37 | 1/37 | 0.050 | 0.050 |
| introgression | 0.02 | 3/38 | 10/38 | 25/38 | 0.769 | 0.562 |
| introgression | 0.03 | 14/42 | 20/42 | 39/42 | 0.940 | 0.696 |

## Held-out neutral error

| Origin | Estimator | False-positive regions | FPR | Wilson 95% interval |
|---|---|---:|---:|---:|
| de_novo | truth | 3/50 | 6.0% | 2.1%–16.2% |
| de_novo | median | 2/50 | 4.0% | 1.1%–13.5% |
| de_novo | mass | 2/50 | 4.0% | 1.1%–13.5% |
| introgression | truth | 0/50 | 0.0% | 0.0%–7.1% |
| introgression | median | 2/50 | 4.0% | 1.1%–13.5% |
| introgression | mass | 2/50 | 4.0% | 1.1%–13.5% |

Wilson intervals, region/event F1 and alternative scans are in performance.tsv. All simulation artifacts remain outside Git. Figures are saved only to PNG and editable vector PDF.
- Fixed target is first 50 accepted replicates in proposal order; actual target is in schedule.
- Deadline-limited strata use a complete proposal prefix; Wilson intervals are exploratory.
- All methods and scan variants have separate null calibration; their union is not a calibrated discovery rule.
- Cumulative 10-kya comparator uses the known simulated onset and is an oracle comparator, not an onset inference method.
- Selection truth [0,10) kya is evaluated on observable [1,10) support: nine one-kyr cells; 0-1 kya is not scanned.
- Region F1 counts joint hits (spatial hits for cumulative baselines), selected misses, and neutral any-call false positives; it depends on cohort mixture.
- Event F1 also penalizes unmatched connected components; each selected region has at most one true event.

## Interpretation and exact all-pairs audit

All 3,505 submitted units completed: 3,167 selected proposals produced 201
segregating survivors, and 338 neutral regions completed. The conservative
elapsed simulation clock was 58 minutes 19 seconds, including setup checks,
replays, and the completion phase. No submitted unit failed or remained
incomplete. Targets of 50 selected survivors per setting were not reached.

Posterior mass improves recent-interval detection for introgressed s=0.02
(25/38 versus median 10/38) and s=0.03 (39/42 versus 20/42). It gives little
observed power for s=0.005–0.01, and no improvement over the median for de novo
s=0.03 (both 4/16). The separately calibrated 1/2/5/10-kyr multiscale scan has
the same posterior-mass detection counts as the primary one-kyr scan.
The 10-kya cumulative comparator detects 22/38 and 38/42 strong introgression
replicates; the 50-kya comparator detects 13/38 and 27/42. These cumulative
comparisons assess detection only, not inferred timing.

Using all 19,900 true pairwise TMRCAs improves the de novo s=0.02 joint result
from 0/15 to 3/15 and s=0.03 from 7/16 to 11/16. Introgression s=0.02 improves
from 3/38 to 16/38, and s=0.03 from 14/42 to 25/42. Full-pair truth still has
zero joint detections at s=0.005–0.01 in both models. Its held-out FPR is 3/50
for de novo and 0/50 for introgression, unchanged from the sampled-pair truth
benchmark. The full-pair diagnostic is **truth only**; GammaSMC still uses
1,000 pairs. See `all_pairs_truth_performance.tsv`.

Sparse carrier representation limits the weak de novo experiment: all 11
s=0.005 survivors and 10/12 s=0.01 survivors have zero carrier–carrier pairs
among the 1,000 decoded pairs. Five and three of those survivors, respectively,
have no sampled carrier at all. Survival was conditioned on whole-population
segregation, exactly as specified, without a sample-frequency threshold.
For comparison, the median surviving whole-population frequency is 0.43%
at de novo s=0.005 and 0.83% at s=0.01. Introgression survivors are more common
(median 17% and 36% for these coefficients), so sparse sampling cannot explain
all weak-selection failures. See `selected_sample_carriers.tsv` and
`selected_frequencies.tsv`.

Detection and onset dating must remain separate. With full-pair truth,
introgression s=0.02 has 35/38 spatial detections but only 16/38 overlapping
the recent interval; s=0.03 has 42/42 spatial detections but only 25/42 recent
overlaps. Posterior mass gives more recent overlaps, which is consistent with
posterior smoothing/redistribution affecting timing. These results do not
demonstrate precise recovery of the 10-kya onset. Among 60 selected replicates
with a posterior-mass call at the exact focal coordinate, none has a wholly
outside-time call component. Broad intervals can pass this overlap criterion;
mean center temporal F1 is only 0.406 and 0.656 for introgression s=0.02 and
s=0.03. The next useful comparison is full-pair posterior inference together
with frequency-stratified power and multiple onset dates.

All primary confidence intervals are in `performance.tsv`; for example, mass
joint power is 65.8% (Wilson 95%: 49.9–78.8%) at introgression s=0.02 and
92.9% (81.0–97.5%) at s=0.03. These are pilot estimates conditional on surviving
segregation. Region F1 uses each selected cohort together with 50 neutral test
regions, so its value depends on that study mixture.

Figures were visually inspected and saved as PNG plus vector PDF with embedded
editable text. Plot rendering was separated from the numerical report after
fixing a floating-point error-bar edge case and spacing the coefficient labels.
`plot_provenance.json` records the final plotting source and figure hashes;
the original numerical-report source is retained under `executed_sources/main`.
