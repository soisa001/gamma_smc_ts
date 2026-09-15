# Removing the ALT/REF penalty and comparing carrier/noncarrier recency

**Current specification:** the [archaic ALT/ALT signal](archaic_alt_alt_signal.md)
uses `fAA * nAA / Npairs` as the primary score. This report records the preceding
method comparisons; its contrast-based labels are historical.

**The ALT/REF penalty can be removed without a meaningful loss in the tested
detector. Raw ALT/ALT and REF/REF fractions, however, behave very differently
from the frequency-weighted score.** At exact archaic sites with T=50 kya,
removing only the penalty leaves every decoded region call unchanged at
p<=.01, .05 and .10. The unweighted ALT/ALT regional maximum has zero power
because nearly every neutral region also reaches a fraction of one.

This analysis uses only EAS selection starting immediately at the 50-kya
introgression pulse: 100 selected regions with s=0.005 and 1,000 neutral
regions, 400 haplotypes, 10,000 fixed pairs, and 2% introgression. Rates, seeds,
and original profiles are unchanged. The later-onset arm is excluded.
All statistics use hard `frac_recent_T` calls; no posterior-mass average is used.

## Scores and comparison design

Let A=fAA, R=fRR, H=fAR, and w=nAA/Npairs. The comparisons are:

| Name | Local score |
|---|---|
| ALT/ALT | A |
| REF/REF | R |
| Unweighted contrast | max(A-R,0) |
| Weighted contrast without ALT/REF | w*max(A-R,0) |
| Original joint | w*max(A-R,0)*(1-H) |
| Carrier mass | w*A = recent ALT/ALT pairs / all sampled pairs |

The first three requested comparisons have no all-pair gate and no carrier
frequency weight. The weighted contrast and original joint are also compared
with exactly the same original median all-pair gate, isolating removal of H.
Ungated versions are retained in the full tables. Empty classes supply no
stand-alone evidence; missing R uses a zero comparator baseline so fixation
is retained. No region is discarded and no AF or minimum pair-count filter
is imposed.

The blind scan takes the maximum score across eligible archaic sites in a
10 Mb region, or across their assigned 10 kb grid positions. A separate
focal-candidate diagnostic uses the known selected SNP, described below.
The selected SNP is not supplied to the blind region scan.

The original five whole-region folds are retained, using 400 fit neutrals,
400 separate calibration neutrals, and 200 test neutrals per fold. Each fold
has 80 selected training and 20 selected test regions. Fixed rules are tested
at T=5,10,20,30,40,50 kya, with no cutoff tuning. Conservative calibration
includes ties. Power is the fraction of selected regions called; neutral call
fraction is the operational regional false-positive rate, not discovery FDR.

## Blind region scan at every archaic SNP

At T=50 kya and regional p<=0.05:

| Rule | True power | True neutral calls | Decoded power | Decoded neutral calls |
|---|---:|---:|---:|---:|
| ALT/ALT alone | 0% | 0% | 0% | 0% |
| REF/REF alone | 48% | 4.9% | 20% | 5.0% |
| Unweighted AA-RR contrast | 15% | 4.7% | 9% | 5.1% |
| Weighted AA-RR, no AR penalty, median gate | 79% | 5.2% | 80% | 4.6% |
| Original joint, median gate | 79% | 5.2% | 80% | 4.6% |
| Carrier mass, median gate | 81% | 5.1% | 84% | 4.9% |
| Archaic AF alone | 85% | 5.2% | 85% | 5.2% |

These percentages represent 100 selected and 1,000 neutral regions. For the
decoded 10 kb grid, the weighted contrast without AR gives 81% power and
4.8% neutral calls, exactly matching the original joint call set at p<=.05.
For grid truth there are four changed calls at p<=.05; the net rates change
from 81% / 5.2% with AR to 80% / 5.3% without AR. At grid p<=.10 there is one
changed decoded call. Thus the penalty is not algebraically redundant at a
nearby grid point, where recombination can separate the marker from the tree
being evaluated, but its observed contribution is small.

At exact sites, true H is identically zero under the eligible single-mutation
model for these cutoffs. The AR term is therefore algebraically redundant
for true TMRCA. The decoded exact-site score can change numerically, but the
T=50-kya call sets were identical at all three tested p thresholds.

Thresholds selected from training regions to aim at 70% power give:

| Decoded rule | 10 kb grid: power / neutral calls | Every site: power / neutral calls |
|---|---:|---:|
| Weighted AA-RR, no AR penalty, median gate | 72% / 2.0% | 72% / 2.4% |
| Original joint, median gate | 72% / 1.9% | 72% / 2.4% |
| Carrier mass, median gate | 71% / 1.3% | 69% / 1.4% |
| Archaic AF alone | 70% / 1.4% | 71% / 1.4% |

One or two changed regions do not establish an advantage for one of these
similar operating points. Carrier mass needs neither REF/REF subtraction nor
an ALT/REF penalty and remains competitive with the more complex scores.

## Why raw ALT/ALT recency fails

For decoded TMRCA at T=50 kya, **970/1,000 neutral regions** contain an eligible
archaic SNP with fAA=1. The same is true of **99/100 selected regions**.
For true TMRCA, the counts are **984/1,000 neutral** and **100/100 selected**.
The other 16 neutral regions have no eligible archaic marker and are retained.

Across regions, the median of the minimum ALT/ALT pair count among positions
attaining the regional maximum is **one pair**. A rare neutral marker with one
sampled ALT/ALT pair can therefore attain exactly the same raw fraction as a
selected allele supported by thousands of pairs. Neutral carrier genealogies
can also be recent: a large within-carrier recent fraction is not by itself
specific to selection.

With conservative upper-tail ranks, a test score of one cannot exceed those
neutral ties. This produces zero calls at p<=.05, not an implementation failure.
Conversely, choosing the raw threshold A>=1 to retain selected regions calls
99% of selected regions and 97% of neutral regions in the decoded site scan.
There is no useful 70% operating point at that tied maximum.

Subtracting R without accounting for carrier support does not solve the
problem: a rare marker can have A=1 and an ordinary R. The weighted score
suppresses that event by multiplying by the small nAA/Npairs.

## Dependence on the time cutoff

For the ungated decoded regional scan at exact archaic sites, p<=.05:

| T (kya) | ALT/ALT power / neutral calls | REF/REF power / neutral calls | Unweighted contrast power / neutral calls |
|---:|---:|---:|---:|
| 5 | 0% / 0% | 68% / 4.4% | 5% / 4.6% |
| 10 | 0% / 0% | 56% / 4.5% | 4% / 4.9% |
| 20 | 0% / 0% | 36% / 5.0% | 10% / 5.0% |
| 30 | 0% / 0% | 26% / 5.1% | 4% / 5.1% |
| 40 | 0% / 0% | 22% / 4.9% | 5% / 4.5% |
| 50 | 0% / 0% | 20% / 5.0% | 9% / 5.1% |

REF/REF does detect some selected regions, especially at short cutoffs. However,
REF at an arbitrary archaic marker need not mean noncarrier of the causal
selected allele: the selected haplotype may carry REF at another marker.
Region-wide REF/REF detection must not be described as direct evidence for
selection on the reference allele or as a causal-allele-specific control.

## Known focal allele versus its noncarriers

To address that distinction, an additional **oracle focal-candidate diagnostic**
uses the known selected SNP in each selected region and the eligible archaic
SNP nearest the same prespecified coordinate in each neutral region. All 100
selected focal coordinates exactly match eligible markers. The 16 marker-free
neutral regions remain with zero evidence. This uses privileged knowledge of
the selected SNP, different ascertainment from a blind scan, and no regional
maximum. Its higher power must not be quoted as genome-wide scanning power.

At T=50 kya and p<=.05, without an all-pair gate:

| Focal-candidate rule | True power | True neutral calls | Decoded power | Decoded neutral calls |
|---|---:|---:|---:|---:|
| ALT/ALT alone | 0% | 0% | 0% | 0% |
| REF/REF alone | 22% | 4.9% | 10% | 5.1% |
| Unweighted AA-RR contrast | 6% | 4.8% | 1% | 5.5% |
| Weighted AA-RR, no AR penalty | 95% | 4.9% | 91% | 5.1% |
| Original joint | 95% | 4.9% | 91% | 5.1% |
| Carrier mass | 99% | 4.8% | 100% | 5.0% |
| Archaic AF alone | 99% | 4.9% | 99% | 4.9% |

Raw A=1 occurs at **242/1,000 decoded neutral focal candidates**, versus none
of the selected focal candidates. In truth, it occurs at 870/1,000 neutrals
and 77/100 selected candidates. Thus removing the regional search does not
rescue the unweighted A statistic at this cutoff.

Class support differs strongly between hypotheses: the median focal ALT/ALT
pair count is 34 in neutral regions and 6,791.5 in selected regions. The median
REF/REF counts are 8,784.5 and 308, respectively. These are not equal-sized
carrier/noncarrier comparisons. At T=5 kya, the focal REF/REF score has 60%
decoded power / 5.0% neutral calls (68% / 4.9% using truth). This is a property
of the remaining noncarrier genealogies and the simulation/ascertainment model;
it does not imply that the reference allele was favored.

## Implication

The simpler defensible choices within this experiment are the weighted
AA-RR contrast without the AR penalty or carrier mass. The experiment has not
shown that the extra AR penalty or RR subtraction improves detection. It has
also not shown a clear improvement over AF alone.

If the intended method must use only within-class recency and avoid frequency
weighting, neutral calibration should account for carrier frequency and class
pair counts. That conditional approach has not been tested here. The current
results concern direct upper-tail scans of A, R and their difference, not every
possible test built from those quantities. Overlapping haplotype pairs are not
independent Bernoulli observations, so naive pair-level binomial p-values would
not provide a substitute for the simulation calibration.

These remain exploratory existing-cohort comparisons with perfect archaic
labels, fixed rates, isolated 10 Mb regions, allele-survival ascertainment, and
the inherited msprime-neutral/scaled-SLiM-selected design. No new simulation
or decoding was needed. A final chosen method requires new-seed confirmation.

## Artifacts and checks

- [Regional metrics](results/allele_class_ablation/metrics.csv),
  [70%-target metrics](results/allele_class_ablation/target70_metrics.csv),
  [saturation diagnostics](results/allele_class_ablation/saturation.csv).
- [Regional site-scan figure, PDF](results/allele_class_ablation/allele_class_archaic_sites.pdf)
  and [PNG](results/allele_class_ablation/allele_class_archaic_sites.png);
  [grid figure, PDF](results/allele_class_ablation/allele_class_grid10kb.pdf)
  and [PNG](results/allele_class_ablation/allele_class_grid10kb.png).
- [Focal diagnostic metrics](results/allele_focal_comparison/metrics.csv),
  [focal allele/pair counts](results/allele_focal_comparison/regions.csv),
  [paired penalty-call comparisons](results/allele_comparison_diagnostics.json).
- [Regional audit](results/allele_class_ablation/audit.json) and
  [focal audit](results/allele_focal_comparison/audit.json).

The regional audit recomputed 1,104,400 prediction/target-threshold rows; the
focal audit directly counted calibration ranks for 283,800 rows. The analysis
matched 550,000 unchanged reference score values. It also checked 66,000
exact-site true-score equalities with and without the AR penalty. Four focused
tests cover penalty isolation, empty/fixed classes, unassigned grid positions,
saturated ties, fold-level calibration, and exact focal-SNP selection.
PNG and editable vector PDF outputs have automated canvas/legend/vector checks
and were not displayed inline. Large intermediate data remain on D.
