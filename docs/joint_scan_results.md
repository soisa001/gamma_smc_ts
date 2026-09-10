# Archaic pair-class scan: regional power and neutral call rates

The decoded pair-class score meets the approximate 70% power objective for
s=0.005 with selection beginning at 50 kya. On the 10 kb grid, a threshold chosen
to target 70% training power calls **72/100 selected regions and 19/1,000 neutral
regions**. Scanning every eligible archaic SNP gives 72/100 and 24/1,000.

However, **archaic allele frequency alone is at least as competitive in this
cohort**: 70/100 selected and 14/1,000 neutral calls on the grid, or 71/100 and
14/1,000 at every site. The large improvement over an all-pair TMRCA scan does
not demonstrate added value from TMRCA beyond carrier frequency. The 10-kya
onset remains difficult: reaching approximately 70% power with the tested
pair-class rules costs roughly 45–49% neutral region calls.

These are exploratory, cross-validated results within the existing fresh EAS
cohort, not prospective validation on new independent simulation seeds.

## Experiment and endpoint

The experiment uses 1,000 neutral regions and 100 selected regions per onset
(50 and 10 kya), s=0.005, dominance 0.5, 400 haplotypes, and the same 10,000
fixed haplotype pairs. All runs retain the 2% archaic pulse at 50 kya, mutation
rate 1.25e-8, recombination rate 1e-8, and 25-year generations. Each original
11 Mb simulation supplies a site-aligned, isolated 10 Mb scored region.

The original fresh simulations, genotypes and pair panel were reused.
Ancestry replays recover oracle archaic labels without replacing the original
genealogies or mutations. The old handoff simulations are not part of this
benchmark. All 100 selected regions per onset remain, including fixation and
insufficient-carrier cases. Power is conditional on the original requirement
that the selected allele survives and is observed in the sample.

A selected region is detected if **any position anywhere in its 10 Mb region
is called**. The neutral endpoint is the fraction of whole neutral regions
with any call. This is the user's operational “FDR,” statistically a regional
false-positive rate (FPR). There is no distance-to-causal-site success window.

Six thresholds are tested: T=5, 10, 20, 30, 40 and 50 kya. Every statistic is
based on hard pairwise `frac_recent_T` calls. For decoded results, each pair's
posterior mean TMRCA is thresholded before calculating the fraction; posterior
probabilities are not averaged.

## Method

At each candidate archaic marker, classify the existing haplotype pairs as
ALT/ALT, REF/REF or ALT/REF. These are allele states of arbitrary haplotypes,
not diploid homozygote classes. The primary confirmation score is

    J_T = wAA * max(fAA_T - fRR_T, 0) * (1 - fAR_T),

where wAA is the proportion of sampled pairs that are ALT/ALT and each f is
the corresponding `frac_recent_T`. Require the all-pair fraction at the same
coordinate and cutoff to exceed the median of neutral grid values estimated
from the fit subset. The primary rule uses T=50 kya and one supported bin.
This is a permissive local screen, rather than a strict 5% regional all-pair
gate, which would cap power at that first stage's detection rate.

Two spatial schemes were evaluated:

- **10 kb grid:** 1,000 regular positions, each assigned its nearest eligible
  archaic marker within 5 kb. A grid position without a marker cannot make a
  joint call. The original all-pair baseline still scans all grid positions.
- **Every archaic site:** TMRCA and pair classes evaluated at the exact position
  of every eligible observed archaic marker. Fixed ALT markers are retained
  when present, to retain complete sweeps.

Candidates obey the same oracle ascertainment in neutral and selected regions.
The causal site is not supplied to the detector. Sixteen neutral regions have
no eligible archaic markers; they remain in the denominator and receive no
joint call. Every selected region has at least one eligible marker.

The region score is the largest supported local score. Alternative rules
require two or three consecutive physical 10 kb bins: several SNPs in the
same bin cannot constitute a run. Different bins may use different markers,
so this does not assert that the same carrier haplotypes support the full run.

The full [method specification](joint_scan_method.md) also defines the simpler
carrier-mass and carrier-enrichment scores, gate/cutoff/run comparisons, empty
class rules, ancestry reconstruction checks, and calibration procedure.

## Held-out calibration

Five whole-replicate folds each use 400 fit neutrals, 400 separate calibration
neutrals, and 200 test neutrals. Each onset supplies 80 selected training regions
and 20 selected test regions per fold. Every region is tested once. Rank
p-values include ties: p=(1 + number of calibration scores at least as large
as the test score)/401. The complete site search is performed before ranking.

The trained joint detector chooses its score family, cutoff, gate and run
using training regions only; the trained all-pair baseline chooses its cutoff
and run the same way. Calibration and test regions do not choose those rules.
The primary J rule was frozen before examining the new joint outcomes. Prior
selected pair-class summaries from this same cohort informed method design,
so cross-validation does not make this an untouched validation set.

## Results at regional p <= 0.05

The primary rule gives nearly identical results with the two position schemes.
Each selected percentage is a count out of 100; each neutral percentage uses
1,000 regions.

| Positions | TMRCA source | Power, onset 50 kya | Power, onset 10 kya | Neutral regions called |
|---|---|---:|---:|---:|
| 10 kb grid | Truth | 81% | 16% | 5.2% |
| 10 kb grid | Decoded | 81% | 17% | 4.8% |
| Every archaic site | Truth | 79% | 15% | 5.2% |
| Every archaic site | Decoded | 80% | 17% | 4.6% |

The decoded comparisons explain where the gain comes from:

| Positions | Rule | Power, onset 50 kya | Power, onset 10 kya | Neutral regions called |
|---|---|---:|---:|---:|
| 10 kb grid | Trained all pairs | 13% | 4% | 5.5% |
| 10 kb grid | Archaic AF alone | 84% | 12% | 4.9% |
| 10 kb grid | Primary J | 81% | 17% | 4.8% |
| 10 kb grid | Carrier mass, T=50 kya | 83% | 17% | 4.8% |
| 10 kb grid | Carrier enrichment, T=50 kya | 83% | 17% | 4.6% |
| 10 kb grid | Trained joint | 82% | 18% | 4.7% |
| Every archaic site | Trained all pairs | 36% | 7% | 4.0% |
| Every archaic site | Archaic AF alone | 85% | 14% | 5.2% |
| Every archaic site | Primary J | 80% | 17% | 4.6% |
| Every archaic site | Carrier mass, T=50 kya | 84% | 17% | 4.9% |
| Every archaic site | Carrier enrichment, T=50 kya | 83% | 17% | 4.7% |
| Every archaic site | Trained joint | 82% | 17% | 4.8% |

Restricting the all-pair scan to archaic sites helps that baseline. This changes
both candidate placement and the set of positions searched; it is not a pure
increase in grid density. Once archaic pair information is included, the dense
site scan provides no clear gain in this cohort. The modest later-onset gains
over AF alone do not approach the desired 70% power at low neutral call rates.

## Thresholds trained to target 70% power

For each fold and onset, the threshold is the lower 30th percentile of the 80
selected training scores. The following are the resulting test-set outcomes,
not thresholds adjusted to force the test set to exactly 70% power.

| Positions | Decoded rule | Onset 50 kya: power / neutral calls | Onset 10 kya: power / neutral calls |
|---|---|---:|---:|
| 10 kb grid | Archaic AF alone | 70% / 1.4% | 74% / 48.1% |
| 10 kb grid | Primary J | 72% / 1.9% | 70% / 49.3% |
| 10 kb grid | Trained joint | 72% / 2.0% | 75% / 46.9% |
| 10 kb grid | Trained all pairs | 71% / 55.6% | 75% / 65.6% |
| Every archaic site | Archaic AF alone | 71% / 1.4% | 72% / 49.1% |
| Every archaic site | Primary J | 72% / 2.4% | 72% / 47.5% |
| Every archaic site | Trained joint | 66% / 1.8% | 70% / 45.5% |
| Every archaic site | Trained all pairs | 71% / 30.3% | 75% / 58.5% |

For the primary grid result, descriptive 95% Wilson intervals are 62.5–79.9%
for power and 1.22–2.95% for the neutral call fraction. For AF alone they are
60.4–78.1% and 0.84–2.34%. These intervals do not include the dependence caused
by shared fitted thresholds. With only 100 selected regions per onset, small
power differences should not be interpreted as a reliable ranking.

The early-onset AF thresholds vary from 0.7375 to 0.755 across grid training
folds. Thus this result is driven by common archaic alleles; it does not
establish detection of low-frequency partial sweeps.

## Time cutoffs and stronger spatial requirements

For the decoded primary family (median gate, one bin), at p<=0.05:

| T (kya) | Grid: onset 50 / 10 kya power | Grid neutral calls | Every site: onset 50 / 10 kya power | Every-site neutral calls |
|---:|---:|---:|---:|---:|
| 5 | 43% / 14% | 5.0% | 46% / 13% | 4.6% |
| 10 | 57% / 12% | 5.0% | 57% / 12% | 4.9% |
| 20 | 73% / 10% | 5.0% | 72% / 10% | 4.6% |
| 30 | 79% / 12% | 5.2% | 81% / 12% | 5.0% |
| 40 | 76% / 16% | 4.6% | 78% / 16% | 4.8% |
| 50 | 81% / 17% | 4.8% | 80% / 17% | 4.6% |

At T=50 kya, increasing required bin support does not produce a consistent gain:

| Positions | Consecutive 10 kb bins | Power, onset 50 kya | Power, onset 10 kya | Neutral calls |
|---|---:|---:|---:|---:|
| Grid | 1 | 81% | 17% | 4.8% |
| Grid | 2 | 79% | 18% | 4.9% |
| Grid | 3 | 77% | 15% | 4.7% |
| Every site | 1 | 80% | 17% | 4.6% |
| Every site | 2 | 82% | 17% | 4.8% |
| Every site | 3 | 78% | 17% | 5.0% |

Each run rule is recalibrated against neutral region scores. Consequently, a
longer-run requirement does not automatically lower the region FPR while
retaining the same nominal p threshold. It must improve the separation between
selected and neutral region scores. These comparisons show little improvement.

## The p=0.001 experiment

The held-out 400-neutral calibration cannot resolve p=0.001; its minimum is
1/401. A separate exploratory analysis calibrates fixed **ungated** rules
against all 1,000 neutrals. A selected region must exceed every neutral score,
giving p=1/1001. At onset 50 kya, ungated J at T=50 kya detects 19% (decoded
grid) and 21% (decoded every-site), compared with 33% and 28% using truth.
AF alone detects 41% and 42%, respectively. All these rules detect 0/100
later-onset selected regions at this threshold.

This is not an independent validation of a 0.1% neutral call rate. Neutral
leave-one-out ranks are constrained by construction, and the table explicitly
marks that limitation. These are ungated scores, not the primary gated rule.

## Conventional FDR and F1

Conventional discovery FDR additionally depends on the fraction pi of tested
regions that are selected. For power P and neutral call fraction Q:

    FDR = (1-pi)*Q / (pi*P + (1-pi)*Q)
    F1  = 2*pi*P / (pi*(1+P) + (1-pi)*Q)

For the decoded grid, early-onset thresholds trained toward 70% power:

| Rule | Assumed selected-region prevalence | Discovery FDR | F1 |
|---|---:|---:|---:|
| Primary J | 10% | 19.2% | 0.762 |
| AF alone | 10% | 15.3% | 0.767 |
| Primary J | 1% | 72.3% | 0.400 |
| AF alone | 1% | 66.4% | 0.454 |

These are illustrative prevalence calculations, not empirical genome-wide FDR
estimates. The full tables also give prevalence 50%, every tested cutoff, and
nominal regional p thresholds of 0.01, 0.05 and 0.10.

## Interpretation and limitations

A 10 kb grid is sufficient for the current pair-class method in this experiment.
The simple primary score offers an interpretable way to require recent
coalescence among archaic-allele carriers at a broadly supported all-pair
position. It meets the early-onset operating target. Nevertheless, the AF-only
comparison is essential: these results do not justify claiming that the
pairwise TMRCA term improves early-onset selection detection over archaic AF.

The later-onset difficulty persists with true TMRCA and oracle archaic labels.
It therefore cannot be attributed solely to Gamma-SMC decoding. The present
frequency-weighted regional scores and null overlap are also limiting. This
does not establish that every possible carrier-based method must fail. A future
alternative could calibrate carrier recency conditional on marker frequency
and carrier-pair count, with the complete marker search repeated in neutral
regions. That would require a frozen rule and new validation seeds.

The benchmark assumes perfect archaic-origin labels, fixed rates, an isolated
10 Mb region, survival/observation ascertainment, and the original
msprime-neutral / scaled-SLiM-selected simulation design. A matched simulator
control remains relevant before attributing every difference exclusively to
selection. The selected allele is aligned to the central grid position in the
original simulations. Candidate-site scans reduce dependence on that alignment.

Because these scores explicitly use archaic carrier frequency, the previous
all-pair statistic's lack of correlation with realized ancestry cannot be used
to dismiss purifying-selection effects on this new detector. No such robustness
claim is established here.

## Artifacts and verification

Large inputs and intermediates remain under
`D:/phase2simselection/sim/eas_joint_scan`. Compact checked results are in
[results/joint_scan](results/joint_scan):

- [All held-out metrics](results/joint_scan/heldout_metrics.csv) and
  [70%-target metrics](results/joint_scan/target70_metrics.csv).
- [10 kb grid figure, PDF](results/joint_scan/power_neutral_calls_grid10kb.pdf)
  and [PNG](results/joint_scan/power_neutral_calls_grid10kb.png).
- [Every-site figure, PDF](results/joint_scan/power_neutral_calls_archaic_sites.pdf)
  and [PNG](results/joint_scan/power_neutral_calls_archaic_sites.png).
- [Primary per-region predictions](results/joint_scan/primary_predictions.csv.gz),
  [training choices](results/joint_scan/training_choices.csv),
  [provenance](results/joint_scan/provenance.json), and
  [independent audit](results/joint_scan/audit.json).

The label and profile phases each completed 1,200/1,200 regions with zero
failures. Original decoded 10 kb counts remain exactly unchanged in every
region after adding output positions. Original true grid counts also match.
The independent audit passed for all 2,145,600 prediction rows, 3,576 metric
rows and 1,891,426 profiled grid/site positions. It recomputed every calibration
rank, verified fit/calibration/test separation, reconstructed the thresholds
trained toward 70% power, checked class partitions and artifact hashes, and
matched the independently calculated AF-only baseline.
Three pilot comparisons additionally verified bit-for-bit identical posterior
alpha/beta values at every original grid position. Native MRCA extraction was
checked against tskit across all marginal trees in a seeded fixture and against
the Python implementation on three full profiles. Focused tests cover pair
classes, empty/fixed cases, physical bins, conservative ranks, fold separation,
and an end-to-end synthetic benchmark with known outcomes.

Figures are saved without inline display, with editable vector PDF text.
Automated checks verify canvas bounds, legend separation from panels, and
absence of raster images in the PDFs.
