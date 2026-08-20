# run4 — results (EPAS1-calibrated)

600/600 replicates simulated and decoded across three arms, zero failures.

## Headline: the EPAS1 calibration works, and the original threshold grid is now correct

run2 concluded that its 1–50 ky cutoff grid sat *below* the admixture floor and
therefore saw almost nothing. run4 tests the prediction that followed: recalibrate
the sweep to EPAS1 and the signal moves back inside that grid. It does.

**CHB archaic arm** (allele fixed in Neanderthal, introgressed, selection on the
introgressed haplotype class from 6,000 years ago at s = 0.05):

| x (years) | tree truth | gamma_smc_frac | gamma_smc_meanp |
|---:|---:|---:|---:|
| 1,000 | 0.04 | 0.00 | 0.52 |
| 4,500 | 0.06 | **0.59** | 0.47 |
| 10,000 | 0.66 | 0.58 | 0.37 |
| **20,000** | **0.73** | 0.27 | 0.23 |
| 30,000 | 0.65 | 0.20 | 0.21 |
| 40,000 | 0.43 | 0.13 | 0.15 |
| 50,000 | 0.23 | 0.10 | 0.13 |

**0.73 power at 20,000 years**, inside the original grid, against 0.20 at best
anywhere in run3. The underlying separation is large: P(TMRCA < 20 ky) is 0.180
in selected replicates against 0.022 in neutral ones, an eight-fold difference.

Final allele frequency across the 100 selected replicates: mean **0.484**, median
**0.502**, interquartile range 0.239–0.735, with **63%** above 0.4 and **32%**
above 0.6. EPAS1's observed 0.63–0.87 sits in the upper part of that
distribution, so the calibration lands in the right regime, if a little low at
the median.

## The de novo "positive controls" are the weakest arms, and that is the result

| arm | final AF (mean) | best power, any estimator |
|---|---|---|
| CHB archaic, single founder | 0.484 | **0.73** |
| CHB de novo | 0.038 | 0.11 |
| EAS de novo | 0.075 | 0.16 |

A single de novo copy at s = 0.05 has only 205 generations to work with, which
takes it from 1/2N to roughly 1.3% deterministically. The simulations agree:
means of 0.038 and 0.075. There is essentially nothing for any statistic to
detect.

This inverts the usual expectation that a hard de novo sweep is the easiest case.
Over a 6,000-year window it is the *hardest*, because it starts from one copy.
**Adaptive introgression is detectable here precisely because it does not start
from one copy** — the archaic haplotype class is already at roughly the archaic
ancestry proportion when selection begins, having drifted neutrally for tens of
thousands of years since admixture. That head start, not the archaic haplotype
background as such, is what makes the signal.

Deterministic check, 205 generations at s = 0.05:

| starting frequency | final |
|---|---|
| single copy (1/2N ≈ 7.5e-5) | 0.013 |
| 0.5% | 0.458 |
| 1.5% (archaic class at onset) | 0.719 |
| 3% | 0.839 |

## Truth versus Gamma-SMC: the crossover is real and reproduces

The archaic arm shows the same pattern run3 found, at the same place:

| x (years) | truth mean | decoded mean | bias | Pearson r | winner |
|---:|---:|---:|---:|---:|---|
| 4,500 | 0.005 | 0.044 | +0.040 | 0.12 | **decoder** (0.59 vs 0.06) |
| 10,000 | 0.049 | 0.166 | +0.117 | 0.42 | truth (0.66 vs 0.58) |
| 20,000 | 0.180 | 0.295 | +0.115 | 0.65 | **truth** (0.73 vs 0.27) |
| 50,000 | 0.406 | 0.435 | +0.028 | 0.91 | truth (0.23 vs 0.10) |

The decoder over-calls recent coalescence, and the bias is largest exactly where
its correlation with the truth is weakest. That bias buys it real power at 4,500
years, where truth is flat (0.59 vs 0.06), and costs it power from 10,000 years
upward once the truth statistic has genuine signal to work with. By 50 ky the
bias is gone and the correlation reaches 0.91.

So the run2 prediction that inference must lose was wrong at the recent end and
right at the deep end, in both run3 and run4.

`gamma_smc_meanp` is again flat near 0.5 at 1,000 years where nothing genuinely
coalesces, which remains a reason not to read it as threshold-specific detection.

## Null calibration

Leave-one-out false-positive rates against the null itself, averaged over
thresholds — all at or below the nominal 0.05:

| arm | tree_truth | gamma_smc_frac | gamma_smc_meanp |
|---|---|---|---|
| CHB archaic | 0.030 | 0.036 | 0.050 |
| CHB de novo | 0.026 | 0.044 | 0.050 |
| EAS de novo | 0.030 | 0.047 | 0.050 |

run4's neutral arms carry no focal allele at all, so there is no conditional
analysis and none of run3's conditional-null starvation.

## What this settles

1. The admixture-floor diagnosis from run2 was correct and actionable: moving the
   sweep onset inside the cutoff grid raised power from ≤0.20 to 0.73.
2. Detectability here is driven by **starting frequency and sweep age**, not by
   the archaic haplotype background. The archaic arm wins because it starts at
   ~1.5%, not because its carriers are divergent.
3. Gamma-SMC's advantage over tree truth is confined to cutoffs more recent than
   the sweep, and comes from a systematic recency bias rather than better
   inference.

## Caveats

- The archaic arm's median final frequency (0.50) is below EPAS1's 0.63–0.87.
  Matching more closely would need a slightly larger `s` or an earlier onset.
- Tree sequences were deleted after decoding to fit the host disk, so re-decoding
  requires re-simulating. The raw per-pair TMRCA matrices are retained, so all
  threshold, window and stratification work remains post-hoc.
- Per-replicate raw data lives only in the WSL working copy and is gitignored;
  the seeds in the run4 config regenerate it.
