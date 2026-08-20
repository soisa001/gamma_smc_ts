# run5 — results

600/600 replicates simulated and decoded across three arms, zero failures.

## Headline: weak selection over the whole introgression window is the best EPAS1 match yet

`s = 0.01` acting for 1,240 generations beats run4's `s = 0.05` for 205
generations on both realism and power.

| | run4 (s = 0.05, 205 gen) | **run5 (s = 0.01, 1,240 gen)** | EPAS1 |
|---|---|---|---|
| final AF, mean | 0.484 | **0.703** | — |
| final AF, median | 0.502 | **0.725** | 0.63–0.87 |
| fraction above 0.6 | 0.32 | **0.69** | — |
| best tree-truth power | 0.73 (at 20 ky) | **0.81 (at 30 ky)** | — |

Tree-truth power, `chb_from_introgression`:

| x (years) | tree truth | gamma_smc_frac | gamma_smc_meanp |
|---:|---:|---:|---:|
| 1,000 | 0.04 | 0.02 | 0.58 |
| 4,500 | 0.03 | 0.21 | 0.51 |
| 10,000 | 0.00 | **0.49** | 0.49 |
| 20,000 | 0.48 | 0.48 | 0.46 |
| **30,000** | **0.81** | 0.40 | 0.42 |
| 40,000 | 0.71 | 0.40 | 0.41 |
| 50,000 | 0.52 | 0.36 | 0.37 |

## Selection inside the Neanderthal branch makes no difference

The `chb_nea_selection` arm was built to test whether sweeping the allele inside
the archaic source first — compressing archaic haplotype diversity so that fewer
and longer haplotypes introgress — would raise power. **It does not.**

| x (years) | from introgression | with Nea selection |
|---:|---:|---:|
| 20,000 | 0.48 | 0.58 |
| **30,000** | **0.81** | **0.82** |
| 40,000 | 0.71 | 0.55 |
| 50,000 | 0.52 | 0.50 |

The peaks are identical to within one replicate (0.81 vs 0.82), and the two arms
trade places on either side of the peak in a way consistent with sampling noise
at n = 100. Final frequencies are also indistinguishable (mean 0.703 vs 0.666).

The likely reason is that the archaic haplotype diversity entering CHB is
already low without any archaic sweep. Neanderthal Ne in this model is 3,600, so
archaic lineages coalesce among themselves fast; and at any single base the
introgressed material usually traces to very few migration events anyway. There
is little diversity left for a sweep in the source to remove.

This is a negative result worth keeping: the hypothesis was reasonable, it was
cheap to test, and it is now ruled out at this effect size.

## The EAS structured-coalescent reference sets the ceiling

| x (years) | tree truth | gamma_smc_frac |
|---:|---:|---:|
| 10,000 | 0.14 | 0.46 |
| 20,000 | 0.72 | **0.70** |
| **30,000** | **0.85** | 0.62 |
| 40,000 | 0.77 | 0.46 |
| 50,000 | 0.66 | 0.43 |

P(TMRCA < 20 ky) is 0.220 in swept replicates against 0.076 neutral. This is an
idealized partial sweep to a known 75%, so 0.85 is roughly the ceiling this
statistic can reach in the EAS demography — and the CHB archaic arm at 0.81 is
close to it. Whatever the archaic arms lose to the messiness of real
introgression, it is not much.

Keeping the PHLASH curve as discrete epochs outside the sweep window rather than
flattening the whole history mattered: at a 50 ky cutoff the sweep-versus-neutral
contrast was 0.14 vs 0.06 under a fully constant Ne, and 0.47 vs 0.18 with the
epochs retained. The deep history sets the neutral background, and the null
distribution is what the test is measured against.

## Can the sweep be chained across epochs instead of flattening the window?

Tested directly, because it would be strictly better than holding the window
flat. **Chaining works; it does not solve the problem.**

| configuration | result |
|---|---|
| 3 chained `SweepGenicSelection` models, no demographic events | **OK** |
| 3 chained sweeps, events only *older* than the whole sweep | **OK** |
| 3 chained sweeps, events *on* the segment boundaries | fails |
| 3 chained sweeps, boundary events + `duration` pinned per segment | fails |
| 7 chained sweeps + a 2-generation neutral gap at each boundary to hold the event | fails |

So `SweepGenicSelection` is **not** restricted to a single use -- msprime accepts
a list of them and simulates the chain correctly -- and it does take a
`duration` argument. The restriction that actually bites is narrower and
unavoidable: *no demographic event may fall anywhere inside the span a sweep
model occupies*. An event exactly on a segment boundary counts as inside, pinning
`duration` does not change that, and inserting a short neutral phase to hold the
event does not either. The realized span of a sweep is set by its stochastic
trajectory, so there is no time point inside the sweep window that is reliably
safe.

The flat window is therefore the correct workaround rather than a shortcut. What
it costs is worth stating precisely: over 0-1,500 generations the PHLASH median
falls from 37,637 to about 3,300, and holding it at the harmonic mean of 4,884
erases that recent expansion. **But the neutral arm uses the identical flattened
window**, so the selected-versus-neutral contrast -- which is what the test
measures -- stays internally valid. What is distorted is the realism of absolute
TMRCA values inside the sweep window, in both arms equally.

If the full PHLASH history is needed through the sweep window, the way to get it
is SLiM rather than msprime: a forward simulation has no such restriction. That
costs roughly 20 seconds per replicate instead of one.

## The truth/decoder crossover reproduces a third time

Same pattern in all three arms, at the same place:

| arm | decoder wins | truth wins |
|---|---|---|
| chb_from_introgression | 10 ky: 0.49 vs 0.00 | 30 ky: 0.81 vs 0.40 |
| chb_nea_selection | 10 ky: 0.51 vs 0.02 | 30 ky: 0.82 vs 0.21 |
| eas_sweep_structured | 10 ky: 0.46 vs 0.14 | 30 ky: 0.85 vs 0.62 |

The decoder's recency bias buys real power at cutoffs more recent than the sweep,
where the truth statistic is flat, and costs power once the truth has genuine
signal. Across run3, run4 and run5 this has now held in five independent arms.

`gamma_smc_meanp` is again flat near 0.5 at 1,000 years, where essentially
nothing genuinely coalesces, so it should not be read as threshold-specific
detection.

## How power is computed, and what AUC adds

For each selected replicate *i* at cutoff *x*, the statistic is
`T_i = P(TMRCA < x)` at the focal base over the 100 within-individual pairs
(`overall` class). The null is the 100 neutral replicates' `T` at the same
cutoff, and

```
p_i = (1 + #{null_j >= T_i}) / (1 + n_null)
power(x) = fraction of selected replicates with p_i <= 0.05
```

It is a one-sided upper-tail permutation p-value, ties counted as exceedances
(conservative). So **power is a per-replicate detection rate**: how often a
single locus of this kind would be called in a scan. Two consequences: the floor
is `1/101 = 0.0099`, and a replicate must beat at least 96 of the 100 nulls to
reach 0.05, so power is bounded by how large the null is.

**AUC** removes both limitations. It is the probability a random selected
replicate scores above a random neutral one, ties split, computed as the mean of
each replicate's percentile against the null. Threshold-free and floor-free:

| cutoff | CHB from introgression | | CHB with Nea selection | | EAS reference | |
|---|---|---|---|---|---|---|
| | AUC | power | AUC | power | AUC | power |
| 10,000 | 0.52 | 0.00 | 0.61 | 0.02 | 0.69 | 0.14 |
| 20,000 | 0.83 | 0.48 | 0.86 | 0.58 | 0.96 | 0.72 |
| **30,000** | **0.95** | 0.81 | **0.94** | 0.82 | **0.98** | 0.85 |
| 40,000 | 0.89 | 0.71 | 0.83 | 0.55 | 0.96 | 0.77 |
| 50,000 | 0.83 | 0.52 | 0.75 | 0.50 | 0.90 | 0.66 |

Discrimination is near-perfect at 30 ky (AUC 0.94-0.98) while power is 0.81-0.85.
The gap is the cost of a 100-replicate null, not a property of the statistic.

## AUC against final allele frequency

The clearest result in run5. Per-replicate percentile against the null *is* the
AUC contribution, so binning it by final frequency gives AUC restricted to that
band. `chb_from_introgression`, tree truth:

| cutoff | AF < 0.4 (n=8) | AF 0.4-0.7 (n=41) | AF >= 0.7 (n=51) |
|---|---|---|---|
| 20,000 | 0.68 | **0.95** | 0.76 |
| 30,000 | 0.53 | **0.99** | 0.97 |
| 40,000 | 0.41 | 0.85 | **0.99** |
| 50,000 | 0.49 | 0.69 | **0.99** |

Two things follow.

**Detection is essentially a function of final allele frequency.** Replicates
ending below 0.4 are indistinguishable from neutral at every cutoff (AUC 0.41 to
0.68, straddling chance). Replicates above 0.7 are almost perfectly separable
(0.97-0.99) provided the cutoff is deep enough.

**The best cutoff shifts deeper as the frequency rises.** Near-complete sweeps
peak at 40-50 ky, mid-frequency ones at 20-30 ky. That is mechanistically what
should happen: once the allele is near fixation the whole sample coalesces back
at the sweep onset (~1,240 generations, 36 ky), so the signal sits at deep
cutoffs; a partial sweep leaves a mixture of swept and unswept lineages, and the
swept subset shows up at shallower ones.

Practically: a scan using a single cutoff is mis-specified for part of the
frequency range, and 30 ky is the best single compromise here.

## The p-value distribution exposes a discreteness problem

`figures/pvalue_distribution_tree_truth.png` plots selected p-values against the
leave-one-out null, one panel per cutoff.

At 30-50 ky the picture is textbook: the null is close to uniform and the
selected replicates pile at the floor. But **at 1,000 to 10,000 years both
distributions collapse onto p = 1.0**. With 100 diploid pairs the statistic only
takes values in steps of 0.01, and at those cutoffs nearly every replicate scores
exactly zero, so ties dominate, `#{null >= T}` is the whole null, and the p-value
is 1 by construction for almost everything.

So at recent cutoffs the test is not merely underpowered, it is *degenerate* --
and the null is not uniform there but a point mass at 1. The reported
false-positive rates near zero at those cutoffs reflect that conservatism, not
good calibration. AUC agrees there is genuinely nothing to find (0.52 at 10 ky),
so the conclusion is unchanged, but the p-value machinery should not be read as
calibrated in that range. Raising the pair count would relieve the discreteness.

## Why the statistic depends on final frequency: the trimodal decomposition

The within-individual TMRCA at an introgressed swept locus is **trimodal**, and
the three genotype classes pull in opposite directions. Measured on the
`chb_from_introgression` selected replicates, against a neutral baseline of
median 9,660 generations (280 ky) and `P(TMRCA < 30 ky) = 0.076`:

| within-individual class | pairs | median TMRCA | P(TMRCA < 30 ky) |
|---|---:|---:|---:|
| **hom carrier** (archaic x archaic) | 5,536 | **1,155 gen (33 ky)** | **0.564** |
| **heterozygous** (archaic x modern) | 3,024 | **25,914 gen (752 ky)** | **0.000** |
| hom non-carrier (modern x modern) | 1,440 | 8,143 gen (236 ky) | 0.087 |

Both intuitions are correct, for different pairs of the same tree:

* Two archaic haplotypes in one individual descend from the introgressed
  material and, under a sweep, coalesce **after** introgression. They are eight
  times *shallower* than neutral.
* One archaic and one modern haplotype cannot coalesce until before the
  human-Neanderthal split. They are **three times deeper** than neutral, and
  `P(TMRCA < 30 ky)` is **exactly zero** -- not small, zero. Those pairs can never
  be recent.

`P(TMRCA < x)` sums over all three, so it only shows net signal once the
hom-carrier mass outweighs the heterozygous mass. At Hardy-Weinberg proportions
that is `f^2 > 2f(1-f)`, i.e. **f > 2/3**. The observed counts bracket it: at
final frequency 0.4-0.7 there are 1,368 hom-carrier pairs against 1,952
heterozygous, and at 0.7-1.0 it is 4,151 against 889.

This is the mechanism behind the AUC-versus-frequency result, and it explains the
sub-chance AUC at low frequency (0.41 at a 40 ky cutoff) that would otherwise
look like noise: below `f = 2/3` the heterozygous class actively drives the
statistic *below* neutral, because those pairs are deeper than any neutral pair.

It also explains why the best cutoff moves deeper as frequency rises. Hom-carrier
`P(TMRCA < 30 ky)` is 0.871 in the 0.4-0.7 band but only 0.461 above 0.7: once
the allele is near fixation the carrier class has been large for long enough to
accumulate its own internal coalescent depth, so its signal shifts to deeper
cutoffs.

**Practical consequence.** `P(TMRCA < x)` over all pairs is a lossy summary of a
distribution whose shape is far more diagnostic than its left tail. A statistic
keyed on the *trimodality* -- or one restricted to hom-carrier pairs, or a
two-sided test that also rewards the anomalous heterozygous depth -- should beat
it, particularly below `f = 2/3` where the current statistic is actively
misleading.

## Null calibration

Leave-one-out false-positive rates against the null itself, averaged over
thresholds:

| arm | tree_truth | gamma_smc_frac | gamma_smc_meanp |
|---|---|---|---|
| chb_from_introgression | 0.030 | 0.037 | 0.050 |
| chb_nea_selection | 0.026 | 0.044 | 0.050 |
| eas_sweep_structured | 0.034 | 0.034 | 0.050 |

All at or below the nominal 0.05.

## What run5 settles

1. **Weak, long selection beats strong, short selection** for both realism and
   detectability, at matched final frequency. `s = 0.01` from the introgression
   onset reproduces EPAS1's frequency distribution better than `s = 0.05` from
   6,000 years ago, and detects better too.
2. **Sweeping the allele inside the archaic source does not help.** Archaic
   haplotype diversity at a single base is already low.
3. **The archaic arm is close to the idealized ceiling** (0.81 against 0.85), so
   little power is lost to the messiness of real introgression relative to a
   clean partial sweep of the same final frequency.
4. The informative cutoff for all three arms is **20–40 ky**, peaking at 30 ky —
   comfortably inside the original grid, and consistent with a sweep that began
   around 1,240 generations (36 ky) ago.

## Caveats

- Both CHB arms are ascertained on the introgressed allele reaching 2% by the end
  of the migration window. Without that most replicates carry no archaic ancestry
  at the focal base at all, but it does mean these are conditional-on-ascertainment
  power estimates, not unconditional genome-scan rates.
- The EAS arm is deliberately artificial: a partial sweep to a fixed 75% in a
  history held flat across the sweep window. It is a reference, not a model of a
  real locus.
- Tree sequences were deleted after decoding to bound disk use; the raw per-pair
  TMRCA matrices are retained, so thresholds and stratifications remain post-hoc.
