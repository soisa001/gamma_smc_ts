# run2 — results

400/400 replicates completed, **zero failures**, 23.2 core-hours (~70 min wall
at 20 workers on WSL2/Ubuntu-22.04, SLiM 4.2.2, stdpopsim 0.3.0).

## Bottom line

**Yes — selection at the introgressed locus is detectable, and the effect is
enormous. But the threshold grid this study was specified with is too recent to
see it, so the headline power numbers understate the signal by a wide margin.**

Within-individual TMRCA at the focal base, median across replicates:

| Arm | selected | neutral | separation |
|---|---|---|---|
| CHB | **2,182 gen (54.6 ky)** | 22,515 gen (562.9 ky) | **10.3×** |
| EAS | **2,158 gen (53.9 ky)** | 10,522 gen (263.1 ky) | **4.9×** |

Restricted to replicates where the allele fixed, the selected median is
2,165 gen (CHB) and 2,148 gen (EAS) — both essentially **the age of the sweep
onset, 2,265 generations**. That is exactly what fixation from standing
variation should do: it compresses the genealogy at the focal base back to the
founding of the sweep and no further. Replicates that lost the allele sit at
22,040 gen (CHB) and 12,905 gen (EAS), indistinguishable from neutral — the
internal control works.

## The problem with the threshold grid

The sweep begins 2,265 generations ago = **56,625 years**. The requested
thresholds stop at **50,000 years = 2,000 generations**, which lands *just below*
the selected distribution's median of ~2,160 generations. The statistic is
therefore measured entirely in the tail where the two distributions have barely
begun to separate.

Power at the requested thresholds:

| x (years) | CHB median obs. | CHB null mean | **CHB power** | EAS median obs. | EAS null mean | **EAS power** |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.000 | 0.003 | 0.03 | 0.000 | 0.000 | 0.04 |
| 4,500 | 0.010 | 0.013 | 0.04 | 0.000 | 0.002 | 0.02 |
| 10,000 | 0.030 | 0.028 | 0.06 | 0.010 | 0.007 | 0.01 |
| 20,000 | 0.060 | 0.055 | 0.08 | 0.040 | 0.041 | 0.04 |
| 30,000 | 0.090 | 0.084 | 0.06 | 0.100 | 0.099 | 0.05 |
| 40,000 | 0.140 | 0.111 | 0.18 | 0.180 | 0.152 | 0.20 |
| 50,000 | 0.235 | 0.132 | **0.67** | 0.365 | 0.200 | **0.59** |

Below 40 ky the selected and neutral distributions are **the same to within
noise** — the median observed value tracks the null mean almost exactly, and
power sits at the nominal false-positive rate. Only at the last threshold does
signal appear, and both power curves are still climbing steeply at the right
edge of the grid (`cross_arm/figures/arm_comparison.png`).

The spatial profiles say the same thing more starkly: at 4,500 y and 20,000 y
there is **no localisation at the focal base at all** — selected and neutral
traces overlap across the whole 10 Mb
(`chb_ancient_eurasia/figures/spatial_profile.png`). The expected "visual spike"
is absent at those thresholds because it does not exist there.

## Selected replicate against the neutral distribution

`<arm>/figures/selected_vs_neutral_band_by_threshold.png` puts one selected
replicate against **all 100 neutral replicates** (median and 2.5-97.5% band),
one panel per cutoff. A single neutral trace is unusable as a reference -- with
100 diploid pairs it moves in steps of 0.01 -- so the null side is the whole
distribution. The selected side stays a single replicate on purpose: averaging
it would blend the sweep with the ~20% of replicates that lost the allele and
are neutral by construction. Every curve including the band edges gets the same
250 kb rolling mean, so the comparison is like for like.

Representatives are the replicates nearest their mode's median focal statistic
at 50 ky, with the selected pick restricted to replicates in which the allele
fixed: CHB rep001, EAS rep002.

- **Below 50 ky, neither arm leaves the neutral band at 5 Mb.** The selected
  trace wanders inside the band across the whole contig. Where it does breach
  the band -- CHB near 1.3 Mb at 10-30 ky, EAS near 0.6 and 2.5 Mb at 20-40 ky --
  it is nowhere near the focal base. Those are the false positives a scan would
  produce at these cutoffs.
- **At 50 ky EAS shows a clean, localised breach at exactly 5 Mb**, rising to
  ~0.40 against a neutral median of ~0.19 and a band top of ~0.38, and decaying
  within a few hundred kb. It is the only decisive excursion on the contig.
- **CHB at 50 ky does not localise.** It reaches ~0.26 at the focal base, only
  touching the band edge, while an unrelated bump near 7.8 Mb (~0.33) breaches
  it more clearly. A single CHB replicate would not be called at 5 Mb.

That CHB is the weaker arm per replicate matters, because CHB is the arm of
interest. It is consistent with the arms' different starting structure: CHB
carriers begin as ~15 complete archaic homozygotes whose own common ancestry is
deep in the Neanderthal lineage, so fixing them leaves a genealogy with several
old internal branches, whereas EAS carriers are drawn from the resident
population and coalesce more tightly. The current data do not settle this --
testing it needs the raw TMRCA distributions, which were not stored.

## The test machinery itself is sound

Leave-one-out false-positive rate of the null against itself
(`results/null_calibration.tsv`):

| x (years) | 1,000 | 4,500 | 10,000 | 20,000 | 30,000 | 40,000 | 50,000 |
|---|---|---|---|---|---|---|---|
| CHB | 0.01 | 0.04 | 0.01 | 0.04 | 0.05 | 0.04 | 0.04 |
| EAS | 0.02 | 0.01 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |

Every value sits at or below the nominal 0.05. The p-value construction is
calibrated and is not over-firing, so the low power is a property of the
statistic's threshold range, not a bug in the test.

## Allele-frequency outcomes

| Arm | lost | segregating | fixed |
|---|---|---|---|
| CHB | 0.16 | 0.08 | 0.76 |
| EAS | 0.13 | 0.19 | 0.68 |

CHB frequency entering Han: mean 0.0904, median 0.0794, and 13% of replicates
had already lost the allele before Han was founded. See `FRAMEWORK.md` §5 — the
placement frequency of 0.0296 is in the Han *ancestor*, not in Han.

The independent Wright-Fisher forecast written before the run predicted CHB
placement AF 0.0296 ± 0.0078 (observed 0.0298 ± 0.0081), entry AF mean 0.0924
(observed 0.0904), loss before Han 14.3% (observed 13.0%), final loss 18.4%
(observed 16.0%) and final fixation 71.7% (observed 76.0%). The agreement is a
two-way cross-check: neither the tick schedule, the migrant-based placement, nor
the fitness callbacks are silently wrong.

## Why the signal sits where it does: the admixture floor

Zhang et al., *Recovering signatures of archaic hominin introgression using
ancestral recombination graphs*, Science (2026), doi:10.1126/science.aef8874,
makes the governing point explicit. Their method (TRACE) detects introgression
because "introgression from a deeply divergent population introduces lineages
that coalesce much further back in time than nonintrogressed lineages... branches
with deep coalescence times that span the interval between the divergence time
of the archaic and modern lineages (T_archaic) and the time of the admixture
event (T_admix)". They run it with `t = 15,000 generations` (420 ky at 28 y/gen)
as the branch-length cutoff.

So the published, validated signature of introgression is **deep** coalescence in
`[T_admix, T_archaic]`. run2's `P(TMRCA < x)` looks for **shallow** coalescence.
Those are opposite ends of the same tree, and in an adaptively introgressed
region they are in direct tension.

Expected TMRCAs under `AncientEurasia_9K19` itself (msprime, 4,000 replicates):

| Pair | mean | median |
|---|---|---|
| archaic at the pulse vs present-day Han | 62,722 gen (1.57 My) | 52,433 gen |
| two archaic lineages at the pulse | 18,533 gen (463 ky) | 2,796 gen (69.9 ky) |
| two present-day Han lineages, neutral | 29,772 gen (744 ky) | 14,584 gen (365 ky) |

The Neanderthal Ne trajectory explains the enormous skew in the second row: it is
192 at the pulse, so archaic lineages coalesce fast if they coalesce at all, but
it balloons to 18,200 by 3,832 generations, so the ones that miss wait hundreds
of thousands of years.

**There is therefore a hard floor.** Every carrier of a fixed introgressed allele
traces back through the pulse, so a pair either coalesces *after* the pulse
(< 2,270 gen = 56.8 ky, only possible because of the sweep) or *before* it, at a
median of ~70 ky and a mean of ~463 ky. **No cutoff below 56.8 ky can capture
anything except sweep-driven, post-pulse coalescence.** All seven requested
cutoffs are below the floor.

That makes the shape of the result inevitable, and the pile-up is extreme:

| Arm | P(TMRCA < 50 ky) | median TMRCA | mass between 50 ky and the median | width of that window |
|---|---|---|---|---|
| CHB | 0.235 | 2,182 gen (54.6 ky) | **0.265** | 182 generations |
| EAS | 0.365 | 2,158 gen (53.9 ky) | **0.135** | 158 generations |

A quarter of all CHB pairs coalesce in the 182 generations between the largest
requested cutoff and the median. Moving the cutoff from 50 ky to ~57 ky — the
floor itself — should take the CHB statistic from 0.235 to above 0.50. That is a
concrete prediction the extended grid would test.

It also explains why CHB detects worse than EAS per replicate despite nearly
identical medians. The two arms differ in their *tails*, not their centres: a CHB
pair that misses the sweep falls back into the archaic lineage (median 70 ky,
mean 463 ky), whereas an EAS pair falls back into the EAS population (neutral
median 10,522 gen vs 22,515 for CHB). Near the floor, EAS has more mass just
above the cutoff, so it crosses first.

## Does fixing the allele in the archaic source distort the truth?

Short answer: **no, and the choice is conservative — but the number of founding
haplotypes does matter, and that is the real modelling commitment.**

1. **Conditional on fixation, the origin time is irrelevant.** If the allele is
   fixed in Neanderthal at the pulse, every archaic haplotype carries it, so the
   carrier genealogy simply *is* the archaic population genealogy. Whether
   fixation happened 3,000 or 20,000 generations ago changes nothing about the
   transmitted material's ancestry. The concern that fixation "could plausibly
   have occurred anytime between the split and introgression" does not bias the
   TMRCA truth.
2. **Fixation is the conservative assumption, not an optimistic one.** Carriers
   of a segregating derived allele form a clade in the gene tree at that locus,
   so their MRCA is at most the population MRCA. Fixation is the limiting case in
   which the clade is the entire population, which gives the *deepest possible*
   carrier genealogy. Assuming a segregating source frequency would push carrier
   TMRCA down and make detection easier, not harder.
3. **The real artifact is that this is a very soft sweep.** Marking every
   introgressing genome starts the sweep from ~30 distinct archaic haplotypes
   whose mutual coalescence is set by the archaic Ne. Real adaptive introgression
   usually means one favoured allele on one archaic haplotype background — a hard
   sweep, in which every carrier coalesces at or after the pulse and TMRCA is
   strictly below the floor. run2 sits at the soft extreme, and a single-founder
   variant would produce a far shallower, far more detectable signal. The two
   bracket the real case, and only the soft end has been simulated.

Nothing here is a bug in the implementation; these are properties of the model
that was specified. But they mean run2 as run measures how well a *soft*
introgressed sweep is detected by a *shallow-coalescence* statistic evaluated
*below the admixture floor* — three choices that each work against detection.

## Would Gamma-SMC have more power than tree truth?

Tested directly, and the answer is no — but the reasoning is worth recording
because the intuition behind the question is sound.

**The intuition is right about the data.** Archaic haplotypes are similar to each
other (Neanderthal Ne is 192 at the pulse), so once the introgressed allele
fixes, heterozygosity at the focal base collapses: 1 het per 18 kb in the
selected arm against 1 per 1.8 kb neutral, a 10x contrast that is glaring in
sequence data. Gamma-SMC will certainly see it.

**But that recency is already in the truth.** The archaic-archaic median TMRCA is
2,796 generations (69.9 ky), which is genuinely recent against the 22,515
generation (563 ky) background. Gamma-SMC would infer roughly that same value.
There is no hidden extra signal for an estimator to recover — truth and estimate
agree that the swept region is ~55-70 ky deep, and the problem is that 55-70 ky
sits above the 50 ky cutoff either way.

**Three quantitative reasons inference should lose, not gain:**

1. *Almost no mutations exist at these depths.* Expected heterozygous sites per
   diploid pair available to infer the TMRCA:

   | TMRCA | in 10 kb | in 100 kb | in 200 kb |
   |---|---|---|---|
   | 1 ky cutoff (40 gen) | 0.01 | 0.10 | 0.20 |
   | 10 ky cutoff (400 gen) | 0.10 | 1.00 | 2.00 |
   | 50 ky cutoff (2,000 gen) | 0.50 | 5.00 | 10.00 |
   | CHB selected median (2,182 gen) | 0.55 | 5.46 | 10.91 |

   At the recent cutoffs there is essentially nothing to infer from. Even at the
   largest cutoff a pair contributes ~10 mutations across the entire sweep block.

2. *The prior will dominate and compress the contrast.* The repo's decoder uses
   `scaled_mutation_rate = 0.000315`, i.e. Ne = 6,300 and a prior mean TMRCA of
   12,600 generations. The selected arm's truth (2,182) sits 5.8x below that and
   the neutral arm's (22,515) 1.8x above, so shrinkage pulls both toward 12,600
   and squeezes the separation from both sides.

3. *Published precedent points the same way.* Zhang et al. report TRACE recall of
   ~80% on true ARGs, ~50% with SINGER-inferred ARGs and <10% with Relate — ARG
   inference retains at most about 60% of the truth's recall on precisely this
   problem.

**The one mechanism that could help — spatial pooling — was tested and is worth
little here.** Gamma-SMC's HMM pools across positions, and the sweep signal is
spatially extended while the sampling noise is local. Averaging the truth
statistic over windows centred on the focal base emulates that:

| window | CHB 30 ky | CHB 40 ky | CHB 50 ky | EAS 30 ky | EAS 40 ky | EAS 50 ky |
|---|---|---|---|---|---|---|
| point | 0.06 | 0.18 | **0.67** | 0.05 | 0.20 | 0.59 |
| +/-50 kb | 0.10 | 0.18 | 0.57 | 0.07 | 0.29 | **0.66** |
| +/-100 kb | 0.13 | **0.25** | 0.64 | 0.09 | **0.32** | 0.65 |
| +/-250 kb | 0.04 | 0.21 | 0.50 | 0.15 | 0.29 | 0.64 |
| +/-1 Mb | 0.07 | 0.17 | 0.43 | 0.11 | 0.12 | 0.46 |

Pooling buys at most ~+0.12 (EAS at 40 ky, 0.20 to 0.32) and actively *hurts* at
the best cutoff for CHB (0.67 down to 0.43-0.64). The optimal window is
+/-50-100 kb, which matches the hard-sweep block half-width
`s / (r * ln(2Ns)) / 2` = 103 kb; wider windows dilute the block with unswept
flanks.

**Conclusion.** Gamma-SMC is worth running because it is what one would apply to
real data, where the true ARG is unavailable. The right question for pass 2 is
"how much of the truth's power does it retain", and the expected answer is
somewhere below 100%, not above it. It is not a substitute for putting the cutoff
above the admixture floor.

## What should change

1. **Extend the threshold grid past the sweep onset.** Add 60, 75, 100 and
   150 ky. At 60 ky the selected median should cross 0.5 while the neutral mean
   stays near 0.15, so power should approach the fixation rate (~0.76 CHB,
   ~0.68 EAS) rather than 0.67/0.59. The current grid cannot show this.
2. **Store the 100 raw focal TMRCA values per replicate.** They were not saved,
   so new thresholds currently require re-simulating all 400 replicates. One
   hundred floats per replicate makes every future threshold a post-hoc
   calculation instead of a 23-core-hour re-run. This is the single change with
   the best cost/benefit.
3. **Re-plot the spatial profile at an informative threshold** (60-75 ky rather
   than 4.5/20 ky), where a localised spike at 5 Mb should actually appear.
4. Only then consider Gamma-SMC decoding — there is no point asking whether the
   decoder recovers a signal at thresholds where the truth shows none.

Items 1-3 are one re-run of ~70 minutes.

## Provenance

- `logs/simulation_status.tsv` — 400 rows, all `completed`
- `validation_report.json` — gates G1, G2 (static), G2b, G4, G5 passing
- Runtime gate G2 confirmed in-simulation: CHB placement used real migrant
  genomes (`inds.migrant`), 30/936 genomes on replicate 0
- Focal-base guard fired on every replicate (2 overlay calls intercepted, 1
  masked — stdpopsim excludes the background DFE itself, run2 masks recapitation)
- Per-replicate intermediates are gitignored; seeds in `config/` regenerate them
