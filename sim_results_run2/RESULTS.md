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
