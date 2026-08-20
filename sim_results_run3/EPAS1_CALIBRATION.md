# EPAS1 as a calibration target

EPAS1 is the best-characterised adaptive-introgression locus in East Asians, so
it is the natural yardstick for whether run2/run3 are simulating a detectable
scenario at all. This note collects the published parameters, compares them to
what run3 actually simulates, and draws out the two consequences.

## Published parameters

| Quantity | Value | Source |
|---|---|---|
| Archaic source | **Denisovan**, not Neanderthal | Huerta-Sánchez et al. 2014 |
| Core introgressed haplotype | **32.7 kb** (extended forms ~40 kb) | Huerta-Sánchez et al. 2014; Hackinger et al. 2016 |
| Frequency, Tibetans | **63–87%** | 63% Hackinger et al. 2016; ~87% Huerta-Sánchez et al. 2014 |
| Frequency, lowland East Asians | **<3%** (~9% Han in some reports) | Hackinger et al. 2016 |
| Admixture time | **~48,700 years ago** (16,000–59,500) | Zhang et al. 2021 PNAS |
| Selection onset | **~9,000 years ago** (2,500–42,000); 2,800–10,000 in an independent estimate | Zhang et al. 2021; Hackinger et al. 2016 |
| Scaled selection | **σ = 2·Ne·s ≈ 118** in Tibetans, ~1 in lowlanders | Hackinger et al. 2016 (SelEstim) |

The locus "remained selectively neutral for a long time in the population before
positive selection occurred", plausibly concurrent with permanent settlement of
the plateau after the Last Glacial Maximum.

## Implied selection coefficient

Taking the allele from ~1% to ~70%:

| Onset | Generations (29 y) | Implied s | Ne implied by σ = 118 |
|---|---|---|---|
| 2,800 ya | 97 | 0.113 | 523 |
| 6,000 ya | 207 | **0.053** | 1,121 |
| 10,000 ya | 345 | 0.032 | 1,869 |

So EPAS1 is a **recent, strong** sweep: `s ≈ 0.03–0.11`, most plausibly ~0.05.

## Tract-length consistency check

Expected surviving tract after introgression is `1/(r·t)`:

- EPAS1: admixture 48,700 ya = 1,679 generations → **59.5 kb** predicted, **32.7 kb**
  observed. Good agreement, which validates using `1/(r·t)` to reason about
  segment length.
- run3 CHB: archaic entry ~940 generations ago → **106 kb** predicted.

run3's archaic tracts are therefore already realistic, and slightly longer than
EPAS1's because the admixture is more recent in that model.

## How run3 compares

| | run3 as built | EPAS1 |
|---|---|---|
| Selection onset | 35,960 ya (1,240 gen) | ~6,000 ya (~207 gen) |
| Selection coefficient | 0.010 | ~0.05 |
| Final frequency | ~0.18 | 0.63–0.87 |
| Founding haplotypes at the locus | many (continuous migration) | one (single introgressed haplotype) |

**run3 simulates a much older, much weaker, much softer sweep than the real
locus it is meant to represent.** All four differences push in the same
direction: less TMRCA compression, less power.

## The consequence that matters

A hard sweep's carriers must coalesce by the sweep onset. So the onset sets a
ceiling on carrier TMRCA:

| Onset | Carrier TMRCA ceiling |
|---|---|
| 6,000 ya | 207 generations |
| 10,000 ya | 345 generations |
| 35,960 ya (run3) | 1,240 generations |

run2 established that the requested 1–50 ky cutoff grid sat *below* the
admixture floor and therefore saw almost nothing. **An EPAS1-calibrated onset
inverts that**: a sweep completing within the last 10,000 years compresses
carrier TMRCA to ≤345 generations, which lands squarely inside the original
1–50 ky grid — precisely the range where run2 found nothing because its sweep
was 36–57 ky old.

Calibrating the onset to EPAS1 is therefore the single highest-value change
available, and it makes the original threshold grid correct again.

## On introducing a focal *haplotype* rather than a focal allele

Genealogy cannot be pasted in. TMRCA is determined by actual descent, so marking
a 50 kb window with mutations on *k* randomly chosen genomes does not make those
genomes share an ancestor: the TMRCA at every base in that window remains the
population TMRCA. (Literally copying *sequence* instead would make Gamma-SMC
infer a recent TMRCA where the true genealogy is old — an artifact, and the
tree-truth analysis would correctly report no signal.)

The operative variable is the **number of founding lineages**, not the segment
length:

- 1 founder (de novo, or a single introgressed haplotype) → hard sweep, all
  carriers coalesce at the origin, maximal TMRCA compression.
- *k* founders → soft sweep, compression diluted roughly by `1/k_eff`, because a
  random carrier pair shares a founder only that often.

This is exactly why the EAS arm is underpowered: it places the allele on ~30–47
*independent* resident genomes, which share no recent ancestry, so amplifying
them compresses nothing.

Segment length controls the **spatial width** of the signal, not its depth — and
it is not a free parameter: it is `1/(r·t)` since introgression. Making the
introduced segment "genome-length" would implicitly assert introgression at
`t ≈ 0`, which is unphysical and would not deepen the signal anyway. EPAS1 is
detectable at 32.7 kb because it descends from **one** founding haplotype swept
**recently and hard**, not because the segment is long.
