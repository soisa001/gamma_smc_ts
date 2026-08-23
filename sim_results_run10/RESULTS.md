# run10 — detection at realistic introgression frequency, against an honest null

Real introgressed regions reach ~30–40% at most. run9's AF-stratified table put
unconditional power at **0.03** in that band. run10 asks whether anything rescues
that, and answers three questions with no new selected simulations — only a new
*null*.

## 1. A bigger panel cannot help

The per-window statistic varies for two reasons: sampling noise, which shrinks
with more haplotypes, and genealogical variance, which does not. Measured on
run9's 300 neutral replicates:

| cutoff | replicate variance ÷ binomial | 50 → 200 haplotypes |
|---|---:|---|
| 10,000 y | — | variance ÷ 10.4 |
| 30,000 y | 75× | variance ÷ 1.65 |
| 50,000 y | 229× | variance ÷ 1.07 |

At 30–50 kya the null's width is almost entirely genealogical. One region has one
genealogy; sampling it harder does not change how much genealogies differ between
regions. **Going from 200 haplotypes to the full 4,000-haplotype panel would buy
almost nothing at the cutoffs that matter.** Only at 10 kya, where the statistic
is small and counts are sparse, does sampling noise dominate.

## 2. Regional averaging cannot help either

Trees decorrelate along the sequence, so averaging over a wide window averages
over quasi-independent genealogies — and the null SD does fall, from 0.061 at
10 kb to 0.016 at 4 Mb (50 kya). But the sweep signal is only ~500 kb wide, so
wider windows dilute signal faster than they shrink noise:

| window | overall AUC, 30 kya, AF 0.25–0.45 |
|---|---:|
| 10 kb | **0.688** |
| 500 kb | 0.617 |
| 4 Mb | 0.541 |

The narrowest window is best. Region merging is a calling convenience, not a
power gain.

## 3. The carrier result was mostly detecting introgression

That left carrier conditioning, which looked decisive — AUC 1.000, power 0.99 at
AF 0.25–0.45. But it was scored against a neutral arm containing **no carriers**,
which cannot separate "this locus was selected" from "this locus is introgressed".

Two archaic haplotypes are already unusual: they either share a post-pulse
ancestor or meet back in the archaic branch. So run10 builds the matched null —
600 neutral replicates of the same grafted model with a census at the pulse,
archaic ancestry called per haplotype, and every stride whose ancestry drifted
into the observed range scored identically. **11,279 loci**, 2,042 of them in
0.25–0.35 and 982 in 0.35–0.45.

Under strict neutrality, the carrier–carrier statistic is already far above the
carrier-free baseline:

| archaic frequency | cutoff | carrier-free | **matched null** | inflation |
|---|---|---:|---:|---:|
| 0.25–0.35 | 10,000 y | 0.0075 | 0.0269 | 3.6× |
| 0.25–0.35 | 30,000 y | 0.100 | **0.411** | 4.1× |
| 0.25–0.35 | 50,000 y | 0.194 | **0.845** | 4.4× |
| 0.35–0.45 | 30,000 y | 0.105 | 0.344 | 3.3× |
| 0.45–0.60 | 30,000 y | 0.108 | 0.270 | 2.5× |

At 50 kya the matched null saturates (q95 = 1.000 in every band), so that cutoff
carries no information about selection at all once introgression is controlled
for.

## 4. What survives

Selected pulse arms from run9, scored against the frequency-matched null on the
identical statistic:

| AF band | AUC 10 ky | AUC 30 ky | AUC 50 ky | power p≤0.05 (30 ky) | power p≤0.01 (30 ky) |
|---|---:|---:|---:|---:|---:|
| 0.15–0.25 | 0.824 | 0.713 | 0.556 | — | 0.000 |
| 0.25–0.35 | 0.736 | 0.773 | 0.616 | 0.238 | 0.000 |
| 0.35–0.45 | 0.731 | **0.819** | 0.661 | 0.294 | 0.000 |
| 0.45–0.60 | 0.737 | 0.793 | 0.663 | 0.250 | 0.038 |
| 0.60–0.71 | 0.652 | 0.743 | 0.611 | 0.231 | 0.077 |

**The signal is real but weak.** AUC ≈ 0.77–0.82 at 30 kya means a selected
region ranks above roughly four fifths of neutrally drifted introgressed regions
at the same frequency — genuinely informative for *ranking* candidates. But
per-locus significance is another matter: **power at p ≤ 0.05 is 0.24–0.29, and
at p ≤ 0.01 it is essentially zero.**

For comparison, the same replicates scored against the carrier-free null gave
AUC 1.000 and power 0.99. Almost all of that was the ascertainment.

## What this means

* **The unconditional statistic does not work at realistic introgression
  frequencies**, and neither more samples nor wider windows change that. Its
  operating range is AF > 0.6, which does occur but is not what you described.
* **Carrier conditioning gives a real signal at 30–40%**, but it is a ranking
  statistic (AUC ≈ 0.8), not a per-locus test. In a genome-wide scan with
  multiple-testing correction, power at p ≤ 0.01 of ~0 will not survive.
* **30 kya is the right cutoff.** 50 kya is destroyed by null saturation once
  introgression is matched; 10 kya is competitive at low frequency.
* Every carrier-conditioned number in runs 6–9 is inflated by 2.5–4.4× from
  ascertainment alone and should be restated against this null.

The honest framing for the empirical analysis is a **ranked candidate list with
calibrated AUC**, not a significance scan — unless the statistic is combined with
independent evidence (frequency outliers, hmmix tract length, functional
annotation) to buy the missing power.

## Caveats

* Selected replicates are matched on final frequency, which within any one arm
  selects atypical replicates — a `pulse_s0p005` replicate landing at AF 0.30 is
  an unlucky sweep. That is the correct conditioning for "given observed
  frequency, is this selected?", but per-band sample sizes are small (34–52).
* The null is msprime and the selected arms are SLiM at Q = 5. The unconditional
  neutral baselines agree (0.097–0.121 versus run9's 0.101), so the engines are
  consistent, but this was not otherwise validated.
* De novo arms are excluded: their carrier class is defined by a modern allele,
  so an archaic-ancestry-matched null is the wrong comparator for them.
* Ancestry here is the model's own census truth, not an hmmix reconstruction;
  real tract-calling error will blur the carrier classes further.

## Reproducing

```bash
python scripts/run10_variance_diagnostic.py
python scripts/run10_regional.py
python scripts/run10_matched_null.py --replicates 600 --workers 20
python scripts/run10_headtohead.py
python scripts/run10_figures.py
```
