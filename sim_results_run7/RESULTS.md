# run7 — selection on an introgressed allele in the grafted EAS model

EAS only, no CHB. The inferred PHLASH EAS history carries the whole human
lineage; a Neanderthal branch is grafted onto it and a single pulse delivers the
focal allele. Archaic parameters are the rounded literature values:

| quantity | value | generations (25 y/gen) | source |
|---|---|---:|---|
| human–Neanderthal split | 700 kya | 28,000 | Prüfer et al. 2014 (550–765 ky) |
| introgression pulse | 55 kya | 2,200 | Sankararaman 2012; Fu 2014 (47–65 ky) |
| archaic Ne | 3,600 | — | Ragsdale & Gravel 2019 |
| admixture | 2.5% | — | Prüfer 2014; Vernot & Akey 2015 |

The focal allele is fixed in the archaic branch just after the split, so what
arrives in EAS rides genuine archaic haplotypes, and it is **required to be at
2.0–3.0% in EAS one tick after the pulse**.

Three families of arms, all sharing one model and one neutral null, so that each
contrast isolates a single variable:

| family | origin | onset | s grid | replicates |
|---|---|---|---|---:|
| introgressed | fixed in archaic branch, 2.5% pulse | 55 kya | 0.002–0.01 | 80 each |
| de novo | single new copy in EAS | 55 kya | 0.001–0.05 | 20 each |
| recent | single new copy in EAS | **5 kya** | 0.05, 0.1 | 20 each |

All conditioned on survival. 620 selected + 80 neutral replicates, zero failures.

## 1. The null is calibrated

Neutral arm scored against its own null:

| cutoff | FPR at α = 0.05 | median p |
|---:|---:|---:|
| 1,000 | 0.025 | 1.000 |
| 4,500 | 0.013 | 1.000 |
| 10,000 | 0.050 | 0.514 |
| 20,000 | 0.050 | 0.661 |
| 30,000 | 0.062 | 0.632 |
| 50,000 | 0.062 | 0.541 |
| 100,000 | 0.050 | 0.516 |

Nominal from 10 ky up. Below ~4,500 y the statistic is zero in almost every
replicate, so p piles at 1.0 — degenerate rather than merely weak.

## 2. Establishment: 2.5% is only ~39 copies

The PHLASH EAS trajectory bottoms out at **Ne = 3,869 around 55 kya**, so a 2.5%
pulse delivers roughly 39 copies. With h = 0.5 the loss probability is about
exp(−2·k·h·s): **14% at s = 0.01, ~68% at s = 0.002**. An unconditioned pilot lost
the allele in 2/12 replicates at s = 0.01, matching. All arms are therefore
conditioned on survival, which is also what the empirical analysis does — it only
ever scans tracts that exist.

## 3. The mechanism: softness, not strength

Power does **not** rise with s. It peaks at intermediate values and falls away,
and the reason is measurable. Hom-carrier pairs at the focal base:

| introgressed | final AF | frac. coalescing after the pulse | median TMRCA |
|---|---:|---:|---:|
| s = 0.002 | 0.54 | 0.894 | 38.9 ky |
| s = 0.003 | 0.65 | **0.933** | 40.1 ky |
| s = 0.005 | 0.87 | 0.834 | 46.6 ky |
| s = 0.010 | 1.00 | 0.803 | **50.5 ky** |
| neutral | — | 0.225 | 289.8 ky |

As s rises, fewer carrier pairs trace to a shared founder and the median TMRCA
moves **older**, creeping toward the 55 kya pulse.

A weakly selected allele lingers at low frequency for many generations, and while
the carrier pool is small its lineages coalesce with each other, so present-day
carriers descend from few founding archaic haplotypes. Strong selection pulls the
allele out of the low-frequency regime quickly, that coalescence never happens,
and many of the ~39 founders survive — a **softer** sweep whose carriers meet only
back in the archaic branch.

So the cutoff must bracket where the mass sits. At s = 0.003 the median is 40 ky
and plenty falls below 30 ky; at s = 0.01 the median is 50.5 ky and almost nothing
does, but the 50 ky cutoff lands on the median. **Fixation relocates the signal to
older cutoffs rather than erasing it**: s = 0.01 has power 0.125 at 30 ky and
0.762 at 50 ky.

The de novo arms show the same gradient (median 38.3 ky at s = 0.005 → 53.0 ky at
s = 0.05), and peak at s = 0.003 (30 ky AUC 0.974, power 0.800).

## 4. A recent sweep is detectable, and at recent cutoffs

Every 55 kya arm finished sweeping long before sampling, which is why nothing
showed signal below ~10 ky. A sweep beginning 5 kya (200 generations) is still
compressing coalescence at the moment of sampling. Scored on a fine grid,
unconditional statistic, against the same neutral null:

| cutoff | recent 5 kya, s = 0.1 | recent 5 kya, s = 0.05 | de novo 55 kya, s = 0.05 | introgressed 55 kya, s = 0.003 |
|---:|---:|---:|---:|---:|
| 1,000 | 0.583 / 0.200 | 0.483 / 0.000 | 0.483 / 0.000 | 0.495 / 0.025 |
| 3,000 | 0.760 / 0.450 | 0.481 / 0.000 | 0.506 / 0.000 | 0.513 / 0.013 |
| **5,000** | **0.932 / 0.850** | 0.464 / 0.050 | 0.561 / 0.050 | 0.510 / 0.025 |
| 10,000 | 0.951 / 0.850 | 0.475 / 0.000 | 0.490 / 0.050 | 0.522 / 0.013 |
| 20,000 | 0.874 / 0.750 | 0.493 / 0.050 | 0.474 / 0.050 | 0.515 / 0.062 |
| 30,000 | 0.877 / 0.650 | 0.454 / 0.100 | 0.453 / 0.050 | 0.671 / 0.188 |
| 50,000 | 0.820 / 0.550 | 0.452 / 0.050 | 0.546 / 0.150 | 0.737 / 0.350 |

(AUC / power.) The mean statistic at 5,000 y is **0.197** for the recent s = 0.1
arm against 0.002–0.004 for every 55 kya arm — roughly a fifty-fold enrichment
where the neutral expectation is near zero.

On hom carriers the recent arm reaches **AUC 1.000 and power 1.000 at every
cutoff from 5,000 y upward**, and still manages 0.823 / 0.667 at 3,000 y.

This is the cleanest confirmation of the timing principle: the detectable window
tracks *when the sweep was happening*, and for a 5 kya sweep it sits at 5–10 ky,
exactly where the 55 kya arms are blind.

### s = 0.05 at 5 kya is not detectable, because the allele never rises

Mean final AF is **0.016**. That is arithmetic, not a failure of the statistic:
additive selection multiplies the odds by exp(s·t/2), which over 200 generations
at s = 0.05 is only ~150-fold, so a single copy in ~8,000 genomes reaches a few
percent even when it survives. Detecting an LCT-like *outcome* requires an
LCT-like *trajectory* — s ≈ 0.1 here, or an older origin, or standing variation.
Published LCT estimates span roughly 0.01–0.19, so s = 0.1 sits inside the
plausible band rather than being an inflated strawman.

## 5. Caveats

* **The hom-carrier panels are computed on a biased subset at low s.** A
  replicate needs ≥10 hom carriers to be scorable, so at introgressed s = 0.002
  only 35 of 80 replicates qualify, and those are the ones that drifted highest
  (mean AF 0.54 among them versus 0.315 across all 80). The de novo low-s arms
  are worse — 7, 13 and 15 of 20 at s = 0.001, 0.002, 0.003. Low-s hom-carrier
  power is therefore conditional on having a scorable carrier class and is
  optimistic. The unconditional statistic has no such problem: all replicates are
  always scorable, and there the introgressed peak sits at s = 0.005 rather than
  at the lowest s.
* **De novo and recent arms carry 20 replicates**, so power has a standard error
  near 0.1; the trend across s is reliable, individual entries are not.
* **Flat maps.** One recombination rate (1e-8) and one mutation rate (1.25e-8)
  throughout, so empirical windows must be matched on map length.
* **Q = 5 rescaling** with selection; stdpopsim warns this is not equivalent under
  selection, and it was not checked against Q = 1.
* **Cutoffs below ~4,500 y are degenerate**, not merely underpowered.

## Figures

* `figures/run7_final_af.png` — final AF by s, against the f = 2/3 line
* `figures/run7_ptmrca.png` — P(TMRCA < x) against the neutral 95% band
* `figures/run7_auc_power.png` — AUC and power against cutoff
* `figures/run7_auc_vs_af.png` — AUC against final allele frequency
* `figures/run7_pvalues.png` — p-value distributions with the neutral calibration
* `figures/run7_scan.png` — spatial extent along the contig. **Top row is the
  unconditional statistic, bottom row is hom carriers**; an earlier version left
  the top row unlabelled because the axis label was set outside the row loop.
* `figures/run7_origin_contrast.png` — de novo against introgressed

## Reproducing

```bash
python scripts/run_run7.py --slim-bin .native-stdpopsim/bin/slim --workers 6 --replicates 80
python scripts/run_run7_denovo.py --slim-bin .native-stdpopsim/bin/slim --workers 6 --replicates 20
python scripts/run_run7_recent.py --slim-bin .native-stdpopsim/bin/slim --workers 6 --replicates 20
python scripts/run7_figures.py
python scripts/run7_softness.py
python scripts/run7_recent_check.py --genotype-class unconditional
python scripts/run7_validate_graft.py --replicates 24
```
