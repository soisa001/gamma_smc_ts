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
arrives in EAS rides genuine archaic haplotypes rather than a random modern
background, and it is **required to be at 2.0–3.0% in EAS one tick after the
pulse**. Selection acts in EAS from the pulse to the present.

A **de novo** counterpart holds everything fixed — same model, same onset, same
shared neutral null — and changes only the origin: one new copy on an ordinary
modern background. Both families are conditioned on survival.

400 introgressed replicates (4 coefficients × 80), 80 de novo (4 × 20), 80
neutral. Zero failures.

## 1. The null is calibrated

The neutral arm scored against its own null, which is the check that makes
everything below interpretable:

| cutoff | FPR at α = 0.05 | median p |
|---:|---:|---:|
| 1,000 | 0.025 | 1.000 |
| 4,500 | 0.013 | 1.000 |
| 10,000 | 0.050 | 0.514 |
| 20,000 | 0.050 | 0.661 |
| 30,000 | 0.062 | 0.632 |
| 50,000 | 0.062 | 0.541 |
| 100,000 | 0.050 | 0.516 |

Nominal at every cutoff from 10 ky up. Below 4,500 y the statistic is zero in
almost every replicate, so p piles at 1.0 and the test is degenerate rather than
merely weak — the same discreteness the earlier runs hit.

## 2. Establishment: 2.5% is only ~39 copies

The PHLASH EAS trajectory bottoms out at **Ne = 3,869 around 55 kya**, so a 2.5%
pulse delivers roughly 39 copies, not a comfortable standing pool. With h = 0.5
the loss probability is about exp(−2·k·h·s): **14% at s = 0.01, rising to ~68% at
s = 0.002**. An unconditioned pilot lost the allele in 2/12 replicates at
s = 0.01, matching. Every arm below is therefore conditioned on survival, which
is also what the empirical analysis does — it only ever scans tracts that exist.

Final frequencies, conditional on establishment:

| arm | s = 0.002 | 0.003 | 0.005 | 0.01 |
|---|---:|---:|---:|---:|
| introgressed, final AF | 0.315 | 0.517 | 0.853 | 0.999 |
| hom-carrier diploids | 15.4 | 33.9 | 75.8 | 99.7 |

| de novo | s = 0.005 | 0.01 | 0.02 | 0.05 |
|---|---:|---:|---:|---:|
| final AF | 0.731 | 0.990 | 1.000 | 1.000 |
| hom-carrier diploids | 59.2 | 98.0 | 100.0 | 100.0 |

## 3. Power peaks at intermediate selection and collapses on fixation

This is the main result, and it runs opposite to the naive expectation that
stronger selection is easier to detect.

**Introgressed**, hom-carrier statistic:

| cutoff | s = 0.002 | 0.003 | 0.005 | 0.01 |
|---:|---:|---:|---:|---:|
| 10,000 | 0.664 / 0.000 | 0.604 / 0.068 | 0.497 / 0.013 | 0.445 / 0.037 |
| 20,000 | 0.847 / 0.457 | 0.772 / 0.339 | 0.703 / 0.269 | 0.484 / 0.013 |
| 30,000 | 0.918 / 0.657 | **0.932 / 0.746** | 0.849 / 0.500 | 0.590 / 0.125 |
| 50,000 | 0.995 / 0.971 | 0.994 / 0.932 | 0.977 / 0.910 | 0.952 / 0.762 |

(AUC / power at α = 0.05.)

**De novo**, hom-carrier statistic:

| cutoff | s = 0.005 | 0.01 | 0.02 | 0.05 |
|---:|---:|---:|---:|---:|
| 10,000 | 0.612 / 0.167 | 0.514 / 0.000 | 0.475 / 0.000 | 0.490 / 0.050 |
| 20,000 | 0.822 / 0.278 | 0.649 / 0.200 | 0.398 / 0.100 | 0.474 / 0.050 |
| 30,000 | **0.950 / 0.889** | 0.781 / 0.350 | 0.450 / 0.000 | 0.453 / 0.050 |
| 50,000 | 0.995 / 0.944 | 0.987 / 0.950 | 0.820 / 0.550 | 0.546 / 0.150 |

The observed statistic against the null at 30 ky makes the mechanism plain:

| arm | observed | null | ratio |
|---|---:|---:|---:|
| introgressed s = 0.002 | 0.302 | 0.100 | 3.00 |
| introgressed s = 0.003 | 0.337 | 0.100 | **3.35** |
| introgressed s = 0.005 | 0.246 | 0.100 | 2.45 |
| introgressed s = 0.01 | 0.115 | 0.100 | 1.15 |
| de novo s = 0.005 | 0.350 | 0.100 | **3.49** |
| de novo s = 0.01 | 0.216 | 0.100 | 2.15 |
| de novo s = 0.02 | 0.092 | 0.100 | 0.91 |
| de novo s = 0.05 | 0.094 | 0.100 | 0.94 |

Both origins peak and then die, and both die once the allele fixes — s = 0.01
for introgression (AF 0.999), s = 0.02 for de novo (AF 1.000). At s = 0.02 and
0.05 the de novo ratio drops slightly *below* 1: a completed sweep leaves a
star-like genealogy rooted at the sweep, so those loci look marginally **older**
than neutral at 30 ky, not younger.

**The reason is timing, not strength.** The statistic sees the epoch during
which coalescence was being compressed. A strong sweep finishes quickly —
at s = 0.05 fixation takes roughly (2/s)·ln(2N) ≈ 385 generations, so it was over
by ~45 kya and 45 ky of ordinary drift has since rebuilt local diversity. A
weaker sweep is still in progress at sampling, so the compression is recent and
visible. Fixation does not erase the signal so much as **relocate it to older
cutoffs**: the introgressed s = 0.01 arm has power 0.125 at 30 ky but 0.762 at
50 ky, because a fixed archaic haplotype's coalescence piles up just below the
55 kya pulse.

The practical consequence is that the cutoff has to bracket the sweep epoch, and
a scan run at a single cutoff is blind to sweeps of the wrong age.

## 4. A prediction of mine that did not hold

I expected a de novo sweep to show signal at cutoffs far more recent than
adaptive introgression could reach, on the grounds that introgression is floored
by the pulse while a hard sweep is not. **It does not.** At 4,500 y and 10,000 y
the best de novo AUC is 0.570 and 0.612, essentially the same as introgression's
0.551 and 0.664. The floor that matters is not the pulse but the *sweep
completion time*, and with an onset 2,200 generations back both origins finish
long before the present. A de novo sweep would need a much more recent onset to
produce recent coalescence.

## 5. The unconditional statistic

The unconditional statistic — every individual, no ascertainment, which is what
run6 established the empirical analysis should use — is weaker but tracks the
same shape:

| cutoff | introgressed s = 0.002 | 0.003 | 0.005 | 0.01 |
|---:|---:|---:|---:|---:|
| 30,000 | 0.582 / 0.037 | 0.671 / 0.188 | 0.781 / 0.362 | 0.588 / 0.125 |
| 50,000 | 0.627 / 0.150 | 0.737 / 0.350 | 0.928 / 0.762 | 0.952 / 0.762 |

Note the peak sits at a *higher* s than for the hom-carrier statistic (0.005
rather than 0.003), which is the f > 2/3 crossover from run5 showing up again:
the unconditional statistic needs the carrier class to be a majority before it
can drag the whole sample's distribution, whereas the hom-carrier statistic is
strongest while carriers are still a distinct minority.

## Figures

* `figures/run7_final_af.png` — final AF by s, against the f = 2/3 line
* `figures/run7_ptmrca.png` — P(TMRCA < x) against the neutral 95% band
* `figures/run7_auc_power.png` — AUC and power against cutoff
* `figures/run7_auc_vs_af.png` — AUC against final allele frequency
* `figures/run7_pvalues.png` — p-value distributions with the neutral calibration
* `figures/run7_scan.png` — spatial extent along the 10 Mb contig
* `figures/run7_origin_contrast.png` — de novo against introgressed

## Caveats

* **De novo arms carry 20 replicates**, so their power estimates have a standard
  error near 0.1 and the individual entries should not be over-read; the trend
  across s is the reliable part.
* **Flat maps.** A single 1e-8 recombination and 1.25e-8 mutation rate
  throughout, so empirical windows must be matched on map length before being
  compared with this null.
* **Q = 5 rescaling** with selection; stdpopsim warns that rescaling is not
  equivalent in the presence of selection, and this was not checked against
  Q = 1.
* **Cutoffs below ~4,500 y are degenerate**, not merely underpowered.

## Reproducing

```bash
python scripts/run_run7.py --slim-bin .native-stdpopsim/bin/slim --workers 6 --replicates 80
python scripts/run_run7_denovo.py --slim-bin .native-stdpopsim/bin/slim --workers 6 --replicates 20
python scripts/run7_figures.py
python scripts/run7_validate_graft.py --replicates 24
```
