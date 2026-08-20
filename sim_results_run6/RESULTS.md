# run6 — the genotype-stratified TMRCA scan, and what the null has to be

run5 found that the within-individual TMRCA at a swept introgressed locus is
**trimodal**: hom-carrier pairs are far younger than neutral, heterozygous pairs
are far older and can never be recent, and hom-non-carrier pairs are
neutral-like. run6 asks the operational follow-up — *at each 10 kb stride, does
the carrier class look younger than a neutrally evolving region?* — and then
asks the question that actually decides whether the answer means anything:
**what is that "neutrally evolving region" allowed to be?**

Two arms are scanned:

| arm | source | model | onset | s | final AF |
|---|---|---|---|---|---|
| `chb_from_introgression` | run5, reused | `OutOfAfricaArchaicAdmixture_5R19` | 1,240 gen (36.0 ky) | 0.01 | mean 0.705, median 0.762 |
| `eas_standing` | **new, run6** | PHLASH EAS median, **full history** | 1,440 gen (36.0 ky) | 0.01 | mean 0.917, median 0.960 |

The CHB side needed no new simulation — run5 retained the raw per-pair TMRCA
matrices and focal genotypes, which is exactly what the scan consumes. The EAS
arm is new and is run in **SLiM rather than msprime**, for two reasons: msprime's
structured-coalescent sweep never marks which lineages carry the beneficial
allele, so it cannot be stratified by genotype at all; and msprime forbids
demographic events inside a sweep, which forced run5 to flatten the PHLASH
history across the sweep window. SLiM has neither limitation, so the
ten-thousand-epoch history is carried exactly. 200/200 replicates completed.

## 1. The scan reproduces across both demographies

AUC of the selected arm against a matched-size neutral null, at the focal base:

**EAS arm** (new; 100 selected + 100 neutral, 25 y/gen)

| cutoff | hom carrier | het | hom non-carrier | all individuals |
|---:|---:|---:|---:|---:|
| 10,000 | 0.548 | 0.422 | 0.590 | 0.528 |
| 20,000 | 0.834 | 0.199 | 0.530 | 0.803 |
| 30,000 | **0.985** | 0.082 | 0.463 | 0.972 |
| 50,000 | **0.996** | 0.187 | 0.647 | 0.987 |
| 100,000 | 0.957 | 0.317 | 0.490 | 0.936 |

**CHB arm** (run5 data, 29 y/gen) peaks at the same place — 0.988 at 30 ky —
with the same het collapse (0.064 at 30 ky, 0.002 at 50 ky). The pattern is a
property of the statistic, not of either demography.

Detection width (hom carrier, AUC >= 0.75): EAS 470 kb at 30 ky, 380 kb at 50 ky,
130 kb at 20 ky, and **zero at <= 10 ky**. CHB: 490 kb at 30 ky. Nothing below
10 ky in either arm.

Figure: `figures/run6_eas_scan.png`.

## 2. The hom-carrier statistic has an ascertainment, and the null must share it

The hom-carrier class is not a random set of pairs. It is selected for *carrying
the same thing on both haplotypes*, and two haplotypes carrying the same derived
allele must coalesce at or after that allele's origin. They are therefore younger
than two random haplotypes **under strict neutrality, with no selection anywhere
in the model**. Scoring an ascertained observation against an unconditional null
credits that ascertainment as though it were selection.

The size of the effect was measured directly, by simulating neutrally, finding a
segregating site of a given frequency near the focal base, and comparing the
homozygous-derived pairs' statistic against the unconditional one at the same
position (`scripts/run6_conditioned_null.py`, 250 replicates per demography).

Null inflation — conditioned / unconditional, at P(TMRCA < 30,000 y):

| class frequency | EAS history | 5R19 history |
|---:|---:|---:|
| 0.25–0.35 | 2.91x | 3.29x |
| 0.45–0.55 | 2.15x | 2.26x |
| 0.65–0.75 | 1.51x | 1.52x |
| 0.85–0.95 | 1.17x | 1.18x |

Two things stand out. The bias is **large** — up to 3x — and it is **almost
completely demography-independent**: the two histories agree to within a few
percent at every frequency. The inflation is a property of the ascertainment,
not of the model, and it decays with frequency exactly as it should, because a
common derived allele is old and constrains the genealogy less.

## 3. Correction: conditioning on a *tract* is stronger than conditioning on a SNP

The observation is not homozygous-*alt*, it is homozygous-*introgressed* — an
hmmix tract call, not a single site. The natural expectation is that a
single-SNP null is the more conservative bar. **It is not**, and the direction
matters.

`scripts/run6_introgressed_null.py` builds the correct comparator: 1,200 neutral
5R19 replicates with migrations recorded, archaic ancestry called per haplotype,
and every stride where at least five homozygous-introgressed individuals exist
scored and pooled by local introgression frequency (13,705 locus x threshold
rows; the drift tail is genuinely populated — 307 loci reach 50–70% and 125
exceed 70% introgression with no selection at all).

Inflation at P(TMRCA < 30,000 y), **both under 5R19, so demography is controlled**:

| class frequency | hom-derived SNP | hom-introgressed tract |
|---:|---:|---:|
| ~0.25–0.35 | 3.29x | 4.75x |
| ~0.45–0.55 | 2.26x | 2.84x |
| ~0.65–0.75 | 1.52x | 1.96x |
| ~0.85–1.00 | 1.18x | 1.39x |

The tract null is **higher at every matched frequency**. Sharing 50 kb of archaic
sequence on both haplotypes is a far stronger statement about recent common
ancestry than sharing one derived allele, so hom-tract pairs coalesce more
recently than hom-SNP pairs do. A hom-alt null therefore sets a bar that is
roughly 1.2–1.3x **too low**, and is mildly *anti*-conservative as a proxy — it
will over-call selection, not under-call it.

Figure: `figures/run6_conditioning_inflation.png`.

## 4. Above the introgression time the statistic is empty

Scoring the same 92 CHB replicates against all four candidate nulls
(`scripts/run6_null_comparison.py`, each replicate matched to the null bin of its
own final frequency):

| cutoff | observed | CHB uncond. | EAS uncond. | SNP-hom | tract-hom | power (CHB / EAS / SNP / tract) |
|---:|---:|---:|---:|---:|---:|---|
| 10,000 | 0.010 | 0.005 | 0.008 | 0.011 | 0.005 | 0.04 / 0.02 / 0.05 / 0.15 |
| 20,000 | 0.156 | 0.022 | 0.046 | 0.066 | 0.031 | 0.67 / 0.44 / 0.38 / 0.58 |
| 30,000 | 0.655 | 0.077 | 0.100 | 0.152 | 0.129 | 0.94 / 0.88 / 0.79 / **0.82** |
| 50,000 | 0.922 | 0.312 | 0.204 | 0.326 | 0.847 | 0.91 / 0.99 / 0.84 / **0.00** |

The last row is the whole argument for conditioning. At 50 ky the unconditional
nulls report power of 0.91–0.99; the correctly ascertained null reports **0.00**.
The introgression pulse sits at ~36 ky, so above it a neutrally introgressed
homozygote has *already* coalesced within the window (null mean 0.847) — there is
nothing left for selection to add. Every apparent 50 ky signal is the
ascertainment, not the sweep.

**The usable window is below the introgression time.** At 30 ky the signal
survives every null: observed 0.655 against a correctly conditioned null whose
mean is 0.129 and whose 99th percentile is 0.215.

Figure: `figures/run6_null_comparison.png`.

## 5. Using the EAS model as the null for a population that has introgression

The EAS history was inferred on the target samples themselves, so it is the right
null *model*; it simply contains no archaic admixture. Substituting it for CHB's
own demography costs less than expected — at 30 ky the null mean moves 0.077 ->
0.100 (1.3x) and power 0.94 -> 0.88 — but the error is not uniform in the cutoff.
At 50 ky the EAS null sits *below* CHB's own (0.204 vs 0.312) and power rises
spuriously to 0.99. A model mismatch that flatters the test at one cutoff and
penalises it at another is not a safe substitution; the cutoff has to be chosen
before the null, on the introgression time, not tuned afterwards.

The deeper problem is that a null model without introgression **cannot produce a
homozygous-introgressed class at all**, so the hom-carrier ascertainment has no
matched counterpart in it. Given that constraint, the defensible route is to drop
the conditioning from *both* sides (`scripts/run6_recommended_design.py`):

Power at alpha = 0.05, empirical statistic vs simulated neutral EAS null:

| cutoff | CHB arm, hom-carrier | CHB arm, unconditional | EAS arm, hom-carrier | EAS arm, unconditional |
|---:|---:|---:|---:|---:|
| 20,000 | 0.44 | 0.16 | 0.41 | 0.33 |
| 30,000 | 0.88 | **0.73** | 0.91 | **0.85** |
| 50,000 | 0.99 | 0.76 | 0.98 | 0.96 |

Dropping the conditioning costs roughly 6–15 points of power at 30 ky and buys
an unbiased test: with no ascertainment on the observation, the null owes no
conditioning, and the unconditional EAS simulation is exactly the right
comparator. This works only because a sweep drags the *whole* window's
distribution once carrier frequency exceeds 2/3 — which is the regime both arms
occupy (CHB median AF 0.76, EAS 0.96) and the regime run5's f > 2/3 crossover
already identified.

## 6. Caveats

* **Local rate heterogeneity.** Every simulation here uses a flat 1e-8
  recombination rate and 1.25e-8 mutation rate. Real windows vary in both, and
  P(TMRCA < x) is sensitive to the local map. Comparing empirical windows to a
  simulated null requires either simulating each window under its own map length
  or restricting to windows whose map length matches the simulation.
* **The 5R19 archaic-ancestry call** uses recorded migrations into the archaic
  deme, which is the model's own ground truth, not an hmmix reconstruction. Real
  tract calls carry error that will blur the genotype classes.
* **Marginal, not genome-wide.** The pooled null is the *marginal* distribution;
  positions within a replicate are linked, so it does not license genome-wide
  multiple-testing claims.
* **Cutoffs below ~10 ky carry no signal** in either arm, at any conditioning.

## Reproducing

```bash
python scripts/run_run6_eas.py --slim-bin .native-stdpopsim/bin/slim --workers 8
python scripts/run6_conditioned_null.py --demography eas --replicates 250
python scripts/run6_conditioned_null.py --demography chb --replicates 250
python scripts/run6_introgressed_null.py --replicates 1200 --workers 6
python scripts/run6_null_comparison.py
python scripts/run6_recommended_design.py
python scripts/run6_figures.py
```
