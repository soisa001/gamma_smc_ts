# run2 — framework sketch

**Status: complete. 400/400 replicates run, zero failures. See [RESULTS.md](RESULTS.md).**
This document is the pre-registration for `sim_results_run2/`; nothing in this
directory yet contains results. See [`README.md`](README.md) for how to run it,
and `validation_report.json` for the gates that already pass.

Sections 2, 5 and 8 were revised once the model parameters were read out of
`AncientEurasia_9K19` rather than estimated; the revisions are marked. §5 in
particular is now backed by 20,000 Wright-Fisher replicates per arm rather than
by closed-form approximations, and it corrects the claim that the CHB sweep
starts at 2.96%.

Question: **can we detect selection at an introgressed locus from the
within-individual TMRCA distribution alone?**

The design deliberately restarts from basics. It fixes one selection
coefficient, one starting frequency, one selection onset, one statistic, and one
null. Everything the existing runs sweep over (nine target allele frequencies,
three selection coefficients, three demographic quantiles, frequency-matched
nulls, de novo vs standing origins, Gamma-SMC decode variants) is out of scope
for this pass.

---

## 1. What is already in the repo, and what run2 changes

| Existing | run2 |
|---|---|
| `introgression_EAS_sim/` — AncientEurasia_9K19 + SLiM, 27 specs (3 s-values × 9 target AFs), each conditioned on a present-day Han AF band | 1 spec: s = 0.01, no AF conditioning |
| Allele drawn as a **single copy** in Neanderthal at 2400 gen, conditioned to reach ≥ 10% in the source before the pulse | Allele is **standing variation at exactly 2.96%** in the recipient at the pulse |
| `no_introgression_EAS_sim/` — PHLASH EAS, **de novo** origin, origin age solved per (s, curve) to hit a 45–55% present-day band | Same EAS demography, but standing variation at 2.96% at the same time point as the CHB pulse |
| `ancient_eurasia_survival_neutral` — null **conditioned on the focal allele still segregating** in the present Han census | Null is plain neutral: no focal allele is drawn, the site may be monomorphic or absent |
| Statistic reported across genotype classes (hom-carrier / het / hom-noncarrier / overall) and across truth + Gamma-SMC | Single statistic: **overall** within-individual `P(TMRCA < x)`, tree truth |

The reason the old runs are hard to read is that the allele-frequency
conditioning is doing three jobs at once: it sets the sweep's starting point, it
sets its endpoint, and it defines the null. run2 separates them. The starting
point is fixed by construction (2.96% at the pulse), the endpoint is a *measured
outcome*, and the null is unconditioned. Note that for the CHB arm the pulse is
in the Han *ancestor*, so the frequency Han itself starts from is measured too —
see §5.

**Inherited inconsistency to fix.** `introgression_EAS_sim/specifications.tsv`
records `realized_han_split_generations = 2010`, while
`ancient_eurasia_survival_neutral.py` records
`REALIZED_HAN_SPLIT_GENERATIONS_Q5 = 2015.0`. Under the stated `floor(t/Q)·Q`
rule with Q = 5, `floor(2016/5)·5 = 2015`. run2 uses 2015 and records nominal
and realized times side by side for every event (gate G1).

---

## 2. Two simulation arms

Both arms share every parameter in the shared table below, except the demography
and the presence of an archaic source.

### Arm A — CHB / Han (`chb_ancient_eurasia/`)

- Demography: `stdpopsim` `HomSap/AncientEurasia_9K19` (Kamm et al. 2019),
  unmodified, `stdpopsim == 0.3.0`.
- Introgression: the catalog's own 2.96% Neanderthal → Loschbour pulse at
  2272 generations ago (realized 2270 at Q = 5).
- Target population: `Han`, split from Loschbour at 2016 gen (realized 2015),
  present-day Ne = 6300.
- **Standing variation is created by fiat, not by rejection.** At the pulse tick
  the focal allele is placed on *every* introgressing genome. This is exactly
  equivalent to the allele being fixed in the archaic source at the moment of
  the pulse, and it buys three things at once:
  - post-pulse recipient allele frequency = **exactly 0.0296** (in the Han
    ancestor; the frequency entering Han itself is measured, §5);
  - carriers are **exactly the archaic haplotypes**, so the divergent
    introgressed background — the thing that makes this question interesting —
    is preserved;
  - **zero rejection sampling**. Reaching the same state by drawing a single
    copy in Neanderthal and conditioning on fixation before the pulse is not
    feasible: fixation probability from one copy is ~1/(2·N_arch), so the accept
    rate would be O(10⁻³) or worse.

### Arm B — EAS custom demography (`eas_phlash/`)

- Demography: single population, piecewise-constant, from the tracked PHLASH
  artifact `no_introgression_EAS_sim/resources/EAS.npz`
  (sha256 `a5bdfd84…db928`), **pointwise median** curve, constant-extrapolated
  from 100 generations to the present.
- No archaic population and no pulse. At **the same time point** (2272 gen
  nominal / 2270 realized) the focal allele is placed on a uniformly random
  2.96% of genomes.
- Consequence, stated plainly: **Arm B is not an introgression model.** It
  matches Arm A on placement frequency, onset time, and selection coefficient,
  and differs in that carriers sit on ordinary EAS haplotypes rather than
  divergent archaic ones. That is the point of the contrast — it isolates how
  much of the detectable signal comes from the sweep itself versus from the
  archaic haplotype background.
- *(Added after implementation.)* The two arms also differ in **initial genotype
  structure**, unavoidably. A pulse migrant in a WF model has both parents in
  the source, so every Arm A carrier starts as a complete archaic homozygote:
  ~1.5% homozygotes, no heterozygotes. Arm B's random-genome placement gives
  Hardy–Weinberg proportions, so almost every carrier is a heterozygote. This is
  not a modelling slip — it is what an introgression pulse actually does — but
  it means the A-vs-B contrast bundles "archaic haplotype background" together
  with "carriers start homozygous", and the two cannot be separated within this
  design.
- Relevant Ne along the trajectory (median curve, diploids):

  | gen ago | 100 | 500 | 1000 | 2016 | 2272 | 5000 | 10000 | 40000 |
  |---|---|---|---|---|---|---|---|---|
  | Ne | 37,637 | 6,222 | 3,327 | 3,645 | 3,947 | 8,243 | 12,121 | 13,138 |

  Note the ~10× expansion inside the last 500 generations; Arm B is the more
  expensive of the two to simulate forward.

### Shared parameters

| Parameter | Value |
|---|---|
| Sequence length | 10 Mb |
| Focal position | 5,000,000 (0-based), reserved from the neutral overlay |
| Mutation rate | 1.25e-8 /bp/gen |
| Recombination rate | 1.0e-8 /bp/gen |
| Generation time | 25 years |
| Selection coefficient | s = 0.01 |
| Dominance | h = 0.5 (additive) |
| Selection window | onset one tick after the pulse — realized **2265** generations ago, the tick at which the introgressed genomes first exist — continuously to the present (Arm A: in Loschbour 2265 → 2015, then in Han 2015 → 0; Arm B: one population throughout) |
| Frequency at placement | 0.0296 — in the Han ancestor for Arm A, in EAS for Arm B. Arm A's frequency *entering Han* is a measured outcome (§5), not 0.0296 |
| Present-day AF conditioning | **none** |
| Replicates | 100 selected + 100 neutral, per arm |
| Sample | 100 diploids drawn uniformly from the present-day target population |
| Engine | SLiM 4.2.2 via `stdpopsim`, Q = 5, burn-in 0.1, recapitation + msprime neutral overlay |

The two-stage pool(500) → panel(100) sampling used by the earlier runs existed
only to support AF matching. run2 drops it and samples 100 diploids directly.

---

## 3. The statistic

For one replicate, at the focal position:

> **T(x) = fraction of the 100 sampled diploids whose two haplotypes coalesce
> less than x years ago**, i.e. TMRCA < x/25 generations, strict inequality.

x ∈ **{1,000; 4,500; 10,000; 20,000; 30,000; 40,000; 50,000}** years.

Two constraints force this exact definition:

1. **Overall, not carrier-stratified.** The null contains no focal allele, so
   there is no genotype class to condition on. Any carrier-stratified statistic
   is uncomparable to the null and cannot produce a p-value. Carrier-stratified
   profiles are still worth plotting for the selected arm — as description, not
   as the test.
2. **Within-individual pairs**, matching
   `eas_sweep_analysis.tree_truth_profiles` and the Gamma-SMC pair convention,
   so a later decode pass is directly comparable.

Granularity is 1/100 = 0.01 per replicate, set by the sample size.

**Pass 1 is tree truth only — no Gamma-SMC decode.** The decoder adds a build
dependency, a posterior-calibration question, and a second source of
disagreement, none of which bear on "is there signal at all". If truth shows a
spike, pass 2 asks whether Gamma-SMC recovers it.

---

## 4. Null and p-value

- 100 neutral replicates per arm: identical demography, identical engine,
  identical Q, identical overlay, **no `DrawMutation`, no fitness event, no
  conditioning of any kind**. The focal base is whatever neutral evolution left
  there, including nothing.
- Same statistic, same position, same sample size.
- One-sided upper tail, per threshold:

  ```
  p_i(x) = (1 + #{ j : T_null,j(x) ≥ T_sel,i(x) }) / (1 + N_null)
  ```

- **Resolution floor: p ≥ 1/101 ≈ 0.0099.** With 100 null replicates no
  replicate can be assigned a smaller p-value. If a genome-wide-scale threshold
  is ever wanted, the null must grow (1,000+ replicates) or its tail must be
  fitted parametrically. Flagged now so the floor is not later mistaken for a
  result.
- Seven thresholds are seven correlated tests. Report per-threshold p-values;
  do not silently take the minimum.
- Headline number is **power**: the fraction of the 100 selected replicates with
  p ≤ 0.05, at each threshold.

### Why the null is run in SLiM and not msprime

Under neutrality an exact coalescent simulation would be faster and, in
principle, distributionally identical. But the selected arm is Q = 5 rescaled
SLiM plus recapitation plus overlay, and a p-value must not quietly absorb an
engine difference. The primary null therefore uses the identical pipeline minus
selection. Gate G3 runs 20 msprime replicates as a cross-check; if they agree,
later expansions of the null can use the cheap engine.

---

## 5. Expected behaviour (pre-registered, so surprises are visible)

The numbers below come from 20,000 Wright-Fisher replicates of each arm's
demography, run before any SLiM simulation. They replace the first version's
Kimura hand-estimates, which were wrong in both arms.

### Where the sweep actually starts, and the Loschbour question

In `AncientEurasia_9K19` the Han split is a proportion-1 mass migration from Han
into Loschbour, so **forwards in time Han is founded out of Loschbour**: between
the pulse (2270) and the split (2015), "Loschbour" *is* the Han ancestor, not a
side branch. Selection there is selection on the Han lineage.

The consequence is that **0.0296 is the frequency at the pulse, in the ancestor
— not the frequency the Han sweep starts from.** Over the intervening 250
generations at Ne = 2,340 the allele grows and drifts, and enters Han at:

| | mean | quartiles | 5-95% | lost before Han exists |
|---|---|---|---|---|
| AF entering Han | 0.092 | 0.021 / 0.072 / 0.142 | [0.000, 0.261] | 14.3% |

This is measured per replicate by a SLiM event registered at the split tick and
reported as `target_entry_af` in `results/final_af.tsv`. **Do not describe the
CHB arm as starting at 2.96% in Han.**

Two alternatives were considered and rejected:

- *Drop the Loschbour fitness callback so selection is Han-only.* This is worse,
  not better: the allele then drifts neutrally at Ne 2,340 for 250 generations
  and **33% of replicates are lost before Han is founded**, 42% overall.
- *Retime the pulse into Han so the AF is pinned at 0.0296.* Clean (2.4% loss),
  but it moves archaic admixture from 56.8 kya to 50.2 kya, leaves the archaic
  tracts unrecombined at sweep onset, and departs from the published model.

The design keeps the published demography and reports the entry frequency
instead. Revisit only if the entry-AF spread turns out to dominate the results.

### Final outcomes

Deterministic additive growth from p = 0.0296 at s = 0.01 needs roughly
`(2/s)·Δlogit ≈ 1,600` generations to reach 99%, and there are 2,265, so
surviving replicates should mostly be at or near fixation.

| Arm | lost | fixed | median AF among survivors |
|---|---|---|---|
| CHB (Ne 2,340 → 6,300) | 18.4% | 71.7% | 1.00 |
| EAS (Ne 3,947 → 37,637) | 11.2% | 65.6% | 1.00 |

So roughly one CHB replicate in five and one EAS replicate in nine carries no
signal at all. **Replicates that lose the allele are kept**, not resampled —
they are the honest left tail of the final-AF distribution and they contribute
genuinely null-like TMRCA profiles, which is exactly why the headline number is
power rather than a single p-value.

### A scaling artefact worth knowing

Q = 5 rescales Loschbour to 468 diploids, so the number of introgressing
individuals is a draw from `Binomial(468, 0.0296)` — about 14, with a standard
deviation of 3.6. The placement frequency is therefore 0.0296 in expectation
with roughly 26% relative spread, about √5 wider than the unscaled model would
give. Gate G2 records the realized value rather than asserting an exact 0.0296.

---

## 6. Deliverables

### Tables

| File | Contents |
|---|---|
| `results/final_af.tsv` | one row per selected replicate: final AF in census and in the 100-diploid sample; lost / segregating / fixed |
| `results/tmrca_profiles.tsv.gz` | long form: replicate × arm × {selected, neutral} × position × threshold → `p_tmrca_lt_threshold` |
| `results/observed_statistics.tsv` | T(x) at the focal position, one row per replicate |
| `results/null_distribution.tsv` | the 100 neutral T(x) values per threshold, with quantiles |
| `results/pvalues.tsv` | per selected replicate per threshold: T, p, and the power summary |
| `results/null_calibration.tsv` | leave-one-out false-positive rate of the null against itself (gate G6) |
| `results/replicate_status.tsv` | one row per attempted replicate: seed, status, wall time, error |
| `validation_report.json` | realized event times, software versions, gate outcomes (written by the `validate` phase) |

Per-replicate intermediates live under `<arm>/<mode>/replicates/` and are
gitignored: they are regenerable from the seeds in `config/`, which are a pure
function of `(arm, mode, index)`. Tree sequences are not written unless
`--save-trees` is passed.

### Figures (PNG + PDF, per arm)

1. **`final_af.png`** — histogram of final allele frequency across the 100
   selected replicates, with the lost fraction called out as its own bar.
2. **`ptmrca_selected.png`** — T(x) against x for the 100 selected replicates:
   per-replicate lines faint, median and IQR ribbon heavy. This is the "visual
   spike".
3. **`ptmrca_neutral.png`** — the same panel for the 100 neutral replicates, on
   identical axes.
4. **`pvalue_panels.png`** — one small panel per threshold: null histogram,
   observed selected values overlaid, p-value annotated.
5. **`spatial_profile.png`** *(cheap bonus)* — T(x) at a fixed x along the whole
   10 Mb on the 10 kb stride grid, selected vs neutral. The focal spike should
   be localised at 5 Mb; a flat profile would be a red flag that something other
   than the sweep is driving the statistic. This costs nothing extra because the
   profiles are already computed per position.
6. **`ptmrca_by_genotype.png`** *(descriptive, selected arm only)* — T(x) split
   by carrier genotype class. Not part of the test.

### Cross-arm

7. **`arm_comparison.png`** — Arm A vs Arm B power against threshold: the direct
   answer to "how much of the signal is the archaic background?"

---

## 7. Directory layout

```
sim_results_run2/
  FRAMEWORK.md                  <- this file
  config/
    run2_chb.json               frozen parameters, Arm A
    run2_eas.json               frozen parameters, Arm B
  chb_ancient_eurasia/
    selected/replicates/rep000..rep099/{simulation.trees, truth_profiles.tsv.gz, endpoint.json}
    neutral/replicates/rep000..rep099/{...}
    results/
    figures/
  eas_phlash/
    (identical structure)
  cross_arm/
    results/
    figures/
  logs/
```

One replicate directory is self-describing: its `endpoint.json` carries the
seed, the realized event times, the final AF, and the hash of the parameters it
was run under. There is no separate ledger / adapter / recovery layer.

---

## 8. Validation gates

Each gate must pass before the next stage runs. G1, G2 (static half), G4 and G5
are already implemented and passing in `validation_report.json`, produced by
`python scripts/run_run2.py validate` — which needs stdpopsim but not the SLiM
binary, and so can gate the design before anything is launched.

- **G0 — environment and timing.** SLiM 4.2.2 and `stdpopsim == 0.3.0` present
  and version-checked. Run 3 replicates per arm and extrapolate wall time before
  committing to 400 runs. *Needs SLiM.*
- **G1 — realized times. PASSING.** The Q = 5 grid puts the pulse at 2270
  generations ago, the standing-variation placement one tick later at 2265, and
  the Han split at 2015, with the two fitness callbacks meeting exactly on the
  split tick so selection has no gap. Nominal and realized values are recorded
  side by side for every event, resolving the 2010-vs-2015 inconsistency in §1.
- **G2 — starting frequency.** Static half **PASSING**: the generated script is
  asserted to select carriers by `inds.migrant` (Arm A) or to place 0.0296
  directly (Arm B), and the placement tick is asserted to be pulse + 1. Dynamic
  half needs SLiM: the patched `add_mut` records the realized post-pulse
  frequency into tree-sequence metadata and raises if no introgressing
  individuals are present, so a wrong `migrant` assumption fails loudly on the
  first replicate rather than silently. Because Q = 5 leaves only ~14 migrants
  (§5), this gate records the realized frequency rather than asserting an exact
  0.0296.
- **G2b — no rejection sampling. PASSING.** `condition_on_allele_frequency` is
  asserted empty in all four generated scripts. The accept rate is 1 by
  construction, which is the single biggest simplification over the earlier runs.
- **G3 — engine agreement.** 20 msprime coalescent neutral replicates per arm;
  compare the null T(x) distribution to the SLiM null (two-sample KS).
  Disagreement means the SLiM null stays primary and the discrepancy gets
  investigated.
- **G4 — focal base reserved. PASSING.** The guard masks `[5000000, 5000001)` in
  every intercepted overlay rate map, is asserted to preserve the integrated
  rate over the other 9,999,999 bases exactly, and is idempotent. stdpopsim
  already excludes a referenced single-site DFE from the *background* overlay,
  but its recapitation DFE spans the whole contig and would otherwise drop
  neutral mutations on the reserved base. Each replicate additionally fails if
  the guard never fired.
- **G5 — sampling. PASSING (static).** All 400 seeds are asserted unique. At
  runtime, every sampled individual must have exactly 2 sample nodes and the
  sample nodes must partition cleanly into the 100 diploids.
- **G6 — null calibration.** Leave-one-out p-values of the null against itself
  must sit near the nominal 0.05. This is the only check that distinguishes "the
  test works" from "the test always fires", and it costs nothing:
  `results/null_calibration.tsv`.

---

## 9. Open decisions

Defaults below are what will be built unless overridden.

| # | Decision | Default | Why it matters |
|---|---|---|---|
| D1 | EAS arm has no archaic source | plain standing variation | Literal reading of the brief. The alternative — attach a ghost archaic population to the PHLASH model with a matched 2.96% pulse and a 27,840-gen divergence — would make the two arms comparable *as introgression models*. Natural pass-2 extension. |
| D2 | Tree truth only, no Gamma-SMC decode | truth only | Removes the decoder build and the calibration question from pass 1. |
| D3 | No conditioning on allele survival | keep lost replicates | Makes the final-AF distribution and the power estimate honest. |
| D4 | Statistic is overall, not carrier-stratified | overall | Forced by the unconditioned null. |
| D5 | Null runs in SLiM, msprime as cross-check | SLiM primary | Keeps the p-value free of engine artefacts. |
| D6 | 100 diploids sampled | 100 | Sets statistic granularity at 0.01. 200 would smooth the curves at modest cost. |
| D7 | 100 null replicates | 100 | Sets the p-value floor at 0.0099. |

---

## 10. Execution environment

Neither SLiM nor `stdpopsim`/`msprime` is installed on the current Windows
workstation, and `bin/gamma_smc` and `.native-stdpopsim/bin/slim` are absent.
The repo's bootstrap scripts (`scripts/bootstrap_eas_sweep_study.sh`) target
Linux. run2 therefore needs a Linux target — WSL, a container built from the
tracked `Dockerfile`, or the AoU workbench — before G0 can run. 24 cores are
available locally.

Total workload is 400 SLiM runs (2 arms × {selected, neutral} × 100). With no
rejection sampling anywhere in the design, wall time is dominated by the forward
simulation of the recent expansion in Arm B. G0 measures it rather than
guessing.

---

## 11. Sequence of work

1. ~~Freeze `config/run2_chb.json` and `config/run2_eas.json`.~~ **Done.**
2. ~~Implement the pulse-tick allele placement and assert the schedule.~~
   **Done** — `validate` passes; the generated scripts are committed under
   `<arm>/validation/` for review.
3. Stand up the Linux environment; run G0 and the dynamic half of G2 on 3
   replicates per arm (`--index 0 --index 1 --index 2`).
4. 100 selected replicates per arm → `final_af.tsv`, figure 1.
5. 100 neutral replicates per arm; G3.
6. Profiles, p-values, figures 2–5.
7. Cross-arm comparison, figure 7.
8. Only then consider pass 2: Gamma-SMC decode, the D1 ghost-archaic EAS
   variant, and a larger null.
