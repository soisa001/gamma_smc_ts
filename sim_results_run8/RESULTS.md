# run8 — three selection scenarios, tree truth and Gamma-SMC decode

One model, four arms, one shared neutral null. The demography is run7's graft,
imported unchanged: the inferred PHLASH EAS history with a Neanderthal branch
splitting at 700 kya (28,000 gen), archaic Ne 3,600, and a 2.5% pulse at 55 kya
(2,200 gen).

| scenario | origin | onset | s | final AF |
|---|---|---|---:|---:|
| `recent_denovo` | de novo, 1 copy | 5 kya | 0.1 | 0.458 |
| `postintro_denovo` | de novo, 1 copy | 55 kya | 0.002 | 0.319 |
| `introgressed_pulse` | 2.5% archaic pulse | 55 kya | 0.002 | 0.364 |
| `introgressed_pulse_weak` | 2.5% archaic pulse | 55 kya | 0.001 | 0.219 |

`postintro_denovo` and `introgressed_pulse` share onset *and* s, so the only
thing separating them is where the allele came from.

**Design.** 100 diploids = 200 haplotypes; the decoder draws **10,000 random
pairs** at a fixed seed and tree truth is computed on exactly those pairs, read
back from the decoder's own manifest — so truth and decode are comparable
pair-for-pair. Simulation at µ = 1.29e-8, r = 1e-8. Decoding at θ = 0.00075,
ρ = 0.0006 with `--recent_call mean`, pinned to the Gamma-SMC paper, so
`frac_recent_<T>` is the *proportion of posterior means below T*. Selected arms
condition on arrival at 2.0–3.0% and on survival; the neutral arm carries no
allele and no conditioning. 100 replicates per arm plus 100 neutral, all seven
thresholds for every arm. **500/500 completed, zero failures.**

## 1. Unconditional: only the recent sweep is detectable

Power, at most one of 100 neutral replicates matching or exceeding:

| cutoff | recent_denovo | postintro_denovo | introgressed_pulse | introgressed_pulse_weak |
|---:|---:|---:|---:|---:|
| 1,000 | 0.010 / **0.850** | 0.010 / 0.000 | 0.000 / 0.010 | 0.010 / 0.000 |
| 5,000 | **0.940 / 0.940** | 0.000 / 0.020 | 0.000 / 0.030 | 0.000 / 0.020 |
| 10,000 | **0.940** / 0.850 | 0.010 / 0.000 | 0.010 / 0.020 | 0.070 / 0.000 |
| 20,000 | 0.860 / 0.610 | 0.010 / 0.010 | 0.050 / 0.010 | 0.040 / 0.000 |
| 30,000 | 0.730 / 0.350 | 0.060 / 0.010 | 0.050 / 0.020 | 0.030 / 0.000 |
| 50,000 | 0.440 / 0.120 | 0.050 / 0.000 | 0.080 / 0.030 | 0.010 / 0.000 |

(truth / decode.) The three 55 kya arms are **not detectable without
conditioning** — power never exceeds 0.09 at any cutoff, and AUC stays in
0.49–0.63.

The mean statistic explains why, and the explanation is not about allele
frequency:

| cutoff | recent_denovo | introgressed_pulse | neutral null |
|---:|---:|---:|---:|
| 5,000 | **0.2303** | 0.0026 | 0.0027 |
| 30,000 | 0.3259 | 0.1147 | 0.1018 |
| 50,000 | 0.3992 | 0.2481 | 0.2134 |

At 5 kya the recent sweep sits **85× above** a near-zero baseline. At 50 kya the
introgressed arm sits 1.16× above a baseline of 0.21 — far too small against
between-replicate scatter. **Detection is driven by contrast against the neutral
baseline at that cutoff, not by carrier frequency.** `recent_denovo` succeeds at
AF 0.458 while `introgressed_pulse` fails at AF 0.364 because the baseline at
5 kya is 0.003 and at 50 kya is 0.21.

This also sharpens run5's f > 2/3 crossover: with random cross-individual pairs,
only ~AF² of pairs are carrier–carrier, so at AF 0.36 roughly 87% of the 10,000
pairs carry no signal and dilute what the carriers contribute.

Figure: `figures/run8_check1_scan.png`.

## 2. Conditioning on carrier pairs is transformative

AUC by carrier class, tree truth:

| cutoff | carrier–carrier | carrier–noncarrier | noncarrier–noncarrier | overall |
|---:|---:|---:|---:|---:|
| **introgressed_pulse** | | | | |
| 10,000 | **0.976** | 0.000 | 0.848 | 0.531 |
| 30,000 | **0.983** | 0.000 | 0.749 | 0.621 |
| 50,000 | **0.992** | 0.000 | 0.665 | 0.620 |
| **postintro_denovo** | | | | |
| 10,000 | **0.994** | 0.000 | 0.841 | 0.522 |
| 30,000 | **0.993** | 0.000 | 0.706 | 0.551 |
| **introgressed_pulse_weak** | | | | |
| 30,000 | **0.998** | 0.000 | 0.712 | 0.538 |

Every scenario that is invisible unconditionally becomes near-perfectly separable
once pairs are stratified — carrier–carrier AUC is 0.98–1.00 at every cutoff from
10 kya up, including the weakest arm.

**Carrier–noncarrier AUC is 0.000**, which is the hard floor predicted from
theory rather than a small number: a pair with one carrier and one non-carrier
cannot coalesce more recently than the allele's origin, so those pairs are
systematically *older* than neutral and rank below the entire null.

**This is an upper bound, not a method.** The null here is unconditioned, and
run6 measured that conditioning on carrier status inflates the statistic 1.2–3×
under strict neutrality. A calibrated version needs a carrier-matched null, and
in the empirical pipeline hmmix enters only as a post-hoc overlap, so this
conditioning is not currently runnable on data at all.

Figure: `figures/run8_check2_classes.png`.

## 3. Decode against truth

Decode tracks truth, with a recency shift that **cancels** because both sides of
the comparison are decoded identically:

* **Below 5 kya the decoder wins outright** — power 0.850 vs 0.010 at 1 kya. Tree
  truth at 1 kya is 0.0006 against a null of 0.0005: almost no pair genuinely
  coalesces inside 1,000 years, so the point statistic is discrete and empty.
  The decoder's posterior is continuous and borrows information from flanking
  sequence, giving 0.0032 against a null of 0.0006 — a real, calibrated 5× excess.
* **At 5 kya they tie** (0.940 both).
* **Above 10 kya truth wins** — 0.730 vs 0.350 at 30 kya, 0.440 vs 0.120 at
  50 kya. The θ = 0.00075 prior implies a constant Ne of ~14,500, which compresses
  deep times and costs discrimination at old cutoffs.

This corrects how I described the same effect in runs 3–5. It is not that
Gamma-SMC "beats the genealogy"; it is that the decoder's prior shifts inferred
times young, which helps where the truth statistic is degenerate and hurts where
it is not. Because the null is decoded under the identical prior, the shift never
threatens validity.

`mean_p_lt` (mean posterior probability) is **uniformly worse** than
`frac_recent` — 0.77/0.69/0.59 against 0.85/0.94/0.85 at 1/5/10 kya — so the
paper's statistic is the right primary.

p-values for the AF-median replicate of each arm confirm the picture: only
`recent_denovo` reaches significance (p = 0.0099, the floor with 100 neutral
replicates), at every cutoff from 5 kya to 30 kya on truth and 1 kya to 20 kya on
decode. No 55 kya arm's median replicate is significant at any cutoff.

Figures: `figures/run8_check3_power.png`, `figures/run8_check3_pvalues.png`.

## 4. A prediction of mine that was wrong

I expected `introgressed_pulse_weak` (s = 0.001) to end near 7% and to have too
few carriers to stratify. It reached **AF 0.219** with ~719 carrier–carrier pairs
per replicate, and stratifies fine. The deterministic logistic I predicted from
ignores that these arms are **conditioned on survival**, which keeps only the
upward-biased tail of trajectories. At weak selection that conditioning does most
of the work, and my closed-form prediction had no term for it.

## 5. Caveats

* **The stratified numbers use an unconditioned null** and are therefore
  optimistic by roughly the 1.2–3× factor run6 measured. They bound headroom.
* **The p-value floor is 1/101 ≈ 0.0099** with 100 neutral replicates, which is
  exactly where `recent_denovo` sits. Finer resolution needs more neutral
  replicates.
* **Decode cannot be stratified** from the summary output, which aggregates over
  the whole pair set. A per-class decode needs three runs against three pair
  files and a re-simulation, since trees are deleted after decoding.
* **Flat maps**: one recombination rate and one mutation rate throughout.
* **Q = 5 rescaling** with selection, unchecked against Q = 1.
* Runs 3–8 used `--recent_call median`; their `frac_recent` values are **not**
  comparable to run8's, because a right-skewed Gamma has median below mean.

## Reproducing

```bash
python scripts/run_run8.py --slim-bin .native-stdpopsim/bin/slim --decoder bin/gamma_smc --workers 6
python scripts/run8_analyse.py
python scripts/run8_figures.py
python scripts/run8_theory_check.py
```

Full parameter provenance, including the implied constant Ne and the significance
rule, is serialised to `config/run8_scenarios.json` (schema v2).
