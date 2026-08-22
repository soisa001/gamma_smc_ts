# run9 — stronger selection, and a delayed EPAS1-like onset

Nine arms on one model, crossing two origins with four selection regimes. The
demography, samples, pairs, rates and decoding are run8's, imported unchanged, so
every difference below is a difference in the selection model alone.

| arm | origin | onset | s | final AF | carrier–carrier pairs |
|---|---|---|---:|---:|---:|
| `pulse_s0p003` | 2.5% pulse | 55 kya | 0.003 | 0.515 | 3,204 |
| `pulse_s0p004` | 2.5% pulse | 55 kya | 0.004 | 0.675 | 5,218 |
| `pulse_s0p005` | 2.5% pulse | 55 kya | 0.005 | 0.803 | 6,819 |
| `pulse_delayed_s0p005` | 2.5% pulse | **10 kya** | 0.005 | 0.280 | 1,066 |
| `denovo_s0p003` | de novo | 55 kya | 0.003 | 0.462 | 2,839 |
| `denovo_s0p004` | de novo | 55 kya | 0.004 | 0.613 | 4,378 |
| `denovo_s0p005` | de novo | 55 kya | 0.005 | 0.736 | 5,897 |
| `denovo_delayed_s0p01` | de novo | **10 kya** | 0.01 | 0.019 | 7 |
| `denovo_delayed_s0p05` | de novo | **10 kya** | 0.05 | 0.789 | 6,579 |

100 replicates per arm plus **300 neutral** (up from run8's 100, dropping the
p-value floor from 1/101 = 0.0099 to 1/301 = 0.0033). 1,200/1,200 completed,
zero failures. Significance: at most one neutral replicate matching or exceeding.

## 1. Raising s works, and the transition is steep

Unconditional power, tree truth, at 50 kya:

| origin | s = 0.003 | s = 0.004 | s = 0.005 |
|---|---:|---:|---:|
| pulse | 0.250 | 0.410 | **0.640** |
| de novo | 0.250 | 0.320 | 0.580 |

run8's s = 0.002 managed 0.05–0.09 at every cutoff. `pulse_s0p005` reaches
AUC **0.946**. So the s grid answers the question directly: an introgressed sweep
becomes detectable without any carrier conditioning somewhere between s = 0.002
and s = 0.005.

**The optimal cutoff is 50 kya, not 20.** Truth AUC for `pulse_s0p005` climbs
monotonically — 0.716 (20 k), 0.818 (30 k), 0.907 (40 k), **0.946 (50 k)** — then
falls to 0.909 at 75 kya. My earlier correction that it was not 20 kya was right;
my replacement guess of 30–40 kya was still too shallow.

## 2. The delayed EPAS1-like onset failed, and instructively

`pulse_delayed_s0p005` never exceeds power 0.09 at any cutoff. Its AUC peaks at
0.696 (20 kya) and then **inverts below 0.5** — 0.533 (40 k), 0.482 (50 k), 0.424
(75 k). At deep cutoffs the locus reads *older* than neutral.

The allele arrives on the pulse, drifts neutrally for 1,800 generations, and only
then experiences 400 generations of selection, reaching AF 0.280. That drift
leaves the carriers as ordinary deep-coalescing archaic haplotypes, and 400
generations of selection cannot compress enough to overcome it. Meanwhile the
locus is enriched for archaic ancestry, which pushes coalescence *deeper* — hence
the inversion.

**A recent onset does not rescue an introgressed sweep; it removes the thing that
made it detectable.** Detectability comes from sustained sweep-driven coalescence
since the pulse, not from recency per se.

## 3. De novo at a recent onset is a different story

`denovo_delayed_s0p05` reaches **AUC 1.000 and power 1.000 at 10 kya**, 0.89 at
5 kya. `denovo_delayed_s0p01` is dead (AF 0.019, 7 carrier–carrier pairs).

The bracket was necessary: at the same 10 kya onset, s = 0.01 goes nowhere while
s = 0.05 sweeps to AF 0.789. Over 400 generations additive selection multiplies
the odds by exp(s·t/2), which is only ≈ 7 at s = 0.01.

Comparing origins at matched onset and coefficient is where the pulse's head
start shows: at 10 kya the pulse arm reaches AF 0.280 from s = 0.005 while de novo
reaches 0.019 from s = 0.01 — twice the coefficient, one fifteenth the frequency.
Starting from 2.5% instead of 1/(2N) is worth more than doubling s.

## 4. The portable result: detection turns on near AF 0.6

Pooling all nine arms and binning on final frequency removes the dependence on
any particular s — selection acts on detection almost entirely *through* the
frequency it produces, since only about AF² of random pairs are carrier–carrier.

Power, tree truth:

| AF bin | 10 kya | 30 kya | 50 kya | n |
|---|---:|---:|---:|---:|
| 0.00–0.05 | 0.000 | 0.000 | 0.009 | 111 |
| 0.05–0.15 | 0.039 | 0.000 | 0.000 | 51 |
| 0.15–0.30 | 0.048 | 0.012 | 0.000 | 83 |
| 0.30–0.45 | 0.032 | 0.032 | 0.011 | 93 |
| 0.45–0.60 | 0.110 | 0.193 | 0.055 | 109 |
| 0.60–0.75 | 0.114 | 0.258 | **0.561** | 132 |
| 0.75–0.90 | 0.184 | 0.313 | **0.718** | 163 |
| 0.90–1.00 | 0.266 | 0.348 | **0.848** | 158 |

Below AF 0.45 there is essentially nothing at any cutoff. Above 0.6 the 50 kya
cutoff works, and improves steadily to 0.85 at near-fixation.

Decode needs a slightly higher frequency (~0.75) but has a **shallower optimum**:
at 30 kya in the top bin it reaches 0.652 against truth's 0.348 — nearly double —
while at 50 kya truth leads 0.848 to 0.741. That is the θ = 0.00075 prior
(implied constant Ne ≈ 14,500) shifting inferred times young, exactly as run8
found, and it means **truth and decode want different cutoffs**: 50 kya for truth,
30 kya for decode.

## 5. The floor prediction, confirmed by where it does *not* hold

Carrier–noncarrier AUC at the 30 kya cutoff:

| arm | onset | carrier–noncarrier AUC |
|---|---|---:|
| all six 55 kya arms | 55 kya | **0.000** |
| `denovo_delayed_s0p01` | 10 kya | 0.416 |
| `denovo_delayed_s0p05` | 10 kya | 0.381 |

A pair with one carrier and one non-carrier cannot coalesce more recently than
the allele's origin. For a 55 kya allele that floor sits above the 30 kya cutoff,
so the statistic is *exactly* zero and those pairs rank below the entire null. For
a 10 kya allele the floor sits below the cutoff, so the constraint does not bite
and the AUC is unremarkable. The prediction is confirmed by both the zero and its
absence.

Carrier–carrier AUC is 0.92–1.00 for every arm, including the ones that are
undetectable unconditionally — but scored against an unconditioned null, so it
bounds headroom rather than measuring calibrated power.

## Figures

* `figures/run9_power.png` — AUC and power by cutoff, all nine arms
* `figures/run9_af_threshold.png` — the AF threshold, pooled across arms
* `figures/run9_quantiles.png` — p-value at each percentile of the statistic
* `figures/run9_classes.png` — carrier-pair stratification by arm
* `figures/run9_scan.png` — spatial profile at each arm's best cutoff

## Caveats

* Stratified classes use an **unconditioned null** and are optimistic by roughly
  the 1.2–3× run6 measured for carrier ascertainment.
* Flat maps: one recombination and one mutation rate throughout.
* Q = 5 rescaling with selection, unchecked against Q = 1.
* The AF-stratified table pools arms with different onsets and origins; within a
  bin those are not interchangeable, so it estimates a marginal relationship
  rather than a causal one.

## Reproducing

```bash
python scripts/run_run9.py --slim-bin .native-stdpopsim/bin/slim --decoder bin/gamma_smc --workers 20
python scripts/run9_analyse.py
python scripts/run9_figures.py
```
