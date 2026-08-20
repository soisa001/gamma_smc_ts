# run3 — results

400/400 replicates simulated and decoded, zero failures.

## Headline: Gamma-SMC beats tree truth at recent cutoffs

This **contradicts the prediction recorded in `sim_results_run2/RESULTS.md`**,
which argued on three grounds that an inferred TMRCA should lose power to the
true genealogy. It does not, at the recent end.

EAS arm, unconditional detection power (`overall` class, 100 selected vs 100
neutral replicates):

| x (years) | tree truth | gamma_smc_frac | gamma_smc_meanp |
|---:|---:|---:|---:|
| 1,000 | 0.00 | 0.00 | **0.51** |
| 4,500 | 0.00 | 0.10 | **0.54** |
| 10,000 | 0.03 | **0.51** | 0.55 |
| 20,000 | 0.32 | **0.60** | 0.54 |
| 30,000 | 0.61 | 0.60 | 0.54 |
| 40,000 | **0.67** | 0.51 | 0.53 |
| 50,000 | 0.65 | 0.48 | 0.52 |
| 100,000 | 0.46 | 0.38 | 0.48 |
| 200,000 | 0.18 | 0.07 | 0.41 |

At 10,000 years the decoder finds signal in **51%** of replicates where tree
truth finds it in **3%**. Truth only overtakes at 30 ky and above.

**Why.** The decoder systematically over-calls recent coalescence, and it does
so more in selected replicates than neutral ones, which manufactures separation
at cutoffs where the truth has none:

| x (years) | truth mean | decoded mean | bias | Pearson r |
|---:|---:|---:|---:|---:|
| 10,000 | 0.008 | 0.139 | **+0.132** | 0.04 |
| 20,000 | 0.077 | 0.366 | **+0.289** | 0.45 |
| 40,000 | 0.538 | 0.555 | +0.018 | 0.96 |
| 50,000 | 0.564 | 0.590 | +0.026 | 0.96 |

The bias is large and the correlation is near zero exactly where the decoder
wins, and the bias vanishes and the correlation reaches 0.96 exactly where truth
wins. So the recent-cutoff gain is **not** the decoder recovering the genealogy
better; it is a systematic shift that happens to be informative. It is real
power — the null calibration below confirms it is not over-firing — but it is
power from a biased estimator, not from better inference.

`gamma_smc_meanp` (the posterior mean probability rather than the point call) is
flat at ~0.5 across every cutoff including 1,000 years, where no pair
genuinely coalesces. That flatness is the signature of a statistic that is
detecting "this replicate looks unusual" rather than anything threshold-specific,
and it should not be read as detection at 1,000 years.

**Null calibration is sound** — leave-one-out false-positive rates against the
null itself, averaged over thresholds:

| arm | tree_truth | gamma_smc_frac | gamma_smc_meanp |
|---|---|---|---|
| CHB | 0.045 | 0.041 | 0.049 |
| EAS | 0.045 | 0.037 | 0.050 |

All at or below the nominal 0.05, so none of the above is inflated type I error.

## The CHB arm is weak, as the design predicted

| | selected | neutral |
|---|---|---|
| allele present at focal base | 0.42 | 0.18 |
| census AF, mean | 0.156 | 0.006 |

Only 42% of selected replicates have any archaic ancestry at the focal base, so
well over half carry no signal at all before selection is even considered.

Power never exceeds 0.20 at any cutoff for any estimator. The reason was
identified before the run: archaic ancestry arrives by continuous low-rate
migration and selection only has 1,240 generations from a frequency that
accumulates from zero, so the allele reaches ~10% rather than sweeping. There is
very little for a within-individual statistic to see.

Carrier-stratified truth power is **0.00 everywhere** for `any_carrier`,
`het` and `hom_carrier`. With the allele at ~10% and 100 sampled diploids, those
classes contain a handful of pairs at most, so the statistic is too noisy to
test. Only `overall` and `hom_noncarrier` are usable.

## A methodological artifact to record: conditional power is not zero, it is undefined

The conditional analysis (restricting both arms to replicates carrying the
allele) reports power of exactly 0.00 for every estimator at every threshold in
both arms. That is **an artifact, not a result**. Conditioning shrinks the null
to 18 replicates (CHB) and 17 (EAS), which puts the p-value floor at
`1/(1+18) = 0.053` and `1/(1+17) = 0.056` — both above the nominal
`alpha = 0.05`, so no replicate can be significant regardless of its statistic.

A conditional test needs at least 19 null replicates to be able to reach
`p <= 0.05` at all. run4 sidesteps this by ascertaining at simulation time, so
every selected replicate carries the allele and the null carries none.

## What this means for run4

1. The recent-cutoff advantage of the decoder is worth carrying forward, and
   run4 keeps all three estimators.
2. run3 confirms the diagnosis that motivated the EPAS1 recalibration: an old,
   weak, soft sweep is close to undetectable, in either estimator, at any cutoff.
3. Conditional-null starvation must be avoided by design, not by post-hoc
   filtering.

## Provenance

- `logs/simulation_status.tsv`, `logs/decode_status.tsv` — 400 rows each, all completed
- `validation_report.json` — schedule and script gates
- Raw per-pair TMRCA matrices and tree sequences are kept in the WSL working
  copy under `<arm>/<mode>/replicates/`, gitignored; seeds in `config/`
  regenerate them
