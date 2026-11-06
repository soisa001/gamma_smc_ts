# ASMC versus Gamma-SMC: small EAS inference demo

ASMC improves temporal inference in this demo. Incorporating EAS demography gives
an additional improvement, but does not recover the narrow true coalescence peak.
The improvement from changing algorithms is larger than the additional improvement
from changing the demographic model within ASMC.

## Design

- Reuse the first 10 selected and first 10 held-out neutral regions in the frozen
  `selection_episode_40_20k_s002_20261106/cohort_manifest.tsv`, in its saved order.
  No regions were chosen by inference performance, and no new simulations were run.
- Selection is additive, s=0.02, active 40–20 kya; origin is a 1% Neanderthal pulse
  at 50 kya. Selected replicates are conditioned on present-day whole-population
  segregation. Neutral regions are unconditioned, with a nearest-segregating-site label.
- Each region has 200 phased haplotypes across 1 Mb. Every estimator uses the same
  saved 1,000 pairs (original pair seed 1729), evaluated at 100 coordinates with
  10-kb stride. The focal comparison is at 500 kb, not at each neutral label.
- Reuse exact true TMRCA, Gamma-SMC posterior-median, and Gamma-SMC posterior-mass
  summaries from verified existing outputs. Independently verify focal truth against
  the tree sequence and haplotype order.
- ASMC 1.4.0, PrepareDecoding 2.2.5. Sequence mode with CSFS disabled: use pairwise
  sequence information, not additional sample allele-frequency information.
  Both ASMC arms have the same inputs and time discretization.
- Mutation rate 1.29e-8/base/generation; recombination rate 1e-8/base/generation;
  generation time 25 years. The constant ASMC model uses diploid Ne=14,534.88372,
  matching the Gamma-SMC reference scaling theta=0.00075.
- ASMC receives *haploid* population sizes. The EAS input approximates the same
  PHLASH pointwise-median diploid history used in simulation. Harmonic interval
  sizes preserve cumulative coalescent hazard at each new boundary, including every
  1-kya boundary through 50 kya. Maximum absolute neutral CDF approximation error
  is 0.0000209 (0.00209 percentage points).
- Time states are 0–1, 1–2, ..., 49–50 kya, followed by logarithmically spaced older
  states to 5 Myr and an unbounded last state. Probabilities are not renormalized
  conditional on coalescence before 50 kya. The one-kyr peak comparison uses bins
  1–2 through 49–50 kya, matching the earlier scan; the five-kyr plot includes 0–5.
- Both ASMC models remain single-population inference models: the EAS arm includes
  population-size history, not an explicit Neanderthal split/pulse or selection.
- Four processes by default. No new RNG draws. Exact source identifiers, file hashes,
  numerical parameters, software versions, and timing receipts are saved alongside
  the summary cache. The separate raw work directory is `/home/mew/asmc-demo-work`.

## Results

Values below are means across the 10 selected regions, at the fixed focal coordinate.
CDF MAE is absolute error against each region's true CDF, averaged over the 50
cutoffs from 1–50 kya and then over regions; units are percentage points (pp).

| Estimator | Below 20 kya | Between 20–40 kya | Below 50 kya | CDF MAE |
|---|---:|---:|---:|---:|
| True TMRCA | 3.88% | 49.67% | 95.21% | 0 |
| Gamma-SMC median | 41.60% | 35.59% | 85.10% | 26.57 pp |
| Gamma-SMC mass | 44.15% | 23.50% | 74.57% | 29.67 pp |
| ASMC constant size | 22.34% | 41.73% | 74.38% | 18.19 pp |
| ASMC EAS history | 15.34% | 41.47% | 69.24% | 15.07 pp |

The selected-minus-neutral mean one-kyr excess peaks at 41.5 kya for truth,
16.5 kya for Gamma medians, 3.5 kya for Gamma mass, 22.5 kya for ASMC constant,
and 24.5 kya for ASMC EAS. These are descriptive maxima, not estimated selection
onset dates. The true peak differs from the earlier 100-replicate figure because
this demo uses only 10 selected replicates. Coalescence peaks need not lie exactly
inside the selection episode.

Relative to Gamma mass, constant-size ASMC reduces mean selected focal CDF error
by about 39%. Adding EAS history to ASMC reduces it by another 17%, for a combined
49% reduction. Both comparisons improve focal CDF MAE in all 10 selected regions.
Across all 100 spatial coordinates, selected CDF MAE is 12.91 pp for Gamma mass,
6.83 pp for ASMC constant, and 5.45 pp for ASMC EAS.

Neutral focal fractions below 20 kya are 4.40% in truth, 12.38% for Gamma mass,
5.38% for ASMC constant, and 3.98% for ASMC EAS. Mean neutral focal CDF MAE is
5.92, 1.50, and 1.35 pp respectively. EAS improves region-wide CDF MAE in every
neutral replicate, but focal CDF MAE in only 4/10 neutral replicates.

The correction is incomplete: ASMC EAS still assigns about four times too much
selected mass below 20 kya, and spreads too much beyond 50 kya. Its main excess
peak remains much younger than the true peak. This demo measures inference accuracy;
it does not establish TPR, FPR, F1, or genome-wide calibrated selection power.

## Figures and data

- [Focal distributions: PDF](01_focal_tmrca.pdf), [PNG](01_focal_tmrca.png).
  Selected and neutral cohort-mean CDFs and one-kyr masses; no uncertainty bands.
- [Selected-minus-neutral excess: PDF](02_selected_excess.pdf), [PNG](02_selected_excess.png).
  One-kyr and five-kyr bins. Shading marks the simulated 20–40-kya episode.
- `comparison.npz`: five estimators × 20 regions × 100 positions × 50 cutoffs;
  includes method names, region IDs, cohort labels, and positions.
- `replicate_metrics.tsv`, `summary.json`: numeric summaries.
- `manifest.json`, `provenance.json`, `decoding_receipts.json`: input identity,
  cache integrity, software, and timing. Large native decoding tables and raw
  haplotypes remain in the work directory, not this report.

## Validation and limitations

All 40 ASMC decodes passed posterior normalization and exact returned-pair identity
checks. Tree truth and source-file hashes were checked. Constant-model state means
match analytic truncated-exponential means; recent state boundaries correspond to
exactly 1 kyr. The effective sequence/CSFS flags and genetic-map units were checked.

ASMC normally outputs at input markers. To obtain exact 10-kb output coordinates,
the input adds observed invariant sites as anchors, including region endpoints.
Both tools retain the same biallelic variants; multiallelic sites are skipped to
match Gamma-SMC's reader. Anchors add no sequence variation, but ASMC's gap
approximations mean their insertion is not numerically invariant. A sensitivity
check on the first region compares common original SNPs with and without anchors:

| Check | ASMC constant | ASMC EAS |
|---|---:|---:|
| Mean absolute CDF difference across original markers | 0.0359 pp | 0.0319 pp |
| Mean difference within 400–600 kb | 0.0195 pp | 0.0191 pp |
| Maximum difference at the nearest original focal marker, 500,247 bp | 0.209 pp | 0.256 pp |
| Maximum anywhere, at 998,310 bp near the boundary | 3.69 pp | 6.76 pp |

The focal sensitivity is small relative to the estimator differences, but boundary
effects remain relevant to whole-region comparisons. Full details are in
`audit.json` and `audit.log`. The inference comparison also includes differences in
transition approximations and time representation; it does not isolate the gamma
posterior approximation alone.

Figures were visually inspected for clipping and overlapping labels. PDFs contain
embedded editable TrueType fonts and no raster image objects; see
`figure_validation.json`. PNGs and vector PDFs are saved, not displayed inline.

## Reproduce saved figures without raw simulations

Run from the repository root. This uses the checked summary cache and works without
ASMC installation or the original local simulation paths.

```bash
uv run --no-project --with numpy==2.4.6 --with matplotlib==3.11.2 \
  python scripts/demo_asmc_gamma.py report --from-cache
```

## Rerun native inference on the existing local simulations

Linux/WSL is required for the ASMC wheels. The source manifest's raw simulation
paths must be accessible. Preparation validates those inputs; it does not simulate
replacement data if they are missing. Use a fresh work directory for a changed
cohort or scientific configuration. Matching caches are verified and reused.

```bash
uv venv --python 3.11 /tmp/asmc-demo-env
uv pip install --python /tmp/asmc-demo-env/bin/python \
  asmc-asmc==1.4.0 asmc-preparedecoding==2.2.5 \
  numpy==2.4.6 matplotlib==3.11.2 tskit==1.0.3
PY=/tmp/asmc-demo-env/bin/python
WORK=/tmp/asmc-demo-work
$PY scripts/demo_asmc_gamma.py prepare --work "$WORK" --replicates 10
$PY scripts/demo_asmc_gamma.py quantities --work "$WORK"
$PY scripts/demo_asmc_gamma.py decode --work "$WORK" --workers 4
$PY scripts/demo_asmc_gamma.py audit --work "$WORK"
$PY scripts/demo_asmc_gamma.py report --work "$WORK"
```

ASMC setup took approximately 54 seconds for constant-size decoding quantities
and 166 seconds for EAS quantities. Individual 1,000-pair region decodes took
roughly 3–5 seconds; native logs and receipts retain the measured times.

References: [ASMC API](https://github.com/PalamaraLab/ASMC/blob/main/docs/asmc_python.md),
[sequence mode and demographic input](https://github.com/PalamaraLab/ASMC/blob/main/docs/asmc.md),
[PrepareDecoding API](https://github.com/PalamaraLab/PrepareDecoding/blob/main/docs/api.md).
