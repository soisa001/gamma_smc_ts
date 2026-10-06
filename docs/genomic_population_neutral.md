# Genomic-region population nulls

This campaign schedules 1,000 neutral, introgressed focal-allele simulations for
each of AFR, AMR, EAS, EUR, MID and SAS (6,000 total). It starts a new campaign;
the earlier median-history campaigns remain paused. There is no selection and
no population-size/time rescaling.

## Sampling and maps

The base seed is 20380101. Region sampling uses a separate SeedSequence domain
(`0x52454749`) from the demographic normals. Propose a uniform integer window start over all valid
11-Mb starts on GRCh38 autosomes 1–22; chromosome probability is proportional
to the number of possible starts. Reject windows with less than 90% usable
coverage. The first 1,000 accepted windows cover every autosome and are shared
across populations for comparison. Windows can overlap. Simulate the padded
11 Mb, choose the nearest eligible focal allele at onset, then retain the
centered 10 Mb, preserving the existing ascertainment rule.

The recombination input is UCSC's sex-averaged deCODE GRCh38 map:
[recombAvg.bw](https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw).
Use its original interval resolution, including hotspots, with cM/Mb multiplied
by 1e-8 to obtain per-base, per-generation rates. Apply the same map during
msprime prehistory and SLiM forward simulation. Unknown map intervals receive
the explicit 1e-8 recombination fallback and are excluded from observed mutations.

No suitable absolute mutation-rate map was present in this checkout or the
inspected PHLASH resources. The requested fallback is **1.29e-8 mutations per
base per generation** on observable sequence. The provisional simulation
observation mask excludes GRCh38 assembly gaps and missing recombination-map
coverage. It is **not the HMMix `strick`/strict mask**, nor a claim to reproduce
sample-level sequencing callability. Mutation overlay omits excluded bases;
Gamma-SMC receives the corresponding included-base BED after focal cropping.
Missing observations are therefore not treated as invariant, called sequence.

The native Gamma-SMC executable currently accepts a scalar recombination rate.
It uses each cropped region's length-weighted mean rate divided by 1.29e-8 for
the r/mu ratio. Thus simulated ancestry includes hotspots, but decoding assumes
the regional mean. This approximation is recorded in every decoder receipt.
Both `exp10=accurate` and `backward_alignment=fixed` remain enabled.

The downloaded resource hashes, intervals, coverage, genomic coordinates and
per-region map hashes are stored in `regions/manifest.json` and `regions.csv`.
Resource contents, parameters, manifests and validated outputs control cache
reuse; a new Git commit alone does not invalidate simulations.

## Demography and conditioning

Generate 1,000 diagnostic draws per population using the archived 10,000-point
log-Ne mean and covariance factor: `log_Ne = mean_log_ne + z @ factor`.
The optional archived jitter is a log-Ne standard deviation. No dense
10,000-by-10,000 covariance matrix, independent time-point sampling, clipping,
or forced quantile matching is used. The factor preserves temporal correlations.
Population seeds are base seed plus population index in the order above.

The default simulation allocation is **100 demographic histories × 10 regions**
per population, using the first 100 of the 1,000 seeded draws. Replicate `r`
uses history `r % 100`. History files preserve all original time points and
hold the nearest endpoint constant outside the archived grid. Generation time
is 25 years. Other assumptions remain: archaic Ne 3,600, archaic split 500 kya,
2% introgression at 50 kya, 200 sampled diploids and 10,000 fixed sampled
haplotype pairs (pair seed 1729).

The focal allele is an archaic-derived, biallelic mutation segregating in the
recipient population immediately after introgression at 50 kya. It is chosen
before neutral forward evolution. Loss or absence from the present sample
causes a logged retry under the **same history and genomic region**. Present
fixation is retained. This remains a conditional focal-site null, not an
unconditional distribution over arbitrary present-day genomic sites.

## CDFs and demographic bands

Report true and decoded TMRCA separately, for all pairs and raw ALT/ALT pairs.
The scored thresholds remain 10, 20, 30, 40 and 50 kya; the requested CDF figures
highlight 10 and 50 kya with empirical 95th/99th percentile score cutoffs.
Finite-null, conservative empirical-p boundaries are also provided and need
not exactly equal those descriptive percentiles, especially with ties.

CDF reporting uses only complete history groups, each weighted equally. Plot
the mean and median of the history-specific CDFs and their pointwise 2.5–97.5%
range. With ten focal scores per history, these bands also contain substantial
Monte Carlo variation. They are not pure demographic uncertainty, simultaneous
bands, or a confidence interval for the mean CDF. The 1% threshold is estimated
from the pooled 1,000 nulls; ten scores cannot resolve a within-history 1% tail.
All completed scores are exported separately, including incomplete groups.
Missing ALT/ALT scores remain in the null denominator below the valid score
support; a missing target score is a no-call.

## Execution

Use Linux/WSL and the existing compiled SLiM 5 and Gamma-SMC binaries. The
existing `.venv-origin` environment includes selection dependencies; install
the additional pinned map reader with:

```bash
uv --no-config pip install --python .venv-origin/bin/python pyBigWig==0.3.25
```

For a new environment, install the project with `[selection,genomic]` extras.
Set `GAMMA_NATIVE_DIR`, `SLIM_BIN`, `ORIGIN_PYTHON` and `GENOMIC_OUTPUT_DIR` when
using a different machine. The phases are restartable:

```bash
bash scripts/run_genomic_neutral.sh prepare
bash scripts/run_genomic_neutral.sh plan
bash scripts/run_genomic_neutral.sh run
# Refresh reports from completed, validated files without starting simulations:
bash scripts/run_genomic_neutral.sh report
```

The runner uses up to 20 workers, 180 GB aggregate worker address-space budget,
and a 200 GB total planning cap that reserves space for the paused campaigns.
Per-task stdout/stderr, retry reasons, exact seeds, input/output hashes, model
scripts and native decoder commands are retained. `run_status.json` reports
live tasks and failures. A lock prevents two runners using the same campaign.

On October 6, after an AFR tree-simplification allocation failure under the
9-GB per-worker limit, the campaign resumed with eight workers and a 160-GB
worker budget (20 GB per worker). This reserves 40 GB of the 200-GB planning cap
for paused campaigns and other overhead. Resource settings do not change the
scientific cache identity, seeds or existing completed results:

```bash
bash scripts/run_genomic_neutral.sh run --workers 8 --worker-memory-budget-gb 160
```

Validation on October 5: 25 existing relevant tests and four new map/history/CDF
tests passed. A separate small end-to-end run exercised variable recombination,
mutation masking, focal ascertainment, cropping, truth scoring and the fixed
decoder. Its tree sequence contained no observed mutations in the excluded
interval. It is a development fixture and is not counted in the 6,000 nulls.

The archived-vs-MVN diagnostic PDFs were rendered and visually inspected.
Across all grid points, typical (median absolute relative) differences in the
median curve are 0.3–1.8%; typical differences in 95% band edges are 0.7–3.6%.
Localized edge differences reach about 22%, especially away from the recent
bottleneck. The Gaussian approximation is therefore close overall, not exact.

The actual first-100-history subset was also plotted and visually checked at
`used_history_diagnostics/`. Its typical median differences are 0.6–1.6%, and
typical band-edge differences are 1.0–3.8%; localized edge differences reach
30%. These are differences across time points, not confidence levels or errors
in the MVN sampler. Both the empirical non-Gaussian fit distribution and finite
draw counts can contribute.
