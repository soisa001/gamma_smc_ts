# AoU within-individual recent-coalescence scan

## Scope and statistic

This branch implements the reduced scan requested for a large phased panel. For
each diploid, Gamma-SMC decodes only its two homologous haplotypes. At every
biallelic segregating SNP in the input VCF, the streaming summary reports

`mean_p_tmrca_lt_threshold = mean_i P(T_i < threshold | sequence data)`.

The default threshold is 4,500 years (150 generations at 30 years/generation).
The Gamma posterior CDF is used, not a hard threshold of posterior mean TMRCA.
The summary also contains mean posterior TMRCA in generations. It has one row
per retained VCF segregating site and does not create a dense base-pair or 200 bp
grid. Indels and SVs are intentionally excluded from inference and can be joined
back by genomic position after candidate loci are defined.

## Install and test

The C++ decoder remains Linux/AVX2 software. Tree-sequence input additionally
needs Python with `tskit`, and `.tsz` needs `tszip`.

```bash
python -m pip install -e '.[test,tszip]'
make clean && make
pytest -q
```

The Dockerfile installs the tree-sequence dependencies. Ordinary `.trees`/`.ts`
files are loaded with `tskit.load`; only `.tsz` files use `tszip.load`. Conversion
uses an argument-safe process launch rather than interpolating paths into a shell
command. For stdin, specify `--input_format`; auto-detection is deliberately not
attempted on a byte stream.

## 1. Empirical within-sample scan

Run one pair for every diploid. `theta` is the Gamma-SMC scaled mutation rate;
`mu` is needed to convert its coalescent time scale back to generations.

```bash
gamma-smc-aou decode \
  --executable bin/gamma_smc \
  --input AFR.phased.bcf --input-format vcf \
  --output AFR.within.tsv \
  --theta 0.0005 --rho-over-theta 0.8 --mutation-rate 1.25e-8
```

The input contract is diploid, phased, biallelic SNP data. Missing alleles are
handled as missing emissions. Unphased heterozygotes fail by default. Use a BED
mask through `--mask`, or per-sample BED files through `--masks-per-sample`.
Use exactly the same callable-region rule for empirical and simulated decoding.

Raw alpha/beta output is optional (`--raw-output`). Omitting it is important at
AoU scale: only five aggregate columns per segregating site are written.

## 2. Neutral simulations

All ancestry simulations explicitly use `msprime.StandardCoalescent()`; DTWF is
not used. A constant demography is the baseline:

```bash
gamma-smc-aou simulate --output-dir sims/AFR/region_001 \
  --replicates 1000 --diploids 2000 --length 1000000 \
  --ne 10000 --mutation-rate 1.25e-8 --recombination-rate 1e-8 \
  --save-trees --seed 1729
```

Mutation and recombination maps are tab-separated files with `position` and
`rate`, beginning at 0 and ending at the simulated length. Piecewise demographic
draws use `draw,time_generations,ne`, as in `examples/phlash_histories.tsv`.
Each replicate cycles through the provided PHLASH MVN draws. This propagates
demographic uncertainty into the null; omit `--histories` for the median/fixed
history analysis. Run both as a sensitivity analysis rather than silently mixing
their null distributions.

The current Gamma-SMC transition kernel accepts one scaled recombination rate per
run, not a native recombination map. For a fine-scale map, split empirical and
simulated data at map intervals (or scientifically chosen regions), decode each
interval with its local rate, and exclude a boundary buffer. Simulating a map but
decoding the whole region with one average rate is a robustness test, not exact
model matching.

The manifest records `model=StandardCoalescent`, history ID, independent ancestry
and mutation seeds, site/tree counts, and elapsed time. The first replicate also
produces a four-panel QC plot with demographic history, SFS, adjacent-site LD,
and the true within-individual first-coalescence curve.

## 3. Decode simulations and calculate p-values

Decode every saved `.trees` file with the same Gamma-SMC settings as the data.
This is naturally a scheduler array: one replicate per task. The full script has
an example. Then calculate one null value per simulation for each empirical site:

```bash
gamma-smc-aou calibrate \
  --observed AFR.within.tsv --sim-glob 'sims/AFR/region_001/decoded/*.tsv' \
  --observed-length 1000000 --simulation-length 1000000 \
  --match nearest --output AFR.region_001.scan.tsv
```

`nearest` matches relative position to the nearest simulated segregating site.
`region_max` compares every observed site with each replicate's maximum and is a
more conservative region-wide scan. Simulated sites are never pooled as if they
were independent replicates. The upper-tail Monte Carlo p-value is

`(1 + number of simulations with S_sim >= S_observed) / (R + 1)`.

With 1,000 simulations the smallest possible p-value is 1/1001. The output also
contains BH q-values, null mean/SD, and number of usable replicates. Treat this as
candidate-locus calibration, not genome-wide significance: LD makes adjacent
sites dependent and 1,000 replicates cannot resolve stringent genome-wide tails.
Cluster candidate SNPs into loci and use region maxima or additional simulations
for final locus-level inference.

## Verification and error summaries

`validate-null` performs an exchangeable leave-one-replicate-out rank test at a
fixed relative position. Its JSON reports mean p-value, KS uniformity p-value,
and the empirical fraction below 0.05; the QQ plot makes tail problems visible.
It reports both the actual conservative `>=` p-value and a randomized-tie rank
used only for diagnosis. A uniform randomized rank plus conservative raw ranks
indicates a discrete statistic, not simulator failure.
This validates the simulator/statistic/calibration machinery. It does **not**
validate Gamma-SMC decoder error. For that, compare each decoded simulation TSV
to the corresponding `truth_summaries` TSV:

```bash
gamma-smc-aou evaluate-decoder --truth-dir sims/region/truth_summaries \
  --decoded-dir sims/region/decoded --output-dir sims/region/decoder_evaluation
```

This reports MAE, RMSE, bias and correlation of `P(T < 4500 years)`, plus the
residual distribution and first true-versus-decoded coalescence curve.

Run calibrations separately for all six populations, using their own PHLASH
history/rates/missingness. A simulation with fewer diploids and all pairwise
comparisons is not an exact substitute for 2,000 disjoint within-individual
pairs: pair dependence and null variance differ. The primary null should match
the empirical 2,000-within-pair statistic; reduced panels are a sensitivity or
pilot analysis.

Within-individual coalescence is a fast locus-screening statistic, not the final
allele-carrier contrast. Use it to reduce loci/samples, then run focused
carrier-vs-carrier and carrier-vs-reference scans for candidate SNPs or overlaid
SVs.

## Minimal and full entry points

- `scripts/run_minimal.ps1` is a CPU-only truth/calibration smoke and scaling run;
  it does not need a GPU or compiled Gamma-SMC.
- `scripts/run_full_hpc.sh` simulates the full null, decodes empirical data, and
  documents the scheduler-array decode/calibration steps.

Run directories are self-contained under `sim_results/`. Their plots and TSV/JSON
summaries are intended to be inspected before a full six-population launch.
