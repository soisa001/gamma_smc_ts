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
  --save-trees --seed 1729 --workers 20
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

For the 2,000-diploid neutral truth comparison:

```bash
gamma-smc-aou plot-truth \
  --summary-dir sim_results/n2000_sensitivity/truth_summaries \
  --sequence-length 200000 --ne 10000 \
  --output-dir sim_results/n2000_sensitivity/truth_plot
```

This compares each replicate's true midpoint recent fraction and mean TMRCA to
the unascertained random-locus standard-coalescent references `2Ne` and
`1 - exp(-threshold_generations / (2Ne))`. Because the plotted point is the
nearest segregating site, those lines are references rather than exact
segregating-site-conditioned expectations.

An optional SLiM 5 hard-sweep validation follows the official conditional-
fixation/tree-sequence recipe. It introduces one strongly beneficial copy,
restores the pre-introduction checkpoint whenever that copy is lost, outputs on
fixation, and recapitates only remaining roots:

```bash
python -m pip install -e '.[selection]'
gamma-smc-aou validate-sweep --output-dir sim_results/hard_sweep \
  --population-size 200 --sequence-length 100000 \
  --selection-coefficient 0.5 --neutral-replicates 39
```

The validation requires SLiM on `PATH` or in `SLIM_BIN`. It plots the local true
TMRCA trough, the corresponding excess of recent within-individual pairs, and
the matched fixed-standard-coalescent null. The automated integration test
requires the selected center to have lower mean TMRCA than its flanks, a
center-to-flank ratio below 0.5, and a one-sided Monte Carlo p-value at most 0.05.

For an age-controlled, paper-scale experiment, use the separate recent-sweep
command. It initializes the full population with an msprime
`StandardCoalescent`, annotates that tree sequence with `pyslim`, introduces one
selected copy at 5 Mb, and runs exactly 180 forward Wright-Fisher generations.
The trajectory is unconditional: the focal allele may be lost, segregating, or
fixed. A lost allele remains at frequency zero, but its 5-Mb coordinate is still
scanned normally. The outcome and realized frequency are written to the results
and must be reported with the p-value. The `s=0` null uses the same 180-generation
SLiM path; every condition starts from an independent msprime
`StandardCoalescent` ancestry, so the p-value does not confound selection with a
different simulator path.

```bash
gamma-smc-aou validate-recent-sweep \
  --output-dir sim_results/recent_sweep_10mb_n2000 \
  --population-size 10000 --sample-diploids 2000 \
  --sequence-length 10000000 --selection-coefficients 0 0.001 0.01 \
  --age-generations 180 --mutation-rate 1.25e-8 \
  --recombination-rate 1e-8 --neutral-replicates 100 \
  --selected-replicates 100 --workers 8
```

The exact-center primary statistic is the fraction of the 2,000 within-diploid
pairs with true TMRCA below 180 generations. The command also evaluates 10-kb
and 100-kb averages, plots the full 10-Mb region so the neutral plateaus are
visible, performs upper-tail Monte Carlo tests, and writes leave-one-out null
calibration, theoretical neutral bias/RMSE, AUC, and power. These are genealogy
truth tests; decoder error still requires `evaluate-decoder` after Gamma-SMC is
run on the mutated tree sequences.

The checked-in 100-neutral/100-per-selected-condition run is an unconditional
power experiment, not a forced-success hard-sweep demonstration. It includes
lost focal alleles in the denominator. This is the relevant design for asking
whether a new selected mutation beginning 180 generations ago creates a
detectable all-sample within-diploid signal; a forced-fixation 50-kb fixture
must not be used to claim that power.

With seed 8675309, the focal allele was lost in 98/100 `s=0.001` replicates and
97/100 `s=0.01` replicates; none fixed. The exact-center statistic had neutral
mean 0.009185, selected means 0.008635 and 0.009225, AUC 0.418 and 0.511, and
power at p < 0.05 of 0.02 and 0.01, respectively. Every exact-center selected
hit at p < 0.05 occurred in a replicate where the focal allele was lost. The
100-kb statistic reached 0.05 and 0.09, but its AUC values were only 0.424 and
0.515, so the latter is compatible with Monte Carlo fluctuation rather than a
reliable sweep signal. Null leave-one-out rejection was 0.04 for p < 0.05 in
the 10-kb and 100-kb tests; the exact-center statistic is more discrete and had
0.01 below 0.05 (0.05 at or below 0.05).

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
