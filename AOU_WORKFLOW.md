# AoU within-individual recent-coalescence scan

## Scope and statistic

This branch implements the reduced scan requested for a large phased panel. For
each diploid, Gamma-SMC decodes only its two homologous haplotypes. At every
biallelic segregating SNP in the input VCF, the streaming summary reports

`mean_p_tmrca_lt_threshold = mean_i P(T_i < threshold | sequence data)`.

The default threshold is 4,500 years, which is 180 generations at the default
25 years/generation.
The Gamma posterior CDF is used, not a hard threshold of posterior mean TMRCA.
The summary also contains mean posterior TMRCA in generations. It has one row
per retained VCF segregating site and does not create a dense base-pair or 200 bp
grid. Indels and SVs are intentionally excluded from inference and can be joined
back by genomic position after candidate loci are defined.

## Install and test

The supported installation is repository-local and locked. It needs network
access once, but does not need root or an existing Python/Conda environment.
On Linux x86_64/AVX2 it installs uv, Python 3.11, SLiM 5.2, the native build
stack, builds Gamma-SMC, audits every executable/import, and runs the tests:

```bash
bash scripts/bootstrap_uv.sh
scripts/aou.sh --help
```

For simulation/calibration/plotting without the C++ decoder, use
`bash scripts/bootstrap_uv.sh --simulation-only`. On Windows, run
`powershell -ExecutionPolicy Bypass -File scripts/bootstrap_uv.ps1`; the C++
decoder itself must run under Linux/WSL2 or in the container. All subsequent
examples may replace `gamma-smc-aou` with the portable wrapper
`scripts/aou.sh` (Linux) or `scripts/aou.ps1` (Windows).

`uv.lock` pins the complete cross-platform Python graph. SLiM is pinned to 5.2
in `.native`; the Linux bootstrap also supplies the compiler, Boost, htslib,
and zstd from conda-forge/bioconda. Ordinary `.trees`/`.ts` files are loaded with
`tskit.load`; only `.tsz` files use `tszip.load`. Conversion uses an
argument-safe process launch rather than interpolating paths into a shell
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

### Official v0.2 container at 1 kb resolution

The locked uv application also wraps the upstream
`docker.io/regevsch/gamma_smc:v0.2` image. It auto-detects Apptainer,
Singularity, or Docker, requests within-individual pairs only, disables output
at every heterozygous site, and uses `--output_at_stride 1000`:

```bash
scripts/aou.sh decode-container \
  --input AFR.phased.vcf.gz --output AFR.within.stride1000.tsv \
  --theta 0.0005 --rho-over-theta 0.8 --mutation-rate 1.25e-8 \
  --generation-time 25 --threshold-years 4500 --output-at-stride 1000
```

The postprocessor reads one Gamma-SMC pair chunk at a time, evaluates the Gamma
posterior CDF at 4,500 years, and averages it across diploids. It then deletes
the large raw posterior by default. Input and output must share a host directory
so one directory can be mounted at `/work`.

The complete retained-sweep validation (selected pseudo-data plus matched nulls,
pointwise Monte Carlo p-values, breakpoints, plots, runtime, and truth error) is:

```bash
scripts/aou.sh run-container-study \
  --source-dir sim_results/two_epoch_growth_s0p05_n2000 \
  --output-dir sim_results/gamma_smc_container_stride1000 \
  --neutral-replicates 100 --output-at-stride 1000 --workers 4
```

The GitHub Actions workflow `Gamma-SMC container stride study` exposes the null
count and stride as manual inputs for machines without a local Linux container
runtime.

## 1b. Whole-genome scan over ~100,000 sampled haplotype pairs

The within-individual scan is capped at one pair per diploid. To reach ~100,000
pairs, sample them uniformly from every haplotype pair in the panel and keep the
per-pair record as one bit per position instead of a Gamma posterior.

```bash
for chrom in $(seq 1 22); do
  scripts/aou.sh decode \
    --executable bin/gamma_smc \
    --input AFR.chr${chrom}.phased.bcf --input-format vcf \
    --output scan/chr${chrom}.tsv \
    --bitmatrix scan/chr${chrom}.bits \
    --theta 0.0005 --rho-over-theta 0.8 --mutation-rate 1.25e-8 \
    --generation-time 25 --threshold-years 4500 10000 \
    --no-output-at-hets --output-at-stride 1000 \
    --n-random-pairs 100000 --pairs-seed 1729 \
    --threads 32 --mask callable.chr${chrom}.bed
done
```

`--n-random-pairs` draws distinct unordered pairs uniformly from all
`C(2N, 2)` haplotype pairs, within-individual pairs included; add
`--exclude-within` to drop them. The same `--pairs-seed` reproduces the same
list, and the list is written into the bit matrix metadata. Use `--pairs-file`
(two 0-based haplotype indices per line) to fix the pairs yourself, which is
what you want if every chromosome must use the same pairs.

**Set `--no-output-at-hets`.** With output at every segregating site, a large
panel produces millions of output positions and the per-thread posterior buffers
grow with it; at stride 1000 chr1 has ~249,000 columns instead of ~3,000,000.

### Per-position summary

`--output` is the small TSV. The first five columns are unchanged from the
within-individual scan, so every existing consumer keeps working, and
`mean_p_tmrca_lt_threshold` is the first threshold's mean probability:

| column | meaning |
|---|---|
| `position_0based`, `position_1based` | output position |
| `n_pairs` | pairs with a usable posterior at this position |
| `mean_p_tmrca_lt_threshold` | alias for `mean_p_lt_4500` |
| `mean_tmrca_generations` | mean posterior TMRCA |
| `n_recent_4500`, `frac_recent_4500` | pairs **called** recent, and the proportion |
| `mean_p_lt_4500` | mean of P(T < 4500 years) across pairs |
| `n_recent_10000`, … | the same block per additional threshold |

There are two families here, and they are not interchangeable:

- **`frac_recent_*`** counts pairs *called* recent. `--recent-call` chooses the
  rule: `median` (default, P(T<t) >= 0.5), `mean` (posterior mean below t), or
  `prob` with `--recent-call-probability`.
- **`mean_p_lt_*`** averages P(T<t) over pairs. It ignores `--recent-call`
  entirely.

**The paper uses neither of the defaults.** It reports "the proportion of
posterior means below a threshold of T", i.e. `frac_recent_*` with
`--recent-call mean`. The `AOU_run` baseline instead averaged the posterior CDF,
which is `mean_p_lt_*`. Both are in every summary, so all three are comparable
from one decode, but pick one deliberately and use it on both observed data and
nulls.

The default here is the median rule, because it is far less sensitive to
posterior width than the soft average: over a CV shift from 0.8 to 1.0 the soft
statistic's inflation over truth moves 2.5x -> 4.6x while the median rule's moves
1.25x -> 1.44x. Since real data and simulations will not share posterior widths
exactly, that difference does not cancel in a null-calibrated p-value.

Do **not** threshold the MAP. For Gamma the mode is `(alpha-1)/beta`, and alpha
sits at its floor of 1 in the recent regime, so the mode collapses to 0 and the
rule becomes a step function of alpha right on the grid boundary — decided by
the clipping, not by the data.

### Which pairs were drawn

Every random draw is written out literally, to a manifest next to the summary:

```
chr1.stride1000.tsv.pairs.tsv
```

No flag is needed — it is derived from the output path whenever
`--n-random-pairs` is used, so it cannot be lost by forgetting one.
`--pairs-manifest` overrides the location, and works for any pair mode.

```
# gamma_smc_pair_manifest_v1
# created_utc	2026-07-29T18:22:41Z
# input	/work/chr1.phased.vcf.gz
# mode	random
# n_pairs	100000
# n_haplotypes	20000
# n_samples	10000
# pairs_seed	1729
# n_random_pairs	100000
# exclude_within	false
# rng	mt19937_64; draw both haplotypes uniformly, reject i==j, deduplicate, sort ascending
# panel_digest	0x3f1a9c04e7b52d18
# pairs_digest	0x9c22b7ff0a41e6d3
# Reuse with --pairs_file to decode exactly these pairs again.
# hap_i	hap_j	haplotype_i	haplotype_j
17	4082	NA12878.1	NA20502.0
...
```

Every header line starts with `#`, which the pairs-file reader already skips,
and the first two columns are the haplotype indices it already parses. So a
manifest is a valid `--pairs-file` with no conversion — decode chr1, then hold
the pair set fixed across the other 21 chromosomes:

```bash
scripts/aou.sh decode --input chr2.vcf.gz --output chr2.tsv --pairs-file chr1.tsv.pairs.tsv --no-output-at-hets --output-at-stride 1000 --threads 0
```

**Why not just keep the seed.** The draw is reproducible from `--pairs-seed`,
but only against the same haplotype ordering. Indices are positions in the
sample list, which depends on the input file, on `--samples`, and on which
records htslib kept. Re-deriving a draw months later from a seed alone, against
a panel that has since been re-exported or re-subset, silently decodes
different pairs.

`panel_digest` is a hash of the ordered sample names. Reusing a manifest
against a panel whose digest differs is refused, not warned about:

```
Error: --pairs_file chr1.tsv.pairs.tsv was drawn against a different panel.
  panel digest: manifest 0x3f1a9c04e7b52d18, this run 0x71ce0b39a2fd4460
  manifest was drawn from: /work/chr1.phased.vcf.gz
  Haplotype indices are relative to the sample order, so these pairs
  would decode different haplotypes than the ones recorded.
  Pass --allow_panel_mismatch to proceed anyway.
```

Range validation alone would not catch this: the indices usually remain inside
the new panel's range, so the run would look perfectly healthy. `pairs_digest`
covers the pair list itself, for checking a manifest has not been edited.

The bit-matrix `.meta` also carries the full `pairs` array and `sample_names`,
so a bit matrix stays interpretable on its own; the manifest is the portable,
reusable form.

### Bit matrix

`--bitmatrix` writes one bit per (pair, position, threshold): ~6 GB per
chromosome before compression for 100,000 pairs at stride 1000 with two
thresholds, against ~200 GB for the equivalent raw alpha/beta. It is a
concatenation of independent zstd frames indexed by a `.meta` sidecar, so a
subset of pairs can be read without touching the rest of the file.

```python
from gamma_smc_aou import bitmatrix

meta = bitmatrix.read_meta("scan/chr2.bits")
counts = bitmatrix.position_counts("scan/chr2.bits")        # (2, n_positions)
carriers = bitmatrix.position_counts("scan/chr2.bits", meta, pair_indices=carrier_pairs)
profiles = bitmatrix.pair_profiles("scan/chr2.bits", [0, 17, 512])
```

`bitmatrix.to_frame` rebuilds the counts table, which is a direct check against
the decoder's own TSV. `scripts/aou.sh bitmatrix-summary` does the same from the
command line and takes `--pairs` to restrict the counts to a subset.

### Cost

The decode is parallel over blocks of pairs; `--threads 0` uses every core.
`--pair-block` sets how many pairs each work unit covers (default 256) and is
also the bit-matrix frame size. Per-thread scratch is roughly
`2 x n_positions x 8 x 4` bytes of posteriors plus `n_segments x 8` bytes of
emission types; the genotype matrix, the callability bitmap and the 490 MB
flow-field cache are shared. The decoder prints its own estimate before starting
and warns when `--output-at-hets` would blow it up.

### Measured

`.github/workflows/decoder-benchmark.yml` builds this branch and `AOU_run` side
by side and runs the same workload through each, as an ablation. On a 4-vCPU
`ubuntu-22.04` runner, 4,000 simulated diploids over 5 Mbp (23,986 segregating
sites), `--only_within`, one threshold, stride 1000, decode time only:

| configuration | decode | step |
|---|---|---|
| `AOU_run` baseline | 18.73 s | — |
| + CPU/memory work, still exact `gamma_p` | 18.30 s | 1.02x |
| + lookup tables | 2.08 s | **8.79x** |
| + 4 threads | 0.85 s | **2.47x** |
| **overall** | | **22.1x** |

Peak RSS 0.77 GB to 0.61 GB.

**Where the time actually went.** Per-phase, single-threaded:

| phase | baseline | + CPU work | + tables |
|---|---|---|---|
| emissions | 0.99 s | 0.29 s | 0.27 s |
| forward | 1.22 s | 1.03 s | 0.82 s |
| backward | 1.21 s | 1.03 s | 0.82 s |
| **statistic** | **15.30 s** | 15.95 s | **0.17 s** |

The baseline spends **82% of its decode inside `boost::math::gamma_p`**,
evaluating P(T<t) once per pair per position. That, not the SMC recursion, was
the bottleneck for this statistic.

This is worth stating plainly, because it bears on whether the rewrite was
warranted at all: **the message-passing core really is close to optimal.**
CPU-level work on it takes 3.42 s to 1.91 s, about 1.8x, and most of that is
the hoisted emission fill (0.99 s to 0.27 s) rather than anything inside the
flow field, whose SIMD kernel was already well tuned. Overall that shows up as
only 1.02x because `gamma_p` drowns it out. The wins here are the statistic
evaluation (90x), threads, and the output volume — not the algorithm.

Total wall time is a poor guide here — it is dominated by fixed startup
(tree-sequence conversion, reading, building the flow-field cache), identical
for both binaries.

**Throughput.** 20,000 random pairs, two thresholds, bit matrix: 22.6 s of CPU
in 5.75 s of wall time, i.e. 3.93x on 4 cores, or **0.058 s per Gbp per pair**.
At the same per-core efficiency, 3.1 Gbp x 100,000 pairs is roughly **35 minutes
on 32 cores**. Treat that as a projection: the flow-field cache is read with an
effectively random access pattern, so memory bandwidth, not arithmetic, is what
will limit 32 threads.

Budget the bit matrix from its raw size, not the compressed size in the
benchmark: a neutral simulation has almost no recent coalescence, so its bits
are nearly all zero and compress ~220x, which real data will not. Raw is
`n_pairs x n_positions x n_thresholds / 8` — 6.2 GB for chr1 and 77.5 GB
genome-wide at 100,000 pairs and two thresholds.

### Output agreement

Run with `--exp10 fast --backward_alignment legacy`, this branch reproduces the
old binary: `mean_p_tmrca_lt_threshold` to 1.3e-6 absolute, `mean_tmrca_generations`
to 4.7e-6 relative. Nothing changed by accident. **The default binary does not
reproduce the baseline** — that is the point of the two corrections below, and
their size is quantified there.

The lookup tables are exact where it matters. Against `--exact_recent_stats`,
which evaluates `boost::math::gamma_p` per element over the same decoded
posteriors, `mean_p_tmrca_lt_threshold` agrees to 1.3e-6 and `n_recent_4500` is
identical at **all 5,000 positions** — the hard-threshold counts, which are the
headline statistic, are not affected at all.

Thread count and `--pair_block` give bit-identical output.

### Two corrections to upstream numerics

Both are **on by default**. `--exp10 fast` and `--backward-alignment legacy`
reproduce the old binary, for comparing against results generated with it.

**`--exp10`** (default `accurate`) controls the `10^x` that converts the message
state into (alpha, beta). Upstream uses Schraudolph's bit trick, which replaces
`2^f` by the straight line `1+f` inside each binade. Measured over the reachable
range that is −3.89%…+2.01% relative error with a mean of +0.03% — a near
zero-mean sawtooth, not a systematic bias. It is nevertheless a deterministic
function of the value rather than noise, so it does not cancel across pairs
whose posteriors land in the same part of the sawtooth, and alpha and beta are
perturbed independently of each other.

It is close to free. `10^x` is evaluated once per output *position*, not once
per segment: on the benchmark panel that is 5,000 calls against ~29,000
flow-field steps, so it is a sixth of the call count and the accurate version's
extra arithmetic (about 11 AVX ops against 3) is diluted accordingly. Measured
single-threaded over three repeats per mode:

| | forward+backward | whole decode |
|---|---|---|
| `--exp10 fast` | 1.954 s | 2.436 s |
| `--exp10 accurate` | 2.009 s | 2.493 s |
| cost | +0.055 s (+2.8%) | +0.057 s (+2.4%) |

Read that as an upper bound rather than a point estimate: run-to-run spread on
a shared runner was ±7% for `fast` and ±11% for `accurate`, wider than the
difference itself. A few percent is the right order; anything larger would have
shown. Note the cost scales with output positions per segment, so it grows if
you drop the stride well below 1000.

**`--backward-alignment`** (default `fixed`) corrects an off-by-one. Because
`output_at_start[k] == output_at_end[k-1]`, the backward pass makes one fewer
write than the forward pass whenever the final segment is itself an output
position; every backward message is then paired with the forward message one
output position to its right, and position 0 receives none at all.

Whether that is reachable depends on the output mode, and this matters for
interpreting old results:

- **`--output-at-hets`** — every segment ending at a segregating site is an
  output position, including the last, so the shift **always** occurs. The
  within-individual scan in section 1 runs in this mode.
- **`--no-output-at-hets --output-at-stride 1000`** — the last segment ends at
  the last segregating site, which is a multiple of the stride only by
  coincidence, so the two settings normally agree. Measured identical over
  5,000 positions on the benchmark panel.

Leaving the final position with only a forward contribution is correct, not a
gap: the backward message there is the Exp(1) prior, and the
`alpha_f + alpha_b − 1` convention makes adding it a no-op.

**Measured effect of each correction**, on the benchmark panel, each against
upstream numerics with everything else held fixed:

| correction | mode | `mean_p` max abs | `n_recent_4500` |
|---|---|---|---|
| accurate `10^x` | stride 1000 | 3.9e-3 | 4,817/5,000 positions differ, max 11 calls, **total −7.0%** |
| backward alignment | stride 1000 | 0 | identical (unreachable, see above) |
| backward alignment | `--output-at-hets` | 1.4e-2 | 1,961/23,984 differ, max 17 calls, total −0.15% |

The `10^x` correction is the one that matters for the stride-1000 scan, and it
is not small: it removes about 7% of the recent calls in aggregate, i.e. the
old binary over-called recent coalescence by roughly that much. That is two
orders of magnitude larger than the lookup-table error, which changes no calls
at all.

**Re-run your nulls.** Observed data and the matched null replicates must be
decoded with the same `--exp10` and `--backward-alignment` settings, or the
p-values are meaningless. Nothing in the pipeline checks this for you.

## Parameter reference

### Required

Only three things, since the rates now have reference defaults:

| flag | notes |
|---|---|
| `--input` / `-i` | vcf, vcf.gz, bcf, `.trees`, `.tsz` |
| one of `--output`, `--recent_summary`, `--recent_bitmatrix` | |
| a pair selection | `--n_random_pairs`, `--only_within` or `--pairs_file`; otherwise exhaustive |

### Defaults

| flag | default | notes |
|---|---|---|
| `--input_format` | `auto` | required for `/dev/stdin` |
| `--allow_unphased` | off | not recommended for haplotype scans |
| `--recent_threshold_years` | `4500` | comma-separated for several |
| `--generation_time` | `25` | 4500 years = 180 generations |
| `--scaled_mutation_rate` / `-m` | `0.00075` | the paper's value; fixed, not per-file |
| `--scaled_recombination_rate` / `-r` | `0.0006` | rho/theta = 0.8 |
| `--unscaled_mutation_rate` | `1.25e-8` | with theta gives 2Ne = 30,000, Ne = 15,000 |
| `--estimate_mutation_rate` | off | estimating theta rescales the time axis per file |
| `--recent_call` | `median` | `median`, `mean`, or `prob` |
| `--recent_call_probability` | `0.5` | only read by `--recent_call prob` |
| `--no_recent_probability` | off | counts only; skips the mean-P accumulation |
| *pair selection* | *exhaustive* | see the footguns below |
| `--only_within` / `-w` | off | one pair per diploid |
| `--n_random_pairs` | `0` (off) | |
| `--pairs_seed` | `1729` | |
| `--exclude_within` | off | drop within-individual pairs when sampling |
| `--pairs_file` | — | |
| `--pairs_manifest` | `<summary>.pairs.tsv` | auto-derived when sampling at random |
| `--allow_panel_mismatch` | off | downgrades the panel-digest check to a warning |
| `--samples` / `-S`, `--samples_against` / `-T` | — | |
| `--output_at_hets` / `-h` | **off unless given** | ~1M positions/chr when on |
| `--output_at_stride` / `-s` | `100000` | the paper's region size; `-1` disables |
| `--mask` / `-a`, `--masks_per_sample` / `-b` | — | mutually exclusive |
| `--cache_size` / `-z` | `1000` bp | no gain above the typical inter-SNP distance |
| `--threads` / `-j` | `0` (all cores) | |
| `--pair_block` | `256` | pairs per work unit and per bit-matrix frame |
| `--exp10` | `accurate` | `fast` reproduces upstream |
| `--backward_alignment` | `fixed` | `legacy` reproduces upstream |
| `--exact_recent_stats` | off | validation only; roughly 90x slower |
| `--flow_field` / `-f` | built-in | |
| `--zstd_compression_level` | `1` | |
| `--only_forward` / `-y`, `--only_backward` / `-d` | `false` | |

### Python CLI

`scripts/aou.sh decode` takes the same values, hyphenated, with four differences:

- `--executable` defaults to `$GAMMA_SMC_BIN`, else `gamma_smc` on PATH.
- `--threshold-years` is space-separated: `--threshold-years 4500 10000`.
- `--no-output-at-hets` is a flag, inverting `--output_at_hets=false`.
- Pair selection has **no exhaustive fallback**: picking none raises.

### Build

| knob | default | notes |
|---|---|---|
| `MARCH` | `native` | use `x86-64-v3` for a portable AVX2 build |
| `OPENMP` | `1` | `OPENMP=0` builds single-threaded |

### Compiled-in constants

Not flags; change them in the source and rebuild.

| constant | value | file |
|---|---|---|
| `parallel_vector_size` | 8 | `common.h`; AVX2 lanes, sets pairs per byte |
| alpha exponent range | `[-16, 16]` | `recent_stats.h` |
| quantile-table mantissa | 8 bits, 8,194 entries | `recent_stats.h` |
| P(T<t) table | 1024 x 1024, 4.0 MB | alpha-coord `[-2,16]`, v `[-8,8]` |
| table accuracy | 4.9e-5 max abs vs Boost | asserted in CI |

### Three defaults that bite

**Pair selection defaults to exhaustive.** With no selector the binary builds
every O(n^2) haplotype pair, which at panel scale allocates before it fails.
Always pass one of `--n_random_pairs`, `--only_within`, `--pairs_file`. The
Python wrapper refuses instead of defaulting.

**Output defaults to every segregating site.** `--output_at_hets` is `true` and
the stride is off, so a whole-genome run without
`--output_at_hets=false --output_at_stride 1000` emits roughly 1M positions per
chromosome rather than 250k, and the bit matrix grows with it.

**`--cache_size` interacts with the stride.** It bounds how far the decoder can
skip in one step; there is no benefit above the typical distance between
segregating sites, and raising it enlarges the flow-field cache.

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

Compatible neutral results can be reused for another selection coefficient
without rerunning `s=0`:

```bash
gamma-smc-aou validate-recent-sweep \
  --output-dir sim_results/recent_sweep_s0p1_n2000 \
  --reuse-null-from sim_results/recent_sweep_10mb_n2000 \
  --population-size 10000 --sample-diploids 2000 \
  --sequence-length 10000000 --selection-coefficients 0.1 \
  --age-generations 180 --mutation-rate 1.25e-8 \
  --recombination-rate 1e-8 --selected-replicates 100 --workers 20 \
  --seed 424242 --no-save-trees
```

The command rejects a reused null if population size, sample size, sequence
length, sweep age, or mutation/recombination rate differs.

With seed 424242, the 20-worker `s=0.1` run finished in 284.2 seconds without
rerunning the neutral simulations. The focal allele was lost in 84/100
replicates, segregating in 16/100, and fixed in none. The unconditional
exact-center mean fraction with TMRCA below 180 generations increased from
0.009185 under neutrality to 0.073995, but the unconditional power was only
0.17 because losses remain in the denominator. All 16 replicates in which the
focal allele was present were significant at the minimum attainable Monte
Carlo p-value, 1/101 = 0.009901; one of the 84 loss replicates was also
significant (false-positive fraction 0.0119). Thus conditional power in this
run was 16/16, while unconditional power was 17/100. The checked-in power
summary reports both estimands explicitly. Future rejection sampling on allele
retention would estimate power conditional on presence in the present-day
sample, not unconditional evolutionary power.

To separate conditional calibration from the within-selected carrier effect,
the first 10 retained `s=0.1` replicate IDs can be analyzed against all 100
saved neutral replicates:

```bash
gamma-smc-aou analyze-retained-sweeps \
  --source-dir sim_results/recent_sweep_s0p1_n2000 \
  --output-dir sim_results/recent_sweep_s0p1_retained10_n2000 \
  --selection-coefficient 0.1 --retained-replicates 10 \
  --workers 20 --seed 424242
```

The source run stripped the selected mutation before overlaying neutral
mutations and did not retain its tree files. The command therefore
deterministically reconstructs only the selected replicate IDs `0, 2, 14, 17,
32, 34, 35, 41, 44, 62` from their original seeds. It fails unless both the
saved population allele frequency and exact-center statistic are reproduced.
This reconstructs the same trajectories; it does not draw 10 replacement
mutations. With 20 workers requested, 10 workers were active because there are
only 10 trajectories, and reconstruction took 22.0 seconds.

In the resulting calibration, the 100 neutral simulations contained 8--29 of
2,000 pairs with TMRCA below 180 generations (mean 18.37). The 10 retained
selected trajectories contained 39--1,703 recent pairs. No neutral replicate
equaled or exceeded any selected result, so each has upper-tail Monte Carlo
`p=(1+0)/(100+1)=0.009901`. The retained trajectories were chosen conditional
on allele survival, so these p-values demonstrate separation for that
conditional set; they are not an estimate of unconditional evolutionary
power.

The within-selected comparison distinguishes three pair classes: neither
haplotype carries the focal mutation, one haplotype carries it, or both carry
it. The primary carrier comparison is two-copy carrier--carrier versus
zero-copy noncarrier--noncarrier; mixed pairs are shown separately rather than
being pooled with carrier pairs. All carrier--carrier pairs had TMRCA below
180 generations, mixed pairs had none, and the noncarrier recent fraction was
0--0.0164. For a single-origin allele introduced 180 generations ago, this is
partly a genealogical identity: two present-day carriers must share the
mutation-bearing ancestral lineage within that age. Consequently, the
carrier--carrier contrast alone is not independent evidence of selection. The
selection test is the excess *all-pair* recent coalescence against the neutral
simulation null; an empirical carrier contrast should additionally use
frequency-, age-, and recombination-matched neutral alleles for calibration.

The same retained trajectories can be plotted as hom-alt versus hom-ref
within-diploid TMRCA profiles across the full chromosome segment and around the
selected site:

```bash
gamma-smc-aou plot-retained-carrier-profiles \
  --source-dir sim_results/recent_sweep_s0p1_n2000 \
  --output-dir \
    sim_results/recent_sweep_s0p1_retained10_n2000/hom_alt_vs_hom_ref_profiles \
  --selection-coefficient 0.1 --retained-replicates 10 \
  --workers 20 --seed 424242 \
  --full-step 50000 --zoom-half-width 500000 --zoom-step 5000
```

Each of the 10 figures contains a 10-Mb panel and a selected-site +/-500-kb
panel. Curves are the mean within-diploid TMRCA among individuals with two
focal alternate haplotypes (hom alt) or two focal reference haplotypes (hom
ref); heterozygous pairs are excluded. Shading is a two-sided 95% log-Wald
confidence interval for the positive mean, using a Student-t critical value.
The log-scale interval avoids the negative lower limits produced by an
additive normal interval for these heavy-tailed TMRCA distributions. Full
profiles are evaluated every 50 kb and the zoom every 5 kb.

At the 5-Mb selected site, mean hom-alt TMRCA ranged from 64.0 to 168.9
generations across the 10 retained replicates, whereas mean hom-ref TMRCA
ranged from 11,680 to 31,492 generations. Every replicate shows the expected
central hom-alt trough, with variable width from recombination and realized
trajectory. Replicate 41 has only 19 hom-alt individuals and therefore visibly
wider intervals; replicate 44 has only 17 hom-ref individuals. The command
records sample sizes in every title, writes all profile points to one TSV, and
also creates a 20-panel overview plus one full-resolution two-panel PNG per
replicate.

## Two-epoch recent-growth validation

The constant-size experiment remains unchanged. A separate experiment models
forward growth from 10,000 to 20,000 individuals 100 generations before the
present, with the selected mutation introduced 180 generations (4,500 years)
ago while the population was still 10,000:

```bash
gamma-smc-aou validate-two-epoch-growth \
  --output-dir sim_results/two_epoch_growth_s0p05_n2000 \
  --ancestral-population-size 10000 --present-population-size 20000 \
  --size-change-generations-ago 100 \
  --sample-diploids 2000 --sequence-length 10000000 \
  --variant-age-generations 180 --generation-time-years 25 \
  --selection-coefficient 0.05 \
  --mutation-rate 1.25e-8 --recombination-rate 1e-8 \
  --neutral-replicates 100 --workers 20 \
  --max-selected-attempts 1000 --minimum-hom-alt-pairs 2 \
  --full-step 50000 --zoom-half-width 500000 --zoom-step 5000 \
  --seed 515151
```

Backward in time, the coalescent population size is therefore 20,000 from 0
to 100 generations and 10,000 earlier. The analytical neutral probability of
coalescing within 180 generations is
`1-exp[-100/(2*20000)-80/(2*10000)]=0.006479`. Across the 100 neutral
simulations, the observed mean was 0.006555 (13.11 of 2,000 pairs), with range
0.0030--0.0115, bias 0.000076, and RMSE 0.001669. This is 28.6% lower than the
0.009185 mean from the separate constant-10,000 null. Leave-one-out rejection
was 0.04 at p <= 0.05. Because the statistic is discrete, the raw-rank KS test
was conservative; randomized-tie ranks had mean p=0.498, KS p=0.960, and
exactly 0.05 at p <= 0.05.

Selected trajectories were evaluated in deterministic seed order and rejected
when the `s=0.05` allele was lost or when fewer than two sampled hom-alt pairs
were available for a confidence interval. Attempts 0--7 were lost; attempt 8
was the first accepted trajectory, so nine seed-ordered attempts were needed.
The 20-worker batch computed attempts 0--19, but later attempts did not affect
which trajectory was accepted. The retained allele had population frequency
0.1021 and sample frequency 0.10275, producing 21 hom-alt, 369 heterozygous,
and 1,610 hom-ref diploids.

The retained selected observation had 37/2,000 recent pairs,
`P(TMRCA<4500 years)=0.0185`. No neutral simulation reached that value, giving
`p=(1+0)/(100+1)=0.009901`. At the selected site, hom-alt mean TMRCA was 130.9
generations (95% CI 117.8--145.4), versus 28,163 generations for hom ref (95%
CI 27,170--29,193). The output includes separate demographic-timing, neutral
calibration, and full/zoom carrier-profile plots. This selected p-value is
conditional on rejection sampling an allele that survives and has an
analyzable hom-alt class; it is not unconditional power for a newly arising
`s=0.05` mutation.

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
