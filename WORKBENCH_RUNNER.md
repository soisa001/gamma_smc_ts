# All of Us Workbench runner

`scripts/run_aou_workbench.sh` stages controlled inputs from their v9 GCS
buckets, samples 100,000 arbitrary within-population haplotype pairs, validates
every chromosome, uploads aggregate and pair-level candidate outputs, and
makes separate plots for each population. It is CPU-only and runs
populations/chromosomes sequentially while using 12 decoder threads by default.

## Fresh notebook cell

This is the complete setup-and-run cell for a fresh Researcher Workbench cloud
environment. Change the final scope to `-chr 1 -pops afr` for the chromosome 1
pilot before launching all 132 population/chromosome jobs.

```bash
%%bash
set -euo pipefail
REPO=/home/jupyter/gamma_smc_ts
REMOTE=git@github.com:soisa001/gamma_smc_ts.git

if [[ -d "$REPO/.git" ]]; then
  git -C "$REPO" remote set-url origin "$REMOTE"
  git -C "$REPO" fetch origin AOU_run_opt
  git -C "$REPO" switch AOU_run_opt
  git -C "$REPO" pull --ff-only origin AOU_run_opt
else
  git clone --branch AOU_run_opt "$REMOTE" "$REPO"
fi

bash "$REPO/scripts/bootstrap_uv.sh" --skip-tests
bash "$REPO/scripts/run_aou_workbench.sh" -chr all -pops all
```

All GCS reads and writes pass `--billing-project "$GOOGLE_PROJECT"` by
default so requester-pays controlled buckets work in Researcher Workbench.
Set `AOU_GAMMA_BILLING_PROJECT` or pass `--billing-project PROJECT_ID` to
override that project explicitly.

The runner intentionally requires both scope arguments. Values are
case-insensitive, and comma lists work too:

```bash
# Pilot
bash scripts/run_aou_workbench.sh -chr 1 -pops afr

# Selected chromosomes and populations
bash scripts/run_aou_workbench.sh -chr 1,2,22 -pops AFR,EUR

# Resolve every path and parameter without accessing GCS
bash scripts/run_aou_workbench.sh -chr all -pops all --dry-run
```

## Suggested Workbench VM

Use an `n2-highmem-16` CPU VM (16 vCPUs, 128 GB RAM), no GPU, with a 1 TB
balanced persistent disk as the starting configuration. This leaves four
vCPUs outside the 12-thread decoder and gives the full-panel genotype data,
plots, and candidate conversion room to coexist. Google currently lists
`n2-highmem-16` as 16 vCPUs and 128 GB RAM in the
[N2 machine-type table](https://docs.cloud.google.com/compute/docs/general-purpose-machines#n2_high-mem).

The controlled BCF object sizes cannot be verified outside Workbench, so 1 TB
is a starting point rather than a guaranteed minimum. Before each download the
runner reads the actual GCS BCF/CSI sizes and requires those bytes plus a 20 GiB
reserve. Before candidate replay it independently requires 16 bytes per
requested-position/pair cell plus 5 GiB. Increase the disk if either preflight
stops the run; 2 TB is the low-intervention choice if disk resizing later is
undesirable.

## Default input contract and paths

The source is the phased lrWGS phase-2 full panel. One BCF is downloaded per
requested chromosome and reused sequentially across the requested populations,
then removed by default. The files may contain SNVs, indels, and SVs; the native
reader retains only biallelic segregating SNPs for Gamma-SMC and rejects an
unphased heterozygous SNP.

The runner reads `ancestry_pred_other`, retaining only `AFR`, `AMR`, `EAS`,
`EUR`, `MID`, and `SAS`. `oth`, missing ancestry rows, and every ID in either
exclusion table are dropped. It obtains the chromosome's actual BCF sample
order with `bcftools query -l`, intersects that order with the filtered
ancestry set, and writes one deterministic sample list per population and
chromosome. Thus an ancestry-table participant absent from the long-read BCF is
also excluded.

The complete default controlled-input and output layout is:

| Artifact | Default path |
|---|---|
| Full-panel phased BCF | `gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr{chr}.bubble.split.bcf` |
| Required BCF index | the full-panel BCF path plus `.csi` |
| Ancestry assignments | `gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv` |
| QC exclusions | `gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv` |
| Relatedness exclusions | `gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv` |
| Exclusion hard mask | `gs://rw-migration-aou-rw-fa99430f/hardmask.hg38.v4.over99.bed` |
| HMMIX strict callable mask | `<WORKSPACE_BUCKET>/hmmix-static/hg38_strick_callability_mask.bed` |
| Staging/results root | `/home/jupyter/gamma_smc_workbench` |
| Chromosome outputs | `gs://rw-migration-aou-rw-fa99430f/gamma_smc/results/{POP}/chromosomes/` |
| Plot outputs | `gs://rw-migration-aou-rw-fa99430f/gamma_smc/results/{POP}/plots/{scope}/` |
| Combined report | `gs://rw-migration-aou-rw-fa99430f/gamma_smc/results/summary/{scope}/` |
| Callable-mask QC | `gs://rw-migration-aou-rw-fa99430f/gamma_smc/results/shared/masks/` |
| Gene labels | GENCODE v50 basic GRCh38 GTF from `ftp.ebi.ac.uk` (staged once locally) |

`WORKSPACE_BUCKET` overrides this default bucket, and `--output-prefix` or
`AOU_GAMMA_OUTPUT_PREFIX` overrides the complete results prefix.

The supplied hard mask describes bases to exclude, whereas Gamma-SMC's
`--mask` accepts bases to include. The runner reads the single-contig name and
length from the BCF/CSI, merges and clips that chromosome's hard-mask
intervals, writes their exact complement using the BCF contig name, and audits
excluded and callable bases. Passing the raw hard mask directly would invert
the intended filter.

The default remains this complemented `hardmask.hg38` analysis under
`/home/jupyter/gamma_smc_workbench` and `gamma_smc/results`. Two mask-mode flags
provide isolated sensitivity analyses:

- `--strict-hardmask` uses the positive HMMIX hg38 strict callable BED without
  complementing it. Local work goes under
  `/home/jupyter/gamma_smc_workbench_strict_hardmask`, and cloud results go
  under `gamma_smc/results_strict_hardmask`.
- `--both-hardmask` runs the default analysis first and the strict analysis
  second, each in the corresponding namespace. It means two separate runs,
  not an intersection of the masks.

If `--local-root` or `--output-prefix` is overridden, strict mode appends
`_strict_hardmask` to that supplied base so the caches cannot collide.

For each chromosome, the output directory receives:

- `chrN.gamma_smc.tsv`: aggregate Gamma-SMC scan;
- `chrN.gamma_smc.tsv.run.json`: exact native command, decoder output, and runtime;
- `chrN.pairs.tsv`: ordered manifest for the 100,000 sampled haplotype pairs;
- `chrN.recent.bits[.meta]`: one deterministic recent/not-recent bit per output
  position and manifest pair;
- `chrN.samples.txt`: BCF-ordered, ancestry/QC/relatedness-filtered sample IDs;
- `chrN.samples.audit.json`: columns, source hashes, and all filter counts;
- `chrN.candidate_regions.tsv`: windows with `frac_recent_4500 > 0.02` by
  default (configurable with `--signal-fraction`), merged
  when the intervening gap is no more than 20 kb;
- `chrN.candidates/`: candidate TMRCA profiles, all-variant score tables,
  representative variants, PNG/PDF plots, and a deterministic artifact
  manifest;
- `chrN.decode.log`: wall-time and peak-memory log;
- `chrN.complete.json`: input fingerprints, provenance commit, exact settings,
  and SHA-256 hashes for required outputs, sample list/audit, pair manifest, and
  callable mask.

Completed chromosome trees are uploaded with checksum-based `gcloud storage
rsync`, so reruns scan but do not recopy unchanged objects. Candidate artifacts
retain their chromosome-qualified `chrN.candidates/` directory. The runner does
not delete unmatched remote objects, and it synchronizes `chrN.complete.json`
only after the chromosome payload succeeds. Plot and combined-report manifests
use the same commit-marker-last rule.

Decode cache compatibility deliberately excludes the raw Git commit. It is
based on the semantic input/settings contract and is then verified against the
exact SHA-256 hashes of the BCF-ordered sample list and pair manifest plus every
required output. Consequently, a plot-only commit reuses legacy chromosome
completions, while a changed cohort, pair draw, input, mask, or decoder setting
does not. The runner also performs one locked uv synchronization up front and
uses `uv run --no-sync` for its per-chromosome subcommands.

Each population gets chromosome PNG/PDF scans, a chromosome summary table, and
a plot manifest. `-chr all` additionally requires all autosomes to validate
before it writes two whole-genome views: the main Manhattan-style PNG/PDF plots
the percentage of pairs whose posterior mean TMRCA is below 4,500 years, with
alternating chromosome colors and the candidate screen marked; a separate
diagnostic PNG/PDF retains the hard-call fraction, mean posterior probability,
and mean posterior TMRCA panels. These are descriptive scans, not
simulation-calibrated p-values.

After the population plots, the runner writes a scope-level report under
`results/summary/{scope}` locally and in GCS. `regions_by_population.tsv` gives
the number of merged candidate regions, chromosomes with regions, signal
windows, covered bases, and strongest peak for each requested population.
`all_candidate_regions.tsv` is recomputed directly from the chromosome scan
summaries at the requested `--signal-fraction`, and
`combined.whole_genome.gamma_smc.{png,pdf}` stacks the population scans on one
shared chromosome axis. Partial chromosome scopes receive the corresponding
`combined.requested_chromosomes` figure. A second
`combined.{scope}.gamma_smc.zoom4pct.{png,pdf}` view fixes the y axis at 4%;
downward triangles mark values clipped at that ceiling. It labels every
candidate-region cluster strictly above 2%. Curated labels are used when a
coordinate-matched override is available; otherwise the nearest protein-coding
gene is shown.

The report keeps the analytical and presentation layers separate:

- `raw_scan_windows.tsv.gz` contains every decoded output position for every
  requested population and chromosome, including its chromosome and cumulative
  genome coordinates. No locus merging is applied.
- `candidate_loci.tsv` contains windows strictly above `--signal-fraction`,
  merged as scan intervals only when the intervening gap is at most
  `--merge-gap` (20 kb by default). This is the candidate-locus definition.
- `plot_loci.tsv` starts from every row in `candidate_loci.tsv` and connects
  adjacent candidate intervals separated by at most `--plot-merge-gap` (1 Mb by
  default). This coarser merging only reduces duplicate labels in the detailed
  whole-genome figure; it never filters or redefines the 20-kb candidate loci.
  `source_region_ids` preserves the exact candidate rows represented by each
  label.
- `ranked_top_windows.tsv` retains the older top-`--top-n` window ranking as a
  diagnostic table. It is not used to decide which candidate loci receive
  labels.
- `report_data_layers.tsv` records the source file, selection rule, merge rule,
  maximum gap, and purpose of all four layers. `all_candidate_regions.tsv` and
  `gene_list.tsv` remain compatibility aliases for the candidate and complete
  labeled plot-locus outputs, respectively.

`gene_list.tsv` has one row per merged candidate-label locus, its contributing
candidate region IDs, peak and GRCh38 coordinates, the nearest protein-coding
gene, and every protein-coding gene overlapping the hit or lying within 500 kb.
The bundled `resources/gamma_smc_2pct_gene_label_overrides.tsv` supplies
literature-curated labels for the current 2% results and records evidence level,
rationale, and a primary reference. Unmatched loci retain the positional nearest
gene and are explicitly marked `positional_only`; neither class is a calibrated
selection call or proof of causality. `--gene-label-overrides` selects a
different override TSV. `--top-n` now controls only
`ranked_top_windows.tsv`; `--plot-merge-gap`, `--gene-context-flank`, and
`--gene-annotation-uri` control the labeled report. The command also prints the
region count for every population and the total.

The candidate pass is population-specific. It plots all decoded-pair TMRCA
quantiles within 500 kb of each peak, queries every BCF record within 100 kb,
and ranks each ALT allele by the fraction of raw TMRCA variance explained by
the ref/ref versus matching-ALT/matching-ALT pair label. These are allele pairs,
not diploid genotypes: for example, haplotype 1 from one person and haplotype 1
from another person can form either class. Mixed pairs and pairs carrying two
different ALT alleles are omitted for that ALT-specific comparison. At the
10 kb scan resolution, each variant is tested against the nearest decoded
TMRCA position; that position is written explicitly in every score row.

Only candidate chromosomes receive the targeted raw-posterior replay. It uses
the literal `chrN.pairs.tsv` as input, verifies the native pair array exactly,
and writes `candidate_tmrca.f32.zst` in position-major order where column `k`
is data row `k` of the pair manifest. The large native alpha/beta stream is
deleted after this validated conversion. No additional samples or pairs are
drawn for candidate regions.

## Analysis defaults

| Setting | Default |
|---|---:|
| decoder threads | 12 |
| posterior call rule | `mean` |
| output stride | 10,000 bp |
| mutation rate | `1.29e-8` |
| transition cache | 1,000 bp |
| scaled mutation rate (`theta`) | `0.00075` |
| recombination/theta ratio | `0.8` |
| recent threshold | 4,500 years |
| generation time | 25 years |
| pair mode | 100,000 distinct unordered haplotype pairs per population |
| pair seed | `1729` |
| signal screen | `frac_recent_4500 > 0.02` (configurable) |
| signal merge gap | 20,000 bp |
| candidate profile | peak +/-500,000 bp |
| variant search | peak +/-100,000 bp |
| minimum class size | 20 ref/ref and 20 matching-alt/matching-alt pairs |
| plot-label loci | all >2% candidate intervals; adjacent candidates with <=1 Mb gap connected |
| ranked-window diagnostic | top 100 hard-call windows/population |
| candidate gene context | protein-coding genes in merged label locus +/-500,000 bp |
| detail plot y ceiling | 4% |
| gene-label threshold | strictly above 2% |

The 1 kb cache is retained because increasing it has linear cache-memory cost
without a demonstrated 10 kb-stride speed benefit. Its steady shared cache is
about 490 MB (plus about 122 MB while constructing it); a 10 kb cache would be
about 4.9 GB (plus about 1.2 GB during construction). Cache size and stride are
independent: the output grid is 10 kb while transition-cache segments stay at
1 kb.

## Overrides and restart behavior

Every path can be changed without editing the script:

```bash
bash scripts/run_aou_workbench.sh -chr 1 -pops AFR \
  --output-prefix gs://my-bucket/gamma_results \
  --bcf-template 'gs://my-bucket/full_panel/chr{chr}.bcf' \
  --index-template '{bcf}.csi' \
  --mask-template 'gs://my-bucket/hardmask.excluded.bed' \
  --ancestry-uri gs://my-bucket/ancestry_preds.tsv \
  --qc-exclusions-uri gs://my-bucket/flagged_samples.tsv \
  --relatedness-exclusions-uri gs://my-bucket/relatedness_flagged_samples.tsv
```

Mask sensitivity examples:

```bash
# Existing/default hardmask only
bash scripts/run_aou_workbench.sh -chr all -pops all

# HMMIX strict callable mask only, in *_strict_hardmask roots
bash scripts/run_aou_workbench.sh -chr all -pops all --strict-hardmask

# Default followed by strict; the completed default run is reused
bash scripts/run_aou_workbench.sh -chr all -pops all --both-hardmask
```

Equivalent `AOU_GAMMA_*` environment variables are documented by
`scripts/run_aou_workbench.sh --help`. Use `--no-mask` only when that is the
intended callable-region policy, and apply the identical policy to the matched
null simulations.

The random draw is over all distinct unordered pairs among the filtered
population haplotypes. A pair can span two people, use either phased haplotype
from either person, or occasionally contain the two haplotypes of one person.
Use `--exclude-within` only if that last category should be removed; it is off
by default.

A chromosome/population is skipped only when the BCF, CSI, ancestry,
QC-exclusion, relatedness-exclusion, and hard-mask GCS fingerprints; sample and
pair manifests; decoder/code settings; summary and bit-matrix structure;
candidate artifact hashes; and output hashes match its completion record. The
completion object is uploaded last. This makes a rerun resume safely after an
interrupted decode, candidate analysis, or upload. Staged chromosome BCFs are
removed only after every requested population for that chromosome succeeds;
`--keep-inputs` retains them, `--force` recomputes, and `--no-upload` keeps the
run local.

Plotting and run-level reporting are independently idempotent. Their manifests
record every input hash, plotting setting, output size, and output hash and are
written atomically after all artifacts. An unchanged rerun reuses the local
plots/report, compares completed files before GCS upload, and uploads the
population/report manifest last. A nonblocking lock on the local root prevents
two notebook cells from writing the same cache and result tree concurrently.
