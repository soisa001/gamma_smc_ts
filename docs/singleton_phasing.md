# Transfer physical singleton phase onto an unchanged SV panel

The implementation looks up the **nearest informative left and right SNVs** by
sample, chromosome, and HiPhase PS.
It does not run SHAPEIT or statistically infer any singleton phase. `phase_rare`
does not retain the target's physical GT/PS linkage, which is the information
needed here. The source audit and pricing references are in
[singleton_phasing_research.md](singleton_phasing_research.md).

## How the lookup works

1. Match panel and HiPhase biallelic SNVs by exact CHROM/POS/REF/ALT. An informative
   anchor must be heterozygous, phased, genotype-concordant, and PASS/unfiltered in
   both files for the singleton carrier. SVs and indels are not used as anchors;
   every original SV/indel remains untouched in the output.
2. A shared `0|1` in the panel and `1|0` in HiPhase means that block is reversed.
   Reverse the singleton's GT when placing it on the panel. For example, a
   HiPhase singleton `1|0` becomes panel `0|1` in that reversed block.
3. Select the nearest qualifying SNV strictly left and strictly right of each
   singleton, within its carrier's original HiPhase PS block. Search outward past
   homozygous and otherwise uninformative markers with no fixed distance cutoff.
   Never borrow an anchor from another sample/PS block. By default both flanks
   must exist and imply the same orientation. Two markers on the same side do
   not substitute for bracketing. If the selected informative flanks conflict,
   reject the singleton; never search past the conflict to find agreement.
   Distant disagreements outside the selected pair do not veto local agreement.
   Tied nearest SNVs at the same coordinate must all agree. An anchor exactly at
   the singleton coordinate counts as neither left nor right.
   `--min-anchors 1` explicitly permits a single-sided result when the other side
   is absent, but still rejects disagreement when both sides exist. The default
   of two is the requested policy, not a calibrated accuracy guarantee.
4. Recompute MAC from the exact cohort's GTs. Add only complete diploid,
   biallelic SNVs with MAC=1 absent from the panel. This includes reference-minor
   singletons (AC=AN-1). Only PASS/unfiltered new SNVs qualify. Missing calls,
   haploid calls, multiallelic sites, non-singletons, and new SVs are counted and
   excluded. No missing genotype is converted to homozygous reference.
5. Exclude unsupported/conflicting singleton candidates from the augmented panel.
   Record their coordinates, carrier, original PS, selected anchor positions,
   orientation counts, and exclusion reason in `excluded_singletons.tsv.gz`.
   There are no newly added unphased records. No statistical fallback overrides
   missing physical evidence, and original panel records are never excluded.
6. Merge added records around the original panel records and independently read
   back the result. Every original record, including SVs, must match the original
   VCF text exactly and remain in its original order. Original header lines and
   sample order must also be retained. The whole compressed file necessarily
   changes when records are added; CRLF/LF line endings are normalized.

The input panel must contain chromosome-wide oriented haplotypes from statistical
phasing/ligation. Independent local PS blocks alone are insufficient. The manifest
requires an explicit `scaffold_phase_scope: chromosome` assertion. Global phase
errors in the existing panel are not corrected by this procedure.

## Inputs and operating limits

- Both inputs must describe the same joint cohort, with exactly the same sample
  IDs/order, genome build, contig names, and variant representations. The manifest
  records the build; the program does not independently prove reference-build
  identity or liftover alleles. It rejects incompatible contig lengths/header
  types, malformed rows, unsorted records, duplicate exact variant keys, unphased
  panel heterozygotes and empty shards. With no shared SNV anchors it logs a
  warning, excludes all candidates, and preserves the original panel.
- Per-chromosome `.vcf`/ordinary `.vcf.gz` inputs work with Python alone. Whole-
  genome inputs require indexed BGZF VCF or indexed BCF plus `bcftools` on PATH.
  Plain multi-contig VCF input is rejected rather than rescanned for every job.
  The BCF/indexed path uses `bcftools view --no-version -Ov -r CHROM
  --regions-overlap 0` and checks the process exit status.
- This does **not** merge 12,000 per-sample HiPhase VCFs or infer reference calls
  from absence in those files. If the 3 TB input is per-sample, a separate joint
  genotyping/merge policy and singleton-discovery stage are needed first.
- Default: four parallel chromosome workers, declared total budget 26 GB, seed
  1729 (the lookup itself is deterministic and uses no RNG). `--memory-gb` records
  the intended allocation; it is not an OS memory cap. Enforce memory in Slurm or
  VM provisioning. A one-chromosome pilot uses one active worker. Use at most 24
  workers locally; a larger authorized VM may use up to 32, bounded by job count.
- Genotypes are streamed by coordinate; informative anchor positions are stored
  in SQLite, with a bounded in-memory update batch. Each singleton uses two
  indexed predecessor/successor lookups on `(sample, PS, position)`. There is no
  full genome-by-sample matrix. The anchor database now retains every informative
  position and can be substantially larger than the former block-count database;
  include it in pilot disk sizing.

## Run on a local VM or HPC

Clone/pull and test synthetic data first. This cell accesses only GitHub code and
synthetic fixtures; it does not stage or access controlled data.

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
cd "$REPO"
[[ "$(git branch --show-current)" == AOU_run_opt ]]
git pull --ff-only origin AOU_run_opt
uv run --no-project --python 3.12 --with pytest \
  python -m pytest tests/test_phase_singletons.py -q
```

Use the following only on the machine where approved inputs already reside.
Replace the local paths and build/contig names. The `{contig}` path placeholder
can be omitted when both inputs are indexed whole-genome cohort files. The script
reads sample IDs from headers and writes them to a local manifest without printing
them. On Workbench, keep the manifest and all outputs within the controlled
environment; this cell contains no cloud data access or Git push.

```bash
%%bash
set -euo pipefail
cd "$HOME/gamma_smc_ts"
uv run --no-project --python 3.12 scripts/make_singleton_manifest.py \
  --scaffold '/approved/local/panel/{contig}.bcf' \
  --hiphase '/approved/local/hiphase/{contig}.vcf.gz' \
  --contigs chr1 \
  --genome-build GRCh38 --assert-chromosome-phase \
  --output /approved/local/singletons/manifest.json
uv run --no-project --python 3.12 scripts/phase_singletons.py \
  --manifest /approved/local/singletons/manifest.json \
  --out /approved/local/singletons/results \
  --workers 4 --memory-gb 26 --min-anchors 2
```

Start with a representative chromosome or pre-extracted representative contig
shard with all 12,000 samples. For a regional pilot, extract **both** inputs with
buffered complete physical blocks; cutting off anchors at interval edges lowers
the apparent transfer rate. Run selected full chromosomes when possible. After
the pilot, expand `--contigs` and reuse the same output directory. Adding unrelated
jobs does not require rescanning the HiPhase file for already prepared contigs.
Do not split a physical block across parallel arbitrary windows.

For Slurm, request `--cpus-per-task=4 --mem=26G` around the same command. If data
staging is needed in a controlled environment, the separate approved staging
operation should use `gcloud storage rsync SOURCE DEST --billing-project
"$GOOGLE_PROJECT"`; neither the implementation nor the run cell accesses buckets.

## Outputs, evidence, and restart behavior

Each contig directory contains:

| File | Meaning |
|---|---|
| `augmented.vcf.gz` | Entire original panel plus supported, phased singleton records |
| `added_singletons.vcf.gz` | Only the supported phased singleton additions |
| `excluded_singletons.tsv.gz` | Rejected candidate coordinates/carriers and reasons, including selected flank positions |
| `anchors.sqlite` | Same/flipped evidence by sample index, PS, and informative SNV position |
| `audit.json` | Counts, record/header hashes, paths, parameters, input identities, elapsed time, process peak memory, preservation check |
| `workflow.log` | Start/completion/failure messages, decode commands, parameters, progress |
| `prepare.done.json`, `done.json` | Atomic completion records with output SHA-256 checksums |

Added records carry `SP_METHOD`, `SP_REASON`, `SP_ANCHORS`, `SP_LEFT_POS`,
`SP_RIGHT_POS`, `SP_SOURCE_PS`, and `SP_SOURCE_GT`. The last two retain the carrier's
original HiPhase evidence. `SP_ANCHORS` counts only observations at the selected
nearest flank coordinates. Exclusion reasons distinguish no physical PS, no
informative anchor, missing left/right anchor, and conflicting nearest anchors.
Old FORMAT/PS and FORMAT/PQ are removed from **added** records because they describe
the old physical coordinate system. Other new-record fields are retained. No
field of an original panel record is changed. `anchors.sqlite` sample indices
refer to the ordered manifest. These artifacts may contain controlled data and
must stay with the input data.

Outputs use ordinary gzip, not BGZF, to keep the core implementation portable.
They are sequentially readable but **not tabix-indexable until reblocked**. On a
Linux VM with `bgzip` and `tabix`, make an indexed copy as a separate phase:

```bash
%%bash
set -euo pipefail
SRC=/approved/local/singletons/results/chr1/augmented.vcf.gz
DST=/approved/local/singletons/results/chr1/augmented.bgz.vcf.gz
# Only create/recreate the indexed copy if it is absent or older than its source.
if [[ ! -s "$DST" || "$SRC" -nt "$DST" ]]; then
  gzip -dc "$SRC" | bgzip -@ 4 -c > "$DST.partial"
  bgzip -t "$DST.partial"
  mv "$DST.partial" "$DST"
fi
tabix -f -p vcf "$DST"
```

Cache decisions use manifest/sample identity, input paths/sizes/mtime, semantic
algorithm version and relevant parameters, plus SHA-256 verification of generated
artifacts. A new repository commit or worker-count change is not a cache key.
Preparation is reusable across downstream failures and unrelated manifest edits.
Algorithm version 2 invalidates version 1's block-consensus caches, since both the
anchor selection and exclusion policy changed.
Inputs are hashed logically during their normal scan, and source identity is
checked again at completion. A same-size input modification with deliberately
preserved mtime cannot be detected on a cached run; use immutable staged inputs or
a new output directory in that situation. Gzip CRC failures and decoder nonzero
exits fail the run. Corrupt generated VCFs are detected by checksums and rebuilt.
A stale `.lock` after a killed process needs removal only after confirming the job
is no longer running. Failed runs never produce a new successful completion marker.

## Local validation, 2026-10-27

The focused suite passes 29 tests covering phase flips, same orientation,
reference-minor singletons, exclusion of SV anchors, expansion past homozygotes,
nearest-flank agreement/conflict, distant disagreements, tied-coordinate conflicts,
missing left/right flanks, sample/PS isolation, unsupported calls, duplicate/sort errors,
sample/header mismatches, truncated gzip, original-record preservation, corrupt-
output repair, CLI/manifest generation, cache reuse, stale completion invalidation,
and randomized orientations.

The synthetic cohort-width benchmark uses 12,000 samples, four chromosome jobs,
160 shared markers bracketing 400 new singleton sites. It checks that all 320
supported phases match truth, all 80 deliberately unphased targets are excluded,
and all 160 original records survive exactly. Updated measured runtime,
per-worker peak memory, and counts are recorded in
[singleton_transfer_benchmark.json](singleton_transfer_benchmark.json).

This is a short parser/correctness benchmark with artificial GT:PS fields, not a
biological accuracy estimate or full-genome timing extrapolation. No user callset
was available locally. Plain/gzip paths were tested on Windows; the external
BCF/indexed-VCF adapter needs a Linux pilot with real `bcftools` (local WSL startup
failed). SHAPEIT is neither required nor invoked.

## Cost for 12,000 samples

Physical transfer removes the statistical-phasing compute stage. Reading and
parsing the 3 TB HiPhase file, scanning the 50 GB panel, writing new cohort-wide
records, and verifying output now dominate. The implementation makes **one full
HiPhase pass**; it rereads only the extracted singleton spool and the much smaller
panel/output for assignment and verification. Each contig is processed separately.

**A runtime forecast cannot be obtained from 3 TB compressed input alone.** Site
count, genotype-field width, compression ratio, shared heterozygosity, anchor count,
and added-output size matter. Synthetic uncompressed end-to-end throughput
cannot be applied to the 3 TB compressed size. A real all-sample pilot is needed
before promising a cost ceiling.

For a provisional budget, GCP Iowa on-demand `n2-highmem-4` (4 vCPU/32 GiB) is
$0.262028/hour. Balanced disk is $0.000136986/GiB-hour. This gives the following
**conditional runtime scenarios**, excluding existing bucket storage, transfer,
platform fees, and taxes:

| Actual elapsed duration | VM + 4 TiB disk | VM + 8 TiB disk |
|---|---:|---:|
| 1 day | $19.75 | $33.22 |
| 3 days | $59.26 | $99.66 |
| 7 days | $138.28 | $232.55 |
| 30 days | $592.65 | $996.64 |

Thus **$100–$235 covers a 3–7 day run with 8 TiB scratch**, if the pilot supports
that duration; it is not evidence the run finishes within a week. A month on that
allocation is about $1,000. Larger VM scenarios and regional-transfer prices are
in the research report. Disk keeps accruing costs while retained, including when
the VM is stopped.

Four TiB is only a price example and may be too small. Peak disk includes staged
inputs plus the complete raw singleton-candidate spool during preparation,
accepted singleton output, the augmented panel, SQLite evidence, excluded-candidate
audit, and any later BGZF copy. A conservative first estimate is
`input_HiPhase + 2*panel + candidate_spool + 2*added_singleton_size + SQLite + audit + headroom`; add another
`panel + added_singleton_size` while retaining a BGZF copy. Do not provision based
only on the input's 3 TB size.

For the pilot, record sample count, shared variants, singleton candidates,
compressed input/output bytes, `prepare_seconds`, total elapsed time, and peak
memory. Extrapolate using **site counts and measured genotype parsing cost** as
well as byte throughput, not genomic span alone. Use chromosomes with varied
variant density and PS fragmentation; retained output size controls disk sizing.
Keep a contingency for I/O contention and uneven chromosome workloads.
