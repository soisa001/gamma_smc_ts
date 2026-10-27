# Chr1 HiPhase singleton transfer in Workbench

The runner selects **chr1 only** by default, stages the following inputs and
their `.csi` indexes, and writes the augmented panel to local disk:

| Input | URI |
|---|---|
| HiPhase cohort | `gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/phasing/hiphase/aou_lr_phase2_v1.hiphase.short_sv.chr1.bcf` |
| Original phased panel | `gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf` |

The supplied paths are wired verbatim. Actual object availability, sizes, sample
ordering, and phase-set contents are checked when running in the authorized VM.
The default build label is GRCh38, consistent with the existing AoU v9 workflow;
override `--genome-build` if appropriate. No coordinates are lifted or transformed.

## Setup and run

Run this cell on the Workbench VM with existing gcloud bucket access. No
`GOOGLE_PROJECT` variable or account/billing override is required by the runner.
It pulls code, installs missing public tools, stages chr1, creates a
local sample manifest, then runs the nearest-flanking-SNV transfer. It does not
upload results or push anything to GitHub.

```bash
%%bash
set -euo pipefail
REPO=/home/jupyter/gamma_smc_ts
ROOT="$HOME/singleton_phased/working"
VCFS="$HOME/singleton_phased/vcfs"

if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt \
    git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
cd "$REPO"
[[ "$(git branch --show-current)" == AOU_run_opt ]]
git pull --ff-only origin AOU_run_opt

bash scripts/bootstrap_singleton_tools.sh
UV="$(command -v uv || true)"
UV="${UV:-$REPO/.tools/uv-bin/uv}"
"$UV" run --no-project --python 3.12 scripts/run_singleton_chr.py \
  --chrom chr1 --root "$ROOT" --vcf-dir "$VCFS" --phase all \
  --workers 4 --chunks-per-chrom 4 --overlap-bp 1000000 \
  --io-threads 4 --memory-gb 26 --min-anchors 2
```

The lightweight bootstrap reuses an available `bcftools`, including the existing
repository `.native/bin/bcftools`. Otherwise it installs bcftools/htslib 1.22 into
`.tools/singleton-native` using a private micromamba environment on Linux x86_64.
It does not build Gamma-SMC or install unrelated simulation tools. Python and the
pinned CRC32C library are managed with `uv`. The runner also accepts
`--bcftools /path/to/bcftools` or `BCFTOOLS_BIN`.

To inspect the resolved paths/commands without contacting GCS or writing files:

```bash
%%bash
set -euo pipefail
cd /home/jupyter/gamma_smc_ts
UV="$(command -v uv || true)"
UV="${UV:-$PWD/.tools/uv-bin/uv}"
"$UV" run --no-project --python 3.12 scripts/run_singleton_chr.py \
  --chrom chr1 --root "$HOME/singleton_phased/working" \
  --vcf-dir "$HOME/singleton_phased/vcfs" --dry-run
```

## Staging and restarts

`--phase stage` performs only downloads/index checks. `--phase prepare` uses
validated local files to create the ordered cohort manifest. `--phase run` uses
validated local files, recreates the small manifest, and runs phasing. `all`
performs the complete sequence. The `prepare` and `run` modes make no GCS calls.

Staging uses `gcloud storage rsync SOURCE_DIRECTORY LOCAL_DIRECTORY` with an exact
filename exclusion filter. It inherits gcloud's existing account and configuration.
An explicit `--billing-project PROJECT_ID` is optional and, when supplied, is
forwarded to rsync and object-metadata requests. The runner does not infer a
billing override from `GOOGLE_PROJECT`. A Requester Pays bucket may still require
a billing project through gcloud configuration or that optional flag; the storage
service enforces that requirement. [Requester Pays documentation](https://docs.cloud.google.com/storage/docs/requester-pays).
The chr1 filter excludes chr10 and all other chromosomes. It never deletes
unmatched destination files. Source size, generation, and CRC32C are read before
staging and checked again afterward. Every newly staged BCF/CSI is verified
against its cloud CRC32C, and index readability is checked with bcftools. A
missing CSI (404) is created locally; permission errors are not treated as a
missing index. [gcloud rsync reference](https://docs.cloud.google.com/sdk/gcloud/reference/storage/rsync),
[object metadata reference](https://docs.cloud.google.com/sdk/gcloud/reference/storage/objects/describe).

A `stage.done.json` manifest records remote generations/checksums and verified
local size/mtime. Matching manifests skip restaging, even after a repository
commit changes. Locally modified files or changed remote generations trigger
checksum-based rsync and verification. The cached fast path assumes staged files
are immutable; a same-size corruption with deliberately preserved mtime is not
detectable without another checksum pass. Remove the relevant `stage.done.json`
and rerun `--phase stage` to force that pass. Never remove an active runner lock.

The disk preflight uses the actual source object sizes plus a **64 GiB minimum
scratch reserve**, adjustable with `--scratch-reserve-gb`. That reserve is a floor,
not a prediction that the full chromosome run will fit. The singleton spool,
SQLite anchor index, audit, and augmented output may need substantially more
space; retain headroom. Staged inputs are retained for restart and future analysis.

## Local outputs

With the root in the cell above, the requested output is:

```text
~/singleton_phased/vcfs/chr1.singleton_phased.vcf.gz
```

Both directories are created as needed. `--root` defaults to
`~/singleton_phased/working`; `--vcf-dir` defaults to `~/singleton_phased/vcfs`.
The runner checks the completed VCF checksum, then creates the final named file
atomically. It uses a hard link on the same filesystem to avoid storing a second
large VCF, or a checksum-verified copy when hard links are unavailable. Treat
generated VCFs as immutable. A matching existing final file is reused on restart.

The working directory `~/singleton_phased/working/results/chr1` contains the
internal `augmented.vcf.gz`, compact `added_singletons.jsonl.gz`,
`excluded_singletons.tsv.gz`, `chunks.json`, `audit.json`, and workflow logs.
Each `chunk_000` through `chunk_003` directory contains its own `anchors.sqlite`,
compact candidate/addition files, extraction/phasing completion markers and log.
Only the final merge expands a sparse carrier into a cohort-wide row. Added
records use FORMAT/GT only: noncarriers are `0|0`, or `1|1` for reference-minor
singletons; the carrier gets the transferred phase. Carrier source GT/PS remain
in INFO. Original panel records retain every field, including FORMAT annotations.
The working root contains `manifest.json`, `staging.json`, `runner.log`, and the staged
inputs in `inputs/hiphase` and `inputs/scaffold`. All are local to the VM. No upload
command is implemented. The output is ordinary gzip VCF; see the
[core workflow guide](singleton_phasing.md) for a later BGZF/indexing step.

## What workers mean

| Setting | Actual behavior |
|---|---|
| `--workers 4 --chunks-per-chrom 4` | Four **process workers within chr1**, each owning approximately one quarter of the coordinate span |
| `--overlap-bp 1000000` | Read a 1 Mb anchor halo on each side; only the disjoint core owns singleton output |
| Extraction | Each worker streams its indexed panel/HiPhase regions into its own SQLite anchor index and sparse candidate file |
| Phase assignment | After all indexes are complete, workers independently resolve their candidates; missing local flanks query other chunks' indexes |
| Final insertion | One streaming merge writes the augmented panel, followed by an independent preservation audit |
| `--chunks-per-chrom 1` | Original whole-chromosome implementation; only one active worker for chr1 and full-format addition intermediates |
| `--io-threads 4` | Up to four gcloud transfer threads in one process; also used to create a missing CSI |
| `--memory-gb 26` | Recorded total application allocation, not an OS-enforced memory limit |

BCF decoding runs as separate bcftools subprocesses feeding each Python worker.
The phasing code uses multiprocessing, not a native multithreaded kernel. Four
workers do not promise fourfold speedup: disk bandwidth, VCF decoding, unequal
variant density and the serial merge/audit limit scaling. Two roughly equal
coordinate intervals can have different workloads. The two input files are staged sequentially, with parallel transfer threads
within each rsync invocation. gcloud concurrency overrides are per subprocess,
not persistent changes to the user's gcloud configuration. [gcloud concurrency
properties](https://docs.cloud.google.com/sdk/gcloud/reference/config/set).

The cross-chunk lookup uses the carrier's same HiPhase PS and the nearest
informative SNV in the required direction. It does not scan the input again or
stop at the configured overlap. Duplicate halo copies of an anchor are never
double-counted. This keeps the whole-chromosome decision rule intact even with
zero overlap. Conflicting informative flanks still exclude the candidate.

If the core manifest later contains multiple chromosomes with chunking enabled,
chromosomes run sequentially, each using the same bounded worker pool size. There
is no nested four-chromosome-by-four-chunk process multiplication. With chunking
disabled, the original parallelism across chromosomes remains available. The
underlying CLI defaults to one chunk for compatibility; this chr1 runner defaults
to four. Plain unindexed VCF is supported for tests but each worker must read it
sequentially; use the staged indexed BCFs for actual parallel region reads.

## Validation status

The local focused suite contains 29 original phasing tests, 9 chunk tests and 13
staging-control tests (51 total). Chunk tests compare serial and parallel
decisions, including exact boundaries, empty cores, zero/large halos, remote
flanks, conflicting/tied anchors, random phase blocks, reference-minor calls and
corrupt intermediate repair. A 12,000-sample synthetic single-chromosome run used
four workers, retained all 80 panel records, correctly phased 320 candidates and
excluded 80 deliberately unphased candidates in 8.84 seconds. All 320 accepted
candidates exercised cross-chunk lookup. See
[the benchmark record](singleton_chunk_benchmark.json); this is not a real-data
runtime forecast or a BCF throughput measurement.
The staging tests mock cloud and bcftools calls and cover exact chr1 selection,
download direction, CRC mismatch detection, restart reuse/repair, missing indexes,
permission errors, changing source generations, disk preflight, and local-only
prepare/run dispatch. These tests do not claim a real cloud download succeeded.

The development host has no gcloud or native bcftools, and WSL startup fails.
Consequently, actual bucket access and the real chr1 analysis must be executed in
Workbench using the cell above. No controlled callset has been downloaded or
analyzed on the development host.
