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

Run this cell on the Workbench VM that has bucket access and `GOOGLE_PROJECT`
configured. It pulls code, installs missing public tools, stages chr1, creates a
local sample manifest, then runs the nearest-flanking-SNV transfer. It does not
upload results or push anything to GitHub.

```bash
%%bash
set -euo pipefail
REPO=/home/jupyter/gamma_smc_ts
ROOT=/home/jupyter/singleton_transfer_chr1
: "${GOOGLE_PROJECT:?GOOGLE_PROJECT must identify the billing project}"

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
  --chrom chr1 --root "$ROOT" --phase all \
  --billing-project "$GOOGLE_PROJECT" \
  --workers 4 --io-threads 4 --memory-gb 26 --min-anchors 2
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
  --chrom chr1 --root /home/jupyter/singleton_transfer_chr1 --dry-run
```

## Staging and restarts

`--phase stage` performs only downloads/index checks. `--phase prepare` uses
validated local files to create the ordered cohort manifest. `--phase run` uses
validated local files, recreates the small manifest, and runs phasing. `all`
performs the complete sequence. The `prepare` and `run` modes make no GCS calls.

Staging uses `gcloud storage rsync SOURCE_DIRECTORY LOCAL_DIRECTORY
--billing-project "$GOOGLE_PROJECT"` with an exact filename exclusion filter.
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
/home/jupyter/singleton_transfer_chr1/results/chr1/augmented.vcf.gz
```

The same directory contains `added_singletons.vcf.gz`,
`excluded_singletons.tsv.gz`, `anchors.sqlite`, `audit.json`, and workflow logs.
The root contains `manifest.json`, `staging.json`, `runner.log`, and the staged
inputs in `inputs/hiphase` and `inputs/scaffold`. All are local to the VM. No upload
command is implemented. The output is ordinary gzip VCF; see the
[core workflow guide](singleton_phasing.md) for a later BGZF/indexing step.

## What workers mean

| Setting | Actual behavior |
|---|---|
| `--workers 4` | Up to four **separate chromosome processes** in the underlying multi-chromosome workflow |
| This chr1 runner | Exactly one chromosome job, so only **one phasing worker** is active |
| Within chr1 | Python parsing, nearest-anchor lookup, and output generation are serial; there is no native multithreaded phasing kernel |
| `--io-threads 4` | Up to four gcloud transfer threads in one process; also used to create a missing CSI |
| `--memory-gb 26` | Recorded total application allocation, not an OS-enforced memory limit |

BCF decoding runs as separate bcftools subprocesses feeding the Python worker;
that is not four-way parallel phasing. Increasing `--workers` will not accelerate
chr1. The two input files are staged sequentially, with parallel transfer threads
within each rsync invocation. gcloud concurrency overrides are per subprocess,
not persistent changes to the user's gcloud configuration. [gcloud concurrency
properties](https://docs.cloud.google.com/sdk/gcloud/reference/config/set).

Later, the core manifest can contain multiple chromosomes, making `--workers 4`
useful. Parallelizing a single chromosome would require a separate implementation
that preserves access to flanking anchors across partition boundaries; this runner
does not split chr1 into arbitrary windows.

## Validation status

The local focused suite contains 29 phasing tests and 10 staging-control tests.
The staging tests mock cloud and bcftools calls and cover exact chr1 selection,
download direction, CRC mismatch detection, restart reuse/repair, missing indexes,
permission errors, changing source generations, disk preflight, and local-only
prepare/run dispatch. These tests do not claim a real cloud download succeeded.

The development host has no gcloud or native bcftools, and WSL startup fails.
Consequently, actual bucket access and the real chr1 analysis must be executed in
Workbench using the cell above. No controlled callset has been downloaded or
analyzed on the development host.
