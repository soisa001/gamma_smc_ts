# chr1–chr21 native Truvari SV merges at 0.7

The final requested output is a **cohort VCF with native Truvari genotype consolidation**, one per chromosome at similarity 0.7. This supersedes the earlier sites-only and two-threshold requests. Four chromosome jobs can run concurrently. Settings are `--keep common --gt all -r 500 -s 50 -S -1`, with both `-p` and `-P` set to 0.7. `--hap`, type mixing and chaining are off.

## Native Truvari behavior

This workflow now calls unmodified `truvari collapse` 5.5.0, as requested. The custom genotype-transfer wrapper has been removed. No installed package code is edited and no GT, phase, or frequency correction is applied after merging.

Sample genotypes and phase fields are supplied to Truvari. Their consolidation and representative selection follow native behavior, including the previously observed 5.5.0 limitations: a donor ALT can remain reference in the representative under `--gt all`, and `--keep common` selected the lower-AC member in our synthetic test. The workflow does not promise complete donor-phase preservation. `--gt all` prevents merging alleles carried by the same sample; `--hap` remains off.

Source: [pinned native implementation](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py).

## Workbench cell

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
  git -C "$REPO" switch AOU_run_opt
  git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
cd "$REPO"
export SV_PHASED_ROOT="$HOME/sv_merge_chr1_21_phased"
export SV_CHROM_JOBS=4
export SV_MEMORY_GB=26
bash scripts/run_sv_phased_batch_workbench.sh all
```

No cloud data or outputs are uploaded, and no plots are displayed. The launcher uses uv and runs synthetic tests before the data merge. Four jobs is the chromosome concurrency limit; each job runs only one Truvari process at a time. The 26 GB total is divided among jobs as a working budget, not an OS-enforced memory limit.

The full-panel template is the repository path already discussed with the user:

```text
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.{chrom}.bubble.split.bcf
```

Staging uses `gcloud storage rsync` restricted to the exact BCF/index for each chromosome and no billing-project flag. Transfers run one at a time under a Linux process lock at `${XDG_CACHE_HOME:-$HOME/.cache}/gamma_smc_ts/gcloud-rsync.lock`, because simultaneous rsync invocations can interfere with gcloud's shared source-list files. This lock covers this workflow's transfers across chromosome workers and batch roots; unrelated manual rsync commands do not honor it. Preparation and merging remain concurrent, with up to four chromosome workers total. Lock waits/acquisitions/releases are logged; the operating system releases the lock if a worker exits. Do not delete the lock file.

Downloaded temporary BCFs are removed only after the prepared input is verified; set `SV_KEEP_STAGED_BCF=1` to keep them. The prepared phased SV input is retained. To use existing local indexed files, set `SV_LOCAL_INPUT_TEMPLATE='/absolute/path/aou_lr_phase2_v1.{chrom}.bubble.split.bcf'`; these external input files are never removed. Failed workers are reported immediately while other chromosomes continue.

## Minimal preparation and direct outputs

Truvari accepts a biallelic, sorted BGZF VCF plus tabix index. An original VCF that already satisfies those requirements can go straight to the native command below. For the source panel BCF, the batch runner uses a streaming bcftools pipeline to convert formats, split multiallelic records with `--multi-overlaps 0`, and select structural alleles (literal net change ≥50 bp, plus symbolic/BND records). It does not impose pure-indel, A/C/G/T-only, PASS-only, SV-type, or upper-size filters. Shared REF/ALT padding is not rewritten. Phased `1|2` splits into `1|0` and `0|1`, with the existing phase-set field retained.

Truvari's `-o` file is already the merged result; `-c` contains removed cluster members. The batch writes both directly as BGZF, checks gzip integrity, indexes them, checks sample order and count conservation, and sorts only if indexing reveals that a tool output is not sorted. It does **not** build `new.unsorted.vcf` or reconstruct INS/DEL partitions.

```bash
# For an already prepared original SV VCF (run inside the pinned uv environment):
truvari collapse \
  -i original.sv.vcf.gz -o chr1.p70.vcf.gz -c chr1.p70.collapsed.vcf.gz \
  -p 0.7 -P 0.7 -r 500 -s 50 -S -1 --keep common --gt all
```

By default only 0.7 is run. Existing 0.9 files are left untouched and are not included in the refreshed batch counts unless 0.9 is explicitly selected.

Final outputs under `$SV_PHASED_ROOT`:

```text
chr1/prepare/input.vcf.gz          # original phased SV baseline retained for audit/resume
chr1/p70/chr1.vcf.gz               # final native merged VCF at 0.7
chr1/p70/chr1.collapsed.vcf.gz     # removed original members, with MatchId links
... through chr21 ...
counts.tsv                        # original, retained and collapsed counts per chr/threshold
```

Every VCF has an index. `CollapseId` and `MatchId` link the retained representative and removed records. `logs/chrN.worker.log`, `chrN/logs/batch.log`, and `chrN/p70/logs/collapse-native.log` contain tool output, paths, full parameters and native Truvari statistics. No summarized controlled data are pushed to GitHub.

Phases are `plan`, `prepare`, `merge`, or `all`. Repeating `all` reuses checksummed preparation and completed threshold outputs. The input snapshot is frozen and retains its source metadata/checksums; use a new output root for a new panel generation. The native-consolidation cache key differs from the former wrapper key, so old wrapped merge outputs are regenerated while verified prepared inputs are reused. Unrelated Git commits do not invalidate results. Use a different root if you want to retain both versions. If interrupted, verify the worker has stopped before removing its stale `.phased_batch.lock`.

## Aggregate INS/DEL comparison

### Chr1 threshold sensitivity comparison

Select `--chromosomes chr1 --similarities 0.7 0.9` on `run_sv_phased_batch.py`, or append these options after the phase on the Workbench launcher. For example, `bash scripts/run_sv_phased_batch_workbench.sh all --chromosomes chr1 --similarities 0.7 0.9 --jobs 1` runs the two thresholds sequentially, independently from the same prepared input. Both sequence and size thresholds change together; all other native Truvari settings stay fixed. Existing verified p70 output is reused. The default remains 0.7 for the whole-genome batch. A subset run writes `counts.chr1.tsv` at the batch root without replacing the full-batch `counts.tsv`.

Then run `uv run --no-project --python 3.12 scripts/compare_sv_native.py --root "$SV_PHASED_ROOT" --chromosomes chr1 --similarities 0.7 0.9 --bcftools "$BCFTOOLS_BIN" --threads 4`. Individual reports are written in `chr1/compare_native_p70` and `chr1/compare_native_p90`; the joined `chr1/compare_p70_p90/threshold_comparison.tsv` gives original, p70, p90, percentage reductions, and p90-minus-p70 values. Input preparation manifests must match before reports can be joined. Each individual report also contains inter-SV distance summaries. The cluster audit accepts `--similarity 0.9` for auditing p90 separately.

Let an existing batch finish before launching another merge against its root. Use a separate code checkout for the comparison if the existing checkout is in use. The source-data and merge cache keys are unchanged by adding a threshold; a new Git commit alone does not trigger rerunning completed p70 merges.

For completed native 0.7 outputs, run `uv run --no-project --python 3.12 scripts/compare_sv_native.py --root "$SV_PHASED_ROOT" --chromosomes chr6 chr8 --bcftools "$BCFTOOLS_BIN" --threads 4`. Omit `--chromosomes` to inspect all completed autosomes. Unfinished chromosomes are skipped. This reads site fields through bcftools, not sample FORMAT fields into Python, and does not rerun merging or touch its files. Up to three site scans run concurrently per chromosome.

Results are written under `chrN/compare_native_p70/`: `changes.tsv` (before/after, absolute and percentage reductions), `before_after.tsv`, `inter_sv_distance_summary.tsv`, and `summary.json`. Only aggregate changes are printed; no plots or individual variants are displayed. Metrics include INS/DEL record and exact unique-allele counts, raw length sums/quantiles, unique insertion-allele bp, and the union of simple deleted reference intervals. Complex replacements and symbolic deletions are excluded from the deletion union; symbolic alleles are excluded from raw length metrics. The baseline is the prepared SV input, not every small variant in the original panel.

Validation reads all three VCFs to EOF and checks exact site-multiset conservation (`before = retained + removed`), preparation counts, and unchanged input identities/manifests. It does not rehash the entire cohort VCFs. Results are cached against input identities and completion manifests. Existing genotype and representative behavior remains native Truvari.

## Cluster audit

For a read-only cluster/representation audit, run `uv run --no-project --python 3.12 scripts/audit_sv_native.py --batch-root "$SV_PHASED_ROOT" --chrom chr8 --bcftools "$BCFTOOLS_BIN"`. Aggregate JSON is printed and saved to `chr8/audit_native_p70/audit.json`; per-cluster positions and membership counts stay on disk in `clusters.tsv`. The audit checks site conservation and cluster annotations, summarizes literal-vs-complex representations, compares INFO type/size to raw REF/ALT lengths, and reports native match scores, cluster sizes, position spans, and identical-POS/REF matches with different ALT sequences. It does not validate genotypes or establish biological equivalence. Results are cached against input file identities and completion manifests.

## Tests

```bash
SV_PILOT_INTEGRATION=1 uv run --no-project --python 3.12 \
  --with truvari==5.5.0 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p test_sv_phased_batch.py -v
```

Tests compare every batch output record, including all sample FORMAT values, against independent unmodified Truvari CLI runs at 0.7 on two synthetic chromosomes. They also check that no 0.9 output is created, phased multiallelic decomposition, same-sample calls staying separate, record conservation, unchanged removed calls, retained symbolic/non-pure alleles, and cache reuse. A four-process staging test simulates gcloud's shared source listing through the actual staging call site and verifies transfers cannot overwrite one another's listing. No controlled panel data are used in local tests.
