# chr1–chr21 phased SV merges at 0.7 and 0.9

The final requested output is a **phased cohort VCF**, one per chromosome at each similarity threshold. This supersedes the earlier sites-only request. Four chromosome jobs can run concurrently; each runs 0.7 then 0.9 independently against the same original input. Settings are `--keep common --gt all -r 500 -s 50 -S -1`, with both `-p` and `-P` set to the requested threshold. `--hap`, type mixing and chaining are off.

## Truvari behavior verified before implementation

The pinned Truvari 5.5.0 code and real command-line tests show:

- `--hap` is for a **single individual**, requires `--keep first`, and merges a best matching pair across haplotypes into a homozygous call. It is incompatible with this cohort/common configuration.
- `--gt all` prevents clustering alleles carried by the same sample, regardless of whether they are `1|0` or `0|1`. It does not compare phase sets or establish haplotype equivalence.
- Native genotype-aware consolidation has a reproducible defect: an ALT present only in a removed record can remain `0|0` in the representative. The source calculates a new GT but does not assign it back in that case. A phased-looking native output is therefore insufficient evidence of correct merged genotypes.
- Native `--keep common` in this release also selected the **lower-AC** record in the synthetic test. Its AC comparator is ascending. We keep the requested native representative choice; the wrapper does not silently change sorting or clustering.

Sources: [collapse documentation](https://github.com/ACEnglish/truvari/wiki/collapse), [pinned implementation](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py). The observed 5.5.0 behavior takes precedence over descriptions that imply reliable genotype consolidation or highest-AC selection.

`scripts/truvari_collapse_phased.py` replaces **only** the genotype consolidation function for this invocation. Truvari's clustering, representative choice and `-o`/`-c` writers remain native. Under `--gt all`, a sample can carry at most one member of a cluster. That member supplies the output sample's exact GT tuple, phased flag and FORMAT fields, including PS when available. The implementation checks the copied GT, phased flag and phase fields, and fails if it encounters two ALT-bearing members in one sample. Unphased and partial calls remain unphased/partial; no phase is invented. Samples with no ALT donor retain the representative's original reference/missing call. AC/AN/AF are recomputed for consolidated clusters.

This preserves the phase information attached to the contributing original allele. It does not establish new phase relationships across distinct phase sets, validate biological equivalence of merged alleles, or make all reference/missing calls into jointly re-genotyped cluster calls. Original input and removed records remain available for audit.

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

Staging uses `gcloud storage rsync` restricted to the exact BCF/index for each chromosome and no billing-project flag. Up to four raw cohort BCFs can be staged concurrently. Downloaded temporary BCFs are removed only after the prepared input is verified; set `SV_KEEP_STAGED_BCF=1` to keep them. The prepared phased SV input is retained. To use existing local indexed files, set `SV_LOCAL_INPUT_TEMPLATE='/absolute/path/aou_lr_phase2_v1.{chrom}.bubble.split.bcf'`; these external input files are never removed.

## Minimal preparation and direct outputs

Truvari accepts a biallelic, sorted BGZF VCF plus tabix index. An original VCF that already satisfies those requirements can go straight to the phase-preserving entry point below. For the source panel BCF, the batch runner uses a streaming bcftools pipeline to convert formats, split multiallelic records with `--multi-overlaps 0`, and select structural alleles (literal net change ≥50 bp, plus symbolic/BND records). It does not impose pure-indel, A/C/G/T-only, PASS-only, SV-type, or upper-size filters. Shared REF/ALT padding is not rewritten. Phased `1|2` splits into `1|0` and `0|1`, with the existing phase-set field retained.

Truvari's `-o` file is already the merged result; `-c` contains removed cluster members. The batch writes both directly as BGZF, checks gzip integrity, indexes them, checks sample order and count conservation, and sorts only if indexing reveals that a tool output is not sorted. It does **not** build `new.unsorted.vcf` or reconstruct INS/DEL partitions.

```bash
# For an already prepared original SV VCF (run inside the pinned uv environment):
python scripts/truvari_collapse_phased.py \
  -i original.sv.vcf.gz -o chr1.p70.vcf.gz -c chr1.p70.collapsed.vcf.gz \
  -p 0.7 -P 0.7 -r 500 -s 50 -S -1 --keep common --gt all
```

Use 0.9 for both similarity arguments for the second independent run. This entry point uses the native interface but includes the explicit phase-consolidation fix described above; it is not an unmodified `truvari collapse` command.

Final outputs under `$SV_PHASED_ROOT`:

```text
chr1/prepare/input.vcf.gz          # original phased SV baseline retained for audit/resume
chr1/p70/chr1.vcf.gz               # final merged phased VCF at 0.7
chr1/p70/chr1.collapsed.vcf.gz     # removed original members, with MatchId links
chr1/p90/chr1.vcf.gz               # final merged phased VCF at 0.9
chr1/p90/chr1.collapsed.vcf.gz
... through chr21 ...
counts.tsv                        # original, retained and collapsed counts per chr/threshold
```

Every VCF has an index. `CollapseId` and `MatchId` link the retained representative and removed records. `logs/chrN.worker.log`, `chrN/logs/batch.log`, and `chrN/p{70,90}/logs/collapse-phased.log` contain tool output, paths, full parameters and phase-copy audit counts. No summarized controlled data are pushed to GitHub.

Phases are `plan`, `prepare`, `merge`, or `all`. Repeating `all` reuses checksummed preparation and completed threshold outputs. The input snapshot is frozen and retains its source metadata/checksums; use a new output root for a new panel generation. A change to the phase-consolidation implementation invalidates merges, because it can change genotypes; unrelated Git commits do not. If interrupted, verify the worker has stopped before removing its stale `.phased_batch.lock`.

## Tests

```bash
SV_PILOT_INTEGRATION=1 uv run --no-project --python 3.12 \
  --with truvari==5.5.0 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p test_sv_phased_batch.py -v
```

Tests run the actual native and wrapped collapse commands, reproduce the native donor-GT defect, verify identical native cluster/representative selection, and verify donor GT/PS/DP and phase flags in output. The two-chromosome batch checks 0.7 versus 0.9 results, opposite-haplotype same-sample calls staying separate under `--gt all`, phased multiallelic decomposition, all original allele keys accounted for, unchanged removed calls, retained symbolic/non-pure alleles, and cache reuse. Additional consolidation tests cover partial, unphased and haploid calls and rejection of unexpected co-carriage.
