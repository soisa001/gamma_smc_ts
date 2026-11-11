# Full chr1–chr21 popped callsets with same-original-position SV merging

This workflow follows the final position/representation policy, which differs
from the earlier same-**trimmed**-position trial. It does not alter those trials.

1. Start from the indexed long-read **popped panel**, not bubble or imputed.
2. Trim identical REF/ALT flanks internally, then classify pure INS/DEL. A
   replacement retaining sequence on both sides is not a merge candidate.
3. Eligible alleles are biallelic simple INS/DEL with absolute net length at
   least 50 bp. SNVs, smaller indels, replacements, symbolic records and
   multiallelic records pass through unchanged. The 50-bp candidate cutoff is
   intentional here to meet the request that smaller indels remain unchanged;
   it differs from the earlier all-size candidate trial.
4. Only alleles with identical **original CHROM, original POS and type** can
   merge. Trimmed positions are not the grouping key. Different original starts
   never merge even when trimming produces identical alleles/positions.
5. Within each group, use native Truvari 5.5.0 `-p 0.7 -P 0.7 --keep common
   --gt all -r 500 -O 0`, no `--chain`, no `--hap`, no reference alignment.
   The 500-bp native check remains an additional constraint on trimmed alleles;
   it never permits crossing original-position groups. Each added allele must
   match its representative, rather than every other member. Genotype
   consolidation and representative ordering are unmodified native behavior.
6. Only representatives that actually absorb another allele receive the trimmed
   POS/REF/ALT and native consolidated sample fields. Updating POS when removing
   shared prefix bases is required for valid VCF coordinates, not left alignment.
   Unmerged candidates retain their entire original records, including untrimmed
   alleles, INFO, FORMAT, phased GTs, and IDs. Removed members are omitted.
7. Write one complete sorted/bgzip VCF and tabix index per chromosome. No global
   SV-only output filter is applied. All source sample columns are retained.

The adapter calls native `collapse_chunk` on one original-position/type group
at a time, in one process per type/chromosome; it does not launch a CLI process
per site. Working inputs are ordered by original position and intentionally
unindexed, since trimming can change coordinate ordering. Native output is
sorted/indexed afterwards. No chromosome-length header guesses are required.

## Workbench cell

The dedicated checkout/runtime/output root keeps this separate from ongoing
chr1 and TLR trials. Up to four chromosomes run concurrently, with INS and DEL
processed sequentially within each chromosome. Native collapse runs on one
thread per chromosome. Streaming decode/sort helpers can add CPU activity;
the 26-GB total budget guides sort buffers and is not an OS-enforced memory cap.

```bash
%%bash
set -euo pipefail
umask 077
REPO="$HOME/popped_full_merge_code/gamma_smc_ts"
mkdir -p "$(dirname "$REPO")"
if [[ ! -d "$REPO/.git" ]]; then
    git clone --single-branch --branch AOU_run_opt \
        git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
git -C "$REPO" pull --ff-only origin AOU_run_opt

export POPPED_FULL_ROOT="$HOME/popped_full_same_original_pos_p70"
export POPPED_FULL_JOBS=4
export POPPED_FULL_MEMORY_GB=26
# Optional, if all indexed popped BCFs are already local:
# export POPPED_INPUT_TEMPLATE='/path/to/aou_lr_phase2_v1.{chrom}.popped.bcf'
# Optional requester billing override; otherwise retain gcloud defaults:
# export POPPED_BILLING_PROJECT="$GOOGLE_PROJECT"
bash "$REPO/scripts/run_popped_full_batch_workbench.sh" --phase all
```

Exact default source template (plus `.csi`):

```
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.{chrom}.popped.bcf
```

Known already-staged chr1/TLR popped inputs can be reused read-only. Otherwise
the existing serialized, checksum-verified `gcloud storage rsync` helper stages
one chromosome/index at a time; simultaneous rsync directory-listing collisions
are avoided. Inputs are retained, not deleted. No controlled data are uploaded.

## Outputs and resume

For each `chrN` under the output root:

- `final/chrN.popped.samepos.p70.vcf.gz` and `.tbi`: full callset.
- `final/old_to_new.jsonl.gz`: exact original and final allele identities for
  every candidate, representative working ID, cluster size and status. Many old
  rows may point to one representative, so the same file supports reverse lookup.
  All records outside the candidate table map to themselves unchanged.
- `final/summary.json`: original/final full-record counts, unchanged records,
  removed INS/DEL counts, merged representative counts and parameters.
- `prepare/source_alleles.jsonl.gz`: original-to-trimmed candidate identities.
- `merge/{INS,DEL}/{retained,collapsed}.vcf.gz`: native working-ID outputs.
- `logs/workflow.log`: provenance, commands, progress every 100,000 records during
  full assembly and error tracebacks. Per-type native logs are alongside it.

At the batch root, `comparison.tsv` reports before/after biallelic net-SV50
counts, numbers removed and percentages by INS/DEL/TOTAL for every chromosome.
`batch_status.json` records completed/failed jobs.

Re-running the same cell validates/reuses completed phases. Cache keys depend
on source file/index identities, sample manifest and scientific parameters;
they do not depend on Git commit IDs. Phase artifacts have SHA256 checksums.
Available phases: `stage`, `prepare`, `merge`, `finalize`, `all`. To rerun only
one chromosome, append `--chromosomes chr1`; completed chromosomes otherwise
remain reusable. Distinct output roots have independent locks.

Writing the complete cohort VCF requires a full chromosome read and sort, so
this step is substantially larger than the previous SV-only outputs. Assembly
streams directly through bcftools sort/compression; there is no giant plaintext
`new.unsorted.vcf`. Only one merged representative's genotype record is loaded
at a time using BGZF virtual offsets. Sorting may spill to chromosome-local
temporary disk. Source BCFs and intermediate/final outputs remain on the VM.

Native genotype consolidation is unchanged, as requested earlier: this workflow
does not implement a custom phase/GT repair or assert that native consolidation
fixes the previously discussed genotype-transfer defect. Unmerged records keep
their original phased genotypes.

## Validation

The native integration test covers two concurrent chromosomes, same original
positions with different trimmed positions, different original positions that
trim to the same position, INS/DEL separation, missing contig lengths, exact
preservation of unmerged padded alleles/SNVs/small indels/replacements and
multiallelic records, a chromosome with no eligible SVs, membership mapping,
count conservation, cache reuse, and repair of a corrupted final output without
rebuilding preparation. All fixtures are synthetic.
