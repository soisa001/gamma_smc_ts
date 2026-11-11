# TLR-locus SV similarity and merge audit

The supplied figure shows hg38 chr4:38,660,000–38,960,000. This matches
`ashg_revised/TLR_CUSTOM_PLOT_CELL.md` and the TLR custom-target table in the
`hmmix_aou` checkout. The user selected the **popped long-read panel** for this
audit. The figure pipeline itself reads precomputed Phase-10 SV links/event
caches; its upstream launcher defaults to the **bubble-split panel**, and
configuration can override those sources. Popped alleles therefore must not be
assumed to correspond one-to-one to the figure's markers.

Observed source chain in `soisa001/hmmix_aou` (local checkout commit `3a0857f`):

- `ashg_revised/shortlist_followup.py` calls `sv_sources(cfg)` for the TLR panel.
- `ashg_revised/exploratory_sv.py` reads `segment_sv_links*.v2.tsv.gz` and
  `event_cache/chr*/events.tsv.gz`, or explicitly configured alternatives.
- `scripts/run_sv_introgression_analysis.sh` defaults `PANEL_TEMPLATE` to
  `panel_bubble_split_vcf/aou_lr_phase2_v1.{contig}.bubble.split.bcf` and records
  the effective template in its run configuration.

The audit defaults to:

`gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.chr4.popped.bcf`

Use `TLR_POPPED_INPUT` for an already-local indexed original panel BCF. The
launcher checks a few known local paths, otherwise stages the full chr4 BCF and
CSI once using the existing serialized, CRC-checked gcloud storage rsync helper.
Initial staging may dominate runtime; it is not a remote range-only download.
After staging, only the target window plus 2 kb flanks is extracted. No cloud
access or controlled data analysis is performed on the development machine.

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts_reports"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
  git -C "$REPO" switch AOU_run_opt
  git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
cd "$REPO"
export TLR_AUDIT_ROOT="$HOME/tlr_popped_sv_audit"
# Optional: avoids staging when the indexed original chr4 popped BCF is present.
# export TLR_POPPED_INPUT='/absolute/path/aou_lr_phase2_v1.chr4.popped.bcf'
bash scripts/run_tlr_sv_audit_workbench.sh --phase all
```

## What is checked

1. Partition and trim simple INS/DEL exactly as in the chr1 trial. Other >=50 bp
   literal SVs pass through. Original records are never changed. Working copies
   receive unique audit IDs to make duplicate original IDs unambiguous; the
   original ID, position and full allele strings remain in the VM-side map.
2. Run both current policies (500 bp matching and same post-trimming CHROM/POS)
   independently at 0.7 and 0.9, with native Truvari 5.5.0 `--keep common --gt all`,
   `-O 0 -s 0 -S -1`. Smaller simple indels in the window/flanks participate.
   No reference-based left-alignment is performed.
3. Compare every pair of >=50 bp literal INS/DEL alleles overlapping the focus
   window: native size and sequence similarity, reciprocal overlap, position
   constraint, simple-allele eligibility, and the native shared-carrier block.
   Whole-input cohort genotypes are used; plotted haplotypes are not the cohort.
4. Report direct pair eligibility separately from actual regional cluster
   membership. Native representative choice and other alleles can change the
   latter. A cluster is not an all-pairs clique. Whole-chromosome clustering can
   differ because the regional subset omits additional outside candidates.

The additional payload score removes shared flanks and reports global
`1 - edit_distance / max(payload_lengths)`. It is explicitly separate from
Truvari's native sequence score, which is the score used for merge eligibility.
Long-sequence pair scores above `--max-sequence-bp` (default 100,000) are marked
unassessed; native local collapse still uses its ordinary rules. The all-pairs
audit stops above `--max-focus-alleles` (default 1,000); narrow `--region` or raise
this explicit compute guard. It does not silently sample or omit pairs.

## Output and interpretation

Only aggregate JSON is printed. Individual records, sequences, positions, IDs,
per-pair carrier counts and native output genotypes remain in VM files:

| File under `$TLR_AUDIT_ROOT` | Contents |
|---|---|
| `audit/alleles.tsv` | Audit label, original ID, original/trimmed POS, type, net length, eligibility |
| `audit/pair_scores.tsv` | Size/sequence scores, exact failed criteria, direct pair pass and actual local cluster result |
| `audit/local_cluster_membership.tsv` | Each focus allele's representative for all four runs |
| `audit/summary.json` | Aggregate counts and interpretation limits |
| `prepare/alleles.jsonl.gz` | Original and trimmed full allele sequences, kept only on VM |
| `current_r500/p70/`, `current_r500/p90/` | Local native retained/collapsed VCFs and final SV subset |
| `same_start/p70/`, `same_start/p90/` | Corresponding same-start outputs |
| `logs/audit.log` | Commands, provenance, parameters and failures |

`pair_passes=True` means those two records satisfy the direct rules in that row.
`same_local_cluster=True` means the actual regional collapse assigned them to
the same representative. `shared_carriers_gt_all` explains why even identical
alleles can remain separate under the requested native genotype constraint.
Representatives below 50 bp are explicitly flagged because final SV filtering
can omit them. Pairs involving replacements report `non_simple_passthrough`.

Phases are `prepare`, `audit`, and `all`. `audit` reuses the existing preparation;
use the same region/flank arguments. Repeating `all` reuses verified regional
outputs. Source identity includes the BCF/index path, size and modification time;
staged objects additionally undergo the existing cloud CRC checks. No full local
chr4 hash is required just to extract a small locus. Keep the source file immutable
during extraction. The synthetic integration test checks size/sequence thresholds,
shifted starts, shared carriers, replacements, native cluster links, small-indel
participation and cache reuse.
