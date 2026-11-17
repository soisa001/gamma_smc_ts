# TLR genotype transfer audit and separate trial

This audit traces original popped genotypes through native Truvari retained records,
the completed final VCF, and saved plot genotype tables at chr4:38,754,377.
It does not modify the completed chromosome callsets or patch installed Truvari.
Run `bash scripts/run_tlr_gt_transfer_audit_workbench.sh` in a separate checkout.
`POPPED_FULL_ROOT` defaults to `$HOME/popped_full_same_original_pos_p70`.
`TLR_GT_AUDIT_ROOT` defaults to `$POPPED_FULL_ROOT/tlr_gt_transfer_audit_38754377`.
Use `--plot-audit-dir PATH` if the focused plot tables are elsewhere.

Only aggregate results are printed. Sample identities, GTs, plot aliases such as
S02.h1, and stage comparisons are saved in `sample_gt_audit.tsv.gz` on the VM.
The source and mapping come from the completed run manifest; the audit checks
source file identities, mapping SHA256, sample order and exact allele identities.
It replays native consolidation independently and checks that input identities
remain unchanged. A run lock prevents concurrent writes. Re-running regenerates
the small audit outputs; `summary.json` is written only on successful completion.

Interpret `summary.json` as follows:

- `samples_losing_known_ALT_dosage_in_native`: sample/cluster comparisons whose
  original union contains more known ALT copies than the native retained GT.
- `native_final_gt_phase_ps_mismatches`: whether final assembly changed native GT,
  phase flag or PS. Zero means these stages agree.
- `native_replay_gt_phase_ps_mismatches`: whether replay of stock 5.5 consolidation
  reproduces the stored native result. Zero supports native consolidation as the
  explanation when ALT loss is already present there.
- `plot_audit.checked_genotype_rows` and `GT_string_mismatches`: whether the saved
  plot GT strings agree with the corresponding source/final VCFs. An empty object
  means no plot table was found; it does not establish plot correctness.

## Confirmed native behavior

Truvari 5.5 `gt_aware_consolidate` computes a replacement GT but omits GT from its
write-back loop. A synthetic representative `0|0` plus donor `1|0` produces `0|0`.
This is independent of `--gt all` preventing clustering of variants sharing an
ALT-carrying sample. If that restriction rejected a pair, both records should
remain; it does not explain removal of an already-clustered donor without transfer.

Simply switching to `--gt off` is insufficient for phased union: the other native
consolidation branch skips a representative that already carries ALT, so `1|0`
plus `0|1` can stay `1|0`. The trial explicitly handles each haplotype instead.

Two indexed locus-only VCFs and old/new mapping TSVs are written:

1. `existing_clusters.gt_union.vcf.gz`: preserve existing cluster membership and
   calculate GT from all original members, separating GT repair from clustering.
2. `intrasample.gt_union.vcf.gz`: re-cluster this position/type with stock Truvari
   matching, size/sequence 0.7, refdist 500, overlap 0, keep common, gt off, no chain,
   native consolidation disabled; then calculate the same explicit GT union.

For a cluster representing one allele class, per-haplotype ALT presence is OR:
`1|0 + 1|0 -> 1|0`; `1|0 + 0|1 -> 1|1`; `0|0 + 1|0 -> 1|0`.
The second example means both chromosomes carry an allele assigned to the class;
it does not claim their original insertion sequences were identical.
GT slots are assumed to have consistent panel orientation at the shared original
position. Conflicting explicit PS values and ambiguous unphased combinations are
flagged and set missing, not assigned guessed phase. Known ALT overrides missing
on that haplotype; otherwise missing remains missing. Singleton GTs are copied.
Merged representatives use trimmed alleles; singletons use their original alleles.
Trial FORMAT fields are GT/PS only, and AC/AN/AF are recomputed. These are diagnostic
locus trials, not replacements for full chromosome VCFs with quality/likelihoods.

## When the GT write-back regression appeared

The introducing commit is
[`9686b61ad0d3a2c4250424c5b37546891a1e8777`](https://github.com/ACEnglish/truvari/commit/9686b61ad0d3a2c4250424c5b37546891a1e8777),
November 2, 2023, titled `adding testing`. The function was called
`super_consolidate` then. Its final loop changed from iterating
`new_fmts.items()` (including GT) to `all_fmts` (which explicitly excludes GT).
It was renamed `gt_aware_consolidate` in commit `9da5ef0` on January 4, 2024.

The first stable tag containing the regression is v4.2.0. We executed the exact
historical helper functions, extracted with Python AST, against the same synthetic
pysam records in one modern environment. This is a function-level historical
reproduction, not an installation and full CLI run of every historical release.

| Source revision | Representative 0\|0 plus donor 1\|0 |
| --- | --- |
| Parent of 9686b61 | 0/1: ALT transferred, phase lost |
| 9686b61 | 0\|0: ALT lost |
| v4.1.0 | This genotype-aware function is absent |
| v4.2.0, v4.2.1, v4.2.2, v4.3.0, v4.3.1 | 0\|0: ALT lost |
| v5.0.0, v5.1.0, v5.1.1, v5.2.0, v5.3.0, v5.4.0, v5.5.0 | 0\|0: ALT lost |

This does not establish that all Truvari modes or all genotypes are affected.
Missing representative GTs follow an additional copy path, and gt off uses a
different consolidation function. Older releases are not a demonstrated phased
union fix. The actual S02 result still requires the VM audit above.

Validation: five focused tests pass, including a full synthetic popped-batch run,
original/native/final audit, lost donor ALT detection, separate intrasample union
yielding 1|1, recomputed allele counts, and unchanged original final-VCF SHA256.
