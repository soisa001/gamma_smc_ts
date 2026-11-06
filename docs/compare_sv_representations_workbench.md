# Compare original bubble, popped and imputed chr1 representations

Run `scripts/compare_sv_representations.py` on the controlled VM. No controlled
data are accessed during development. The default sources, explicitly supplied
or previously confirmed by the user, are:

```text
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.chr1.popped.bcf
gs://long-reads-phase-2-imputation-vwb/batch-001/aou2_50k_init.batch-001.chr1.prod-05fea60.glimpse2.popped.bcf
```

The staging phase uses the existing exact-file `gcloud storage rsync` helper,
without a billing-project flag, under the shared workflow rsync lock. It checks
cloud checksums, source generations and indexes. It downloads the corresponding
`.csi` when present or builds a local index. There are no uploads. Use `--bubble`,
`--popped`, or `--imputed` to reuse local indexed originals. Prefer the original
BCFs rather than prepared filtered inputs when investigating record packing.

Phases: `stage`, `compare`, or `all` (default). Default root is explicitly
provided by the caller. Default threads=4 allows up to three concurrent input
scans; staging is serialized. Memory budget=26 GB is recorded, not enforced.
SQLite stores hashes and counts on disk; REF/ALT strings are retained only while
scanning records or extracting selected examples. Keep disk available for three
BCFs and their SQLite indexes. No full cohort VCF is rewritten.

## What is compared

The default scope is literal net INS/DEL alleles with
`abs(len(ALT)-len(REF)) >= 50`. `--min-net-bp` changes that reporting scope.
Each ALT is handled independently, without changing input files, GTs, IDs or
annotations. No AF, AC, FILTER or genotype eligibility filter is applied.
Symbolic/BND, spanning `*`, equal-length replacements and SNPs are outside this
SV scope. Consequently this comparison cannot establish equality of whole BCFs.

For every dataset pair and INS/DEL type, report both:

- **exact:** same CHROM/POS/REF/ALT, compared through SHA256 keys. One multiallelic
  row and several biallelic rows can have identical allele sets. Unique-set and
  record-multiplicity equality are reported separately.
- **canonical:** uppercase REF/ALT, remove common prefix then common suffix,
  and advance POS by the prefix length. Empty residual REF/ALT is permitted in
  this internal key. This recognizes shared-flank representation differences;
  it is NOT reference-guided repeat left alignment, atomization, or a proof of
  general haplotype equivalence. It assumes a shared reference genome.

An imputed set may be a subset of a panel. Check both directional shared
percentages and only-in-one-file counts before interpreting a lower allele
count as a representation change. These statistics do not assess imputation
accuracy, phase accuracy, or genotype concordance.

## VM-only examples

Examples are selected deterministically in SHA256 order, default three per
category. Shared simple/complex alleles and file-specific alleles are included.
They are illustrative, not random prevalence estimates. Each example retains
source ID, ALT index/count, full REF/ALT, trimmed residual strings and lengths.
Edlib aligns the residual ALT (query) to residual REF (target), reporting one
minimum-edit CIGAR and counts: `=` matches, `X` substitutions, `I` inserted ALT
bases, `D` deleted REF bases. Repeats can admit other equally optimal placements.

Pairs of distinct ALTs from selected complex POS/REF groups have native Truvari
5.5.0 sequence and size scores, alongside direct trimmed REF/ALT scores. Source
SVTYPE/SVLEN are used when present. No GTs are read for these scores. Passing
sequence/size thresholds is NOT a predicted actual merge: genotype conflicts
and other criteria are not evaluated. These are not scores from normalized
Truvari calls. Alignment/score work is skipped above `--max-align-bp` (10,000
default); the example strings and skip reason are still saved.

Outputs under `--root`:

- `record_structure.tsv`: original rows, multiallelic rows, ALT entries,
  selected alleles, and sample-column counts.
- `pairwise_summary.tsv`: aggregate three-way exact/trimmed comparisons.
- `allele_examples.txt`: readable previews and edit explanations.
- `allele_examples.jsonl`: full individual allele strings and metadata.
- `within_site_pairs.jsonl`: selected pair scores, linked by example hashes.
- `comparison_metadata.json`: source paths and definitions.
- Per-input directories: local headers, logs, hash/count SQLite database,
  counts, and checksummed cache manifest.

Only the two aggregate TSVs and aggregate progress are printed. Individual
alleles and headers remain on the VM. No plots are displayed. Completion checks
include bcftools decoding, sortedness, SQLite integrity, input/index identities
unchanged, recovered example hashes, and checksummed output manifests. Caches
depend on input identities/settings/schema rather than Git commit. RNG is not
used; repository seed 1729 is recorded. An OS lock prevents concurrent runs at
the same root and releases on exit; do not delete the persistent lock file.

Synthetic tests exercise row packing, duplicate multiplicities, changed site
sets, common-prefix trimming with a POS shift, actual BCF indexed queries,
native complex-DEL scoring, empty residual edit alignments, and cache reuse.
