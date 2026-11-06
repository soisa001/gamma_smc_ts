# Inspect popped chr1 before choosing a merge representation

The user-specified source is
`gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.chr1.popped.bcf`.
The corresponding `.csi` is staged when available; otherwise bcftools builds a
local index. No controlled data are accessed during local development.

`scripts/inspect_sv_representation.py` stages with the existing exact-file
`gcloud storage rsync` implementation and shared rsync process lock, validates
cloud size/checksums and source generation, and retains the original BCF.
There is no billing-project flag, upload operation, or merging. A 4 GiB free
scratch reserve is required in addition to staging capacity. Default staging
threads: 4; memory budget recorded as 26 GiB (not an OS-enforced cap).

Phases are `stage`, `inspect`, and `all`. Use `--input-bcf` to inspect an already
local indexed BCF/VCF instead of staging. Use `--compare-input` for an optional
local, indexed bubble-split baseline. Reports for the two inputs scan concurrently
(up to two scans with default four threads); native `bcftools query` streams
site fields, and Python does not decode cohort genotypes. Analysis stores
length histograms and same-position allele hashes rather than all allele strings.
The lock is released by the OS on process exit; do not delete its persistent file.

## Metrics

Original rows are counted before any splitting, including multiallelic row
counts and ALT-count quantiles. Each literal ALT is then inspected independently
in memory, without modifying the input. Scopes are `all_literal` and `net_ge50`
(absolute ALT-minus-REF length >=50 bp, applied only to reporting).

- Net INS/DEL, equal-length replacement, SNP, and reference-equal counts.
- Simple versus complex alleles: trim maximal identical prefix, then suffix
  without overlap; both residual sequences nonempty means complex. This is an
  operational representation category, not a biological event decomposition.
- Shared flank lengths, REF/ALT lengths, net lengths, and residual lengths.
  One anchor base is ordinary for indels; `shared_flank_gt1` distinguishes more
  than one shared base. Internal shared sequence is not measured by trimming.
- Exact distinct alleles/duplicates and distinct POS/REF groups within each
  scope/type. These groups do not establish original graph-site identity.
- Complex alleles with <50 bp net change but >=50 bp residual REF or ALT span.
- Symbolic/BND and missing/spanning alleles counted separately and excluded from
  literal length summaries. Header INFO identifiers are stored in JSON; the full
  header, including sample identifiers, stays in `header.txt` on the VM only.

Compare `net_ge50` allele metrics when the existing bubble baseline was previously
size-filtered. Its raw multiallelic-row fraction is not comparable if that file
was already split. Both files should use the same genome build and cohort; the
script does not prove allele/haplotype equivalence between different callsets.
It does not validate genotypes or calculate Truvari pairwise similarity.

Evidence of a simpler representation would be lower complex-allele fractions,
less substantial shared-flank sequence, and shorter REF/ALT lengths for similar
net sizes. Fewer records alone are insufficient: popping may select a different
set of intervals or decompose sites. The filename is not evidence that alleles
are minimal or safe to merge. Actual results require running on the VM.

## Outputs and reproducibility

Under `--root`, `net_ge50_comparison.tsv` joins both datasets. Each label directory
contains `summary.json`, `allele_summary.tsv`, `length_distributions.tsv`,
`header.txt`, completion hashes, and tool logs. Only aggregate counts/tables are
printed; no plots or individual sites. `inspection.log` records parameters,
paths, Python/script provenance, progress, and the fixed repository seed 1729.
Report caches depend on input identity, bcftools version, and analysis schema;
unrelated Git commits do not invalidate them. A fresh scan checks bcftools exit
status, sorted chr1 records and unchanged input identity. Cloud staging separately
validates full object checksums. Local comparison inputs are never changed.

Synthetic tests cover multiallelic rows, duplicate records, complex replacements,
trimmable simple indels, small-net complex alleles, sortedness, actual indexed
BCF decoding, cache reuse, and preservation of the comparison input:

```bash
SV_PILOT_INTEGRATION=1 BCFTOOLS_BIN=/path/to/bcftools \
  uv run --no-project --python 3.12 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p test_inspect_sv_representation.py -v
```
