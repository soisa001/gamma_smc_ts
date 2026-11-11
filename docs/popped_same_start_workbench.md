# Chr1 same-start sensitivity trial

This trial reuses the checksummed preparation from `chr1_popped_simple_merge`.
It runs 0.7 and 0.9 independently, only allowing comparisons within the same
**trimmed VCF CHROM/POS**, separately for simple INS and DEL. The equality rule
uses the valid anchored VCF POS, not the first changed base inferred from an
alignment; right-anchored records retain their existing anchoring convention.
No identical-REF restriction is imposed: same-start DELs of different lengths
remain candidates. Replacements pass through unchanged. Sub-50 bp simple
indels still participate, with net >=50 bp selection applied afterward.

The adapter `collapse_same_start.py` imports pinned Truvari 5.5.0, inserts
same-position groups before native size/distance subclustering, and calls native
matching, representative selection, genotype consolidation and VCF writers.
It does not modify the installed package or monkey-patch its functions. This is
an explicit candidate-grouping adapter, not an assertion that `-r 0` enforces
same POS. Native arguments stay `-p T -P T -r 500 -O 0 -s 0 -S -1 --keep common
--gt all`. Each emitted cluster is checked for equal CHROM/POS. Native
`CollapseId`/`MatchId` remain local to threshold/type; no cross-run ID comparison.

Original preparation is linked into the new run, not copied or rewritten.
The existing root is shared-locked while in use and is never used as the new
output root. Missing contig lengths use the existing all-record BED workaround.
All outputs stay on the VM. There is no staging, upload, or genotype repair.
The source preparation must remain available because the new root links to it.

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
export SV_CURRENT_ROOT="$HOME/chr1_popped_simple_merge"
export SV_SAME_START_ROOT="$HOME/chr1_popped_same_start"
bash scripts/run_popped_same_start_workbench.sh all
```

Phases: `merge` runs only the new merges; `compare` scans completed results;
`all` does both. Defaults are four threads and 26 GB, with two concurrent native
INS/DEL processes and sequential thresholds. Site comparisons use up to four
scan workers. Checksums are verified on resume; a new code commit alone does not
invalidate completed scientific results. Full integrity checks may take time
before progress messages appear. A failed new merge can be restarted with the
same cell without repeating preparation. Comparison failure can be retried with
the launcher `compare` phase alone.

New final outputs:

```text
$HOME/chr1_popped_same_start/p70/chr1.sv50.vcf.gz
$HOME/chr1_popped_same_start/p90/chr1.sv50.vcf.gz
```

Both contain the original cohort sample columns, native merged genotypes and
unchanged replacement SVs. Per-type retained/collapsed VCFs remain for audit.
Comparison files under `$HOME/chr1_popped_same_start/compare_current/`:

| File | Contents |
|---|---|
| `comparison.tsv` | Printed aggregate counts and unique bp before/after for all four policy/threshold combinations |
| `metrics.tsv` | Record counts, exact unique alleles, raw bp sums, duplicate counts, DEL reference union and INS unique allele bp |
| `same_start_vs_current.tsv` | Current and same-start result side by side, plus same-start minus current |
| `definitions.json` | Scope, definitions, exclusions, input manifests and provenance |

All metrics use literal net >=50 bp INS/DEL alleles and include unchanged
replacements in the allele-count denominator. Symbolic records are outside
this metric scope. Baseline uses the shared prepared representation, which
preserves record counts/net lengths and pure-deletion coverage while trimming
can create exact duplicates. Unique allele metrics explicitly deduplicate exact
CHROM/POS/REF/ALT after preparation. Final site multisets must be subsets of
that baseline; counts/bp must also agree with the preparation manifest.

For **DEL**, `unique_bp_*` is the union of deleted reference intervals, with
shared anchor/flank bases excluded. Each reference base counts once regardless
of how many DELs cover it. `unique_bp_removed` is the size of the before-union
minus the after-union, and the comparison checks that final deletion coverage
does not gain new bases. Complex/ambiguous DELs are excluded from this union
and their count is reported in `metrics.tsv`; the supplied popped diagnostic
had none. Both left- and right-anchored simple deletions are handled.

For **INS**, `unique_bp_*` sums net insertion length once per exact allele.
It is not a union of unique inserted sequence. TOTAL rows report allele counts
but leave combined unique bp blank because the two definitions differ.

Tests use only synthetic inputs. They exercise same-start vs shifted overlapping
deletions, padded insertions whose POS changes during trimming, matching at both
thresholds, small indels before final selection, unchanged replacement/sample
fields, native cluster links, correct deletion-union loss, duplicate alleles,
reuse of original files, and cached reruns. Run:

```bash
SV_PILOT_INTEGRATION=1 uv run --no-project --python 3.12 \
  --with truvari==5.5.0 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p 'test_popped*.py' -v
```
