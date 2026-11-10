# Chr1 popped-panel trial: merge pure indels only

Run `scripts/run_popped_simple_sv_workbench.sh` for independent 0.7 and 0.9 runs.
The default source is the confirmed full popped panel:

`gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_popped_vcf/aou_lr_phase2_v1.chr1.popped.bcf`

The launcher first reuses the indexed popped input recorded in
`$HOME/chr1_sv_representation_comparison/comparison_metadata.json`, then checks
the known popped inspection locations. Set `SV_POPPED_INPUT` to override this.
Otherwise it stages the exact BCF and index with the existing serialized,
checksummed `gcloud storage rsync` helper. No billing override is required;
`--billing-project` remains available if needed. No data or results are uploaded.

## Selection and native merging

- Classify each biallelic allele after maximal nonoverlapping shared-prefix and
  shared-suffix trimming. Only unambiguous A/C/G/T alleles with exactly one empty
  residual sequence enter collapse. A one-base REF alone is not sufficient.
- Keep a shared anchor when writing the trimmed VCF; adjust POS and an existing
  END when trimming. Do not left-align, atomize, split records, change GT/PS, or
  repair native genotype consolidation. An inconsistent existing SVTYPE/SVLEN
  causes an error rather than silently influencing native matching.
- Run INS and DEL separately with native Truvari 5.5.0:
  `-p T -P T -r 500 -O 0 -s 0 -S -1 --keep common --gt all`.
  Both size and sequence similarity change with T. No `--hap` or `--chain`.
  With the default four-thread budget, INS and DEL run concurrently; the two
  thresholds run sequentially, independently from identical prepared inputs.
- **No minimum input size:** all simple nonzero-length indels can participate.
  Equal-length literal changes cannot be pure indels or meet the final net-length
  criterion, so bcftools omits them from the extraction stream. Other sub-50 bp
  replacements are not merged and are omitted from the final SV subset.
- Final literal INS/DEL selection is `abs(len(ALT)-len(REF)) >= 50`, applied after
  merging. Replacements and ambiguous alleles satisfying that cutoff are copied
  unchanged, including sample fields. Symbolic SVs/BNDs, if present, pass through
  unchanged because literal net size is undefined; they are outside INS/DEL
  summary denominators. Equal-length replacements/SNPs remain in the original
  source, which is never changed or deleted.

The final file is an **SV subset**, not a rewritten 10-million-record SNP panel.
There is no giant `new.unsorted.vcf`. Partition files are BGZF-compressed while
streaming. Only partitions that need sorting are sorted, using bounded memory.
No sample fields or individual variants are printed to the notebook. Native
genotype consolidation is used as requested; matching the native output is not
a claim that native GT consolidation fixes the previously discussed limitations.

## Workbench cell

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
export SV_SIMPLE_ROOT="$HOME/chr1_popped_simple_merge"
# Optional, if the automatic lookup does not find your staged popped BCF:
# export SV_POPPED_INPUT='/absolute/path/aou_lr_phase2_v1.chr1.popped.bcf'
bash scripts/run_popped_simple_sv_workbench.sh --phase all --similarities 0.7 0.9
```

Phases are `prepare`, `merge`, and `all`. Repeat the same cell to resume verified
preparation and completed partition/threshold results. Cache keys use source
checksums and analysis settings, not the Git commit. Keep this new root separate
from earlier bubble/all-allele experiments. A Linux process lock prevents two
runs in the same root and releases automatically on process exit. Prepared inputs
are frozen during a merge-only phase.

Outputs under `$SV_SIMPLE_ROOT`:

| Path | Contents |
|---|---|
| `comparison.tsv` | INS, DEL, total counts/bp before/after and percentages, both thresholds |
| `prepare/summary.json` | Original SV50 counts, simple/passthrough counts, all-size partition counts |
| `prepare/INS.vcf.gz`, `prepare/DEL.vcf.gz` | Trimmed simple indels, all sizes, with genotypes |
| `prepare/passthrough.sv50.vcf.gz` | Untouched replacement/ambiguous SVs and symbolic SVs/BNDs |
| `prepare/trimmed_alleles.tsv.gz` | Original-to-trimmed allele mapping, only on VM |
| `p70/chr1.sv50.vcf.gz` | Final combined 0.7 VCF, indexed, with native sample genotypes |
| `p90/chr1.sv50.vcf.gz` | Final combined 0.9 VCF, indexed, with native sample genotypes |
| `p70/INS/retained.vcf.gz`, `collapsed.vcf.gz` | Unmodified native partition outputs; also DEL and p90 |
| `logs/workflow.log` | Provenance, progress, parameters and commands |

Native `CollapseId`/`MatchId` annotations are **local to each INS or DEL run**;
use the type-specific files for membership lookup. They are not renumbered in
the combined VCF. Earlier `compare_sv_native.py`/`audit_sv_native.py` batch-root
layouts do not apply to this new partitioned workflow.

Read the `all_INS_DEL` rows in `comparison.tsv` for the proportion removed out
of all original literal SVs, including unchanged replacements. `simple_only`
uses only mergeable simple alleles as its denominator. `unchanged_passthrough`
must show zero reduction. All these rows use the same final net >=50 bp rule on
both sides. `merged_out = before_sv50 - after_sv50`; raw bp sums are net allele
length sums, not deleted-reference union coverage or genotype-weighted burden.
The confirmed input is expected to give **86,733** literal SVs: **66,369 INS**
and **20,364 DEL**, of which **1,174 INS replacements** are left unchanged.
Actual merge counts are determined on the VM, not inferred from the earlier run.

Synthetic tests cover prefix/suffix trimming, one-base replacements, unchanged
passthrough sample fields, exact agreement with independent native calls, small
indels participating before final filtering, empty partitions, checksums and
resuming, and repair of a deliberately corrupted final file:

```bash
SV_PILOT_INTEGRATION=1 uv run --no-project --python 3.12 \
  --with truvari==5.5.0 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p test_popped_simple_sv.py -v
```
