# Compare the popped TLR locus against the completed full p70 merge

This report reads the existing original popped BCF and completed chr4 final VCF.
It does not stage data, run Truvari, modify either input, or regenerate figures.
Only aggregate counts are printed. Individual allele and plot crosswalks remain
on the VM. A separate code checkout/output lock keeps it independent of jobs.

The plotted window and title describe different intervals. Both are reported:

| Scope | Region (1-based inclusive) |
|---|---|
| `plot_window` | chr4:38,660,000–38,960,000 |
| `TLR_interval` | chr4:38,780,001–38,890,001 |

The latter converts the plot target's BED interval [38,780,000, 38,890,001) to
VCF coordinates. The plot-window interval follows the previous TLR audit.
Count literal INS/DEL alleles with absolute `len(ALT)-len(REF)` at least 50 bp,
including net-length replacements, and require REF-span overlap with the window
(including the VCF anchor). These are whole-panel allele counts, not counts of
markers or carriers among the ten displayed haplotypes. Multiallelic ALT entries
are counted individually; the confirmed popped panel was biallelic.

```bash
%%bash
set -euo pipefail
umask 077
REPO="$HOME/tlr_comparison_code/gamma_smc_ts"
mkdir -p "$(dirname "$REPO")"
if [[ ! -d "$REPO/.git" ]]; then
    git clone --single-branch --branch AOU_run_opt \
        git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
git -C "$REPO" pull --ff-only origin AOU_run_opt
export POPPED_FULL_ROOT="$HOME/popped_full_same_original_pos_p70"
bash "$REPO/scripts/run_tlr_popped_comparison_workbench.sh"
```

Files under `$POPPED_FULL_ROOT/tlr_before_after/`:

- `counts.tsv`: before/after, removed count, percentage, simple allele counts
  and distinct exact allele counts, separately for INS/DEL/TOTAL and each window.
- `base_pairs.tsv`: total net allele bp and unique deleted reference bp for
  simple deletions. The deletion union excludes shared flanks/anchor, counts
  overlapping deleted bases once and is clipped to the selected window. The
  raw net length sum is a sum of whole selected allele lengths, not clipped.
- `allele_changes.tsv`: original and representative IDs, positions, lengths,
  cluster size and status for SVs in the full plotted window.
- `plotted_allele_changes.tsv`: matches the figure's `SV1`, etc. and marker colors
  to original and final allele identities, with cluster size and status.
- `summary.json`: input provenance, exact window definitions, policy and aggregate
  plot crosswalk results. This crosswalk describes identities; it does not infer
  post-merge phased GT carrier assignments or regenerate genotype-based markers.

Source paths come from the completed batch manifest. Before reporting, the tool
checks source/index stat provenance, the lineage file SHA256, successful indexed
regional reads, and equality between the regional final VCF and representatives
predicted by the source-to-final mapping. It does not rehash the entire large
chr4 VCF just to inspect this locus. Trimming-induced boundary exits are accounted
for when validating the mapping; before-minus-after includes such exits if any.
Completed plot tables are optional; absent tables are reported explicitly.
Individual ID/sequence/genotype values are not printed in notebook output.

The report is intentionally separate from the earlier TLR trial audit: that
trial used local candidate windows and same-trimmed-position grouping. Here the
comparison uses the actual full-chromosome, same-original-position p70 output.
