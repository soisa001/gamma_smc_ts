# Compare the popped TLR locus against the completed full p70 merge

This report reads the existing original popped BCF and completed chr4 final VCF.
It does not stage data, run Truvari, or modify either input. The optional
haplotype plot reads the actual before/after genotypes and saves figures locally.
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
bash "$REPO/scripts/run_tlr_popped_comparison_workbench.sh" \
    --plot-haplotypes --haplotype-window plot_window --max-haplotypes 8 --threads 4
```

Files under `$POPPED_FULL_ROOT/tlr_before_after/`:

- `counts.tsv`: before/after, removed count, percentage, simple allele counts
  and distinct exact allele counts, separately for INS/DEL/TOTAL and each window.
- `counts_by_class.tsv`: separate simple INS/DEL, complex replacements and
  ambiguous-sequence alleles, with their own before/after denominators. These
  are the counts printed in the notebook. Complex means both REF and ALT retain
  sequence after removing their common prefix and suffix; it takes precedence
  over the ambiguous category if both apply. The report verifies that no complex
  or ambiguous allele entered the simple-only merge mapping. The original merge
  already excluded them; this change separates their reporting.
- `base_pairs.tsv`: total net allele bp and unique deleted reference bp for
  simple deletions. The deletion union excludes shared flanks/anchor, counts
  overlapping deleted bases once and is clipped to the selected window. The
  raw net length sum is a sum of whole selected allele lengths, not clipped.
- `allele_changes.tsv`: original and representative IDs, positions, lengths,
  cluster size and status for SVs in the full plotted window.
- `plotted_allele_changes.tsv`: matches the figure's `SV1`, etc. and marker colors
  to original and final allele identities, with cluster size and status.
- `summary.json`: input provenance, exact window definitions, policy and aggregate
  crosswalk and optional haplotype-plot results.

With `--plot-haplotypes`, these additional files are saved in the same folder:

- `tlr_haplotypes_before_after.png` and `.pdf`: unmerged on top, merged on bottom,
  same haplotypes and order, covering **chr4:38,660,000–38,960,000** by default,
  including the merged insertion at **38,754,377**. The PDF
  uses editable text and vector marks. No notebook display is invoked.
- `selected_haplotypes.tsv`: figure labels, actual sample IDs, phase slots and
  documented introgression status. Sample labels such as S01 refer to the same
  individual in both panels; h1/h2 are its two VCF GT slots.
- `haplotype_genotypes.tsv.gz`: actual GT, PS, allele identity, class and plotted
  state for each selected haplotype and before/after allele.
- `haplotype_placement_audit.tsv`: compares the union of source cluster ALT
  placements with actual merged ALT placement, per cluster and haplotype.
- `haplotype_plot_summary.json`: selection parameters, phase assumptions, input
  fingerprint/SHA256, aggregate placement changes, software versions and output
  checksums. It is written only after successful plotting.

The candidate pool is the samples in the existing plot's
`popped_selected_lanes.tsv.gz`, not the entire cohort. Called tracts in that table
supply the shading. Selection first favors introgressed haplotypes carrying
merged-cluster members, then other merged-cluster carriers, then context lanes;
within these groups it favors distinct original alleles. Ties are deterministic.
If a sample's ALT placement changes, its partner phase slot is included within
the requested lane limit so a phase-slot move is visible. If no merged cluster
is carried by these saved samples, the aggregate output reports zero clusters
shown; the workflow does not invent examples or resample the whole cohort.

Use `--haplotype-window TLR_interval` to reproduce the earlier narrow extent
(chr4:38,780,001–38,890,001). The wider default selects lanes using merged
clusters throughout the wider region, so it can prioritize insertion carriers;
the selection can differ from the previous narrow figure, but both panels
always use the same selected haplotypes. To preserve the narrow outputs, set
`TLR_COMPARE_ROOT="$POPPED_FULL_ROOT/tlr_before_after_with_insertion"` before
running the wider comparison. Genotypes are reread from both original and
merged files, and existing merge outputs are never changed.

Cluster labels such as `C1:54` mean cluster C1, 54 bp; the representative keeps
the same cluster label after merging. Simple INS/DEL and complex replacements
use different marker shapes. Only merged-cluster alleles receive text labels;
other regional SVs remain visible as markers. Both panels show the actual VCF
genotypes, **not** source carriers copied onto the representative. Unphased
heterozygotes, missing phase slots and non-diploid calls get an unknown marker;
homozygotes have phase-independent placement. Known ALT in a partial phased GT
is retained on its known slot. A blank tract background does not establish
absence of introgression because the saved lane table is a selected subset.

The introgression decoder and original VCF are assumed to use the same hap1/hap2
orientation, as in the original popped plot. This does not rephase the merged
VCF or prove consistency across phase sets. Differences in ALT presence or
phase-slot placement are shown and counted, not silently repaired.
PS values are saved to support follow-up checks.

## Show different source members of the insertion merge

When the saved plot samples carry only the representative, choose
`--carrier-pool cohort --focus-merge-pos 38754377`. This screens source GTs across
all panel samples for members of the merge at that original position. It keeps
at most two candidate haplotypes per source allele and their partner slots,
then reads regional original/merged GTs only for the small resulting sample
pool. Selection prioritizes coverage of distinct target source alleles before
documented introgression. It does not select based on successful genotype
transfer or conceal losses/phase changes. Both panels use the same haplotypes.

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
export TLR_COMPARE_ROOT="$POPPED_FULL_ROOT/tlr_merge_carrier_examples"
bash "$REPO/scripts/run_tlr_popped_comparison_workbench.sh" \
    --plot-haplotypes --haplotype-window plot_window \
    --carrier-pool cohort --focus-merge-pos 38754377 \
    --max-haplotypes 8 --threads 4
```

Alongside the full regional figure, `tlr_focus_merge_before_after.png` and `.pdf`
zoom to the target merge with 2-kb flanks. `carrier_discovery.tsv` records the
number of known ALT carriers for each original member. Aggregate summary fields
`target_source_alleles` and `target_source_alleles_shown` show whether every
member is illustrated; members without any known phased/phase-independent ALT
carriers are reported, never fabricated. `selected_haplotypes.tsv` identifies
cohort-discovered samples and marks ancestry as unassessed when no saved tract
is available. An unshaded newly selected haplotype is **not** thereby confirmed
unintrogressed. Existing called tracts are reused where available; this option
does not regenerate introgression calls. All carrier identities stay on the VM.

Set `TLR_PLOT_DIR` if the existing private popped plot tables are elsewhere.
The haplotype option requires the selected-lanes table; the older identity-only
comparison can still run without it. Rerunning this small regional report
rebuilds its outputs, repairing corrupt figures, without restaging data or
rerunning any merge. Scope remains literal alleles with absolute net length
at least 50 bp; equal-length replacements and replacements with smaller net
length are outside this report even if their changed spans are long.

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
