**Paused EAS array: preliminary evaluation, 17 September 2026.**

The user subsequently specified a positional/gene-level FPR. The
[matched-position evaluation](positional_eas_evaluation.md) is now the primary
analysis. This document retains the earlier whole-region-maxima endpoint.

The user requested a pause before completion. No simulation or decoder process
was active at inspection. A validation-only pass finished at 23:50:58 UTC and
verified 1,947 saved regions: 1,000 neutral and 947 selected. There are 53
unfinished selected simulations. No simulation or decoding was started for
this evaluation. The archive, partial attempts, and deterministic seeds remain
on D for a later explicit resumption.

**Results at T=50 kya and regional p<=0.05.** All detection columns below
score the known selected site against neutral whole-region maxima. They are
lower bounds on full regional power, not completed new regional scans.

| s | Saved selected regions | Mean selected AF | AF-only detection | ALT/ALT carrier-mass detection | All-pair detection |
|---|---:|---:|---:|---:|---:|
| 0.001 | 92 | 17.9% | 0.0% | 1.1% | 0.0% |
| 0.002 | 93 | 30.8% | 8.6% | 7.5% | 0.0% |
| 0.003 | 95 | 43.9% | 25.3% | 26.3% | 5.3% |
| 0.004 | 94 | 59.5% | 50.0% | 47.9% | 8.5% |
| 0.005 | 100 | 78.4% | 79.0% | 77.0% | 33.0% |
| 0.006 | 94 | 81.4% | 83.0% | 81.9% | 47.9% |
| 0.007 | 94 | 90.7% | 92.6% | 92.6% | 54.3% |
| 0.008 | 95 | 95.1% | 97.9% | 94.7% | 71.6% |
| 0.009 | 95 | 98.5% | 100.0% | 94.7% | 61.1% |
| 0.010 | 95 | 99.2% | 100.0% | 98.9% | 73.7% |

The corresponding held-out neutral region call fractions are 4.9% for AF,
5.0% for carrier mass, and 5.2% for all-pair recency. At p<=0.01, s=0.005
selected-site detection is 61% for AF and 57% for carrier mass, with neutral
region call fractions of 1.2% and 0.9%, respectively. All-pair detection is
21% with 0.7% neutral calls. The 77/100 carrier result at p<=0.05 has a 95%
Wilson interval of 67.8% to 84.2%. These focal results do not replace the
earlier full-scan s=0.005 results (83% truth carrier-mass power and 84% AF
power), which can also detect signals away from the selected site.

Across arms there is no clear incremental advantage over AF in this focal
evaluation. At s=0.003, carrier mass adds two regions missed by AF and loses
one AF hit; at s=0.005, it adds three and loses five. The paired exact tests
at T=50 kya have p>=0.0625 across all ten arms; these exploratory comparisons
do not establish equivalence. Using the 91 replicate IDs completed in every
arm changes the primary detection fractions by at most 1.26 percentage
points. The saved tables include all six T cutoffs and p<=0.01/p<=0.05.
Shorter cutoffs do not show a consistent large improvement: at s=0.003,
carrier detection is 28.4% at 10 kya versus 26.3% at 50 kya; at s=0.004 it
is 42.6% versus 47.9%; and at s=0.005 both give 77%.

The lineage diagnostic explains why AF and carrier mass are so similar in
the weaker arms. At s=0.001, 0.002 and 0.003, respectively, 85/92 (92.4%),
83/93 (89.2%) and 76/95 (80.0%) regions have one sampled-carrier ancestral
lineage at 50 kya at the selected site. Median ALT/ALT `frac_recent_50000`
is 1 in these arms. When all carriers share a lineage, carrier mass equals
the ALT/ALT pair weight, approximately AF squared, and supplies no independent
recency ranking among those regions. The AF/carrier-mass Spearman correlations
in the three arms are 0.995, 0.991 and 0.981. Thus, varying s in this model
does not provide a balanced test of the proposed advantage across multiple
contributing lineages. A future founder-diversity comparison would need to
measure or vary that feature explicitly and retain matching neutral controls.
No such simulations have been launched.

The [arm summaries](results/eas_h400_paused_evaluation/arm_summary.csv),
[detection estimates and intervals](results/eas_h400_paused_evaluation/selected_site_metrics.csv),
[neutral region call fractions](results/eas_h400_paused_evaluation/neutral_region_metrics.csv),
and [paired comparisons](results/eas_h400_paused_evaluation/af_mass_paired_comparison.csv)
are exported to the repository. Full per-region diagnostics and predictions
remain on D. The [audit](results/eas_h400_paused_evaluation/audit.json) passed:
1,900 old s=0.005 focal values matched exactly, neutral call counts reproduced
the previous analysis, and the lineage partition agreed with direct pairwise
TMRCA for every evaluated pair.

The archive is `/mnt/d/phase2simselection/sim/eas_q02_h400`; the evaluation is
`/mnt/d/phase2simselection/sim/eas_h400_pause_evaluation_20260917`. The root
`simulation_status.json` now records the validation result rather than the
stale running status. `pause_requested.json` records the user's pause, and
`simulation_status_at_pause.json` preserves the last coordinator report.
`pause_state.json` records the 53 missing identities and that resumption has
not been authorized.
`simulation_inventory.csv` indexes all 1,947 verified inputs. The validation
checked scientific identities, file hashes, readable raw/cropped trees,
sample counts, focal coordinates, and genotype/carrier-mask agreement.

**What is evaluated.** For each saved selected region, extract the true local
TMRCA of the same 10,000 pairs at the selected allele at 5 Mb. At each of
5, 10, 20, 30, 40 and 50 kya, calculate the all-pair `frac_recent_T` and
carrier mass `n_recent_alt_alt / n_total_pairs`. The AF comparison is the
selected allele's sampled frequency. There is no REF/REF subtraction,
ALT/REF penalty, additional all-pair gate, or terminal-frequency filter.

Compare each selected-site score to the corresponding pre-existing neutral
**whole-region maximum** on the 10 kb grid. These references use the same
1,000 neutral regions, pair panel, demographic model and fixed rates. Their
input tree hashes are checked against the newly validated archive. Retain
the original five folds: 400 neutral calibration regions per fold, with
200 different neutral regions held out. New selection arms inherit the
existing selected replicate-to-fold assignments. The rank is
`p = (1 + count(neutral_maximum >= selected_site_score)) / 401`.
Report p <= 0.05 and p <= 0.01 separately at each T; do not select the best T
and treat its p-value as calibrated. This design cannot resolve p=0.001.

This endpoint is a **selected-site detection rate and a lower bound on
full regional power**, because a full scan can also find a signal away from
the selected site. It is not a completed blind regional scan of the new
arms. The selected site is an eligible archaic marker under this model and
lies exactly on the grid; fixed markers are retained by the existing marker
definition. The lower-bound interpretation assumes that the subsequent scan
uses this same marker definition, pair panel, grid and statistic. Neutral
call fractions still refer to the complete 10 Mb search. They are regional
false-positive rates, not a prevalence-independent discovery FDR.

Recompute the old neutral held-out call counts and require exact agreement
with the earlier analysis. Also require newly extracted s=0.005 focal truth
scores to agree with the saved earlier focal analysis. Wilson intervals
describe sampling uncertainty in the selected-site detection proportion;
they do not make the incomplete cohort or earlier method development
prospective validation. Report both all saved regions and the common set of
replicate IDs completed across every arm. Completion can depend on simulation
duration, so the paused array remains preliminary.

**Carrier ancestry diagnostic.** At the focal tree, partition carrier samples
by their ancestral branches immediately before the 50-kya cutoff. Two samples
share a group exactly when their TMRCA is strictly below 50 kya. Verify this
equivalence against every sampled pair. Save the number of groups, their
effective number `1 / sum(descendant_fraction^2)`, the largest group's share,
and the exact all-carrier-pairs recent fraction. These are local ancestral
lineages among sampled carriers, not a direct reconstruction of the number
of introgressing founders or distinct haplotypes over 10 Mb. Their relation
to carrier coalescence is algebraic; a correlation is not independent evidence
for a founder-diversity benefit.

The evaluation script writes per-region diagnostics, per-arm AF summaries,
selected-site predictions and intervals, neutral regional call fractions,
paired AF/carrier-mass comparisons, an audit and hash/provenance manifests.
Intermediate focal calculations are cached by input tree/carrier hashes,
pair manifest, relevant parameters and analysis-source content, independently
of unrelated Git commits. No new randomization is introduced.

**Reproduce the evaluation on the configured WSL host.** This command reads
existing simulation files and analysis references. The `simulation-status`
phase validates existing inputs without generating the 53 unfinished regions.

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
cd "$repo"
git pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
bash scripts/launch_eas_h400_array.sh simulation-status
bash scripts/evaluate_paused_eas_array.sh
```

Both launchers use the configured `uv` environment. Only validation and
evaluation run; the simulation array remains paused.
