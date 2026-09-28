# Evidence audit for the introgression and selection manuscript

Reviewed 28 September 2026. Companion to `introgression_selection_draft.md`.

## Scope and source hierarchy

The workspace root is not itself a Git repository. It contains two scientific codebases, multiple checkouts, local supplied-data reviews, public reference annotations, test outputs, and extensive selection-simulation outputs. The review used repository-wide file inventories and targeted reading of the empirical preparation, inference, annotation, analysis, and reporting code and documentation. It did not read every generated simulation replicate, dependency file, or binary byte as manuscript evidence. Selection simulations, synthetic fixtures, and performance benchmarks were excluded from biological Results.

| Checkout | Repository / branch | Scientific source revision reviewed |
|---|---|---|
| `code` | `soisa001/gamma_smc_ts`, `AOU_run_opt` | `eba1c82` |
| `hmmix_aou_postdecode` | `soisa001/hmmix_aou`, `refactor/postdecode-analysis` | `aa5c684` |
| `hmmix_aou_plan_inspection` | `soisa001/hmmix_aou`, `optimize/fastcore-prep_hardmask` | `09cbc34` |
| `hmmix_aou_linuxtest` | Linux validation checkout of the HMMIX repository | `aa5c684` |
| `ci-native-fix` | Gamma-SMC build/CI checkout | `37d27cb` |

The last two checkouts provide implementation/validation context, not independent empirical cohorts. The source revision of this draft is separate from the revision used to generate a historical controlled-data output. A present-day default cannot establish a historical run parameter.

Evidence was prioritized as follows: supplied aggregate tables with checkable records and checksums; their local annotation/review records; current executable code; current run documentation; historical narrative or template reports. Missing production records were not filled with benchmark or test values. No controlled cloud data were accessed, no analysis was run on Workbench, and no figures were displayed.

## Results that can be supported locally

| Draft result | Local source, relative to workspace root | Verification and limitations |
|---|---|---|
| 38,219 exact alleles; INS/DEL counts, singleton fraction, lengths, annotation-field counts | `hmmix_aou_plan_inspection/resources/loci/ad_hoc_sv_20260921.tsv.gz` | Full gzip read, schema/row validation, unique event IDs, direct tally. Compact export has no called-allele denominators or complete run provenance. Annotation fields were counted, not independently reannotated for all 38,219 alleles. |
| 4,706 ranking rows; 1,096 alleles; 626 original record intervals; 48 strata | `sv_exclusive_review_20260915/summary.json` and the supplied attachment path recorded in that JSON | Original input checksum and row counts checked when accessible. This is a selected top-ranking set, not a complete catalogue. |
| Functional classes, 601 callable anchors, 212 callable/unflagged alleles, candidate carrier counts | `sv_exclusive_review_20260915/annotated_unique_events.tsv`, `priority_candidates.tsv`, `overlap_details.json`, `refseq_crosscheck.json` | Annotated table SHA-256 checked against the existing review manifest. Feature totals recomputed. Candidate consequences follow the existing GENCODE v47/RefSeq review; full alleles were not reconstructed here. |
| Intersection of the two SV exports: 601 exact event IDs | The two SV tables above | Exact-ID set intersection. This numerical equality with the callable subset does not by itself establish identical inclusion criteria or a common production run. |
| LCT EM results in 12,238 individuals | `lct_dprime_supplied_counts.tsv`; `hmmix_aou_plan_inspection/docs/lct_panel_ld.md` | Population counts summed; D-prime and r-squared values read directly; haplotype-frequency normalization and LD formula consistency checked. Original joint genotype counts are not retained in this TSV, so the EM fit was not rerun. This is optional contextual evidence. |
| Historical selection leads only | `hmmix_aou_postdecode/resources/selection/gamma_smc_selected_regions.placeholder.tsv` | 17 rows, all ranked by `frac_recent_4500`; 1 AFR, 4 AMR, 3 EAS, 5 EUR, 3 MID, 1 SAS. Explicit placeholder and display-cluster output, excluded from quantitative 50,000-year Results. |

`verified_summary.json` records recalculated counts, input hashes, integrity checks, and source-revision information. The full local file inventory and audit-generation script are kept in `manuscript_review_20260928/` outside the Git repository. Only the aggregate evidence and manuscript files are published with this draft.

## Methods evidence map

Paths below are relative to the named checkout. They point to implementation, not proof of completed execution.

| Topic | Principal source files |
|---|---|
| Observation preparation and ancestry labels | HMMIX `src/aou_hmmix/prep/cli.py`, `prep/fastparse.py`, `prep/manifest.py`, `observations.py`, `analysis_gate.py` |
| Cohort and denominator | HMMIX `docs/postdecode.md`, `src/aou_hmmix/postdecode_io.py`, `postdecode_tasks.py`, `analysis/cohort.py`, `analysis/population_contract.py` |
| Training, HMM grid and alpha | HMMIX `bench/run_end_to_end_12k.sh`, `src/aou_hmmix/allsamples.py`, `fused_tasks.py`, `fastcore/grid.py`, `artemis.py`, `policy.py` |
| AFR/OTH model treatment and exceptions | HMMIX `src/aou_hmmix/decode_african.py`, `oth_strategy.py`, `donor_overrides.py`, `decoded_snapshot.py` |
| Posterior and lineage-affinity rules | HMMIX `src/aou_hmmix/introgression.py`, `segment_annotate.py`, `supported_segments.py` |
| Masks, burden, saturation, bins and deserts | HMMIX `src/aou_hmmix/masking.py`, `analysis/callable_introgression.py`, `analysis/regional_introgression.py`, `analysis/phase12b_landscape.py` |
| Genic and functional tests | HMMIX `src/aou_hmmix/analysis/permutation.py`, `frequency.py`, `constraint.py`, `high_introgression_enrichment.py` |
| Exact-allele SV semantics and statistics | HMMIX `src/aou_hmmix/analysis/sv_panel.py`, `sv_events.py`, `sv.py`, `variant_introgression.py`, `sv_source_classification.py`, `sv_genotypes.py`, `phase13.py` |
| Candidate annotations and local review | HMMIX inspection checkout `docs/ad_hoc_sv_loci.md`, `docs/phase13.md`; workspace `sv_exclusive_review_20260915/methods_review.md`, `results_review.md` |
| Empirical selection input/defaults | Gamma-SMC `scripts/run_aou_workbench.sh`, `WORKBENCH_RUNNER.md`, `python/gamma_smc_aou/decoder.py`, `workbench.py`, `workbench_candidates.py` |
| Pair sampling, time units, and posterior summary | Gamma-SMC `src/pair_sampling.h`, `pair_manifest.h`, `gamma_smc.cpp`, `gamma_smc.h`, `recent_stats.h` |
| Integrated selection/introgression test | HMMIX `src/aou_hmmix/analysis/selected_region_enrichment.py` |
| LCT two-locus estimation | HMMIX inspection checkout `docs/lct_panel_ld.md`, `scripts/lct_panel_ld.py` |

Auxiliary implemented analyses include AFR-only ADMIXTURE, PCA projection, gene-constraint correlation, tract-length/affinity mixtures, SHARED/LOW_EVIDENCE diagnostics, and curated loci including FOXP2, GHR, BNC2, KRT, and TMEM106B. These were considered in the review, but no validated full-cohort results for them were found locally. They were not turned into completed biological findings. Mixture components were not equated with introgression pulses. Public LCT sequence notes were not treated as panel-specific results. The TMEM106B and ad-hoc plotting directories contain synthetic QA examples and cannot establish allele frequencies or cohort haplotype structure.

## Resolved discrepancies and remaining decisions

1. **Requested selection statistic.** The draft uses mean posterior P(TMRCA < 50,000 years), not a hard threshold applied to posterior mean TMRCA. Under the verified defaults, 50,000 years = 2,000 generations and 0.0688 coalescent units. Both the posterior mean and the posterior CDF depend on the gamma parameters, but thresholding one is not equivalent to integrating the other.
2. **Runner defaults do not implement the requested candidate screen automatically.** The Workbench runner currently defaults to 4,500 years, `recent_call=mean`, and candidates with `frac_recent_4500 > 0.02`. Supplying `--threshold-years 50000` changes the threshold but does not convert the candidate selector into a posterior-probability selector. The manuscript request did not authorize changing the analysis implementation; no decoder or candidate-selection code was altered. A final 50k CDF-based rule and matching outputs remain necessary.
3. **Direction of evidence.** A recent-coalescence excess corresponds to larger P(TMRCA < 50,000), not smaller values. The HMMIX selected-region module and old report generator contain misleading narrative examples referring to low P(TMRCA < 4500). The actual overlay consumes intervals and does not establish that direction. The draft corrects the scientific wording without claiming a production rerun.
4. **Pair design.** The current empirical runner samples 100,000 distinct unordered within-population haplotype pairs, including same-individual pairs by default. The initial `AOU_WORKFLOW.md` introduction describes the older within-individual-only design. These cannot be conflated.
5. **HMMIX alpha policy.** The production launcher calls the all-sample path with sample-specific Artemis alpha, not the older 100-calibrator population-median design. The README retains both descriptions. The actual exported methods/parameter records must verify which policy generated the final snapshot.
6. **Population labels and exclusions.** `ancestry_pred_other` supplies manifest labels; `ancestry_pred` is the forced canonical prediction used for OTH training safety. Some report-generator wording incorrectly describes canonical predictions as the general cohort labels. Introgression QC/relatedness exclusions are optional in the original launcher and preserved by the snapshot. Gamma-SMC excludes both lists by default. Do not call the introgression cohort unrelated or claim matching cohorts without the manifests.
7. **Primary mask.** Current introgression code uses the positive strict callable mask as primary. The old August PDF described the hardmask complement as primary. The Gamma-SMC runner still defaults to the hardmask complement and offers strict mode separately. The manuscript distinguishes these; the final integrated analysis must document compatibility.
8. **OTH migration.** The latest described snapshot has mixed self-training/EUR-donor treatment and one recorded AMR donor override. Older documents describe donor decoding for every OTH individual. The single override must remain in provenance but participant identifiers need not appear in the manuscript.
9. **Historical SV denominators.** The selected ranking's 24,476 haplotypes do not prove all genotypes were observed at every event. Complete-genotype sidecars introduced later provide more appropriate AC/AN estimates but were not present with that historical ranking. Do not retrofit them into the reported counts.
10. **Genic permutation unit.** The code's `collapse` option removes duplicate coordinate intervals; its default is off. No production genic P values are claimed until the input unit, settings, and results are available. The same applies to secondary enrichment, mixture, and constraint analyses.
11. **Reference annotations.** The selected-SV review used GENCODE v47; regional maps default to a pinned UCSC RefSeq resource; the Gamma-SMC runner stages GENCODE v50 basic. Different gene labels across products need not reflect a biological disagreement. Annotation overlap alone does not identify an SV's effect or the target of selection.
12. **Selection evidence status.** Posterior probabilities are not frequentist P values. Neither the placeholder locus table nor a gene's known selection history validates a 50,000-year empirical result. Adaptive-introgression claims require evidence connecting the selected haplotype/allele to archaic ancestry beyond regional overlap.

## Production fields needed to finish the Results

These are missing manuscript inputs, not requests to rerun completed work. Prefer already-computed, approved aggregate exports with the run identifier and settings.

| Required information | Expected artifact(s) or equivalent aggregate export |
|---|---|
| Final numbers included, excluded, failed, fallback and zero-tract, by population and model treatment | Snapshot `cohort_counts.tsv`, `cohort.tsv` aggregate summaries, `methods.json`, `parameters.tsv`, decoder counts |
| Tract counts, source partition, callable genome, burden mean/median/IQR, union and saturation | `callable_introgression_summary.tsv`, `regional_introgression_summary.json`, haplotype-summary aggregates and saturation endpoint |
| High-region/desert counts, principal loci and exact thresholds | `high_introgression_regions.tsv`, desert tables, regional summary and invocation parameters |
| Functional/genetic enrichment effect sizes and corrections | Functional-track manifest, high-introgression, genic, block-permutation and constraint summaries with tested families |
| Final SV totals and exact AF distributions | Complete event/catalogue summaries, source-classification summary and complete called-allele denominator sidecars |
| The 50,000-year scan | Scan summary explicitly recording threshold 50,000, generation time, mutation/time scale, pair/sample identities, mask, posterior-probability column, and valid pair counts |
| The 50,000-year candidate definition | Exact threshold/ranking and eligibility rule, unmerged scan windows, candidate intervals, and mask sensitivity summary |
| Integrated result | Matching population/callability definitions, verified candidate input, selected-region enrichment summary and source-specific allele follow-up |

The oldest PDFs under `hmmix_aou_postdecode/results/paper/` were inspected as text. `scientific_results.pdf` contains six pages and `detailed_methods.pdf` eight; both identify themselves as TEMPLATE and were generated on 5 August 2026. Their manifest marks the expected production artifacts absent. They are not numerical evidence and contain superseded method descriptions. No replacement PDF was generated, because editable Markdown preserves equations and explicitly unresolved fields without making the draft look final.

## Suggested use

Use the main Methods as a detailed working section, subject to the recorded run values. The two SV Results subsections are supported by local supplied summaries, with their export-specific scope preserved. The introgression-landscape and 50k-selection Results paragraphs are structured fill-in text. The LCT paragraph is optional supporting material. Before journal submission, remove editorial notes only after replacing the corresponding fields and reconciling the run-specific choices above.
