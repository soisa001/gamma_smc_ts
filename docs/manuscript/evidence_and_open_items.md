# Evidence audit for the demography, introgression, and selection manuscript

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

## PHLASH addition: sources and independently verified results

The demographic review extends beyond the original workspace inventory to two local checkouts of `git@github.com:soisa001/phlash_ld.git`:

- `C:/Users/Lenovo/OneDrive/Documents/phlash_ld_opt_audit`, reviewed at `70e177ea0e77fe38d408eb8042e13a44c94d15c9`, contains the six published `artifacts/aou-mvn/<POP>.npz` files and adjacent JSON provenance. They were introduced by commit `7130ff5`. This is the primary artifact and implementation source for the addition.
- `C:/Users/Lenovo/OneDrive/Documents/phlash_run`, reviewed at `ec38cb256e185a16d03acf0512eecb5dd5f23795`, provides a second implementation checkout and newer input/cache documentation. Its defaults must not be substituted for those of a historical fit.

The method paper is [Terhorst (2025), Nature Genetics](https://doi.org/10.1038/s41588-025-02323-x). Study-specific LD extension and preprocessing details come from the local source, principally `scripts/aou_pipeline.py`, `scripts/run_aou_ld_bootstrap.sh`, `scripts/aou_mvn.py`, `scripts/fit_aou_mvn.py`, `src/phlash/mcmc.py`, `src/phlash/data.py`, and `src/phlash/ld/data.py`. The declared package version is 1.0.6; the repository revision also matters because this is a modified implementation.

`summarize_phlash.py` directly re-read all six NPZ/JSON pairs using `allow_pickle=False`, checked ZIP CRCs and the sidecar's SHA-256/size, validated schema/population/dimensions and all finite positive histories, verified the exact shared float32 geometric grid, checked unique `fit000`--`fit099` source names, and reconstructed the saved log-space mean and covariance factor. It then calculated the pointwise median and linear 2.5th/97.5th percentiles across the 100 fit-level posterior medians. It did not deserialize original fit pickles, fit a model, access controlled buckets, or sample new demographies. The four-thread run used NumPy 2.4.6 and wrote `phlash_summary.json`; source NPZ and sidecar checksums are recorded there. No participant identifiers are included.

The resulting Results table uses the nearest stored grid point to 2,000 generations, at 1,999.40087890625 generations, and rounds N_e to whole individuals. Other point summaries are retained at the nearest grid points to 100, 500, 1,000, 10,000, and 40,000 generations. Curve extrema are taken over the entire archived grid and do not estimate uncertainty in event timing. The main table is checked against this aggregate JSON by `verify_draft.py`; the portable verifier does not itself reopen the external source artifacts.

The existing `artifacts/aou-mvn/plots` figures instead summarize 1,000 log-MVN draws produced with a downstream simulation seed schedule. Their bands are not the empirical percentiles across 100 fit-level histories. Neither those plots nor the simulated validation reports supply the demographic numbers newly added to this manuscript. No selection-simulation methods or results were added.

### What the aggregate handoff establishes, and what remains unresolved

| Item | Evidence / manuscript treatment |
|---|---|
| Six populations; 100 fit-level histories each; 10,000 points over 100--40,000 generations | Verified directly in all six artifacts and sidecars. |
| Log-space mean/covariance representation; zero jitter; covariance rank at most 99 | Verified from arrays and reconstructed from the archived curves. |
| Original input fits | Sidecars record `phlash_bootstrap_trial_1mb_2000i_200p/bootstrap/fits/<POP>/fit000.pkl` through `fit099.pkl`, with sizes/timestamps. They do not store fit-file hashes, original per-fit metadata, or sample manifests; source names alone do not establish run parameters. |
| Mutation-rate scaling | Optimized launcher explicitly defaults to 1.29e-8; generic caller defaults to 1.25e-8. The actual value used for the archived fits remains unconfirmed. No rescaling was applied to the archived curves. |
| Mask choice | The artifact checkout's executable generic caller defaults `--hardmask` to enabled and its launcher explicitly passes `--hardmask`. Its prose documentation describes disabled-by-default `--masking`, as does the newer checkout. Executable source takes precedence for the reviewed launcher, but the historical invocation is still required. |
| Sample selection and seed | Code implements 100 unique samples per replicate, 100 replicates/population, repeated without-replacement subsampling with pool reset. Launcher seed defaults to 20260724. Historical sample counts, eligible pool sizes, and persisted seeds are absent from the handoff. This is not a claim of 10,000 distinct participants per population. |
| HMM/AFS/LD fit settings | Draft records reviewed launcher settings, including LD weight 20, 200 particles, at most 2,000 iterations, 100-bp windows, and 1,000-window chunks plus warm-up. Confirm historical LD enablement, weights, HMM geometry, priors/parameterization, map version, and convergence from per-fit metadata and diagnostics. |
| Sequence filtering and denominators | Code defines coding-gene/flank exclusions, chromosome trimming, mask rules and contig length thresholds. Retrieve `run_sanity.json`, `viable_contigs.tsv/json`, annotation/map hashes and retained-base totals for the actual fits. |
| Uncertainty | Reported bands are between-fit quantiles of posterior medians. Within-fit posterior dispersion, sample-overlap dependence, and optimizer variability are not separately recoverable from the handoff. No calibrated coverage claim or simultaneous confidence band is made. |
| Calendar time | Archive time is in generations. The draft gives the conditional conversion at 25 years/generation to relate 2,000 generations to the Gamma-SMC horizon. It does not establish a historical PHLASH plotting convention. |
| Empirical model fit | Original AFS, LD, objective, held-out predictive-density and stopping diagnostics are needed to assess convergence and adequacy. Checksum/covariance checks are integrity tests, not demographic validation. |

The original introgression/selection audit above retains its earlier source scope; this addition does not claim to refresh those analyses from later repository outputs.
