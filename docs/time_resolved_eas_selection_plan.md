# Time-resolved all-pairs selection scan: EAS simulation plan

Status: full benchmark remains proposed; a bounded no-introgression pilot completed on 2026-11-02. See [measured pilot performance](results/time_resolved_eas_20261102/README.md) and [execution details](time_resolved_pilot_execution.md). The design below specifies 1-50 kya, 10-kb spatial scoring, and posterior-mass versus median scoring. The initial objective is to test whether an excess-coalescence interval recovers a selection episode from 30 to 20 kya. All times below are years before present unless stated otherwise.

## 1. Estimand and statistic

At every position x on the fixed 10-kb grid, average the posterior CDF over the fixed, genotype-unfiltered haplotype-pair panel:

`F_x(t) = sum_p w_p Pr(T_p(x) < t | sequence) / sum_p w_p`.

Use equal pair weights initially. For tree-sequence truth, replace the posterior probability with the indicator `T_p(x) < t`. The existing genomic null design uses 200 diploids and 10,000 fixed sampled haplotype pairs (pair seed 1729). Preserve this panel convention initially: "all pairs" means no focal-genotype filtering within that sampled panel, rather than claiming that all 79,800 unordered pairs of 400 haplotypes were decoded. Check exact enumeration against pair subsampling in a small truth-only sensitivity experiment.

Calculate CDF values at 1, 2, ..., 50 kya inclusive: 50 cutoffs and 49 one-kya intervals over [1,50). Also retain F(1 kya) as the below-1-kya diagnostic; this mass is not discarded or renormalized into the search interval. Probability older than 50 kya remains outside the searched intervals and in the pair denominator. A sliding cumulative cutoff alone cannot localize an episode: extra mass at 25 kya also raises F(30), F(40), and F(50).

For an interval I=[a,b), use interval mass `M_x(I) = F_x(b) - F_x(a)`. Compare this with the same quantity in matched neutral simulations. Report excess mass in percentage points, and optionally divide by interval width to display excess mass per kyr. Do not difference p-values or threshold posterior mean TMRCAs.

Primary discovery family: all 49 one-kya intervals with endpoints in [1,50] kya. Evaluate 2-, 5-, and 10-kyr windows stepping by 1 kyr as a prespecified secondary multiscale comparison. If their calls are combined with the primary calls, calibrate the combined family (182 time windows per position). This is a fine-grid time-resolved scan, not a claim of one-kyr biological resolution; posterior smoothing may make narrow intervals uninformative.

Score the entire region from the first pilot at zero-based positions x=0,10,000,... strictly below sequence length. A 10-Mb region has 1,000 positions, including the selected coordinate at 5 Mb: 50,000 position-by-cutoff summaries and 49,000 primary position-by-time-interval scores per method. The 10-kb stride specifies evaluation spacing, not a 10-kb averaging window. Any subsequent spatial smoothing must be separately frozen and calibrated. Report both the known-center temporal diagnostic and the blind spatial-by-time scan; only the latter establishes scan performance. Use exactly the same positions, pairs, validity masks, and simulated regions for every scoring method. Coordinates exported by native tools must be reconciled explicitly with this zero-based grid.

### Paired posterior-mass versus median comparison

For each pair p, cutoff t, and position x, retain `q_p,x(t)=Pr(T_p(x)<t | sequence)`. The mass arm contributes q itself to the pair sum. The median arm contributes `1[median(T_p(x)|sequence)<t]`; for the continuous Gamma posterior this equals `1[q_p,x(t)>0.5]`, with strict-cutoff ties handled consistently. The truth arm contributes `1[T_true,p(x)<t]`.

Thus a pair with q(50 kya)=0.10 contributes 0.10 to the cumulative pair sum under posterior mass, versus zero under the median rule. For an interval [20,30), a pair with q(20)=0.10 and q(30)=0.35 contributes 0.25, even if its posterior median lies outside the interval. Divide by the number of valid pairs, not by the sum of q values. Interval contributions are differences of the same pair's CDF, never differences between independently estimated medians.

This preserves each marginal posterior's probability mass in the statistic. The averaged mass remains a posterior expected fraction, not a complete joint posterior over the fraction or a confidence interval: dependencies between pairs and decoder/model misspecification still require whole-region simulation calibration. Use the same pairs valid at every cutoff for both arms. Median estimates must be actual medians, not posterior means mislabeled as medians. Do not reuse a median-calibrated neutral threshold for the mass arm.

## 2. Null calibration and interval calls

Pairs share genealogies; positions and time windows are dependent. The independent statistical unit is a complete simulated region/panel. Do not calculate binomial standard errors from the pair count or bootstrap individual pairs as independent observations.

Separate simulated regions into development, calibration, and final test cohorts. Keep all outputs of one replicate together. Prevent shared underlying trajectories across cohorts. When testing robustness to demographic uncertainty, also reserve unseen coherent history draws; repeated replicates within a history require clustered uncertainty summaries.

1. On development nulls, fit a position/interval-specific or prespecified matched-stratum center and scale for M_x(I), giving a one-sided excess score Z_x(I). Diagnose tied/zero-variance intervals and prespecify a rank-based alternative or exclude unidentifiable intervals; do not silently add a small denominator.
2. For each independent calibration null r, calculate `S_r = max_(x,I) Z_r(x,I)` using the same frozen scoring rule.
3. For a target's maximum `S_target`, calculate `p_global = (1 + count(S_r >= S_target)) / (B + 1)`.
4. An interval can be called using its simultaneous scan-adjusted value `p_adj(x,I) = (1 + count(S_r >= Z_target(x,I))) / (B + 1)` and a primary threshold of 0.05.
5. At each genomic position, form temporal support as the union of significant intervals, merging only overlap or adjacency. For spatial event reporting, group neighboring significant grid cells into components using a frozen adjacency rule; retain position-specific temporal support, contributing windows, scores, corrected p-values, and the unmerged cells for audit.

This controls the chance of any call in a matched entirely neutral region under the calibrated model. It does not make the union a confidence interval for selection onset/cessation, nor establish strong time-specific error control inside a selected genealogy. Peak location is the maximum standardized excess window's center; retain the raw excess profile so normalization effects are visible.

For the primary blind scan, use each neutral region's maximum over the full position-by-time family from the outset. Calibrate posterior-mass and median scoring separately at the same region-level alpha=0.05; any combined "either method" discovery rule requires calibration across methods too. The known-center diagnostic uses its own time-only correction and is labeled separately. Genome-wide discovery requires additional region-level multiplicity control or a genome-level calibration; per-region 0.05 is not genome-wide 0.05. Null matching must include sample/pair panel, decoder settings, masks, ascertainment, region/map design, mutation rate, and demographic model.

## 3. EAS simulation design

Reuse validated EAS PHLASH resources and existing forward-simulation/tree-sequence infrastructure. Distinguish the EAS Ne model from the separate AncientEurasia Han introgression model. Use 25 years/generation: 30 to 20 kya is 1,200 to 800 generations ago. Selection is zero outside this interval. Record the simulator's fitness convention for s and h, requested event times, realized rounded times, and the complete allele-frequency trajectory.

The first complete benchmark is a general selection scan under the EAS Ne history with no introgression or archaic-allele ascertainment in either selected or neutral simulations. Use a standing allele from a valid neutral prehistory, ascertained at selection onset, with onset frequencies near 1% and 5%; record realized frequencies and mutation age. Do not construct standing variation by arbitrarily changing present genomes. Add a de novo-at-onset sensitivity arm, since mutation age and selection onset can otherwise be confounded.

| Parameter | Proposed initial specification |
| --- | --- |
| EAS history | Median history for diagnostic pilot; coherent correlated PHLASH history draws for main robustness study |
| Selection episode | 30 to 20 kya |
| Time cutoffs | 1 through 50 kya inclusive, step 1 kya |
| Spatial scoring | Every 10,000 bp across the region, zero-based |
| Scoring arms | True-tree indicators, posterior-median indicators, posterior probability mass |
| Introgression | None in the initial benchmark |
| Selection coefficients | 0, 0.001, 0.005, 0.01; dominance h=0.5 under the documented simulator convention |
| Standing frequency at onset | Approximately 0.01 and 0.05, with a predefined ascertainment band |
| Sample/pair panel | 200 diploids; 10,000 genotype-unfiltered sampled pairs |
| Sequence | Existing 10-Mb region design, selected coordinate at center |
| Pilot maps/rates | Fixed recombination rate 1e-8 and mutation rate 1.29e-8 per bp/generation; main extension uses existing mapped regions and masks |
| Scaling | Unscaled timing-validation pilot; any accelerated simulation must pass a comparison before production |
| Randomness | Preserve genomic campaign base seed 20380101 and pair seed 1729; derive separate deterministic study/stage/replicate domains to avoid reusing old random streams |

Choose the allele and panel before observing the resulting signal. Retain losses and fixation in unconditional power; additionally report conditional power among alleles observed today. Do not impose a present-day frequency target in the main experiment. Stratify results afterward by final frequency and realized frequency rise. Neutral controls follow the identical onset ascertainment and observation rules.

Only after the no-introgression benchmark, consider an optional extension of selected and null designs in the existing EAS introgression context (2% pulse at 50 kya), selecting an eligible introgressed allele during 30-20 kya. Match any survival/present-observation conditioning in both classes and label the resulting power conditional. Existing introgressed, present-observed focal-site nulls are not automatically valid for an arbitrary-locus or no-introgression scan.

Temporal validation: hold out episodes such as 40-30 and 20-10 kya, plus a shorter 30-25 kya episode. Intervals at the 1/50-kya boundaries should be flagged as censored by the scan range. Add duration/strength combinations with similar cumulative allele-frequency change to test whether apparent duration is identifiable. Explore a wider diagnostic time range if signals fall beyond the primary boundaries; recalibrate any expanded discovery search.

## 4. Separate detection, overlap, and localization metrics

For replicate r, let truth `A_r=[20,30)` kya and inferred support be `D_r`, a union of called intervals. "Overlap" requires positive duration; merely touching an endpoint is not a hit.

| Metric | Definition and purpose |
| --- | --- |
| Detection power | Fraction of selected replicates with any call, irrespective of timing |
| Overlap TPR | Fraction of selected replicates with positive-duration `D_r intersect A_r`; this is the requested recovery criterion |
| Neutral FPR | Fraction of held-out neutral replicates with any call anywhere in the searched genomic-by-time grid |
| Any off-target call rate | Fraction of selected replicates with at least one reported component wholly outside A_r, even if another component overlaps |
| Entirely mistimed replicate rate | Fraction of selected replicates with nonempty D_r and no overlap with A_r |
| Event precision/recall/F1 | Match calls one-to-one to true episodes using positive overlap. Extra unmatched calls are FP, missed episodes are FN. Precision=TP/(TP+FP), recall=TP/(TP+FN), F1=2TP/(2TP+FP+FN). Report the selected:null composition used. |
| Temporal precision | Duration(D_r intersect A_r) / duration(D_r) |
| Temporal recall | Duration(D_r intersect A_r) / duration(A_r) |
| Temporal F1 | `2 duration(D_r intersect A_r) / (duration(D_r) + duration(A_r))` |
| Temporal FPR | Duration(D_r outside A_r) / duration(search range outside A_r) |
| Localization | Intersection-over-union, peak-center error, onset/cessation error, inferred width, and boundary censoring |

The user's outside-interval false positive is reported explicitly as an off-target/mistimed rate. It must not replace neutral FPR or be used as the false-discovery proportion in an F1 formula. A broad [1,50) call overlaps [20,30), but its temporal precision is 10/49 (about 0.204) and temporal F1 is 20/59 (about 0.339); overlap TPR alone would reward this poor localization.

For blind-scan event matching, require spatial as well as temporal recovery. Prespecify a spatial tolerance of 100 kb around the selected coordinate for the initial benchmark and report sensitivity to that tolerance. Match a component only if it contains a called cell within this tolerance whose time interval overlaps the true episode; do not combine a spatial hit at one time with an unrelated temporal hit elsewhere. Count distant unmatched components separately as spatial false calls. For the temporal metrics above, use support at the exact selected-coordinate grid point as the primary diagnostic; do not project calls across the entire 10-Mb region into one temporal interval. Report blind joint spatial-and-temporal TPR and exact-center temporal TPR separately.

For selected replicates with no call, temporal recall and temporal F1 are zero; precision and peak/boundary errors are undefined. Count misses separately rather than averaging them out. For neutral replicates, truth is empty; evaluate neutral FPR and called duration, not temporal recall. Report pooled event counts and macro-averaged temporal metrics separately. Bootstrap independent regions, or histories then regions for hierarchical designs; give binomial intervals for simple independent replicate rates.

## 5. Phased experiment and decision points

1. **Truth-only feasibility:** approximately 50 neutral and 50 selected simulations for the central episode, initially s=0.01. Verify selected-event scheduling, allele trajectories, and whether excess true TMRCA mass lies inside the episode. These are development diagnostics, not final error estimates.
2. **Decoded pilot:** decode the same no-introgression regions once at every 10 kb, retaining posterior parameters. Derive both posterior-mass and posterior-median scores at 1-50 kya without rerunning the HMM. Compare truth, median, and mass across the identical spatial/time grid. Establish temporal smoothing and whether secondary 2/5/10-kyr widths help the primary one-kyr scan. Compare with fixed F(50 kya) and max-over-cumulative-cutoffs baselines at matched empirical neutral FPR; calibrate each searched family separately. Report paired differences in detection, overlap TPR, neutral FPR, temporal F1, and timing error, with replicate-level uncertainty; do not assume posterior mass will necessarily improve performance.
3. **Frozen benchmark:** provisionally 500 development nulls, 1,000 calibration nulls, 1,000 independent test nulls, and 200 selected replicates per strength-by-onset-frequency cell. Existing compatible complete outputs may contribute only after contract and split checks. Budget from pilot timings before expanding. At B=1,000 the smallest empirical p is 1/1,001; this is not a genome-wide tail calibration. At true neutral FPR 0.05, 1,000 independent tests give standard error about 0.007. At power 0.5, 200 selected replicates give standard error about 0.035.
4. **Generalization:** unseen episode times/durations, coherent EAS histories, realistic region maps/masks, and sample-size sensitivity, initially all without introgression. Introgression is an optional later extension. Publish all planned cells and failures; never select only the first completed/favorable histories.
5. **Biological timing decision:** if detection is strong but coalescence timing systematically misses the selection episode, report an excess-coalescence epoch. A separate forward-simulation inference model is then needed to infer selection onset/cessation. Train that model on disjoint trajectories and validate interval coverage on held-out parameters before naming its output a selection-time estimate.

## 6. Implementation and execution contract

Extend the general all-pairs posterior reducer to evaluate the fine CDF grid and derive both median and probability-mass arms from each posterior once, at all 10-kb positions. Persist aggregate CDFs, interval masses, valid-pair counts, and method labels; retain raw posterior parameters until replay and integrity checks pass. Stream pair/position chunks and accumulate sums and valid counts rather than allocating a time-by-position-by-pair tensor. The existing `focused_selection_decode.py` demonstrates Gamma posterior-CDF evaluation, but its within-diploid panel and full probability tensor should not be copied as the new all-pairs execution design.

Proposed phases: plan -> simulate -> decode -> calibrate -> evaluate -> report. Proposed artifacts: manifest, seed ledger, truth/decoded time profiles, neutral calibration scores, called intervals, per-replicate metrics, aggregate metrics with uncertainty, and figure manifest. These phases/outputs describe future implementation; they are not currently runnable commands.

Use uv and default to four total workers/threads with a 26-GB memory budget; budget simulation workers times decoder threads to avoid oversubscription. Local development may use up to 24 total threads. Preserve the validated native simulator/decoder versions and Linux requirements from existing workflows. Add a SLURM wrapper only when targeting the cluster.

Cache by sample/pair manifest, scientific parameters, demographic/map/input hashes, relevant source contracts, and verified output integrity. Record the Git revision for provenance, but do not invalidate completed work solely because it changes. Fine-grid CDFs can be recomputed from retained posterior parameters; coarse CDF summaries alone cannot reconstruct missing thresholds, so missing raw posteriors may require re-decoding, not re-simulation.

Validate monotone CDFs, [0,1] bounds, nonnegative interval mass, consistent denominators/masks across times, exact pair identities, time-unit conversion, manifest splits, checksums, completed compressed streams, and realized selection timing. Use common validity masks across cutoffs. Do not silently clip substantive inconsistencies. Relevant tests should cover known synthetic truth/CDFs, scan calibration logic, interval merging/matching, and no-call/extra-call metric cases; no full-suite run is needed for this planning document.

Write plots to PNG and editable vector PDF only: large letter-page-readable text, horizontal legends above panels, no overlap, no notebook/inline plot display. Set `os.environ['MPLBACKEND']='Agg'` before Matplotlib imports; set PDF/PS fonttype 42 and prevent rasterization. Inspect rendered files locally. Suggested figures: allele trajectory plus truth/decoded excess profiles; replicate-by-time support heatmap; power and false-call rates by strength/frequency; temporal F1 and timing error by true episode.

## 7. Relation to Akbari et al.

The likely reference is [Akbari et al., Ancient DNA reveals pervasive directional selection across West Eurasia (Nature, 2026)](https://doi.org/10.1038/s41586-026-10358-1). Extended Data Fig. 9 uses 2,000-year windows sliding by 100 years to estimate time-varying polygenic selection intensity. Its observations are dated ancient allele frequencies with population-structure modeling. The analogy is the time-window strategy; it does not supply calibration or a timing interpretation for modern all-pairs TMRCA excess. See [the accompanying primary-source research note](akbari_time_resolved_selection_research.md).

Repository starting points: [genomic null contract](genomic_population_neutral.md), [existing EAS sweep methods](../no_introgression_EAS_sim/METHODS.md), and [posterior reducer](../python/gamma_smc_aou/focused_selection_decode.py).
