# EAS lab meeting figures - 21 September 2026

15 figures in slide-sized 16:9 layout. Each has a 3199 x 1800 PNG and a vector PDF with editable text. The combined PDF has one figure per page; full captions and suggested speaking points are below.

The figures use the completed 1,000 neutral + 1,000 selected array (100 selected per s). Gamma-SMC decoding is available only for s=0.005 and the neutral cohort. Simulation generation and new Gamma-SMC decoding were not started for this figure request.

Suggested main narrative: 01 -> 02 -> 03 -> 04 -> 05 -> 07 -> 08 -> 09 -> 10. Figure 06 verifies FPR calibration; 11-15 are backup figures.

All primary error rates are positional FPR, not FDR among discoveries. Tests at a pre-specified gene require the same statistic in a matching neutral interval. No maximum across unrelated positions or T cutoffs enters the primary p-values. The 100-kb iHS window is a separate endpoint.

## 01. Design and carrier score (main)

[PNG](figures/01_design_and_score.png) | [Vector PDF](figures/01_design_and_score.pdf)

**Takeaway:** Require the same haplotype pairs to carry the archaic allele and coalesce recently.

**Caption:** EAS; generation time 25 years; mutation rate 1.25e-8 and recombination rate 1e-8 per base per generation; dominance 0.5. Selected alleles must survive and be observed; fixation is retained. ALT/ALT denotes two haplotypes, not a diploid homozygote. The score uses frac_recent_T and no ALT/REF penalty or REF/REF subtraction.

## 02. Positional versus regional null (main)

[PNG](figures/02_matched_position_null.png) | [Vector PDF](figures/02_matched_position_null.pdf)

**Takeaway:** Only 4.3% of neutral focal positions reach 18% AF; unrelated peaks elsewhere do not enter the local test.

**Caption:** The focal score assigns the nearest eligible archaic marker within 5 kb of the fixed focal stride. All 1,000 neutral positions are retained; 169 have an assigned marker. Marker-free positions score zero. The whole-10-Mb maximum is shown only to explain the earlier calibration mismatch. A pre-specified gene requires calibration of the same summary in a matching neutral gene window.

## 03. Allele frequency by selection strength (main)

[PNG](figures/03_allele_frequency_distribution.png) | [Vector PDF](figures/03_allele_frequency_distribution.pdf)

**Takeaway:** Weak selection produces a broad AF distribution; its mean does not imply every replicate is detectable.

**Caption:** Boxes show interquartile range and median; whiskers show the 5th and 95th percentiles; diamonds show means. Selected data are conditional on survival and sample observation. Neutral data are unconditional pre-specified positions, with absent nearby markers scored zero. These are not frequency-matched cohorts.

## 04. Power across selection coefficients (main)

[PNG](figures/04_positional_power.png) | [Vector PDF](figures/04_positional_power.pdf)

**Takeaway:** Carrier-based local tests are much more sensitive than all-pair recency at weak selection.

**Caption:** All TMRCA curves use true TMRCA. Each position is compared with 400 neutral calibration positions in one of five folds; 200 neutral positions per fold are held out. All 1,000 selected and 1,000 neutrals are included. iHS uses the nearest scorable core within 5 kb. Wilson bands do not include uncertainty from estimating the neutral reference or cross-fold dependence.

## 05. Power across TMRCA cutoffs (main)

[PNG](figures/05_time_cutoff_power.png) | [Vector PDF](figures/05_time_cutoff_power.pdf)

**Takeaway:** Shorter cutoffs can add discrimination beyond frequency, particularly at s=0.001.

**Caption:** Cutoffs are 5, 10, 20, 30, 40 and 50 kya; there is no per-replicate maximum over cutoffs. Each cell uses 100 selected replicates and the matching neutral statistic. Comparing cutoffs is exploratory, not an independently validated tuning procedure.

## 06. Neutral FPR calibration (main)

[PNG](figures/06_positional_false_positive_rate.png) | [Vector PDF](figures/06_positional_false_positive_rate.pdf)

**Takeaway:** Matched nulls control local false calls close to the nominal level for both truth and decoding.

**Caption:** Cell labels are percentages, with separate color scales for the 5% and 1% panels. Neutral marker-free positions contribute zero carrier mass. The false discovery rate among selected discoveries also depends on the prevalence of selection and is not estimated by this figure. No chromosome-maximum error rate is substituted for the positional rate.

## 07. Truth versus decoded TMRCA (main)

[PNG](figures/07_truth_vs_decoding.png) | [Vector PDF](figures/07_truth_vs_decoding.pdf)

**Takeaway:** At T=50 kya and p<=0.05, all-pair power falls from 86% to 61%; carrier mass remains at 100%.

**Caption:** This comparison is available only for the already decoded s=0.005 cohort. Both sources use the same 10,000 sampled pairs and local marker assignments. Decoded frac_recent_T thresholds posterior-mean TMRCA per pair; it does not average posterior mass. Weak-selection arms have not been decoded. Matching nulls calibrates the test but does not remove decoding-related loss of discrimination.

## 08. Incremental gain from coalescence (main)

[PNG](figures/08_gain_beyond_allele_frequency.png) | [Vector PDF](figures/08_gain_beyond_allele_frequency.pdf)

**Takeaway:** The added coalescence signal is clearest at a shorter cutoff in the weakest selection arm.

**Caption:** Left: power for pre-specified AF and carrier-mass scores. Right: within-replicate gains and losses relative to AF at s=0.001; negative bars are losses, not negative probabilities. Null FPR is approximately 5% for each method but not identical. Paired p-values and all cutoffs are in paired_af_comparison.csv; they are unadjusted for examining several cutoffs and arms.

## 09. iHS power and target-site eligibility (main)

[PNG](figures/09_ihs_and_sweep_completion.png) | [Vector PDF](figures/09_ihs_and_sweep_completion.pdf)

**Takeaway:** Near fixation, the selected SNP is often unscorable by iHS; nearby SNPs may still carry signal.

**Caption:** The two iHS curves represent separately calibrated endpoints. Positional iHS uses the nearest finite-scoring core with MAF>=5% within 5 kb; the fixed 100-kb endpoint uses the fraction of finite-scoring cores with |standardized iHS|>2, requiring at least 20 cores. Right: exact focal-SNP eligibility, not eligibility of the fallback positional statistic. Fixation and missing/undefined scores never remove replicates from the denominator.

## 10. Detection and FPR across the region (main)

[PNG](figures/10_spatial_decay.png) | [Vector PDF](figures/10_spatial_decay.pdf)

**Takeaway:** Strong local carrier signal extends into linked sequence while neutral per-position FPR stays near 5%.

**Caption:** Each point is the fraction of 100 selected or 1,000 neutral replicates called at that coordinate. No smoothing or maximum over positions is used. Detection at linked positions describes the extent of the signal, not a separate causal target or a localization false discovery. Curves are restricted to +/-500 kb for readability; full 10-Mb pointwise results are in the source table.

## 11. Example neutral and selected regions (backup)

[PNG](figures/11_example_regions.png) | [Vector PDF](figures/11_example_regions.pdf)

**Takeaway:** The gene or position defines the scope of the null; peaks outside it are irrelevant.

**Caption:** Neutral example neutral/rep0388: choose focal AF<5% and whole-region maximum AF>=40%, then the maximum nearest the median of eligible examples (ties by ID). Selected example 0014: AF nearest the s=0.005 median (ties by replicate). Gold shading marks a fixed 100-kb interval centered at 5 Mb. AF and carrier mass have different scales of interpretation despite both lying in [0,1].

## 12. Power at stricter thresholds (backup)

[PNG](figures/12_threshold_stringency.png) | [Vector PDF](figures/12_threshold_stringency.pdf)

**Takeaway:** The s=0.005 signal remains strong at stringent ranks; very weak selection loses substantial power.

**Caption:** p<=0.05 and 0.01 use 400 neutral calibration positions per fold. The separate p<=0.001 sensitivity uses all 1,000 neutrals for selected ranks and leave-one-out ranks for neutrals. One extreme neutral rank gives 1/1,000 calls for the shown scores; this is not precise independent validation of a 0.1% tail. Comparisons remain positional, not family-wise over a genome.

## 13. Why AF and carrier mass can be similar (backup)

[PNG](figures/13_carrier_genealogies.png) | [Vector PDF](figures/13_carrier_genealogies.pdf)

**Takeaway:** When most carriers share recent ancestry by T=50 kya, weighting by coalescence adds less beyond AF.

**Caption:** Carrier mass equals the sampled ALT/ALT pair fraction times within-carrier frac_recent_T. AF squared approximates the first factor for random distinct pairs, with finite-panel and pair-sampling differences; it is not an exact upper bound. Branches are the sampled-carrier ancestral groups immediately below the 50-kya cutoff, verified against pairwise TMRCA, not a census of introgressing founders. The model has a strong archaic bottleneck.

## 14. Demography and decoder assumptions (backup)

[PNG](figures/14_demography_and_assumptions.png) | [Vector PDF](figures/14_demography_and_assumptions.pdf)

**Takeaway:** Matched nulls address calibration under the model; they do not establish robustness to every model misspecification.

**Caption:** The EAS trajectory is read from the same PHLASH-derived model used by the archived simulations. The Ne history is displayed only to 50 kya; simulation genealogies also contain the older archaic branch. Decoder Ne is theta/(4*mu)=15,000. Neutral and selected simulations share the archaic bottleneck and rate assumptions. No variable-Ne Gamma-SMC flow field or empirical recombination-map benchmark is claimed here.

## 15. Pair-class coalescence diagnostic (backup)

[PNG](figures/15_pair_class_coalescence.png) | [Vector PDF](figures/15_pair_class_coalescence.pdf)

**Takeaway:** The carrier score combines this within-ALT/ALT recency with the fraction of the complete pair panel that is ALT/ALT.

**Caption:** ALT/ALT and REF/REF are haplotype-pair classes at the focal archaic marker. Only regions with an assigned marker contribute to this within-class diagnostic: 100 selected and 169 neutral. Counts are pooled across pairs, so regions with larger classes contribute more weight; this is not a replicate-average or the unconditional neutral denominator used for FPR. REF/REF is a diagnostic and is not subtracted from the carrier score; ALT/REF is not scored.

## Interpretation limits

Selected power is conditional on allele survival and observation. The unconditional positional null includes marker-free positions; it is not a null conditional on first selecting an observed archaic SNP. The archaic-marker oracle is more restrictive than simply lying on an introgressed segment and is not a literal finite-reference Neanderthal/Denisovan/African ascertainment experiment.

This array uses a strong shared archaic bottleneck, uniform mutation/recombination rates and isolated 10-Mb regions. It does not demonstrate genome-wide FDR control, causal localization, or a gain from multiple independent introgressing haplotypes. Cutoff comparisons are exploratory. Wilson intervals are conditional on the estimated null and do not include its Monte Carlo uncertainty.

## Methods references

Voight et al. (2006), iHS: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072

scikit-allel implementation: https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/selection.html
