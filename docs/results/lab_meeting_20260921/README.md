# EAS lab meeting figures - 21 September 2026

15 figures in slide-sized 16:9 layout. Each has a 3199 x 1800 PNG and a vector PDF with editable text. The combined PDF has one figure per page; full captions and suggested speaking points are below.

All primary TMRCA plots use Gamma-SMC decoded frac_recent_T. Decoding is available for 1,000 neutrals and 100 selected regions at s=0.005, so decoded power is restricted to that coefficient. AF and iHS use observed simulated genotypes and retain all 1,000 selected regions (100 per s). Figure 07 alone includes true TMRCA as an explicitly labeled validation comparison. No simulations or new decoding were started for this revision.

Suggested main narrative: 01 -> 02 -> 03 -> 04 -> 05 -> 06 -> 08 -> 09 -> 10. Figure 07 is truth-versus-decoding validation; 11-15 are backup figures. Figures 03_a and 03_b remain deferred because actual trajectories were not recorded and replays were declined.

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

## 04. Decoded power at s=0.005 (main)

[PNG](figures/04_positional_power.png) | [Vector PDF](figures/04_positional_power.pdf)

**Takeaway:** At nominal 5%, decoded carrier power is 100%, all-pair power 61%, and positional iHS power 67%.

**Caption:** Carrier mass and all-pair recency use Gamma-SMC decoded posterior-mean TMRCA hard calls. AF and iHS use the observed simulated genotypes. All four methods are restricted to the same 100 s=0.005 selected regions; the other selection arms have not been decoded. Each test position uses 400 neutral calibration positions in its fold, with 200 held-out neutral positions. Wilson intervals do not include fitted-null uncertainty or cross-fold dependence.

## 05. Decoded power across TMRCA cutoffs (main)

[PNG](figures/05_time_cutoff_power.png) | [Vector PDF](figures/05_time_cutoff_power.pdf)

**Takeaway:** Decoded carrier mass has 100% power at nominal 5% across the tested cutoffs in the s=0.005 cohort.

**Caption:** All TMRCA scores are decoded. Cutoffs are 5, 10, 20, 30, 40 and 50 kya, each calibrated separately against the matching decoded neutral statistic. No per-replicate maximum over cutoffs is used. Only s=0.005 has selected decoding; the full s grid is not represented. Cutoff comparisons are exploratory.

## 06. Neutral FPR calibration (main)

[PNG](figures/06_positional_false_positive_rate.png) | [Vector PDF](figures/06_positional_false_positive_rate.pdf)

**Takeaway:** Matched decoded nulls keep local false calls close to the nominal level.

**Caption:** Both rows use Gamma-SMC decoded statistics. Cell labels are percentages, with separate color scales for the 5% and 1% panels. Neutral marker-free positions contribute zero carrier mass. Discovery FDR also depends on the prevalence of selection and is not estimated here. No chromosome-maximum error rate is substituted for positional FPR.

## 07. Truth versus decoded TMRCA (validation)

[PNG](figures/07_truth_vs_decoding.png) | [Vector PDF](figures/07_truth_vs_decoding.pdf)

**Takeaway:** At T=50 kya and p<=0.05, all-pair power falls from 86% to 61%; carrier mass remains at 100%.

**Caption:** This comparison is available only for the already decoded s=0.005 cohort. Both sources use the same 10,000 sampled pairs and local marker assignments. Decoded frac_recent_T thresholds posterior-mean TMRCA per pair; it does not average posterior mass. Weak-selection arms have not been decoded. Matching nulls calibrates the test but does not remove decoding-related loss of discrimination.

## 08. Decoded carrier mass versus allele frequency (main)

[PNG](figures/08_gain_beyond_allele_frequency.png) | [Vector PDF](figures/08_gain_beyond_allele_frequency.pdf)

**Takeaway:** The decoded s=0.005 cohort is near the power ceiling for both AF and carrier mass; weaker-selection gains remain untested.

**Caption:** The left panel has a zoomed 90-100% power axis. At nominal 1%, AF detects 99/100 and decoded carrier mass detects 97-99/100 across cutoffs. Right: paired gains and losses, with negative bars representing AF-only detections. At nominal 5%, both methods detect all 100 regions at every cutoff. These near-ceiling results do not establish a gain from decoded coalescence at weaker selection, which has not been decoded. Comparisons across cutoffs are exploratory.

## 09. iHS power and target-site eligibility (main)

[PNG](figures/09_ihs_and_sweep_completion.png) | [Vector PDF](figures/09_ihs_and_sweep_completion.pdf)

**Takeaway:** Near fixation, the selected SNP is often unscorable by iHS; nearby SNPs may still carry signal.

**Caption:** The two iHS curves represent separately calibrated endpoints. Positional iHS uses the nearest finite-scoring core with MAF>=5% within 5 kb; the fixed 100-kb endpoint uses the fraction of finite-scoring cores with |standardized iHS|>2, requiring at least 20 cores. Right: exact focal-SNP eligibility, not eligibility of the fallback positional statistic. Fixation and missing/undefined scores never remove replicates from the denominator.

## 10. Detection and FPR across the region (main)

[PNG](figures/10_spatial_decay.png) | [Vector PDF](figures/10_spatial_decay.pdf)

**Takeaway:** Decoded carrier signal extends into linked sequence while neutral per-position FPR stays near its nominal level.

**Caption:** TMRCA curves use Gamma-SMC decoded scores at both nominal thresholds; AF uses genotypes. Each point is the fraction of 100 selected or 1,000 neutral replicates called at that coordinate. No smoothing or maximum over positions is used. Detection at linked positions describes signal extent, not a separate causal target or a localization false discovery. The FPR panels have different scales matching the nominal thresholds. Full 10-Mb pointwise results accompany the plotted +/-500-kb interval.

## 11. Example neutral and selected regions (backup)

[PNG](figures/11_example_regions.png) | [Vector PDF](figures/11_example_regions.pdf)

**Takeaway:** The gene or position defines the scope of the null; peaks outside it are irrelevant.

**Caption:** Neutral example neutral/rep0388: choose focal AF<5% and whole-region maximum AF>=40%, then the maximum nearest the median of eligible examples (ties by ID). Selected example 0014: AF nearest the s=0.005 median (ties by replicate). Gold shading marks a fixed 100-kb interval centered at 5 Mb. AF and carrier mass have different scales of interpretation despite both lying in [0,1].

## 12. Decoded power at stricter thresholds (backup)

[PNG](figures/12_threshold_stringency.png) | [Vector PDF](figures/12_threshold_stringency.pdf)

**Takeaway:** At the pooled p<=0.001 threshold, decoded carrier mass detects 92%, AF 90%, and decoded all-pair recency 23%.

**Caption:** Carrier and all-pair scores use decoded TMRCA in the s=0.005 cohort. p<=0.05 and 0.01 use 400 neutral calibration positions per fold. The separate p<=0.001 sensitivity uses all 1,000 neutrals for selected ranks and leave-one-out ranks for neutrals. One extreme neutral rank gives 1/1,000 calls for these scores; this is not precise independent validation of a 0.1% tail. These are positional tests, not genome-wide family-wise tests.

## 13. Decoded carrier-score components (backup)

[PNG](figures/13_carrier_genealogies.png) | [Vector PDF](figures/13_carrier_genealogies.pdf)

**Takeaway:** Decoded carrier mass weights allele-pair abundance by the fraction of those same pairs called recent.

**Caption:** Both panels use Gamma-SMC decoded TMRCA at the focal selected allele in the 100 s=0.005 regions. Left: carrier mass; right: within-ALT/ALT frac_recent_T. The score is exactly the sampled ALT/ALT pair fraction times the right-panel quantity. AF squared approximates the pair fraction for random distinct pairs, with finite-panel and pair-sampling differences; it is not an exact bound.

## 14. Demography and decoder assumptions (backup)

[PNG](figures/14_demography_and_assumptions.png) | [Vector PDF](figures/14_demography_and_assumptions.pdf)

**Takeaway:** Matched nulls address calibration under the model; they do not establish robustness to every model misspecification.

**Caption:** The EAS trajectory is read from the same PHLASH-derived model used by the archived simulations. The Ne history is displayed only to 50 kya; simulation genealogies also contain the older archaic branch. Decoder Ne is theta/(4*mu)=15,000. Neutral and selected simulations share the archaic bottleneck and rate assumptions. No variable-Ne Gamma-SMC flow field or empirical recombination-map benchmark is claimed here.

## 15. Pair-class coalescence diagnostic (backup)

[PNG](figures/15_pair_class_coalescence.png) | [Vector PDF](figures/15_pair_class_coalescence.pdf)

**Takeaway:** The carrier score combines this within-ALT/ALT recency with the fraction of the complete pair panel that is ALT/ALT.

**Caption:** All recency calls are Gamma-SMC decoded. ALT/ALT and REF/REF are haplotype-pair classes at the focal archaic marker. Only regions with an assigned marker contribute: 100 selected and 169 neutral. Counts are pooled across pairs, so regions with larger classes contribute more weight; this is not a replicate-average or the unconditional neutral denominator used for FPR. REF/REF is diagnostic and is not subtracted from the score; ALT/REF is not scored.

## Interpretation limits

Selected power is conditional on allele survival and observation. The unconditional positional null includes marker-free positions; it is not a null conditional on first selecting an observed archaic SNP. The archaic-marker oracle is more restrictive than simply lying on an introgressed segment and is not a literal finite-reference Neanderthal/Denisovan/African ascertainment experiment.

This array uses a strong shared archaic bottleneck, uniform mutation/recombination rates and isolated 10-Mb regions. It does not demonstrate genome-wide FDR control, causal localization, or a gain from multiple independent introgressing haplotypes. Cutoff comparisons are exploratory. Wilson intervals are conditional on the estimated null and do not include its Monte Carlo uncertainty.

## Methods references

Voight et al. (2006), iHS: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072

scikit-allel implementation: https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/selection.html
