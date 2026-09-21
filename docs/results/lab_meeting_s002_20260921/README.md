# EAS s=0.002: focal carrier comparison

15 figures: PNG and vector PDF with editable text, plus one combined PDF.

Primary TMRCA figures use Gamma-SMC decoded frac_recent_T; only figure 07 plots truth as an explicit validation reference. 100 saved selected regions, 1,000 existing decoded neutrals, 400 haplotypes and 10,000 identical sampled pairs. No simulations, ancestry recovery, or trajectory replays were run.

All-pair spatial scores are available across the 10 Mb region. Archaic carrier classes are known only at the focal selected allele in this arm. Spatial carrier plots are omitted, and missing labels are never treated as observed absence of archaic ancestry.

Error rates are matched-position FPR, not discovery FDR or genome-wide family-wise error. Selected power is conditional on survival and observation; fixed alleles stay in the denominator. 5-fold evaluation uses 400 neutral calibration and 200 held-out positions per fold. The p=0.001 pooled-rank sensitivity has only one extreme neutral rank. Cutoff comparisons are exploratory. Wilson intervals omit uncertainty in the estimated null.

## 01. Design and score

[PNG](figures/01_design_and_score.png) | [Vector PDF](figures/01_design_and_score.pdf)

**Takeaway:** Compare the same decoded statistic at the same coordinate in selected and neutral regions.

**Caption:** Same saved 100 s=0.002 selected simulations; no new simulations or ancestry replays. Selection begins at the 2% pulse, 50 kya. Fixed mutation/recombination rates, 25-year generations, h=0.5. Selected power is conditional on survival and observation; all fixed replicates are retained. Carrier information is available at the selected allele only. The carrier score has no ALT/REF or REF/REF penalty. Spatial plots show all-pair recency only.

## 02. Positional versus regional null

[PNG](figures/02_matched_position_null.png) | [Vector PDF](figures/02_matched_position_null.pdf)

**Takeaway:** Only 4.3% of neutral focal positions reach 18% AF; unrelated peaks elsewhere do not enter the local test.

**Caption:** The focal score assigns the nearest eligible archaic marker within 5 kb of the fixed focal stride. All 1,000 neutral positions are retained; 169 have an assigned marker. Marker-free positions score zero. The whole-10-Mb maximum is shown only to explain the earlier calibration mismatch. A pre-specified gene requires calibration of the same summary in a matching neutral gene window.

## 03. Allele frequency by selection strength

[PNG](figures/03_allele_frequency_distribution.png) | [Vector PDF](figures/03_allele_frequency_distribution.pdf)

**Takeaway:** Weak selection produces a broad AF distribution; its mean does not imply every replicate is detectable.

**Caption:** Boxes show interquartile range and median; whiskers show the 5th and 95th percentiles; diamonds show means. Selected data are conditional on survival and sample observation. Neutral data are unconditional pre-specified positions, with absent nearby markers scored zero. These are not frequency-matched cohorts.

## 04. Decoded focal power

[PNG](figures/04_positional_power.png) | [Vector PDF](figures/04_positional_power.pdf)

**Takeaway:** Decoded all-pair power at T=50 kya is 5% at nominal 5%.

**Caption:** Each selected replicate is ranked against 400 neutral values at the same focal coordinate in its assigned fold; 200 held-out neutral positions per fold measure FPR. Five folds cover 1,000 neutrals. Wilson intervals describe binomial uncertainty conditional on the estimated null, not its Monte Carlo uncertainty. Carrier and all-pair scores use Gamma-SMC frac_recent_T; AF and iHS use genotypes. iHS is a separately calibrated endpoint.

## 05. Decoded power across time cutoffs

[PNG](figures/05_time_cutoff_power.png) | [Vector PDF](figures/05_time_cutoff_power.pdf)

**Takeaway:** Calibration and power use the same decoded score and genomic coordinate.

**Caption:** All scores are Gamma-SMC decoded and separately calibrated at each cutoff (5, 10, 20, 30, 40, 50 kya). No maximum across cutoffs enters the test. Cells use all 100 selected replicates. Cutoff comparisons are exploratory, not per-replicate tuning.

## 06. Positional neutral FPR

[PNG](figures/06_positional_false_positive_rate.png) | [Vector PDF](figures/06_positional_false_positive_rate.pdf)

**Takeaway:** Calibration and power use the same decoded score and genomic coordinate.

**Caption:** All scores are Gamma-SMC decoded and separately calibrated at each cutoff (5, 10, 20, 30, 40, 50 kya). No maximum across cutoffs enters the test. FPR is the fraction of neutral positions called, not discovery FDR. Each neutral is held out from its own calibration set; all 1,000 positions remain in the denominator. Color scales differ between nominal thresholds.

## 07. Truth versus decoded validation

[PNG](figures/07_truth_vs_decoding.png) | [Vector PDF](figures/07_truth_vs_decoding.pdf)

**Takeaway:** At 50 kya and nominal 5%, all-pair power is 19% from truth and 5% after decoding.

**Caption:** This is the only figure that uses true TMRCA as a plotted reference. Each source has its own matching neutral calibration. The decoded rule thresholds posterior-mean TMRCA per pair; it does not average posterior probability mass. Matching the null calibrates errors but need not restore discrimination lost in decoding.

## 08. Carrier mass versus AF

[PNG](figures/08_gain_beyond_allele_frequency.png) | [Vector PDF](figures/08_gain_beyond_allele_frequency.pdf)

**Takeaway:** At nominal 1%, AF detects 41%; decoded carrier mass detects 33%-38% across cutoffs.

**Caption:** Both panels use decoded carrier mass at the focal selected allele. Negative bars are AF-only detections, positive bars carrier-only detections. Exact paired tests in the accompanying table are exploratory and unadjusted across cutoffs. No tuning of T per replicate is used.

## 09. iHS power and target-site eligibility

[PNG](figures/09_ihs_and_sweep_completion.png) | [Vector PDF](figures/09_ihs_and_sweep_completion.pdf)

**Takeaway:** Near fixation, the selected SNP is often unscorable by iHS; nearby SNPs may still carry signal.

**Caption:** The two iHS curves represent separately calibrated endpoints. Positional iHS uses the nearest finite-scoring core with MAF>=5% within 5 kb; the fixed 100-kb endpoint uses the fraction of finite-scoring cores with |standardized iHS|>2, requiring at least 20 cores. Right: exact focal-SNP eligibility, not eligibility of the fallback positional statistic. Fixation and missing/undefined scores never remove replicates from the denominator.

## 10. All-pair detection across the region

[PNG](figures/10_spatial_decay.png) | [Vector PDF](figures/10_spatial_decay.pdf)

**Takeaway:** The positional null is matched separately at every stride; unrelated peaks never set the cutoff.

**Caption:** Only all-pair scores are used. Each point is the fraction of 100 selected or 1,000 neutral regions called at that stride. Full 10-Mb results are saved; +/-500 kb is displayed without smoothing. Linked signal is not counted as a new causal target or localization error. Spatial carrier scores are unavailable from the saved labels and are omitted.

## 11. Example all-pair profiles

[PNG](figures/11_example_regions.png) | [Vector PDF](figures/11_example_regions.pdf)

**Takeaway:** All-pair profiles are available across the entire saved 10-Mb region.

**Caption:** T=50 kya. Neutral neutral/rep0286 has focal decoded all-pair score closest to the neutral median. Selected replicate 0006 has selected AF closest to the s=0.002 median; ties use IDs. Gold shading marks a fixed 100-kb interval around 5 Mb. These examples were not selected for significant peaks and do not estimate power.

## 12. Stringency sensitivity

[PNG](figures/12_threshold_stringency.png) | [Vector PDF](figures/12_threshold_stringency.pdf)

**Takeaway:** All-pair power at pooled p<=0.001 is 0%.

**Caption:** Nominal 5% and 1% use 400 calibration neutrals per fold; 0.1% uses all 1,000 neutrals to rank selected regions and leave-one-out ranks for neutrals. This is a different calibration size, and only one extreme neutral rank is available. It is not precise independent validation of a 0.1% tail or genome-wide error control.

## 13. Focal decoded carrier-score components

[PNG](figures/13_carrier_genealogies.png) | [Vector PDF](figures/13_carrier_genealogies.pdf)

**Takeaway:** The same pairs must both carry the allele and have recent decoded TMRCA.

**Caption:** Carrier mass is exactly within-ALT/ALT frac_recent_T times the ALT/ALT fraction of the 10,000-pair panel. AF squared only approximates this sampled-pair fraction. Replicates with no sampled ALT/ALT pair have mass zero and undefined within-class recency; they remain in the power denominator. All labels are from the saved selected allele; no spatial archaic labels are inferred.

## 14. Demography and decoder assumptions

[PNG](figures/14_demography_and_assumptions.png) | [Vector PDF](figures/14_demography_and_assumptions.pdf)

**Takeaway:** Matched nulls address calibration under the model; they do not establish robustness to every model misspecification.

**Caption:** The EAS trajectory is read from the same PHLASH-derived model used by the archived simulations. The Ne history is displayed only to 50 kya; simulation genealogies also contain the older archaic branch. Decoder Ne is theta/(4*mu)=15,000. Neutral and selected simulations share the archaic bottleneck and rate assumptions. No variable-Ne Gamma-SMC flow field or empirical recombination-map benchmark is claimed here.

## 15. Pair-class coalescence diagnostic

[PNG](figures/15_pair_class_coalescence.png) | [Vector PDF](figures/15_pair_class_coalescence.pdf)

**Takeaway:** The carrier score combines this within-ALT/ALT recency with the fraction of the complete pair panel that is ALT/ALT.

**Caption:** All recency calls are Gamma-SMC decoded. ALT/ALT and REF/REF are haplotype-pair classes at the focal archaic marker. Only regions with an assigned marker contribute: 100 selected and 169 neutral. Counts are pooled across pairs, so regions with larger classes contribute more weight; this is not a replicate-average or the unconditional neutral denominator used for FPR. REF/REF is diagnostic and is not subtracted from the score; ALT/REF is not scored.

## Model limits

Uniform mutation/recombination rates; isolated 10-Mb regions; shared strong archaic bottleneck; fixed-Ne decoder matched in neutral and selected sims. These comparisons do not test variable-Ne flow fields, empirical maps, genome-wide FDR, or independent multiple-haplotype introgression. Carrier labels use an oracle, not a finite archaic-reference/outgroup experiment. Trajectory figures 03_a/03_b remain deferred.
