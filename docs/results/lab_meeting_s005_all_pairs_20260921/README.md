# EAS s=0.005: all pairs only

10 figures in the same layout as the s=0.002 all-pair pack. Each has a PNG and a vector PDF with editable text; the combined PDF has one figure per page.

All primary TMRCA figures use saved Gamma-SMC decoded frac_recent_T across all 10,000 sampled pairs, without carrier conditioning. Figure 07 alone plots true TMRCA as an explicit validation reference. The cohort contains 100 selected and 1,000 neutral regions, with 400 haplotypes per region. No new simulations, decoding or ancestry replays were run.

Selection starts at the 2% pulse, 50 kya. The primary null is matched to the same coordinate. Five-fold evaluation uses 400 calibration and 200 held-out neutrals per fold; all fixed selected replicates remain included. Power is conditional on original allele survival/observation. Positional FPR is not discovery FDR or genome-wide family-wise error.

## 01. Design and score

[PNG](figures/01_design_and_score.png) | [Vector PDF](figures/01_design_and_score.pdf)

**Takeaway:** Compare the same decoded statistic at the same coordinate in selected and neutral regions.

**Caption:** Same saved 100 s=0.005 selected simulations; no new simulations or ancestry replays. Selection begins at the 2% pulse, 50 kya. Fixed mutation/recombination rates, 25-year generations, h=0.5. Selected power is conditional on survival and observation; all fixed replicates are retained. All-pair scoring uses every sampled pair regardless of archaic allele status.

## 02. Focal all-pair score distributions

[PNG](figures/02_matched_position_null.png) | [Vector PDF](figures/02_matched_position_null.pdf)

**Takeaway:** Separation of the score distributions determines positional power.

**Caption:** The score uses all 10,000 sampled pairs at the same 5-Mb coordinate in each region. There is no requirement for an archaic marker and no maximum across the region. Selected simulations were conditioned on allele survival and observation when generated; no further selected or neutral filtering is performed here.

## 03. Decoded focal power

[PNG](figures/04_positional_power.png) | [Vector PDF](figures/04_positional_power.pdf)

**Takeaway:** Decoded all-pair power at T=50 kya is 61% at nominal 5%.

**Caption:** Each selected replicate is ranked against 400 neutral values at the same focal coordinate in its assigned fold; 200 held-out neutral positions per fold measure FPR. Five folds cover 1,000 neutrals. Wilson intervals describe binomial uncertainty conditional on the estimated null, not its Monte Carlo uncertainty. The all-pair score uses all 10,000 pairs.

## 04. Decoded power across time cutoffs

[PNG](figures/05_time_cutoff_power.png) | [Vector PDF](figures/05_time_cutoff_power.pdf)

**Takeaway:** Calibration and power use the same decoded score and genomic coordinate.

**Caption:** All scores are Gamma-SMC decoded and separately calibrated at each cutoff (5, 10, 20, 30, 40, 50 kya). No maximum across cutoffs enters the test. Cells use all 100 selected replicates. Cutoff comparisons are exploratory, not per-replicate tuning.

## 05. Positional neutral FPR

[PNG](figures/06_positional_false_positive_rate.png) | [Vector PDF](figures/06_positional_false_positive_rate.pdf)

**Takeaway:** Calibration and power use the same decoded score and genomic coordinate.

**Caption:** All scores are Gamma-SMC decoded and separately calibrated at each cutoff (5, 10, 20, 30, 40, 50 kya). No maximum across cutoffs enters the test. FPR is the fraction of neutral positions called, not discovery FDR. Each neutral is held out from its own calibration set; all 1,000 positions remain in the denominator. Color scales differ between nominal thresholds.

## 06. Truth versus decoded validation

[PNG](figures/07_truth_vs_decoding.png) | [Vector PDF](figures/07_truth_vs_decoding.pdf)

**Takeaway:** At 50 kya and nominal 5%, all-pair power is 86% from truth and 61% after decoding.

**Caption:** This is the only figure that uses true TMRCA as a plotted reference. Each source has its own matching neutral calibration. The decoded rule thresholds posterior-mean TMRCA per pair; it does not average posterior probability mass. Matching the null calibrates errors but need not restore discrimination lost in decoding.

## 07. All-pair detection across the region

[PNG](figures/10_spatial_decay.png) | [Vector PDF](figures/10_spatial_decay.pdf)

**Takeaway:** The positional null is matched separately at every stride; unrelated peaks never set the cutoff.

**Caption:** Only all-pair scores are used. Each point is the fraction of 100 selected or 1,000 neutral regions called at that stride. Full 10-Mb results are saved; +/-500 kb is displayed without smoothing. Linked signal is not counted as a new causal target or localization error. All 10,000 sampled pairs contribute at each stride.

## 08. Example all-pair profiles

[PNG](figures/11_example_regions.png) | [Vector PDF](figures/11_example_regions.pdf)

**Takeaway:** All-pair profiles are available across the entire saved 10-Mb region.

**Caption:** T=50 kya. Neutral neutral/rep0286 has focal decoded all-pair score closest to the neutral median. Selected replicate 0014 has selected AF closest to the s=0.005 median; ties use IDs. Gold shading marks a fixed 100-kb interval around 5 Mb. These examples were not selected for significant peaks and do not estimate power.

## 09. Stringency sensitivity

[PNG](figures/12_threshold_stringency.png) | [Vector PDF](figures/12_threshold_stringency.pdf)

**Takeaway:** All-pair power at pooled p<=0.001 is 23%.

**Caption:** Nominal 5% and 1% use 400 calibration neutrals per fold; 0.1% uses all 1,000 neutrals to rank selected regions and leave-one-out ranks for neutrals. This is a different calibration size, and only one extreme neutral rank is available. It is not precise independent validation of a 0.1% tail or genome-wide error control.

## 10. Demography and decoder assumptions

[PNG](figures/14_demography_and_assumptions.png) | [Vector PDF](figures/14_demography_and_assumptions.pdf)

**Takeaway:** Matched nulls address calibration under the model; they do not establish robustness to every model misspecification.

**Caption:** The EAS trajectory is read from the same PHLASH-derived model used by the archived simulations. The Ne history is displayed only to 50 kya; simulation genealogies also contain the older archaic branch. Decoder Ne is theta/(4*mu)=15,000. Neutral and selected simulations share the archaic bottleneck and rate assumptions. No variable-Ne Gamma-SMC flow field or empirical recombination-map benchmark is claimed here.

## Interpretation limits

Uniform mutation/recombination, isolated 10-Mb regions, a shared strong archaic bottleneck and a fixed-Ne decoder. Cutoffs are compared separately, without per-replicate tuning. Wilson intervals omit uncertainty in the estimated neutral null. Pooled p=0.001 is a sensitivity analysis with only one extreme neutral rank. The data do not establish empirical genome-wide FDR control or causal localization.
