# Chronological decision log

## 1. Input review

The task provides error-free haploid parent/offspring assemblies from one transmission, a new parent assembly, and one-based query coordinates. Only single-nucleotide changes occur. No external scientific sources, other projects, hidden outcomes, or other agents are consulted. The six supplied inputs will remain unchanged. The working hypothesis is a stochastic substitution process whose rates may depend on the parental sequence context. Alternatives to examine include base-only rates, neighboring sequence context, reverse-complement symmetry, longer motifs, and spatial or contig variation. All random sampling will use seed 0. Computation defaults to four threads with a 26 GB memory ceiling.

The first analysis will check sequence correspondence, nucleotide validity, contig lengths, query coverage, and empirical substitution patterns before choosing a model. Model selection will use held-out calibration sites/contigs, never the unavailable new offspring.

## 2. Initial diagnostics

All calibration sequences are valid A/C/G/T, with matched identifiers and lengths: 80 contigs of 50,000 bases (4 million sites). The new parent has 20 contigs of 50,000 bases. The 999,880 queries cover positions 4 through 49,997. Mutation rates are approximately 0.5% overall. The most obvious context effect is an elevated C-to-T rate when the right neighbor is G and the complementary G-to-A rate when the left neighbor is C (approximately 0.9–1.0%). Longer context tables contain sparse, extreme empirical rates that could be sampling noise. Complementary substitutions have similar empirical frequencies, suggesting a reverse-complement pooling model worth testing, without assuming it a priori.

Next decision: compare base-only, CpG-specific, trinucleotide, pentanucleotide, and heptanucleotide models using five held-out groups of whole contigs. Use the same inner three-base margin for comparable scoring. Test both strand-specific and reverse-complement-pooled tables, with shrinkage to simpler models to avoid zero probabilities and sparse-context overfitting. Preserve out-of-fold metrics and uncertainty in model comparisons.

## 3. First held-out results and revised hypothesis

The CpG model improves total held-out negative log likelihood by about 2,718 relative to base-only rates. Trinucleotide and pentanucleotide detail does not improve the CpG model. However, a heptanucleotide table shrunk with 10,000 prior observations improves by 56.7 NLL over the reverse-complement CpG model (paired contig standard error 6.4). This is a reproducible signal that should not be dismissed as sparse noise. The best oriented heptanucleotide table beats its symmetric counterpart by 9.0 NLL, but this is only about 1.4 paired standard errors. The next analysis will test whether the long-context signal reflects a few hotspots, mutation-rate heterogeneity, alternative-specific effects, or simpler long-range sequence features. These are exploratory tests motivated by calibration evidence and require further held-out comparison.

## 4. Discovery of a compact long-range rule

Aggregate residual tests over flank subsets identify one very strong pair: offsets -3 and +2 (Pearson statistic approximately 3,686 with 60 nominal degrees of freedom). Inspection shows that equality of these two parental bases predicts elevated A-to-C/G and T-to-C/G rates, approximately fourfold. A-to-T/T-to-A and C/G-centered changes do not show the same increase. Ordinary reverse-complementing of the context moves this oriented equality test to different offsets, explaining why fully symmetric heptanucleotide pooling weakened the signal. A compact model can preserve the observed assembly orientation for the equality flag while pooling complementary substitution outcomes.

Next: compare binary equality, equality nucleotide identity, full flank-pair identity, and complement pooling while preserving the oriented equality indicator. Check whether the same flank pair is selected separately inside each training fold. Because exploratory feature tests used all calibration data, initial cross-validation is selection-informed; reporting will explicitly distinguish it from truly untouched validation. The strength and reproducibility of the discovered rule will be assessed with fold-specific rediscovery and paired contig metrics.

## 5. Validation of compact rules and secondary signal

All five training folds independently select the (-3,+2) equality rule by a large margin. The compact two-rule model reaches held-out NLL 140,740.4 versus 142,637.8 for CpG alone. Equality nucleotide identity yields only a 1.35 NLL improvement (paired standard error 5.50), so the simpler pooled rate model is preferred. Pair identity, unrestricted strand effects, and adding equality dependence to all C/G rates do not show convincing improvement.

A substitution-specific residual scan discovers a second equality rule: when parental positions -2 and -1 match, A-to-T/T-to-A rises from approximately 0.00079 to 0.00158. This effect was diluted in total-mutation residuals because it affects only one relatively rare outcome. Each training fold independently rediscovers this pair as the strongest remaining effect. Adding this third rule improves held-out NLL by 105.76 (paired contig standard error 14.61), to 140,634.66. After fitting the three rules, the largest residual statistic across the 90 flank-equality/substitution comparisons is only 5.98; no further equality rule is supported after multiple comparisons.

The selected final model has nine pooled substitution-rate parameters: four for A/T changes to C/G stratified by (-3,+2) equality, two for A/T interchange stratified by (-2,-1) equality, two for C/G transitions stratified by CpG, and one common C/G transversion rate. All triggers are evaluated in the given parental assembly orientation. The complement pooling is on event rates, not a reverse-complement transformation of the oriented equality triggers. The probabilities of the three possible changes are estimated from pooled event counts with a small Jeffreys-style half-count; the parental probability is their complement. Contig-level mutation dispersion is modest (approximate chi-square 98.91 on 79 degrees of freedom, p=0.064), insufficient to justify an unobservable new-contig random effect.

## 6. Final output plan

Fit the nine-parameter model to all 3,999,520 calibration sites having full three-base flanks, predict all 999,880 supplied queries in their exact order, and compute the expected changed count by summing one minus the parental probability. Independently re-read the output to verify coordinates, probability finiteness/nonnegativity, row sums, row count, and expected-count consistency. Recheck SHA-256 hashes of all six supplied files. Report parameter and outcome uncertainty separately; neither the realized offspring sequence nor exact changed count can be inferred from these data.

## 7. Prediction result and covariate-shift diagnostic

The output expectation is 6,974.537599 changed alleles across 999,880 queries. This is higher than a simple rescaling of the calibration burden for a concrete reason: among A/T query sites, 399,403 of 498,750 have the (-3,+2) equality, about 80.08%, versus about 25.08% in calibration. These are common, well-supported conditional states in calibration, but their mixture weights change in the new parent. This confirms the practical importance of modeling the long-range conditional rates. The new parent's feature frequencies were inspected only after fitting and choosing the final model.

The saved output passed a complete row-by-row coordinate and probability check. All 999,880 rows are finite and nonnegative, maximum row-sum error is 4.44e-16, and the expectation reconstructed from the saved probabilities agrees to floating-point precision. The nine-rate exposure-based reconstruction independently agrees. Approximate uncertainty scales are 61.03 changes for fitted-mean parameter uncertainty and 83.14 for conditional independent-site outcome variation; these calculations do not account for adaptive-selection uncertainty or unknown site dependence.

## 8. Completion checks

An independently written explicit string-based implementation reproduced all 16,384 possible seven-base-context probability vectors to absolute tolerance 1e-15. All six input checksums and the saved prediction checksum matched. The prediction file has 999,881 lines including its header. `REPORT.md` was generated and inspected; it documents alternatives, cross-validation, adaptive-selection limits, final rate formulas, new-parent covariate shift, uncertainty, provenance, and phased uv reproduction commands. Required outputs are complete. No source files were changed, no external information was used, and no plots were displayed.
