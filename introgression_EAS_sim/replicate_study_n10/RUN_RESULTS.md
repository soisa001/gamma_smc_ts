# Replicated EAS introgressed-sweep results

## Bottom line

This bounded campaign requested 270 independently seeded trajectory slots across 27 biological selection-by-AF-floor cells. 270 slots were attempted, 155 produced accepted 100-diploid Han panels, and 155 were decoded with Gamma-SMC. 24 biological cells had at least one accepted trajectory.

The headline cell-equal conditional recovery estimates were Normalized integrated TMRCA-CDF area (alt/alt - ref/ref; expected > 0): 91.9% across 21 contributing cells; Mean TMRCA in generations (alt/alt - ref/ref; expected < 0): 95.1% across 22 contributing cells. These are qualitative recovery checks in accepted, tree-truth-positive trajectories, not estimates of power.

## Execution and coverage

| Selection coefficient | Requested slots | Attempted | Accepted | Decoded | Cells represented |
|---|---|---|---|---|---|
| 0.01 | 90 | 90 | 23 | 23 | 8/9 |
| 0.005 | 90 | 90 | 84 | 84 | 9/9 |
| 0.001 | 90 | 90 | 48 | 48 | 7/9 |

### Accepted/requested slots by nested AF floor

| Selection | AF≥10% | AF≥20% | AF≥30% | AF≥40% | AF≥50% | AF≥60% | AF≥70% | AF≥80% | AF≥90% |
|---|---|---|---|---|---|---|---|---|---|
| s=0.01 | 3/10 (d=3) | 3/10 (d=3) | 5/10 (d=5) | 4/10 (d=4) | 1/10 (d=1) | 3/10 (d=3) | 2/10 (d=2) | 2/10 (d=2) | 0/10 (d=0) |
| s=0.005 | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 4/10 (d=4) |
| s=0.001 | 10/10 (d=10) | 10/10 (d=10) | 10/10 (d=10) | 8/10 (d=8) | 5/10 (d=5) | 4/10 (d=4) | 0/10 (d=0) | 1/10 (d=1) | 0/10 (d=0) |

Among accepted sample panels, achieved selected-allele AF was 12.0%–94.5% (median 73.0%).
Observed per-panel genotype-count ranges were ref/ref 2–78, heterozygous 7–64, alt/alt 2–91.

The AF targets are nested acceptance floors `[f, 0.95]`, not disjoint bins. An achieved AF can satisfy every lower floor beneath it, although each cell has its own seed stream. Cross-floor differences therefore do not estimate a response across mutually exclusive AF strata.

All effects are conditional on trajectories that returned within the bounded time/draw budgets and passed population AF, sample AF, and genotype-comparator checks. Failures are nonrandom. Intervals describe variation among accepted trajectories; they do not account for conditioning selection or missing cells. Tree-truth summaries use accepted trajectories, whereas Gamma and paired summaries additionally require a valid decode.

## Truth-versus-Gamma recovery

### Headline: biological-cell-equal macro summary

| Metric | Eligible cells | Paired trajectories | Truth expected sign | Gamma expected sign | Raw sign agreement | Gamma detection given truth | Bias | MAE | RMSE |
|---|---|---|---|---|---|---|---|---|---|
| Normalized integrated TMRCA-CDF area (alt/alt - ref/ref; expected > 0) | 22 | 153 | 73.6% | 73.6% | 87.7% | 91.9% | 0.0204 | 0.1126 | 0.1338 |
| Mean TMRCA in generations (alt/alt - ref/ref; expected < 0) | 22 | 153 | 88.0% | 87.0% | 91.4% | 95.1% | 8,943.0 | 12,466.8 | 16,289.1 |

The macro summary gives each eligible biological cell equal weight. Eligibility requires at least two paired trajectories; conditional detection also requires at least one truth-positive trajectory in the cell. No confidence interval is computed for this macro mean.

### Secondary: accepted-trajectory-weighted micro summary

| Metric | Paired trajectories | Truth expected sign (95% CI) | Gamma expected sign (95% CI) | Raw sign agreement (95% CI) | Gamma detection given truth (95% CI) | Bias | MAE | RMSE |
|---|---|---|---|---|---|---|---|---|
| Normalized integrated TMRCA-CDF area (alt/alt - ref/ref; expected > 0) | 155 | 76.1% (68.8%–82.2%) | 74.8% (67.5%–81.0%) | 88.4% (82.4%–92.5%) | 91.5% (85.1%–95.3%) | 0.0097 | 0.1038 | 0.1367 |
| Mean TMRCA in generations (alt/alt - ref/ref; expected < 0) | 155 | 91.0% (85.4%–94.5%) | 87.7% (81.6%–92.0%) | 92.9% (87.7%–96.0%) | 94.3% (89.2%–97.1%) | 9,681.4 | 12,835.3 | 18,736.3 |

The micro summary gives each accepted paired trajectory equal weight and therefore overweights easier, higher-coverage cells. It is a different estimand from the headline cell-equal macro summary.

Raw sign agreement can count a jointly contrary truth/Gamma pair as agreement. Conditional detection does not: it includes only trajectories with the expected tree-truth sign and asks whether Gamma-SMC also has that sign.

No matched neutral/null trajectories were simulated. These rates are descriptive checks, not p-values, sensitivity, specificity, false-positive rate, power, or ROC AUC. The normalized integrated CDF area is an area under a TMRCA CDF, not classifier AUC.

## TMRCA probability-profile accuracy: overall micro summary

This threshold-error summary pools accepted, decoded paired trajectories across biological cells with equal weight per trajectory. It is therefore a secondary overall micro/trajectory-weighted estimand: easier, higher-coverage cells contribute more heavily and can dominate the pooled error.

At the largest decoded cutoff (50 kya):

| Genotype class | Paired | Bias | MAE | RMSE |
|---|---|---|---|---|
| Overall | 155 | -0.1160 | 0.1468 | 0.1956 |
| Ref/ref | 155 | -0.0412 | 0.0876 | 0.1420 |
| Alt/alt | 155 | -0.2572 | 0.3035 | 0.3701 |
| Heterozygous | 155 | 0.0005 | 0.0005 | 0.0049 |

Cell-stratified overall, ref/ref, and alt/alt curves and their trajectory-bootstrap bands are in `results/aggregate_cdf_profiles.tsv` and the PNG/PDF figure set. They are not pooled into an AF-response curve.

## Interpretation and provenance limits

- These results characterize the configured `Q=10`, `slim_burn_in=0.1` process. Q rescaling can materially reduce small simulated populations and collapse short demographic phases; in this campaign, the scaled source check and introgression pulse both occur at 2,270.0 generations ago, so their requested separation is not retained. Q/burn-in sensitivity runs are required before biological calibration.
- Archaic specificity/no ILS is bound to the validated `serialized_forward_event_contract`: mutation mode `single archaic-specific de novo origin`, declared origin population `Neanderthal`, origin after the human-Neanderthal split and before the pulse, and Han entry route `introgression_pulse_only`. The population attached to the retained mutation node is contractually not origin evidence; donor AF, migrant-copy count, and the full trajectory are not recovered from the returned tree.
- Artifact validity does not establish semantic equivalence across implementation source epochs. Any format-only transition that is not cryptographically bound per task remains an operator attestation; a frozen-source rerun is required to remove that provenance limitation.

## Machine-readable outputs

- `results/cell_coverage.tsv`: requested and observed coverage for all biological cells.
- `results/replicate_acceptance.tsv`: achieved sample AF, genotype counts, seeds, and effort for accepted panels.
- `results/cell_effect_summary.tsv`: exact cell/source focal-effect summaries and trajectory-bootstrap intervals.
- `results/truth_gamma_accuracy.tsv`: explicitly labeled macro and micro recovery estimands.
- `results/aggregate_cdf_profiles.tsv` and `results/threshold_accuracy.tsv`: TMRCA-CDF profiles and paired truth/Gamma errors; the overall threshold rows are explicitly trajectory-weighted micro summaries.
- `results/replicate_results_manifest.json`: artifact-relative paths, SHA-256 checksums, contracts, decoder hashes, and environment provenance. This `RUN_RESULTS.md` file is itself checksum-covered by that manifest.
