# Carrier mass versus REF/REF subtraction

Mass: w_AA * frac_recent_T_AA. Contrast: w_AA * max(frac_recent_T_AA - frac_recent_T_RR, 0).
Both omit ALT/REF and the all-pair gate. T=50 kya; immediate-onset EAS s=0.005.
The unit of comparison is an independently simulated 10 Mb region, not a haplotype pair.

| Operating point | Positions | TMRCA | Cohort | Mass calls | Contrast calls | Mass only | Contrast only | Nominal paired p | Holm p (16) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| regional_p_0.05 | grid10kb | decoded | selected | 83/100 | 81/100 | 3 | 1 | 0.625 | 1 |
| regional_p_0.05 | grid10kb | decoded | neutral | 48/1000 | 48/1000 | 6 | 6 | 1 | 1 |
| regional_p_0.05 | grid10kb | truth | selected | 83/100 | 80/100 | 6 | 3 | 0.507812 | 1 |
| regional_p_0.05 | grid10kb | truth | neutral | 50/1000 | 53/1000 | 9 | 12 | 0.663624 | 1 |
| regional_p_0.05 | archaic_sites | decoded | selected | 84/100 | 80/100 | 4 | 0 | 0.125 | 1 |
| regional_p_0.05 | archaic_sites | decoded | neutral | 49/1000 | 46/1000 | 7 | 4 | 0.548828 | 1 |
| regional_p_0.05 | archaic_sites | truth | selected | 81/100 | 79/100 | 5 | 3 | 0.726562 | 1 |
| regional_p_0.05 | archaic_sites | truth | neutral | 51/1000 | 52/1000 | 10 | 11 | 1 | 1 |
| training_power_target_0.70 | grid10kb | decoded | selected | 71/100 | 72/100 | 6 | 7 | 1 | 1 |
| training_power_target_0.70 | grid10kb | decoded | neutral | 13/1000 | 20/1000 | 0 | 7 | 0.015625 | 0.234375 |
| training_power_target_0.70 | grid10kb | truth | selected | 69/100 | 70/100 | 1 | 2 | 1 | 1 |
| training_power_target_0.70 | grid10kb | truth | neutral | 18/1000 | 20/1000 | 4 | 6 | 0.753906 | 1 |
| training_power_target_0.70 | archaic_sites | decoded | selected | 69/100 | 72/100 | 4 | 7 | 0.548828 | 1 |
| training_power_target_0.70 | archaic_sites | decoded | neutral | 14/1000 | 24/1000 | 0 | 10 | 0.00195312 | 0.03125 |
| training_power_target_0.70 | archaic_sites | truth | selected | 70/100 | 69/100 | 2 | 1 | 1 | 1 |
| training_power_target_0.70 | archaic_sites | truth | neutral | 18/1000 | 21/1000 | 3 | 6 | 0.507812 | 1 |

Two-sided exact McNemar calculation: binomial test of mass-only versus contrast-only calls.
These are conditional, exploratory p-values on the saved fold decisions. Shared training/calibration
regions create dependence among cross-validated predictions; the calculation does not propagate
threshold-fitting uncertainty. Holm adjustment covers only these 16 reported comparisons, not
earlier score/cutoff development. Independent new-seed validation is needed for confirmation.

The 70% operating point trains each score's threshold separately toward 70% training power;
held-out power is not forced to match. It is not a p=0.70 significance threshold.
Neutral call fraction is regional false-positive rate, not prevalence-dependent FDR.
