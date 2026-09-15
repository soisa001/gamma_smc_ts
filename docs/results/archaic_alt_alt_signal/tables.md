# Archaic ALT/ALT signal: existing results

Score: `frac_recent_T_among_ALT_ALT * n_ALT_ALT / n_total_pairs`.
The saved method identifiers retain the historical `mass_g0_r1_` prefix.
Only immediate selection at the 50-kya pulse is included.

## Regional p <= 0.05, T = 50 kya

| Positions | TMRCA source | Selected called / 100 | Neutral called / 1,000 |
|---|---|---:|---:|
| 10 kb grid | Truth | 83 (83%) | 50 (5.0%) |
| 10 kb grid | Decoded | 83 (83%) | 48 (4.8%) |
| Every archaic site | Truth | 81 (81%) | 51 (5.1%) |
| Every archaic site | Decoded | 84 (84%) | 49 (4.9%) |

## Thresholds trained toward 70% power

| Positions | TMRCA source | Selected called / 100 | Neutral called / 1,000 |
|---|---|---:|---:|
| 10 kb grid | Truth | 69 (69%) | 18 (1.8%) |
| 10 kb grid | Decoded | 71 (71%) | 13 (1.3%) |
| Every archaic site | Truth | 70 (70%) | 18 (1.8%) |
| Every archaic site | Decoded | 69 (69%) | 14 (1.4%) |

Thresholds use the selected training regions; reported outcomes use held-out regions.
These are exploratory existing-cohort estimates, not new-seed validation.
