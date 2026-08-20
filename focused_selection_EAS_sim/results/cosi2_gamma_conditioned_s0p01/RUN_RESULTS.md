# CoSi2 Gamma-SMC conditioned-neutral results

EAS neutral loci match the selected present population count exactly, 
but their mutation ages come from a fixed global importance-resampled 
reverse-time approximation. Han neutral loci match the 1.1%-1.3% CHB 
gate at generation 645 and survive to the present, without matching the 
selected present-day AF. Neither null conditions on sampled AF/carriers.

## Primary Gamma-SMC pointwise results

| Demography | Threshold years | Selected score | Upper p | Two-sided p | Descriptive AUC |
|---|---:|---:|---:|---:|---:|
| EAS | 1000 | 0.00113516 | 0.0891089 | 0.178218 | 0.92 |
| EAS | 4500 | 0.00390572 | 0.485149 | 0.970297 | 0.52 |
| EAS | 10000 | 0.0243482 | 0.435644 | 0.871287 | 0.57 |
| EAS | 20000 | 0.08 | 0.643564 | 0.732673 | 0.36 |
| EAS | 30000 | 0.129687 | 0.722772 | 0.574257 | 0.28 |
| EAS | 40000 | 0.169351 | 0.831683 | 0.356436 | 0.17 |
| EAS | 50000 | 0.200981 | 0.891089 | 0.237624 | 0.11 |
| Han | 1000 | 0.00020202 | 0.871287 | 0.633663 | 0.22 |
| Han | 4500 | 0.002886 | 0.405941 | 0.811881 | 0.6 |
| Han | 10000 | 0.0254161 | 0.227723 | 0.455446 | 0.78 |
| Han | 20000 | 0.103155 | 0.554455 | 0.910891 | 0.45 |
| Han | 30000 | 0.186474 | 0.633663 | 0.752475 | 0.37 |
| Han | 40000 | 0.267754 | 0.554455 | 0.910891 | 0.45 |
| Han | 50000 | 0.340404 | 0.485149 | 0.970297 | 0.52 |

## Correlation-aware seven-threshold minP

| Demography | Upper omnibus p | Two-sided omnibus p |
|---|---:|---:|
| EAS | 0.287129 | 0.514851 |
| Han | 0.623762 | 0.910891 |

## Conditioning diagnostics

| Demography | Statistic | Tail | Selected | Extreme / 100 | Plus-one p |
|---|---|---|---:|---:|---:|
| EAS | birth_generation | lower | 1475 | 4 | 0.049505 |
| Han | present_chb_population_af | upper | 0.195139 | 0 | 0.00990099 |

With one selected locus per demography, AUC is its neutral-bank 
midrank, not a stable ROC power estimate. Gamma TMRCA is a genealogy 
quantity and is not the same as the EAS mutation-age diagnostic.
