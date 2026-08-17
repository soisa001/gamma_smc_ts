# CoSi2 Gamma-SMC conditioned focal pair results

REF/REF and ALT/ALT denote unordered pairs from the allele-blind fixed 
100-haplotype all-pairs panel; they are pseudo-diploid pair classes, 
not the original biological diploid pairing. Probabilities were reduced 
post hoc from retained float32 Gamma shape/rate streams with SciPy.

The primary ALT/ALT-minus-REF/REF contrast uses the shared neutral bank 
in which both homozygous pair classes exist. Direct REF/REF and ALT/ALT 
tests use their maximal eligible banks and are secondary comparisons.

## Pointwise results plotted

| Role | Demography | Statistic | Best upper-tail threshold | Selected score | Upper p | Two-sided p | Neutral n |
|---|---|---|---:|---:|---:|---:|---:|
| primary_shared_homozygote_bank | EAS | p_hom_alt_minus_hom_ref | 1000 | 0.0330859 | 0.164948 | 0.329897 | 96 |
| primary_shared_homozygote_bank | Han | p_hom_alt_minus_hom_ref | 1000 | 0.0177119 | 0.84375 | 0.34375 | 63 |
| secondary_class_specific_bank | EAS | p_hom_alt | 1000 | 0.0410107 | 0.28866 | 0.57732 | 96 |
| secondary_class_specific_bank | EAS | p_hom_ref | 1000 | 0.00792481 | 0.958763 | 0.103093 | 96 |
| secondary_class_specific_bank | Han | p_hom_alt | 1000 | 0.0393645 | 0.765625 | 0.5 | 63 |
| secondary_class_specific_bank | Han | p_hom_ref | 40000 | 0.612489 | 0.03 | 0.06 | 99 |

## Whole-vector seven-threshold minP

| Demography | Statistic | Upper omnibus p | Two-sided omnibus p |
|---|---|---:|---:|
| EAS | p_hom_alt_minus_hom_ref | 0.28866 | 0.525773 |
| Han | p_hom_alt_minus_hom_ref | 0.953125 | 0.296875 |

Pointwise p-values use conservative ties and the plus-one correction. 
The minP test retains each neutral unit's complete seven-threshold vector.
