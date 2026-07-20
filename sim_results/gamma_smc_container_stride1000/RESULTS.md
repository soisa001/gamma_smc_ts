# Gamma-SMC v0.2 1 kb decoded validation

## Design

- Selected pseudo-empirical panel: the retained two-epoch `s=0.05` replicate,
  10 Mb long, selected variant added at 5.000 Mb 180 generations ago.
- Demography: present `N=20,000`, changing to ancestral `N=10,000` at 100
  generations before present.
- Samples: 2,000 diploids = 4,000 haplotypes = 2,000 within-individual pairs.
- Null: 100 independent matched neutral replicates, decoded identically.
- Rates: `mu=1.25e-8`, `r=1e-8`, `theta=0.0005`, `rho/theta=0.8`.
- Decoder: official `regevsch/gamma_smc:v0.2` container, digest
  `sha256:739f5ce1b06a4de75da1b9865747cc0d6375d9443390ca969b2e4320a4714274`.
- Output: `--only_within --output_at_hets=false --output_at_stride 1000`.
- Statistic: across-pair mean posterior `P(TMRCA < 4,500 years)` at each 1 kb
  output position.

Gamma-SMC emitted 10,000 selected positions. One neutral replicate lacked the
final 9.999 Mb position, so pointwise calibration uses the 9,999 positions
present in every replicate. No simulation or internal position was pooled as an
independent null replicate.

## Practical performance

On a four-logical-CPU AMD EPYC 7763 GitHub runner, the selected 10 Mb panel took
15.61 seconds inside the container and 1.80 seconds to stream/reduce. Gamma-SMC
reported 0.791 GB peak memory, 52,127 segregating sites, 62,080 segments, and a
125,467,705-byte raw posterior. The independent one-null benchmark took 15.75
seconds for the neutral decode.

The complete 100-null job used four concurrent simulation/decode workers. It
finished every simulation and decode in 3,290.3 seconds (54.84 minutes), then
encountered the single missing edge coordinate during calibration. The recovered
profiles were finalized locally without rerunning. Thus the measured bottleneck
is repeated SLiM/tree-sequence reconstruction and 4,000-haplotype VCF writing,
not Gamma-SMC on an already prepared VCF.

A simple prepared-input extrapolation is about 78 minutes for 3 Gb at the same
2,000 within-individual pairs and 1 kb stride, before chromosome-boundary and I/O
effects. Raw posteriors would total about 37.6 GB if retained; this workflow
reduces and deletes them replicate by replicate.

## Exact selected base

At 5.000 Mb:

- selected decoded mean probability: 0.0305024;
- neutral mean: 0.0234616;
- neutral 95% interval: 0.0153090 to 0.0370474;
- null exceedances: 13 of 100;
- corrected upper-tail Monte Carlo p-value: `14/101 = 0.1386`;
- BH q-value: 0.5729.

Therefore the exact selected base is not significant in this decoded run.

## Spatial signal

There are 872 raw pointwise `p<0.05` windows across 10 Mb, merged into 51
tracts. Within +/-500 kb of the selected site there are 309 positive windows.
No individual window survives BH at q<0.05, so the raw highlighted tracts should
not be read as independently genome-wide-significant loci.

The density of positive windows is itself unusual when each neutral replicate
is scanned leave-one-out:

- full 10 Mb: selected 872; neutral mean 399.96, median 369;
  regional-density p=`3/101=0.0297`;
- +/-500 kb: selected 309; neutral mean 40.04, median 29;
  regional-density p=`1/101=0.0099`.

The strongest nearby pointwise tracts are displaced to the left of the selected
base: 4.6945-4.7815 Mb, 4.8135-4.8405 Mb, 4.8485-4.8685 Mb, and
4.9485-4.9785 Mb. This run supports a regional enrichment signal but does not
accurately localize selection to the causal 5.000 Mb coordinate.

## Decoder error against tree-sequence truth

Across the selected profile, decoded probability minus the true fraction of
within-individual pairs coalescing before 180 generations has bias 0.01730,
MAE 0.01730, RMSE 0.01806, and correlation 0.358. Gamma-SMC overestimates the
absolute recent-coalescence probability in this setup. Matched decoding of the
neutral simulations is therefore required; an uncalibrated probability cutoff
would be misleading.

## Interpretation

The actual decoder works and is fast for prepared VCFs. In this difficult recent
sweep example it has regional power when evidence is aggregated across linked
windows, but inadequate point localization at the causal base and no
BH-significant individual windows. For empirical AoU scans, use the regional
density statistic (or a prespecified tract statistic) with matched demographic,
mask, mutation-rate, and recombination-rate simulations. Do not report isolated
raw pointwise p-values as final discoveries.
