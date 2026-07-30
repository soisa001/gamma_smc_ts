# Why the selected sweep was recovered before

## Conclusion

The previous recovery was not caused by the 1 kb output stride, and the failed
low-mutation test was not a Monte Carlo power failure. The decisive change was
the matched mutation-rate configuration:

- At `mu=1.29e-9`, the selected center is below the decoded-null mean at both
  1 kb and 10 kb. The pointwise p-values are 0.901 and 0.891.
- At `mu=1.25e-8`, the same selected genealogy is recovered at both strides.
  At 1 kb, the selected and null-mean statistics are `0.071757` and
  `0.021938`; at 10 kb they are `0.071845` and `0.021981`. Both independently
  calibrated cells have 0/100 null exceedances and `p=1/101=0.009901`.

The low and high overlays contain 5,483 and 53,528 segregating sites,
respectively. The mutation rate changes by 9.69-fold and the observed site
count by 9.76-fold. The high-mutation posterior center is 3,066-fold larger
than the low-mutation center at the same 1 kb stride.

## Matched test

All new cells use the same selected genealogy: 39,515 trees, 4,000 sampled
haplotypes, selected population AF 0.3216, sampled AF 0.31125, additive
`s=0.05`, allele age 180 generations, and selected attempt 1,468 from base
seed 930777. Sites and mutations were cleared before the high-rate binary
mutation overlay, and the mutation-free table collections are exactly equal.
Decoder rates were matched to each mutation rate:

| Mutation rate | theta | rho/theta | Selected sites |
|---:|---:|---:|---:|
| `1.29e-9` | `5.16e-5` | `7.7519379845` | 5,483 |
| `1.25e-8` | `5.0e-4` | `0.8` | 53,528 |

The primary statistic in every comparison is the mean across pairs of the
posterior probability `P(TMRCA < 180 generations)`. It is not the fraction of
pairs called recent under a posterior-mean or posterior-median hard rule.
In the matched low-mutation 10 kb comparison, both hard-call fractions were
exactly zero for selected and null profiles, giving `p=1.0` under both rules.
Changing the hard summary from median to mean therefore does not recover the
sweep; the more informative soft-probability statistic still gives `p=0.891`.

| Configuration | Selected center | Null mean | Null 95% range | Exceedances | p |
|---|---:|---:|---:|---:|---:|
| low mu, 1 kb | `0.00002340` | `0.00006246` | `0.00001925-0.00016328` | 90/100 | `0.90099` |
| low mu, 10 kb | `0.00002720` | `0.00006595` | `0.00002118-0.00016630` | 89/100 | `0.89109` |
| high mu, 1 kb | `0.07175696` | `0.02193752` | `0.01478478-0.03198568` | 0/100 | `0.00990` |
| high mu, 10 kb | `0.07184464` | `0.02198136` | `0.01481146-0.03204353` | 0/100 | `0.00990` |

Changing the stride does alter the approximate Gamma-SMC posterior slightly
because output positions split flow-field segments. The effect is modest and
cannot explain recovery:

- Low mu: center changes from `2.340e-5` at 1 kb to `2.720e-5` at 10 kb; both
  fail with p approximately 0.9.
- High mu: center changes from `0.071757` at 1 kb to `0.071845` at 10 kb.
- Across the 1,000 shared low-mutation selected positions, Pearson correlation
  is 0.958 and the largest absolute probability difference is `7.45e-5`.
- Across the 100,000 shared low-mutation null rows, probabilities correlate
  0.958 between strides.
- The corresponding high-mutation correlations are 0.999997 for the selected
  profile and 0.999990 across 100,000 null rows. The largest absolute
  probability differences are `2.72e-4` and `4.09e-4`, respectively.

The earlier successful official-v0.2 result also used `mu=1.25e-8`; it gave a
selected center of `0.09169`, null mean `0.02346`, and `p=1/101`. The new
native, new-seed result reproduces that qualitative result, so neither the old
container backend nor its 1 kb stride is required for recovery.

## Carrier conditioning

The visible carrier-class truth contrast is real:

| Pair class | n | True recent fraction | True mean TMRCA (generations) |
|---|---:|---:|---:|
| hom-alt | 194 | `1.0000` | `159.9` |
| heterozygous | 857 | `0.0000` | `51,287.1` |
| hom-ref | 949 | `0.002107` | `51,272.7` |
| all pairs | 2,000 | `0.0980` | `46,320.9` |

However, the selected hom-ref class is not equivalent to an unconditioned
neutral replicate. Across the 100 truth nulls, mean center TMRCA is 19,376.7
generations and the empirical 95% range is 8,252.9-34,834.5. Conditioning on
the ancestral allele at the selected site instead identifies the deep
background around the sweep.

Carrier-stratified decoding explains the aggregate result:

| Mutation rate | Hom-alt decoded P(recent) | Heterozygous | Hom-ref |
|---:|---:|---:|---:|
| `1.29e-9` | `0.0001370` | `0.00001874` | `0.00001239` |
| `1.25e-8` | `0.63331` | `0.01162` | `0.01145` |

At high mutation rate Gamma-SMC recovers the hom-alt carrier signal, and its
9.7% sample weight lifts the all-pair statistic to 0.07176. At low mutation
rate the hom-alt posterior collapses to 0.000137, so averaging is not the
cause of failure; the carrier-specific information has already been lost by
the decoder.

## P-value power

One hundred nulls are sufficient for the pointwise and regional question. The
smallest attainable Monte Carlo p-value is 1/101=0.00990, and the high-mutation
cells at both strides reach it. In the low-mutation 10 kb cell, 89/100 nulls
exceed the selected center. The exact 95% interval for the underlying
exceedance probability is 0.812-0.944, so adding more nulls would refine a
large p-value; it would not make this center significant.

One hundred nulls are not sufficient for BH correction across 10,000
pointwise windows. Even the recovered high-mutation scan has 0 BH-significant
windows. Its regional-density tests are significant: 2,467/10,000 positive
windows globally and 863/1,001 within +/-500 kb, both `p=1/101`, versus null
means of 400 and 40.04. The independent 10 kb scan agrees: 250/1,000 positive
windows globally and 86/101 locally, versus null means of 40 and 4.04; both
regional p-values are again `1/101`.

## Optimizations found during the test

- Null preparation now uses a process pool and decoding a thread pool in
  bounded 12-replicate batches. This avoids the Python/GIL bottleneck in
  mutation overlay and VCF writing without retaining all temporary inputs.
- Transient gzipped VCFs use compression level 1. On a controlled 1 Mb,
  high-mutation slice, level 1 took 0.602 s versus 7.232 s for level 9
  (12.01-fold faster); the file was 2.04-fold larger (2.145 MB versus
  1.052 MB).
- Pointwise null calibration now pivots once to a position-by-replicate
  matrix. The complete 1,000,000-row, 10,000-position low-mutation calibration
  runs in 0.293 s instead of repeatedly rescanning the full null table.
- The optimized 100-null high-mutation 1 kb study completed in 341.0 s with
  12 workers, and the complete high-mutation 10 kb study completed in 308.4 s.
  These are observed workflow times, not a controlled before/after benchmark.
- The 1 kb flow-field cache was retained. The separate cache benchmark found
  larger caches slower and more memory-intensive.

## Main files

- `analysis/selection_recovery_hypothesis_tests.png`: four-panel summary.
- `analysis/stride_mutation_factorial_center.tsv`: matched center results.
- `analysis/carrier_class_truth_center.tsv`: carrier-conditioned tree truth.
- `analysis/carrier_class_decoded_center.tsv`: carrier-conditioned decodes.
- `low_mu_stride1kb/`, `high_mu_stride1kb/`, and `high_mu_stride10kb/`:
  complete selected and 100-null profiles, pointwise scans, density tests,
  plots, and metrics.
- `run_manifest.json`: seeds, hashes, commands, assertions, and dimensions.

The high-versus-low comparison changes observed site density and the matched
decoder rate parameters together, as required for a scientifically matched
mutation-rate analysis. It does not separately identify which internal
Gamma-SMC approximation fails in the sparse low-rate regime. No controlled
All of Us Workbench data were used. Validation completed with 59 passing
Windows tests (30 platform skips), 89/89 passing WSL/native-decoder tests, and
an independent profile/p-value audit.
