# No-introgression EAS sweep: bounded local run

## Execution outcome

This report freezes the exploratory local campaign completed on 2026-08-11.
All nine requested selection-by-EAS-demography specifications produced one
accepted 100-diploid sample panel, and all nine were decoded with Gamma-SMC.
Accepted sampled selected-allele frequencies ranged from 0.455 to 0.535 under
the fixed inclusive 0.45-0.55 contract.

| Selection coefficient | `q025` origin | Median origin | `q975` origin | Coverage |
|---:|---:|---:|---:|---:|
| 0.01 | 1,275 generations | 1,289 generations | 1,308 generations | 3/3 |
| 0.005 | 2,635 generations | 2,725 generations | 2,781 generations | 3/3 |
| 0.001 | 15,523 generations | 15,728 generations | 16,306 generations | 3/3 |

The two oldest weak-selection cells required long conditioning runs. The
`q025, s=0.001` cell accepted at sample AF 0.505 after two completed external
sample-panel draws and one interrupted internal trajectory; it completed
inside the declared 45-minute retry window. Attempts were never counted until
stdpopsim returned a complete sample panel, and no AF or genotype threshold was
relaxed.

## Descriptive focal-site results

The expected sweep direction is positive `alt/alt - ref/ref` normalized
integrated TMRCA-CDF area and negative mean-TMRCA difference. The CDF-area sign
agreed in 9/9 cells for both tree-sequence truth and Gamma-SMC. Both signs
occurred together in 8/9 cells for each source. The exception was the oldest
`q025, s=0.001` cell: its CDF-area contrasts remained positive (+0.130 truth,
+0.078 Gamma-SMC), while its mean-TMRCA contrasts were slightly positive
(+738.76 and +66.73 generations, respectively).

For the prespecified representative cell
`no_intro_s0p010_median_age1289g`:

| Source | Normalized integrated CDF area, alt-ref | Mean TMRCA, alt-ref |
|---|---:|---:|
| Tree-sequence truth | +0.548984 | -12,206.18 generations (-305,154.53 years) |
| Gamma-SMC | +0.410077 | -7,247.42 generations (-181,185.53 years) |

These are internal genotype-class contrasts from a selected simulation. The
area is the normalized integral of the TMRCA CDF over the requested 1-50-kya
cutoffs; it is **not** ROC AUC. There is no matched neutral null in this pilot,
so these values do not define a p-value, false-positive rate, or detection
power.

## Runtime and integrity scope

The accepted simulations were generated on Windows with Python 3.11.15,
stdpopsim 0.3.0, tskit 1.0.3, and SLiM 4.2.2. Gamma-SMC decoding ran under WSL
with the existing Linux AVX2 binary whose SHA-256 is
`e57182c105e216c5a3c4afad991d7f51343dbf21d4c1889b38af97a3c52e6d40`.
The bootstrap rebuilds and verifies a decoder for a fresh checkout; this local
result does not claim that the checked binary was freshly compiled during this
campaign.

The final result manifest contains nine plotted specifications, nine complete
decoder rows, artifact-relative paths, and SHA-256 checksums. Bulk tree
sequences and decoder intermediates remain restartable under the ignored
`work/` directory.

All manifest-declared artifacts were re-read and SHA-256-verified after final
aggregation; compressed result tables also passed gzip CRC validation.

## Pilot limitations

- The `q025`, median, and `q975` Ne histories are pointwise quantiles over 100
  bootstrap curves, not simultaneous confidence bands or coherent posterior
  demographic draws.
- Origin ages are deterministic trajectory-derived inputs. At `Q=10`, the
  one-copy starting frequency is approximately `Q/(2Ne)`; the ages are not
  direct mutation-age estimates and are not the paper's complete fixed-age
  grid.
- `Q=10` and `slim_burn_in=0.1` are speed-oriented approximations. Results need
  sensitivity checks with less rescaling and longer forward burn-in before
  biological interpretation.
- The 50% target remains a placeholder. No qualifying genome-wide distribution
  of archaic-specific, no-ILS Han SNV frequencies was identified for replacing
  it in this run.
- One trajectory per cell is suitable for method inspection only. Replicates,
  matched neutral simulations, and null-calibrated p-values remain follow-up
  work.

The machine-readable source of truth is
`results/results_manifest.json`; methods and the reusable run contract are in
[`METHODS.md`](METHODS.md).
