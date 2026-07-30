# High-frequency selected validation at 10 kb and mu=1.29e-9

## Result

The requested configuration does **not** recover the selected center with the
native decoder. Tree-sequence truth at 5 Mb is
`P(TMRCA < 4,500 years)=0.138`, but the decoded value is
`0.0000918843`, compared with a decoded-neutral mean of `0.0000706178`.
Twenty-four of 100 matched decoded nulls are at least as large, giving
`p=(24+1)/(100+1)=0.247525`.

This negative result is not a 10 kb-stride or native-decoder regression. With
the same selected genealogy, current native binary, 10 kb stride, and 1 kb
cache, the historical `mu=1.25e-8` mutation overlay gives a center decoded value
of `0.0880513`. The requested `mu=1.29e-9` overlay has only 5,436 retained
segregating sites instead of 53,186 and gives `0.0000918843`. The controlled
comparison is in `mutation_rate_diagnostic.tsv`.

## Design

- Two-epoch demography: ancestral `Ne=10,000`, present `Ne=20,000`, with the
  size change 100 generations ago.
- Region and panel: 10 Mb, 2,000 diploids, one within-individual haplotype pair
  per diploid.
- Selected allele: single-copy origin at 5 Mb, 180 generations (4,500 years)
  ago, additive `s=0.05` and `h=0.5`.
- Mutation and recombination rates: `1.29e-9` and `1e-8` per bp per generation.
- Matched decoder rates: theta `4*10,000*1.29e-9 = 5.16e-5`,
  rho/theta `1e-8/1.29e-9 = 7.7519379845`; scaled rho remains `0.0004`.
- Decoder grid/cache: exactly 1,000 positions at 10 kb stride, no segregating-
  site output, 1 kb flow-field cache, one decoder thread per replicate and four
  concurrent null workers.
- Calibration: 100 newly simulated and decoded neutral replicates.

The selected genealogy is the deterministic zero-based attempt 2,514 from the
earlier high-frequency study. It was re-materialized and validated rather than
re-screening 2,515 trajectories. After removing sites and mutations, its tskit
tables are exactly equal to the historical tree sequence: 40,196 trees and
4,000 sampled haplotypes. The mutation overlay was regenerated at the requested
rate.

## Tree-sequence truth

- Population AF is `0.379525`; sampled AF is `0.38275`.
- The sample has 274 hom-alt, 983 heterozygous, and 743 hom-ref diploids.
- At 5 Mb, 276/2,000 pairs coalesce within 180 generations, so the true recent
  fraction is `0.138`. No truth-null replicate is as large (`p=1/101`).
- Hom-alt mean TMRCA is 161.4 generations (95% CI 159.0-163.9), versus
  19,008.0 generations (95% CI 18,380.8-19,656.6) for hom-ref diploids.
- The separate 100-replicate truth null has mean `0.006555`, versus the
  two-epoch theoretical value `0.006479`; bias is `0.0000761` and RMSE is
  `0.001669`.

## Native decode and spatial calibration

- The selected VCF contains 5,436 retained SNPs and 13,550 decoder segments.
  The decoder used 8.121 seconds internally (8.210 seconds including the
  wrapper) and 0.585 GB peak RSS.
- The 100 null decodes average 8.383 seconds; the sum across selected plus all
  null decodes is 846.55 seconds. Four-way simulation/decode plus
  postprocessing completed in 621.59 seconds.
- Across the 1,000-position profile, decoded-minus-truth bias is `-0.0133235`,
  MAE `0.0133235`, RMSE `0.0230084`, and Pearson correlation `0.0397`.
- There are 49/1,000 raw pointwise `p<0.05` positions, merged into 37
  descriptive tracts, but none has BH `q<0.05`.
- The selected site itself is not pointwise significant (`p=0.247525`).
  Global positive-window density is 49 versus a null mean of 40
  (`p=0.178218`); within +/-500 kb it is 4/101 windows versus a null mean of
  4.04 (`p=0.534653`). The raw tracts therefore do not support a selected
  region in this run.

The older 1 kb result is not a stride-only comparator: it used an approximately
9.69-fold higher mutation rate and the official v0.2 container. The diagnostic
above holds backend, stride, cache, and genealogy fixed, isolating the matched
mutation-rate configuration.

## Main outputs

- `decoded_study_metrics.json`: machine-readable design, timing, error, center,
  spatial, and density-calibration metrics.
- `run_manifest.json`: commands, versions, hashes, dimensions, and validation
  assertions.
- `neutral_decoded_recent_probability_profiles.tsv.gz`: 100 x 1,000 complete
  decoded-null rows.
- `decoded_spatial_pointwise_calibration.tsv`: selected/null pointwise
  calibration on the 10 kb grid.
- `selected_decoded_vs_truth.tsv.gz`: decoded and tree-truth profiles.
- `selected_decoded_vs_neutral_pvalues.png`: full-profile and local
  calibration plot.
- `decoded_scan_density_calibration.png`: global and local density nulls.
- `mutation_rate_diagnostic.tsv`: controlled historical-versus-requested
  mutation-overlay diagnostic.

The code was validated locally on Windows and WSL2. No controlled All of Us
Workbench dataset was accessed, so data-path, quota, and production-scale
runtime remain to be verified in Workbench.
