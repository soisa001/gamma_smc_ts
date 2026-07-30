# Posterior mean versus median at 10 kb and mu=1.29e-9

## Answer

No: using the posterior median was not the reason selection was missed.

The previously reported primary statistic was already the soft mean posterior
probability, `mean_i P(TMRCA_i < 4,500 years)`. The decoder's
`--recent-call` option changes only the hard per-pair calls summarized by
`n_recent` and `frac_recent`; it does not change that soft probability.

This study decodes the same selected tree and each of 100 newly simulated nulls
twice, once with a posterior-mean hard call and once with a posterior-median
hard call:

- Tree-sequence truth at 5 Mb is `196/2,000 = 0.098`, and no truth-null
  replicate is as large (`p=1/101`).
- The decoded soft probability is `2.72040768e-5` at 5 Mb, below the
  decoded-null mean `6.59522290e-5`. Eighty-nine of 100 nulls are at least as
  large, so `p=(89+1)/(100+1)=0.891089`.
- The posterior-mean hard-call fraction is zero at every one of the 1,000
  selected positions and every position in all 100 null profiles. Its center
  test is therefore `p=1`.
- The posterior-median hard-call fraction is also identically zero in the
  selected and null profiles, also giving `p=1`.
- The soft columns are bit-for-bit identical between the paired call-rule
  decodes, as expected. The complete mean-call and median-call tables are
  numerically identical here because neither hard rule calls a recent pair.

Thus changing median to mean does not rescue the sweep. The failure is upstream
of that summary choice: at the requested low mutation rate, the decoded
posterior assigns far too little recent-coalescence signal even though the
genealogy contains a strong, localized truth peak.

The soft pointwise scan has 47/1,000 raw `p<0.05` positions, merged into 28
descriptive tracts, but no position has BH `q<0.05`. These scattered raw hits
do not support the selected center and are compatible with the expected number
of pointwise false positives.

## Fresh simulation design

- Demography: ancestral `Ne=10,000`, present `Ne=20,000`, with growth 100
  generations ago.
- Region and panel: 10 Mb and 2,000 diploids; one within-individual haplotype
  pair per diploid.
- Selected allele: single-copy origin at 5 Mb, 180 generations (4,500 years)
  ago, additive `s=0.05`, `h=0.5`, and retained only at population AF >=0.30.
- Rates: mutation `1.29e-9` and recombination `1e-8` per bp per generation.
- Decoder: 10 kb stride (exactly 1,000 positions), 1 kb cache, no
  segregating-site output, and one decoder thread per replicate.
- Calibration: 100 new truth-null simulations and a separate set of 100 new
  simulated-and-decoded nulls.

The high-frequency rejection search used base seed `930777` and accepted
zero-based attempt 1,468, whose derived SLiM seed is `945457`. Population AF is
`0.3216`, sample AF is `0.31125`, and the sample contains 194 hom-alt, 857
heterozygous, and 949 hom-ref diploids.

The truth-null base seed is `525252`. These seed series are disjoint from the
historical selected run. A preliminary run based on seed `920241` was excluded
because its accepted attempt derived the already-used historical SLiM seed
`935381`; none of that preliminary run's outputs are included here. The final
selected tree SHA-256 is
`75954f2dca55cfcc226785e1d9cacd9770ad96ae14e6fb7c80c604ccbd99b82b`.

## Runtime and workers

`run-native-study` was invoked without `--workers` and recorded
`neutral_workers=12`, directly validating the new 12-worker default. The
matched two-rule native study completed in 496.215 seconds:

- selected mean-call decode: 7.874 seconds;
- selected median-call decode: 7.752 seconds;
- null mean-call decodes: 9.475 seconds on average;
- null median-call decodes: 9.469 seconds on average.

Live spot checks found exactly 12 decoder processes. Their aggregate RSS varied
from about 5.3 to 7.2 GiB across sampled batches; the Python driver and retained
result frames added further overhead. These are spot observations, not a
formal peak-memory measurement. The 1 kb cache remains the default: the
separate cache benchmark found no decode speedup above 1 kb and nearly linear
memory/startup growth.

## Main outputs

- `posterior_mean_median_rule_comparison.png`: matched soft, posterior-mean,
  posterior-median, and tree-truth comparison.
- `posterior_summary_rule_comparison_metrics.json`: machine-readable center,
  error, pointwise, and region metrics for all three statistics.
- `posterior_summary_rule_center_comparison.tsv`: compact three-row comparison.
- `neutral_decoded_mean_call_profiles.tsv.gz` and
  `neutral_decoded_median_call_profiles.tsv.gz`: 100 x 1,000 matched null rows.
- `comparison_*_spatial_calibration.tsv`: separate soft, mean-call, and
  median-call calibration tables.
- `decoded_study_metrics.json`: design, timing, worker count, selected-center,
  error, and spatial metrics.
- `run_manifest.json`: commands, versions, hashes, dimensions, seed audit, and
  validation assertions.

The code and synthetic simulations were validated locally on Windows and
WSL2. No controlled All of Us Workbench dataset was accessed, so controlled
data paths, quotas, and production-scale runtime remain unverified.
