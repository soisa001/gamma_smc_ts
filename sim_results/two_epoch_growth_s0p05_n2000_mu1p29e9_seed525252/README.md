# Fresh two-epoch truth/null simulation at mu=1.29e-9

This directory is the independent truth/null input for the matched posterior
mean-versus-median native-decoder study.

- Design: 10 Mb, 2,000 diploids, ancestral `Ne=10,000`, present
  `Ne=20,000`, and the size change 100 generations ago.
- Rates and selected-variant criteria: `mu=1.29e-9`, recombination rate
  `1e-8`, age 180 generations (4,500 years), and additive `s=0.05`.
- Reproducibility: base seed `525252`, 100 neutral replicates, 20 workers, and
  228.357 seconds elapsed.
- Null truth: theoretical recent fraction `0.00647892`; observed mean
  `0.00665000`, bias `0.000171079`, and RMSE `0.00207528`.
- The validator's internal lower-frequency selected check accepted
  deterministic attempt 39 (population AF `0.107825`) and had a center truth
  fraction of `0.0125`, with `p=1/101`. This selected check was not used as the
  high-frequency pseudo-data in the decoder comparison.

`neutral_statistics.tsv` has 100 complete rows. `metrics.json` records the
design, calibration, sample counts, seeds, worker count, and elapsed time.
