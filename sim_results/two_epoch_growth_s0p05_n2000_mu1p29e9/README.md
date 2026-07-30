# Two-epoch truth simulation at mu=1.29e-9

This directory is the deterministic truth/null input for the 10 kb
high-frequency native-decoder study.

- 10 Mb, 2,000 diploids, ancestral `Ne=10,000`, present `Ne=20,000`, size
  change 100 generations ago.
- `mu=1.29e-9`, recombination rate `1e-8`, selected-variant age 180
  generations, `s=0.05`, base seed 515151.
- 100 neutral replicates: theoretical recent fraction `0.00647892`, observed
  mean `0.006555`, bias `0.0000760793`, RMSE `0.00166891`.
- The lower-frequency selected validation accepted deterministic attempt 8
  (population AF `0.102125`) and has truth center fraction `0.0185`,
  `p=1/101`.

`neutral_statistics.tsv` has 100 complete rows. `metrics.json` records the
design, calibration, sample counts, seeds, worker count, and elapsed time.
