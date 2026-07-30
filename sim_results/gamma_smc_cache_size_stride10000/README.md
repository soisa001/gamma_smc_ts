# Gamma-SMC cache-size benchmark at 10 kb stride

## Decision

Keep `--cache_size 1000`.

Increasing the cache did not speed the decode phase, while cache construction
time and peak memory rose nearly linearly. Relative to 1 kb, 2 kb increased
wall time by 41.8% (29.5% lower throughput) and used 0.570 GB more peak RSS;
5 kb increased wall time by 159.8% (61.5% lower throughput) and used 2.280 GB
more. A 500 bp cache cut wall time by 25.8% (1.348x throughput) and used
0.285 GB less, but it changed the segmentation and 13 of 1,000
`n_recent_4500` values by one call.
Cache size is therefore part of the numerical configuration, not just a
performance knob. Observed data and nulls must use the same value.

## Workload

- Code: `AOU_run_opt` commit
  `6139b7ee217f6650a376358d873d2bfdd9524a50`.
- Input:
  `../gamma_smc_container_stride1000_s0p05_af30/selected_s0p05_af30.vcf.gz`
  (10,702,259 bytes; SHA-256
  `4850f36fcd8e43a57c8426a917cadb0bb722ac81ac21c7ea2d505b5715259adc`).
- Panel: 2,000 diploids, one within-individual pair each.
- Region: 10 Mb; 53,186 retained segregating sites in this optimized reader.
- Output: 10 kb stride, 1,000 positions.
- Rates: theta `0.0005`, rho/theta `0.8`, mutation rate `1.25e-8`;
  threshold 4,500 years at 25 years/generation.
- Decoder: one thread, corrected `exp10` and backward-alignment defaults.
- Host: WSL2 Linux 6.18.33.2 on an AMD Ryzen Threadripper PRO 5945WX
  (12 cores/24 logical CPUs).
- Peak RSS and elapsed time: `/usr/bin/time -v`; decoder phase times: the
  decoder's own summary.

The mutation rate above belongs to the retained historical input and affects
time unscaling, not cache construction. The requested new simulation uses
`1.29e-9`; cache-size conclusions depend on the site spacing, panel, and flow
field rather than that reporting conversion.

Representative invocation:

```bash
/usr/bin/time -v bin/gamma_smc \
  --input sim_results/gamma_smc_container_stride1000_s0p05_af30/selected_s0p05_af30.vcf.gz \
  --input_format vcf --only_within \
  --scaled_mutation_rate 0.0005 \
  --recombination_to_mutation_ratio 0.8 \
  --unscaled_mutation_rate 1.25e-8 \
  --recent_threshold_years 4500 --generation_time 25 \
  --recent_summary cache1000.tsv \
  --output_at_hets=false --output_at_stride 10000 \
  --cache_size 1000 --threads 1
```

`benchmark.tsv` contains the timing/memory measurements.
`numerical_differences.tsv` compares each alternative with 1 kb.
