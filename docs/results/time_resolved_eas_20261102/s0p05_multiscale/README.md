# Bounded no-introgression EAS time-scan pilot

Completed 2026-11-02. Independent whole-region calibration at alpha=0.05; time-window widths [1, 2, 5, 10] kyr.

Cohort counts: {'calibration': 99, 'development': 40, 'neutral_test': 40, 'selected': 20}

| Method | Neutral FPR | Any-call power | Joint space/time TPR | Center overlap TPR | Mean temporal F1 | Event F1 |
|---|---:|---:|---:|---:|---:|---:|
| truth | 0/40 (0.0%) | 20/20 (100.0%) | 16/20 (80.0%) | 14/20 (70.0%) | 0.290 | 0.348 |
| median | 4/40 (10.0%) | 20/20 (100.0%) | 0/20 (0.0%) | 0/20 (0.0%) | 0.000 | 0.000 |
| mass | 0/40 (0.0%) | 19/20 (95.0%) | 18/20 (90.0%) | 16/20 (80.0%) | 0.263 | 0.480 |

## Design and limitations

- Independent calibration and neutral test regions; no time/pair pseudoreplication.
- Test is maximum over all genomic positions and 182 time intervals; widths=[1, 2, 5, 10] kyr.
- Event matching requires a cell within 100 kb of the selected center and in [20,30) kya.
- Temporal F1 is computed at the selected coordinate, with no calls scored zero.
- This is a 1-Mb bounded pilot, not the planned 10-Mb/10,000-pair benchmark.
- Standing allele chosen on a valid onset genealogy near 5% frequency; no final-frequency conditioning.
- Analytic Gamma mass is a marginal posterior expectation; full-region calibration handles dependence.

Settings: {'schema': 'time-resolved-pilot/v1', 'scaling': 1, 'length': 1000000, 'diploids': 100, 'pairs': 1000, 'selection': 0.05, 'slim': '/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native/bin/slim', 'decoder': '/mnt/d/phase2simselection/code/bin/gamma_smc', 'base_seed': 20380101, 'pair_seed': 1729, 'stride': 10000, 'generation_time': 25, 'cutoffs': [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000, 40000, 41000, 42000, 43000, 44000, 45000, 46000, 47000, 48000, 49000, 50000], 'source_sha256': 'e422a72f33de498a79fa3e1cd4d0bd9ce436a5fbb0511cba76a27581be2037cc', 'slim_sha256': '7af149c29dc65d4aceb0dea19b2745b385eee4c254a71e18ebfff8f76735751a', 'decoder_sha256': '05bc30adeb6a93f08a4a51d5d80eb1e26b8b1ce1c728273dde7a1b28702dbd54'}

Rates, Wilson 95% intervals, paired differences and uncertainty are in `performance.json`. Raw simulation trees/posteriors remain in the local campaign directory; only aggregate simulation results and provenance are versioned.
