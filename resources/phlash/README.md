# Aggregate PHLASH demographic histories

AFR, AMR, EUR, MID, and SAS NPZ files are the published aggregate
`phlash.aou.log-ne-mvn/v1` resources from
`git@github.com:soisa001/phlash_ld.git`, `artifacts/aou-mvn/`, reviewed at
70e177ea0e77fe38d408eb8042e13a44c94d15c9. They contain 100 fit-level population
history curves and aggregate log-Ne covariance information, without individual
genotypes. EAS is already tracked at `no_introgression_EAS_sim/resources/EAS.npz`
and is byte-identical to the EAS artifact from that source.

Expected SHA-256 values are in `configs/all_population_neutral.json`. All six
resources were validated with `docs/manuscript/summarize_phlash.py`, including
ZIP integrity, sidecar hashes/sizes, population identity, positive finite Ne,
time grids, and reconstruction of log-Ne mean/covariance factors. Historical
fit settings are discussed in `docs/manuscript/evidence_and_open_items.md`.
