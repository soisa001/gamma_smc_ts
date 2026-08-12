# EAS demographic resource provenance

`EAS.npz` is an exact byte-for-byte copy of the tracked artifact
`mvn/EAS.npz` in `soisa001/SimulatePhase2`.

- SHA256: `a5bdfd843629d48050cd56a84da4d70e356d0339017d3c1f67be646d727db928`
- Size: 4,110,207 bytes
- Artifact schema: `phlash.aou.log-ne-mvn/v1`
- Population: `EAS`
- Source repository: `git@github.com:soisa001/SimulatePhase2.git`
- Source path: `mvn/EAS.npz`
- Content-introducing commit: `4dd5f1afe66efcb06ecb131f2b339339c6ac34bd` (`Install validated PHLASH MVNs and summary plots`, 2026-08-05)
- Local verification date: 2026-08-11

The NPZ contains:

| Array | Shape | Dtype | Meaning |
|---|---:|---|---|
| `schema` | scalar | Unicode | `phlash.aou.log-ne-mvn/v1` |
| `population` | scalar | Unicode | `EAS` |
| `time` | `(10000,)` | `float32` | geometric grid, 100-40,000 generations ago |
| `mean_log_ne` | `(10000,)` | `float32` | empirical mean log diploid Ne |
| `covariance_factor` | `(100, 10000)` | `float32` | low-rank factor for empirical log-Ne covariance |
| `bootstrap_ne` | `(100, 10000)` | `float32` | 100 retained posterior-median bootstrap curves in diploid Ne units |
| `jitter` | scalar | `float32` | independent log-Ne noise SD; exactly 0 |

The covariance rank is at most 99 because it is estimated from 100 centered
curves. This study uses empirical pointwise 2.5%, 50%, and 97.5% quantiles of
`bootstrap_ne`; it does not reinterpret those pointwise envelopes as a
simultaneous confidence band.

Machine-readable provenance is in [`EAS.provenance.json`](EAS.provenance.json).
