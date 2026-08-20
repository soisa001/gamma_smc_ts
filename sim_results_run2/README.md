# run2

Can selection at an introgressed locus be detected from the within-individual
TMRCA distribution alone?

**Status: code written and pre-flight validated. No simulations have been run.**
[`FRAMEWORK.md`](FRAMEWORK.md) is the pre-registration; this file is how to run it.

## Design in one paragraph

Two arms. **CHB** is `stdpopsim` `AncientEurasia_9K19`: the focal allele is placed
on every genome arriving in the catalog's own 2.96% Neanderthal → Loschbour
pulse, so it starts at exactly the pulse proportion on genuinely archaic
haplotypes. **EAS** is the tracked PHLASH median curve with no archaic source at
all: the allele is placed on a uniformly random 2.96% of genomes at the same
tick. Both then get `s = 0.01`, `h = 0.5` continuously to the present, with **no
allele-frequency conditioning anywhere** — the final frequency is a measured
outcome, and replicates that lose the allele are kept. The null is 100 plain
neutral replicates per arm: no mutation drawn, no fitness event, nothing
conditioned, so the focal base may be empty. The statistic is the **overall**
within-individual `P(TMRCA < x)` at the focal base over 100 sampled diploids,
and the p-value is its one-sided upper tail against the null.

## Running it

`validate` needs `stdpopsim == 0.3.0` but **not** the SLiM binary. It generates
the exact SLiM script each arm and mode would run and asserts the run2
invariants against it. Run this first and read the output.

```bash
python scripts/run_run2.py validate
```

Then, with SLiM 4.2.2 available:

```bash
python scripts/run_run2.py simulate --slim-bin /path/to/slim --workers 12
```

```bash
python scripts/run_run2.py analyze
```

```bash
python scripts/run_run2.py plot
```

A three-replicate smoke test before committing to all 400 runs:

```bash
python scripts/run_run2.py simulate --slim-bin /path/to/slim --index 0 --index 1 --index 2
```

## What validate asserts

| Check | Gate |
|---|---|
| `condition_on_allele_frequency` is empty in every generated script — zero rejection sampling | — |
| The focal draw is scheduled exactly one tick after the pulse (offset 453 vs 454; 2265 vs 2270 generations ago) | G1, G2 |
| The Han split tick is shared by the two fitness callbacks, so selection has no gap | G1 |
| One drawn mutation and one fitness callback per carrying population in `selected`; none of either in `neutral` | — |
| The CHB arm selects carriers by `inds.migrant`; the EAS arm places `0.0296` directly | G2 |
| The reserved-base guard removes exactly one base of mutational mass and is idempotent | G4 |
| All 400 replicate seeds are unique | G5 |

Two further gates need SLiM and therefore the smoke run: **G0** (environment and
timing) and the in-simulation half of **G2** — the SLiM script records the
realized post-pulse frequency into tree-sequence metadata and errors out if no
introgressing individuals are present, so a wrong `migrant` assumption fails
loudly on replicate 0 rather than silently.

## Layout

```
sim_results_run2/
  FRAMEWORK.md                 pre-registration
  README.md                    this file
  config/run2_{chb,eas}.json   frozen parameters and the full seed list
  validation_report.json       the validate phase output
  <arm>/validation/            the generated SLiM scripts, as reviewed
  <arm>/{selected,neutral}/replicates/   per-replicate intermediates (gitignored)
  <arm>/results/               the study result tables
  <arm>/figures/               figures 1-6
  cross_arm/                   figure 7 and the power comparison
```

Per-replicate directories are gitignored: they are regenerable from the seeds in
`config/`, which are a pure function of `(arm, mode, index)`.

## Code

| Module | Role |
|---|---|
| [`run2_config.py`](../python/gamma_smc_aou/run2_config.py) | every frozen parameter and the seed rule |
| [`run2_models.py`](../python/gamma_smc_aou/run2_models.py) | demography, contig, events, tick arithmetic, the two SLiM patches |
| [`run2_simulate.py`](../python/gamma_smc_aou/run2_simulate.py) | one replicate end to end, and the tree-truth statistics |
| [`run2_analysis.py`](../python/gamma_smc_aou/run2_analysis.py) | null distribution, p-values, power, null calibration |
| [`run2_plots.py`](../python/gamma_smc_aou/run2_plots.py) | figures 1-7 |
| [`run2_workflow.py`](../python/gamma_smc_aou/run2_workflow.py) | phases and the validation gates |

Tests: `tests/test_run2_models.py`, `tests/test_run2_analysis.py`,
`tests/test_run2_simulate.py` (46 tests, none requiring SLiM).
