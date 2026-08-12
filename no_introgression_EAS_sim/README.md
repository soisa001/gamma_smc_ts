# No-introgression EAS sweep study

This directory runs the de novo EAS sweep design described in
[`METHODS.md`](METHODS.md). Raw tree sequences and decoder intermediates are
restartable under `work/`; compact summaries and cross-specification figures
are written under `results/`. The checked execution outcome and its limitations
are summarized in [`RUN_RESULTS.md`](RUN_RESULTS.md).

## Bootstrap

From a fresh or existing checkout of the run branch:

```bash
if [[ ! -d gamma_smc_ts/.git ]]; then
  git clone git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git fetch origin AOU_run_opt
git switch AOU_run_opt
git pull --ff-only origin AOU_run_opt
bash scripts/bootstrap_eas_sweep_study.sh --skip-tests
```

The bootstrap keeps the general repository environment at `.native` with SLiM
5.2, builds Gamma-SMC, and installs SLiM 4.2.2 separately at
`.native-stdpopsim`. The study runner uses the latter for stdpopsim-generated
SLiM scripts. A local 10-Mb diagnostic found that unmodified stdpopsim 0.3 code
fails under SLiM 5.2 at the removed `Subpopulation.genomes` API. A narrow
`genomes`-to-`haplosomes` transform ran but had still not accepted a trajectory
after 200 seconds; SLiM 4.2.2 accepted the same design/seed in 84.9 seconds.
Thus 4.2.2 is the validated study default; this is a workload-specific result,
not a general SLiM performance claim.

## Run the phases

Run phases explicitly from the repository root. Defaults are four concurrent
specifications; every active Gamma-SMC decoder uses four threads (up to 16
decoder threads when four specifications are ready).

```bash
bash no_introgression_EAS_sim/run.sh plan
bash no_introgression_EAS_sim/run.sh simulate
bash no_introgression_EAS_sim/run.sh decode
bash no_introgression_EAS_sim/run.sh plot
```

`decode` requires Linux x86_64 with AVX2. `all` runs the four phases in order,
but separate invocations are easier to restart and audit. To filter a pilot or
resume a subset, repeat `--spec` with a substring from `specifications.tsv`:

```bash
bash no_introgression_EAS_sim/run.sh simulate --spec s0p010_median
```

Override parallelism for a run:

```bash
EAS_SWEEP_WORKERS=8 EAS_SWEEP_THREADS=8 \
  bash no_introgression_EAS_sim/run.sh simulate
```

`work/` and per-specification bulk figure directories are ignored by Git.
`specifications.tsv`, study design/provenance records, compact result tables,
result manifests, and `results/figures/cross_specification/` remain trackable.
The restart ledger binds each launched seed to the simulation contract and
advances past a timed-out internal stdpopsim trajectory; it does not mislabel
that interrupted trajectory as a completed external sample-panel draw.
