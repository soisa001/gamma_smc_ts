# Introgressed EAS sweep study

This directory runs the AncientEurasia Han introgression design described in
[`METHODS.md`](METHODS.md). Raw tree sequences and decoder intermediates are
restartable under `work/`; compact summaries and cross-specification figures
are written under `results/`. The checked execution outcome and its limitations
are summarized in [`RUN_RESULTS.md`](RUN_RESULTS.md).

The completed exploratory replicate extension is documented in
[`replicate_study_n10/README.md`](replicate_study_n10/README.md), with methods
in [`replicate_study_n10/METHODS.md`](replicate_study_n10/METHODS.md) and final
results in
[`replicate_study_n10/RUN_RESULTS.md`](replicate_study_n10/RUN_RESULTS.md). It
requested ten independently seeded slots per biological cell rather than
guaranteeing ten accepted replicates: 270 slots were requested, 155 were
accepted and decoded across 24/27 represented cells (22 with at least two
paired trajectories for the macro headline), 104 exhausted their search budget,
11 ended timed out, and 0 failed.

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
SLiM scripts. Direct SLiM 5.2 execution fails because stdpopsim 0.3 emits the
removed `Subpopulation.genomes` API. A compatibility-transformed 10-Mb test was
slower than the accepted 4.2.2 comparison and did not finish within 200
seconds, so 4.2.2 remains the validated study default.

## Run the phases

Run phases explicitly from the repository root. Defaults are four concurrent
specifications; every active Gamma-SMC decoder uses four threads (up to 16
decoder threads when four specifications are ready).

```bash
bash introgression_EAS_sim/run.sh plan
bash introgression_EAS_sim/run.sh simulate
bash introgression_EAS_sim/run.sh decode
bash introgression_EAS_sim/run.sh plot
```

`decode` requires Linux x86_64 with AVX2. `all` runs the four phases in order,
but separate invocations are easier to restart and audit. To filter a pilot or
resume a subset, repeat `--spec` with a substring from `specifications.tsv`:

```bash
bash introgression_EAS_sim/run.sh simulate --spec s0p010_af50
```

For a bounded exploratory pass over rare conditioned cells, declare the limit
explicitly; timed-out cells are killed as process groups and written to the
simulation status table rather than treated as successful output:

```bash
bash introgression_EAS_sim/run.sh simulate --spec-timeout-minutes 30
```

Override parallelism for a run:

```bash
EAS_SWEEP_WORKERS=8 EAS_SWEEP_THREADS=8 \
  bash introgression_EAS_sim/run.sh simulate
```

`work/` and per-specification bulk figure directories are ignored by Git.
`specifications.tsv`, study design/provenance records, compact result tables,
result manifests, and `results/figures/cross_specification/` remain trackable.
Repeated bounded invocations advance deterministic seeds for an interrupted
internal stdpopsim trajectory while preserving every completed sample-panel
draw; a timeout is never counted as accepted data.
