# Replicated EAS introgressed-sweep study

This directory contains the completed exploratory extension of the
AncientEurasia Han introgressed-sweep pilot in the parent directory. The
biological model is unchanged: one archaic-specific focal allele enters the
recipient lineage only through the Neanderthal pulse, selection begins at the
pulse, and present-day Han panels are conditioned on nested allele-frequency
floors from 10% through 90%.

The campaign contains 27 biological cells (`s` = 0.01, 0.005, or 0.001 by nine
frequency floors) and requests ten independently seeded trajectory slots per
cell, for 270 requested slots. Ten is a search target, not a guarantee of ten
accepted trajectories. The completed campaign accepted and decoded 155 slots
across 24/27 represented cells; 22 cells had the at-least-two paired
trajectories required for the macro headline. The remaining terminal statuses
were 104 `exhausted`, 11 `timed_out`, and 0 `failed`. `specifications.tsv` and `study_design.json` are
immutable once created. Bulk tree sequences and decoder intermediates live in
the ignored `work/` directory; compact reports and letter-page figures in both
PNG and PDF are in `results/`.

See [`METHODS.md`](METHODS.md) for the replicate and stopping rules, the parent
[`../METHODS.md`](../METHODS.md) for the full demographic/selection model, and
[`RUN_RESULTS.md`](RUN_RESULTS.md) for the completed local run.

## Reproduce or resume

Bootstrap the repository first:

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

Create or validate the immutable plan:

```bash
uv run --frozen --all-extras python \
  scripts/run_eas_introgression_replicates.py plan \
  --repo-root . --campaign replicate_study_n10 --replicates 10
```

On Windows, run or resume one bounded simulation pass with the separately
installed stdpopsim-compatible SLiM 4.2.2 executable:

```powershell
uv run --frozen --all-extras python `
  scripts/run_eas_introgression_replicates.py simulate `
  --repo-root . --campaign replicate_study_n10 --replicates 10 `
  --slim-bin .native-stdpopsim\Library\bin\slim.exe `
  --workers 12 --threads 4 --max-external-draws 20 `
  --max-launched-draws 40 --spec-timeout-minutes 5 `
  --cumulative-spec-timeout-minutes 30
```

Decode under Linux/WSL because the bundled AVX2 Gamma-SMC binary is a Linux
executable. After a read-only throughput audit, the completed campaign used 24
concurrent workers with one decoder thread per worker:

```bash
.venv-wsl/bin/python scripts/run_eas_introgression_replicates.py decode \
  --repo-root . --campaign replicate_study_n10 --replicates 10 \
  --decoder-bin bin/gamma_smc --workers 24 --threads 1 \
  --max-external-draws 20 --max-launched-draws 40 \
  --spec-timeout-minutes 5 --cumulative-spec-timeout-minutes 30
```

Worker and thread counts are execution settings recorded as provenance; the
validated Gamma-SMC output contract treats thread count as provenance-only
under its output-invariance contract rather than as scientific cache identity.

Aggregate all 270 requested statuses together with every available validated
truth trajectory and every available accepted-and-decoded Gamma-SMC trajectory:

```bash
uv run --frozen --all-extras python \
  scripts/run_eas_introgression_replicates.py plot \
  --repo-root . --campaign replicate_study_n10 --replicates 10
```

Each invocation appends immutable started/completed/failed records under
`phase_invocations/`. For an incomplete slot, the effective timeout is the
smaller of five minutes and the time remaining in its 30-minute cumulative
budget. Re-running a phase validates checksum-current complete artifacts and is
a cache/no-op for those slots; it does not consume a new seed or repeat the
simulation/decode. Plotting deliberately rejects `--spec`, preventing a
filtered run from overwriting campaign-wide results.

SLiM 5.2 is not used here: stdpopsim 0.3 emits a model using an API removed in
SLiM 5.2. The runner version-checks and requires SLiM 4.2.2.

Some early invocation records predate per-task source-hash recording. Current
validation establishes artifact integrity and parseability, but it cannot by
itself prove semantic equivalence across those source epochs; removing that
evidence limitation requires a frozen-source rerun.
