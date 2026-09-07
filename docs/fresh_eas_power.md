# Fresh EAS decoded-power study

This study starts from new RNG draws and reuses no previous simulation outputs.
The PHLASH EAS demographic resource is a model input, not a simulated dataset.
The canonical specification is `configs/eas_q02_50k.json`.

- EAS, 2% introgression at 50 kya; selection begins at that pulse.
- Ten coefficients, 0.001 through 0.010 in steps of 0.001; 100 accepted selected replicates per coefficient, plus 1,000 neutral replicates shared across coefficients.
- 100 diploid samples, 10,000 unique sampled haplotype pairs, 10 kb output stride.
- Constant mutation rate 1.25e-8 and recombination rate 1e-8 per bp per generation; 25 years per generation.
- The matched archaic branch retains the prior split at 700 kya, Ne=3,600, and Ne=10 bottleneck for 100 generations immediately before the pulse. Both selected and neutral simulations use this branch. No conditioning on realised ancestry or a pulse AF band is imposed.
- Selected alleles are fixed in the archaic source, required to survive in the recipient population, and required to appear in the sampled panel. Sample AF=1 is retained. Failed sample ascertainment is logged and retried with a deterministic attempt seed.
- Both engines simulate 11 Mb. Selected crops center the selected allele; neutral crops center the nearest segregating site, with no MAF floor. The decoder processes the resulting 10 Mb, with exact one-based VCF conversion and a common pair/position manifest.
- SLiM 4.2.2 with scaling factor 5 and burn-in 0.1 supplies selected simulations. msprime 1.4.2 supplies neutral simulations. stdpopsim is pinned to 0.3.0.

The primary statistic is `frac_recent_T`, the fraction of pairs whose posterior mean TMRCA is below T. Cutoffs are 5, 10, 20, 30, 40, and 50 kya; each is reported separately. Power at alpha=0.05 is the fraction of selected focal-stride statistics with `(1 + number of neutral statistics >= selected statistic)/(1000 + 1) < 0.05`. Neutral reference values come from one aligned focal stride per independent replicate. There is no consecutive-window requirement and no claim of genome-wide error control.

## Phases and local invocation

Run inside WSL. The checkout is `/mnt/d/phase2simselection/code`; outputs are `/mnt/d/phase2simselection/sim/eas_q02`. The launcher uses 20 single-threaded workers, as requested. It uses the local pinned uv environment and resolves native library/SLiM locations from `GAMMA_NATIVE_DIR` and `SLIM_BIN`, with defaults for this machine. `UV_BIN` defaults to `$HOME/.local/bin/uv` so a background WSL launch does not depend on login-shell PATH setup.

The permitted RAM budget is 300 GB, subject to memory exposed by WSL/the operating system. It is an allowance, not a reservation. The current workflow bounds concurrent work by worker count and does not preallocate the memory budget. Resource changes requested during a run are recorded in `resource_allowance.json`; they do not change simulation seeds or invalidate completed outputs.

```bash
bash scripts/launch_fresh_eas.sh smoke
bash scripts/launch_fresh_eas.sh run
```

The smoke phase checks three fresh neutral replicates and one fresh selected replicate each at s=0.001, 0.005, and 0.010. These become the first replicates of the production study. The run phase validates and reuses only complete replicates from this fresh study and performs final analysis after all 2,000 are complete. A file lock prevents overlapping runners.

An interrupted replicate resumes from its hashed `simulated.json` receipt when the simulation phase completed. Completed replicates require all saved hashes and decoded profile checks to pass; missing or corrupt completed files stop the run for investigation. Adding coefficients or replicate counts preserves existing seeds. Cache eligibility uses the scientific configuration, executable identity, software versions, exact pair manifest, and saved outputs, not the raw Git commit.

All logs and temporary files are on D. `run_status.json` records progress; `sample_manifest.json` records planned replicate identities and seeds. A storage check stops new submissions before reaching the 3 TB aggregate limit, retaining 20 GB of free/reserved space. Restarting WSL still stops a running process; rerun the same command to resume.

Final outputs are `analysis/focal_statistics.csv`, `selected_pvalues.csv`, `power.csv`, `power.png`, and `power.pdf`. Power intervals are 95% Wilson intervals conditional on the estimated null; they do not include null Monte Carlo uncertainty. Full simulation and cropped decoder-input trees are retained, along with compact true-statistic profiles for checks.

## Reproducible environment

Install with uv in Linux/WSL, using Python 3.11.15. `requirements/fresh_eas_runtime.txt` pins the tested environment. The repository's older uv version bound is bypassed with `--no-config` for this standalone workflow.

```bash
uv --no-config venv --python 3.11.15 .venv
uv --no-config pip install --python .venv/bin/python -r requirements/fresh_eas_runtime.txt
uv --no-config pip install --python .venv/bin/python --no-deps -e .
```

Build `bin/gamma_smc` with the repository Makefile and available Boost/HTSlib/zstd development libraries. Point `GAMMA_NATIVE_DIR` to their prefix and `SLIM_BIN` to SLiM 4.2.2. Native binary hashes are captured in every study manifest.
