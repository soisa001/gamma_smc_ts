# Bounded EAS time-resolved pilot

This implementation runs the no-introgression feasibility experiment from
`time_resolved_eas_selection_plan.md`. It is a smaller benchmark designed to
complete within a 30-minute interactive session; it does not replace the full
planned power study.

The completed 2026-11-02 results are in
[`results/time_resolved_eas_20261102/README.md`](results/time_resolved_eas_20261102/README.md).
All 219 main runs and 20 additional strong-selection controls completed;
simulation and scoring took 21.52 minutes including setup.

## Scientific contract

- Use the checked EAS PHLASH median history and 25 years per generation.
- Initialize the complete diploid population at 30 kya using neutral ancestry.
  Use DTWF for the first 100 biological generations of prehistory, then Hudson
  ancestry deeper in the same EAS history.
- At the fixed central coordinate, choose an existing ancestral branch whose
  descendants are nearest 5% of the onset population, restricted to 2-10%.
  Place the focal standing mutation on that branch. This is explicit
  conditioning on a standing allele's onset frequency, applied identically to
  neutral and selected replicates; it is not an unconditional mutation sample.
- Run an unscaled SLiM 5.2 Wright-Fisher simulation from 30 kya to the present.
  In the selected arm, homozygote advantage is s=0.01 and dominance is h=0.5
  (heterozygote advantage 0.005), active for the 400 reproductive generations
  spanning 30-20 kya, then zero. The neutral arm has s=0 throughout.
- Do not reject loss, low final frequency, or fixation. Record every trajectory.
  A missing eligible initial branch is an explicit failure, not an unlogged retry.
- Sample 100 diploids uniformly at present, then use the same seeded panel of
  1,000 unordered haplotype pairs without genotype filtering in every replicate.
- Use 1-Mb regions, recombination rate 1e-8 and neutral mutation rate 1.29e-8 per
  bp per generation. Remove the deliberately placed focal marker before adding
  independent neutral sequence mutations; decoder evidence therefore comes from
  linked neutral sequence, without using the selected genotype as a label.
- Score positions 0,10,000,...,990,000 bp and all 50 cutoffs 1,2,...,50 kya.
  The center (500 kb) is on this grid. Reconcile native VCF coordinates explicitly.
- Decode once with Gamma-SMC, retain its Gamma posterior parameters, and reduce
  to both exact posterior mass and median indicators. Evaluate true-tree
  indicators at the identical positions and pairs. A median call is implemented
  as Gamma CDF >0.5 with a strict below-cutoff convention.

The runner exposes a scaling option for future diagnostics, but the reported
initial campaign uses Q=1. Results from different Q values must not be pooled.

## Independent cohorts and performance

The initial schedule has 40 development nulls, 99 calibration nulls, 40 held-out
neutral tests, and 40 selected tests. Seeds derive deterministically from base
20380101, a study/stage namespace, and the replicate ID. Pair seed is 1729.

Fit interval-specific nuisance means and standard deviations on development
nulls, pooling their genomic positions. This pooling defines a scoring rule;
it does not treat positions as independent tests. Each calibration null then
contributes one maximum over the entire genomic-by-time grid. Conservative
empirical p=(1+#calibration maxima >= score)/(B+1) includes ties. Zero-variance
training intervals are excluded and counted explicitly.

Primary windows are 49 one-kyr intervals. The separately calibrated secondary
scan tests 182 intervals with widths 1,2,5,10 kyr, all stepping by one kyr.
These are separate method comparisons, not an uncorrected combined discovery
rule. The report retains each scoring method's separate null calibration.

Neutral FPR means any call anywhere in a held-out neutral region. Joint TPR
requires a significant cell within 100 kb of the known selected coordinate and
overlapping [20,30) kya. Temporal overlap and temporal F1 are also measured at
the exact selected coordinate. Connected components use orthogonal grid
adjacency. Event F1 allows one true event per selected region; unmatched
components count as false events, including all components in neutral tests.
Consequently, event F1 depends on the selected-to-neutral test mixture.

Rate uncertainty uses Wilson 95% intervals. Paired mass-minus-median differences
use seeded bootstrap resampling of whole test regions. No-call temporal F1 is
zero. Peak errors are unconditional descriptions, not confidence bounds on a
selection date. Neither this pilot nor an all-pairs marginal posterior average
establishes a full joint posterior over selection histories.

## Execution and restart

Use a Linux environment with the repository's uv-managed dependencies, SLiM 5.2,
and built Gamma-SMC. The local run uses eight worker processes, one decoder
thread each, and a 1,200-second campaign deadline plus 180-second per-unit
timeouts. The controller terminates a timed-out worker's complete process group.
Its sequentially issued task manifest records every completed, failed, timed-out,
or unstarted unit. No incomplete subset should be described as a full cohort.

The completed run also launched a separate four-worker, 20-replicate s=0.05
positive-control batch, for at most 12 simultaneous workers. This exploratory
batch was added after observing weak preliminary recovery at s=0.01. It reused
the first 20 selected-ancestry seeds to provide a paired strength diagnostic.
Its nulls were taken from the main run only after checking that all parameters
other than the unused selected-arm coefficient matched. The two strengths are
reported separately and must not be pooled as independent replicates.

```bash
uv run --frozen --all-extras python scripts/run_time_resolved_pilot.py run \
  --out /path/to/time_scan_pilot \
  --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --scaling 1 --length 1000000 --diploids 100 --pairs 1000 \
  --selection 0.01 --workers 8 --seconds 1200 --unit-timeout 180

uv run --frozen --all-extras python scripts/report_time_resolved_pilot.py \
  --root /path/to/time_scan_pilot --out /path/to/report/one_kyr
uv run --frozen --all-extras python scripts/report_time_resolved_pilot.py \
  --root /path/to/time_scan_pilot --out /path/to/report/multiscale --widths 1 2 5 10
```

Re-running with identical scientific settings and source contracts reuses
completed units; Git commit changes alone do not invalidate them. Larger cohort
targets may be scheduled under the same manifest. Raw trees/posteriors stay in
the simulation directory; aggregate simulation results and scripts can be pushed.
The local pilot uses an existing uv-managed Linux environment, recording its
actual package versions in the manifest; it does not claim a fresh lockfile sync.

Every unit checks carrier identity at initialization, population-size trajectory,
pair/position ordering, finite valid posteriors, monotone bounded CDFs, and complete
posterior streams. Reports verify score and tree hashes and cohort seed uniqueness.
PNG and vector PDF figures are written to disk only, with headless plotting and
editable PDF fonts. Relevant focused tests are in `tests/test_time_resolved_pilot.py`.

## Ongoing selection extension (2026-11-03)

**Superseded initial design:** the unconditional run described immediately below
was interrupted when the user required present-day segregation. Its outcomes
are proposal diagnostics only, not the requested conditional power benchmark.
The corrected design and runnable phases follow at the end of this section.

`scripts/run_ongoing_selection_pilot.py` runs a separate experiment with selection
enabled from 50 kya to the present, s in 0.01, 0.02, 0.03, 0.04, 0.05, h=0.5,
and no scaling. The de novo allele starts as one copy. The introgressed allele
is fixed among Neanderthal donors and absent among EAS founders. At 50 kya,
the nearest integer to 1% of the 2N founder haplotypes comes from Neanderthals;
founders are randomly paired into diploids. Their actual genealogies are
simulated using a 500 kya human–Neanderthal divergence and Neanderthal Ne=3,600.
These demographic assumptions follow the existing repository configuration.
The selected marker tags donor membership at onset; it is not an additional
simulation of the allele's earlier fixation within Neanderthals.

The complete population at onset is simulated with DTWF for the first 100
generations backward, then Hudson ancestry. EAS size changes use the same checked
PHLASH median history as the earlier pilot. Population-size schedules, founder
identities, realized pulse proportion, loss/fixation times, and all 2,001
trajectory time points are retained. The focal marker is removed before adding
independent neutral mutations for decoding. No failed allele establishment is
retried or discarded, and there is no present-day allele-frequency filter.

The local campaign requests 20 development, 99 calibration, 20 held-out neutral,
and 20 selected replicates per coefficient **for each origin**: 478 regions.
It uses 20 workers, one decoder thread each, and a 1,200-second campaign deadline.
CLI defaults remain four workers. The two initial smoke regions are excluded.
Four additional workers process the final 40 scheduled calibration units from
the same manifest while the main controller handles earlier units, for at most
24 simultaneous workers. Their schedule and status are recorded separately;
the main controller verifies and reuses their receipts when it reaches them.
The pipeline samples 1,000 unfiltered unordered haplotype pairs from 100 diploids
per region; it does not enumerate all 19,900 possible pairs.

The report uses independent, origin-matched neutral calibration. The primary
1-kyr interval scan, secondary multiscale scan, and cumulative 50-kya baseline
are evaluated separately, each corrected over its entire genomic/time search.
Spatial TPR means a significant cell within 100 kb of the selected site; neutral
FPR means any significant cell in a held-out neutral region. Conditional TPR
among survivors and among alleles that ever reached 10% is also reported, with
the actual denominator. These are descriptive strata, not replacement cohorts.

Because selection covers the entire 1–50 kya scan grid, temporal overlap does
not validate dating and temporal false positives outside selection cannot be
measured. Temporal F1 here is coverage of [1,50) kya. The omitted 0–1 kya interval
is explicitly outside the grid. Fixation also removes segregating variation
even though the selection coefficient stays enabled. Comparisons with the
earlier 30–20 kya experiment change initial allele frequency and origin as well
as duration, so they do not isolate duration's causal effect.

Run in phases from a Linux environment with SLiM 5.2 and the native Gamma-SMC
binary installed. The repository bootstrap can provision the native stack.
No controlled data are accessed by this simulation-only command.

```bash
%%bash
set -euo pipefail
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
OUT="$PWD/results/ongoing_eas_20261103"
export MPLBACKEND=Agg
export LD_LIBRARY_PATH="$PWD/.native/lib:${LD_LIBRARY_PATH:-}"
uv run --frozen --all-extras python scripts/run_ongoing_selection_pilot.py run \
  --out "$OUT" --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --development 20 --calibration 99 --test 20 --selected 20 \
  --selection 0.01 0.02 0.03 0.04 0.05 --workers 4 --seconds 1200
uv run --frozen --all-extras python scripts/report_ongoing_selection_pilot.py \
  --root "$OUT" --out "$OUT/report"
```

Use `--workers 20` to reproduce the local concurrency; four workers may leave
units incomplete at the same deadline. Re-running the simulation phase with
identical scientific settings checks and reuses completed files. The report
lists every missing scheduled unit and verifies tree, posterior, and score
hashes before computing performance. Figures are saved only to PNG and vector
PDF with editable text. Focused tests are in
`tests/test_ongoing_selection_pilot.py` and `tests/test_time_resolved_pilot.py`.

### Final design: segregating selected allele and relabeled neutral focal site

Selected focal alleles must be segregating at present: 0 < whole-population
frequency < 1. Fixation does not qualify as survival. Following the user's final
clarification, neutral genealogies are unconditioned, and their focal coordinate
is relabeled to the nearest variant segregating in the present-day sample.
No particular neutral mutation introduced at 50 kya must survive. The 10 kb
grid and whole-region calibration/FPR remain unchanged. The revised starting
coefficient grid is 0.005, 0.01, 0.02,
0.03. Selection from 50 kya is the current pilot; stronger coefficients with
recent onset around 5 kya are reserved for a separate comparison. In a recent
introgression-selection experiment, the donor pulse must remain at 50 kya and
selection begin later. No 5-kya result has been simulated or claimed here.

`run_conditioned_ongoing_pilot.py` performs exact rejection sampling with early
termination at absorbing loss/fixation. This saves computation without changing
the accepted distribution. Selected simulations use the strict segregation gate
and keep their preassigned seeds. Neutral rejection is an optional experimental
mode (`--condition-neutral`), disabled by default and excluded from final power.
Old proposals can be reused from a source campaign, with source-receipt hashes.
The source campaign can also provide the unconditioned neutral reference.
The 0.005 coefficient has its own `s0.005` directory and cannot collide with 0.01.

`report_conditioned_ongoing_pilot.py` reports proposal acceptance and checks
whether valid power estimation is possible. It uses only each stratum's fully
completed deterministic proposal prefix, avoiding bias toward quickly rejected
simulations at a deadline. The minimum independent null sizes are 5 development,
19 calibration (the minimum for empirical p<=0.05), and 1 held-out neutral test.
Those minima permit a calculation but do not imply precise power or FPR estimates.
Insufficient selected acceptances produce **not estimable**, never a zero-power
result. Larger cohorts remain necessary for a useful study. The selected allele
need not be observed as segregating in the random sample; its condition is on
the whole simulated population. Neutral focal labels are polymorphic within
the sample, which also establishes population segregation.

```bash
%%bash
set -euo pipefail
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
export MPLBACKEND=Agg
export LD_LIBRARY_PATH="$PWD/.native/lib:${LD_LIBRARY_PATH:-}"
SOURCE="$PWD/results/ongoing_neutral_reference"
OUT="$PWD/results/ongoing_segregating"
# Phase 1: independent, unconditioned neutral regions for each origin.
uv run --frozen --all-extras python scripts/run_ongoing_selection_pilot.py run \
  --out "$SOURCE" --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --development 10 --calibration 20 --test 10 --selected 0 --workers 4 --seconds 1200
# Phase 2: bounded rejection sampling; defaults use four workers.
uv run --frozen --all-extras python scripts/run_conditioned_ongoing_pilot.py \
  --source "$SOURCE" --out "$OUT" --workers 4 --seconds 1200 \
  --attempts 100 --calibration-attempts 200 --selection 0.005 0.01 0.02 0.03
# Phase 3: integrity checks, acceptance counts, and performance when estimable.
uv run --frozen --all-extras python scripts/report_conditioned_ongoing_pilot.py \
  --root "$OUT" --neutral-reference "$SOURCE" --out "$OUT/report"
```

The local corrected run uses 24 persistent workers and a six-minute final batch
after the earlier design corrections. Changing the attempt budget on restart
preserves existing scientific outputs and seeds. All figures stay on disk.
The executed runner snapshots preserve the exploratory neutral-conditioning
phase; the current runner defaults to the user's final selected-only condition.
The final report is `docs/results/ongoing_eas_segregating_20261103/README.md`.

### Recent onset: 10 kya to present

The next experiment starts selection 10,000 years ago (400 generations), at
s = 0.005, 0.01, 0.02, and 0.03, with h = 0.5. A de novo mutation enters as
one copy at that date. For introgression, the 1% Neanderthal pulse stays at
50 kya and the donor allele evolves neutrally until 10 kya. The split remains
500 kya and donor Ne remains 3,600. Present-day segregation is required only
for selected alleles; neutral regions retain the nearest segregating variant
as their focal label. The marker used for conditioning is removed before the
neutral mutation overlay and GammaSMC inference.

The runner rejects absorbing loss/fixation early, checks exact allele counts,
and verifies that selection acts for exactly 400 reproductive transitions.
The time scan retains 1–50 kya at 1-kyr resolution and a 10-kb genomic stride.
The pilot samples 1,000 unfiltered haplotype pairs from 100 diploids, using
fixed pair seed 1729, to approximate the all-pairs statistic within the budget.
It does not evaluate all 19,900 possible pairs. Each region is 1 Mb.

The report separates spatial detection within 100 kb of the simulated site
from joint detection requiring temporal overlap with [1,10) kya. The interval
[0,1) is outside the requested grid. Center temporal F1 measures overlap of
called one-kyr cells with those nine observable selected cells. A call wholly
outside the selected interval is recorded as mistiming; any call in an
independent neutral test region contributes to neutral FPR. Max-statistic
calibration covers the entire region and all scanned intervals, separately
for each estimator and scan variant. This is region-level error control,
not genome-wide family-wise error control.

Targets are 20 development nulls, 99 calibration nulls, 50 held-out nulls per
origin, and the first 50 accepted selected proposals per coefficient and
origin. A wall-clock cap can leave smaller cohorts. Reporting uses complete
proposal prefixes, preserves seed order, and records all losses, fixations,
incomplete units, and failures. Smaller samples are reported with counts and
Wilson intervals rather than treated as definitive power estimates.

```bash
%%bash
set -euo pipefail
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
export MPLBACKEND=Agg
export LD_LIBRARY_PATH="$PWD/.native/lib:${LD_LIBRARY_PATH:-}"
OUT="$PWD/results/recent10k_eas"
# One-time setup on a new Linux machine (does not run the full test suite).
if [ ! -x .native/bin/slim ] || [ ! -x bin/gamma_smc ]; then
  bash scripts/bootstrap_uv.sh --skip-tests
fi
export PATH="$PWD/.tools/uv-bin:$PATH"
# Phase 1: simulate, infer pairwise posteriors, and score; bounded to 55 minutes.
uv run --frozen --all-extras python scripts/run_recent_selection_pilot.py \
  --out "$OUT" --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --onset-years 10000 --selection 0.005 0.01 0.02 0.03 \
  --development 20 --calibration 99 --test 50 --accepted 50 \
  --workers 4 --seconds 3300
# Finish only already-submitted units, with a further 2.5-minute hard cap.
uv run --frozen --all-extras python scripts/finish_recent_selection_pilot.py \
  --root "$OUT" --workers 4 --seconds 150
# Phase 2: validate artifact hashes and trajectories, then report performance.
uv run --frozen --all-extras python scripts/report_recent_selection_pilot.py \
  --root "$OUT" --out "$OUT/report"
# Phase 3: exact all-pairs genealogy benchmark, using the same saved simulations.
uv run --frozen --all-extras python scripts/audit_recent_all_pairs_truth.py \
  --root "$OUT" --report "$OUT/report" --workers 4
```

The setup phase is outside the simulation time cap. The local experiment used 24 workers; four remains
the CLI default. Re-running with the same scientific manifest reuses verified
completed units. Figures are saved to PNG and editable vector PDF on disk.
The all-pairs audit compares all 19,900 true TMRCA pairs against the 1,000-pair
benchmark, with independent calibration for the full-pair statistic. This
diagnoses subsampling effects; it does not infer 19,900 GammaSMC posteriors.

## Three selection episodes and expanded coefficient grid

Completed local results: [three-episode EAS report](results/selection_arrays_eas_20261104/README.md).

`run_selection_arrays.py` compares 10 kya to present, 50 kya to present, and
40–20 kya, with s = 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05. The target is
100 present-day segregating selected replicates per coefficient, origin, and
episode: 42 treatments. De novo alleles enter at selection onset. The 1%
Neanderthal pulse always occurs at 50 kya. In the 40–20 kya episode, selection
acts for exactly 800 reproductive transitions and then stops; drift continues
for another 800 generations. Strong sweeps that fix are excluded by the
requested segregation condition. A stratum with no surviving replicates has
unestimable conditional power, not zero power.

The scheduler shares one wall-clock deadline across treatments, attempts them
round-robin, and stops dispatching 150 seconds before the deadline to finish
submitted work. Any remaining worker process groups, including SLiM and
GammaSMC children, are killed at the deadline. Targets can remain incomplete.
Analysis uses uninterrupted proposal prefixes to avoid selecting only the
fastest surviving simulations. The scientific manifest, source snapshots,
exact trajectories, seeds, artifact hashes, losses, fixations, and incomplete
tasks are retained. Re-running the same output directory validates and reuses
completed work. Scheduler-only source changes are accepted only if scientific
parameters and the parsed simulation/timing functions are unchanged.

Optional `--reuse-recent /path/to/recent10k_eas_fast_20261103` references verified
prior raw simulations without copying them. Those source files must remain
available. Omit this option for a fresh experiment. Introgression nulls are
shared across episodes because their selection coefficient is always zero;
de novo null pools match their episode-specific forward duration.

The primary detection grid covers all 49 adjacent intervals [1,2), …, [49,50)
kya and 100 positions separated by 10 kb, irrespective of simulated selection
timing. Truth is used only for evaluation. Independent null maxima calibrate
the whole-region scan separately for posterior median, posterior mass, and
true trees. The posterior pilot uses the same 1,000 unfiltered pairs as earlier
experiments. An independent all-19,900-pair true-tree audit diagnoses pair
subsampling effects. Neither calibration is genome-wide family-wise control.
As a secondary analysis, widths of 1, 2, 5, and 10 kyr are scanned with 1-kyr
starts and jointly calibrated across all widths and positions. This tests
whether narrow temporal bins lose power. Primary and secondary scans each
have their own calibration; their union is not a calibrated discovery rule.

```bash
%%bash
set -euo pipefail
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
export MPLBACKEND=Agg
export LD_LIBRARY_PATH="$PWD/.native/lib:${LD_LIBRARY_PATH:-}"
# Phase 0: one-time native/runtime setup, outside the simulation budget.
if [ ! -x .native/bin/slim ] || [ ! -x bin/gamma_smc ]; then
  bash scripts/bootstrap_uv.sh --skip-tests
fi
export PATH="$PWD/.tools/uv-bin:$PATH"
OUT="$PWD/../sim/selection_arrays"
# Phase 1: one shared 55-minute cap, including inference and drain time.
uv run --frozen --all-extras python scripts/run_selection_arrays.py \
  --out "$OUT" --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --selection 0.001 0.005 0.01 0.02 0.03 0.04 0.05 --accepted 100 \
  --development 20 --calibration 99 --test 50 \
  --workers 4 --seconds 3300 --drain-seconds 150
# Optional phase 1b: finish only submitted units, with a further 150-second cap.
uv run --frozen --all-extras python scripts/finish_selection_arrays.py \
  --root "$OUT" --workers 4 --seconds 150
# Phase 2: validate saved simulations and report posterior performance.
uv run --frozen --all-extras python scripts/report_selection_arrays.py \
  --root "$OUT" --out "$OUT/report"
# Phase 3: exact all-pair truth benchmark; no new simulations.
uv run --frozen --all-extras python scripts/audit_selection_arrays.py \
  --root "$OUT" --out "$OUT/report" --workers 4
```

The local November 4 experiment used 24 workers and reused the earlier 10-kya
cohort. Four workers remains the portable default. Each invocation has its
own time cap; resuming after a cap explicitly spends additional compute.
Figures are written as PNG and vector PDF with editable text, without inline
display. The temporal heatmaps use vector cells in PDF, not raster images.
