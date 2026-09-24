# Introgressed standing-allele experiment

Authorized 2026-09-24. This supersedes the mixed-origin pilot. The user specified
constant Neanderthal Ne=3,600, a 500-kya split, a 2% pulse into EAS at 50 kya,
an existing archaic-derived allele segregating in EAS when selection starts,
100 selected simulations per specification, and no de novo simulations for now.

Configuration: `configs/segregating_introgression.json`. The generation time is
25 years; population rescaling is disabled (Q=1). The EAS PHLASH history, mutation
rate 1.25e-8, recombination rate 1e-8, h=0.5, 400 sampled haplotypes, 10,000
sampled pairs, seed base 20380101, and pair seed 1729 are retained.
The expanded 28-cell grid is in `segregating_introgression_grid.csv`.

The donor's constant Ne is inspired by the
[Ragsdale–Gravel stdpopsim implementation](https://raw.githubusercontent.com/popsim-consortium/stdpopsim/0.3.0/stdpopsim/catalog/HomSap/demographic_models.py).
Our rounded split date, EAS history, and discrete introgression pulse make this
a hybrid specification, not an exact reproduction of that published model.
There is no added donor bottleneck. The donor Ne is an effective diploid size.

## Cohorts and focal-allele ascertainment

| Cohort | Focal choice / selection onset | Number |
|---|---|---:|
| Selected I50 | Immediately after the 50-kya pulse | 100 per s |
| Selected I10 | At 10 kya, after 40 kyr of neutral drift | 100 per s |
| Calibration null I50 | Same 50-kya choice, s=0 | 1,000 |
| Calibration null I10 | Same 10-kya choice, s=0 | 1,000 |
| Independent neutral targets | Same choice rules, s=0 | 100 per onset |

The 12 positive coefficients are 0.001, 0.002, 0.003, 0.004, 0.005, 0.006,
0.007, 0.008, 0.009, 0.010, 0.020, and 0.030. Total: 2,400 selected targets,
2,000 calibration nulls, and 200 independent neutral targets = 4,600 retained
regions. Nulls are shared across coefficients, cutoffs, pair classes, and TMRCA
sources within an onset. Onsets need separate calibration distributions because
the focal is chosen at different times. No de novo tasks are generated.

At the specified onset, choose the closest eligible variant to 5.5 Mb of the
11-Mb simulation. Eligibility requires exactly one mutation at a biallelic site,
an origin in the Neanderthal branch between the split and pulse, and frequency
strictly between 0 and 1 in the full EAS census at onset. The allele may have
been fixed or segregating in the donor. Trace migration along the mutation's
branch rather than assuming its lower node identifies its population of origin.

Eligible coordinates must lie in [5,000,000, 6,000,000] so the full 10-Mb window
fits. Resolve equal distances using the lower coordinate. Shift the retained
window to place the chosen focal at exactly 5 Mb. Choose the I10 focal at 10 kya;
do not follow a previously chosen I50 focal or reset its frequency to 2%.

An attempt with no eligible variant is logged and regenerated. After choosing
the allele, retain the existing ascertainment requirement of present-day sample
observation. Loss or absence from the sampled panel is logged and retried using
a new deterministic attempt seed and a newly generated prehistory. Fixation
after onset is retained. Never choose a different allele retrospectively after
observing its survival or selection response. Thus detection is conditional on
the observed focal allele, not unconditional power per mutation introduced.
Rejected-attempt receipts preserve progress across restarts. The operational
ceiling is 10,000 attempts per retained region; reaching it stops the queue for
inspection rather than silently reducing a cohort or changing ascertainment.

## Simulation phases and biological checks

1. Simulate the full EAS census's neutral genealogy at onset with msprime's
   discrete-time Wright–Fisher model through the archaic split. Use Hudson
   ancestry only deeper in the shared ancestral human branch. No population or
   time rescaling is applied. At the immediate pulse, construct the post-pulse
   cohort explicitly by assigning each of its 2N haplotypes a donor source with
   probability 0.02. This avoids msprime's same-time event/sample activation
   ordering, without moving the 50-kya date. At 10 kya, simulate backward through
   the earlier pulse normally.
2. Realize neutral mutations before onset and choose the focal using only that
   state. Save `ascertainment.trees`, including absolute ages and migration
   records, so eligibility and the nearest-site decision can be independently
   reproduced. Randomly pair the sampled chromosomes into EAS diploids.
3. Load that state into SLiM 5.2 using pyslim 1.1.1. Preserve the existing
   mutation IDs and genotype states. Assign fitness to the chosen mutation and
   run exactly 2,000 (I50) or 400 (I10) subsequent WF generations under the
   recorded EAS size schedule. Null and selected runs use the same code.
4. Sample 200 diploids without replacement. Overlay neutral mutations only in
   the forward interval; never redraw the pre-onset variation used for allele
   ascertainment. Protect the focal coordinate from recurrent overlay. Convert
   mutation states to nucleotides, crop, and verify focal carrier identity in
   both the original and shifted tree sequences.
5. Hash/read audit, Gamma-SMC decoding, and analysis follow as separate phases.

The hybrid initialization follows the general
[pyslim initialization workflow](https://tskit.dev/pyslim/docs/latest/tutorial.html),
with the explicit DTWF ancestry and onset ascertainment above. Runtime versions,
binary hashes, configuration, source hashes, seeds, attempt outcomes, chosen
mutation age/origin/position, crop offset, and allele-frequency trajectories are
recorded. Scientific cache compatibility is based on inputs and runtime semantics,
not the Git commit. Earlier predefined-focal and de novo artifacts are incompatible.

## Analysis

Evaluate the prespecified focal coordinate at cutoffs 10, 20, 30, 40, and 50 kya.
Use separate truth and Gamma-SMC results, all-pairs and raw within-ALT/ALT
fractions, empirical upper-tail p=(1+#null >= target)/1001, and p<=0.05.
The decoded fraction uses pairs' posterior mean TMRCA classifications, as in
the existing pipeline; it is not a posterior-probability average.

Report power and neutral detection/FPR with counts and Wilson intervals, AF
distributions, accepted AF trajectories, raw focal scores, and representative
spatial profiles chosen by proximity to the cohort's median final AF. FPR is
the neutral-only error endpoint requested as FDR; no mixed-prevalence discovery
FDR is inferred from the target:null ratio. All FPR heatmaps share one scale.
Missing ALT/ALT pair classes are unavailable/no-call outcomes and coverage is
reported. Figures are saved as PNG and vector PDF, with editable large text,
without notebook or inline display.

## Run or resume

Default output: `/mnt/d/phase2simselection/sim/segregating_introgression_h400`.
Default resources: four workers, 26 GB total worker address-space limit, 20-GB
disk reserve, and the existing 3-TB volume-use ceiling. A per-study lock prevents
concurrent runners. SLiM 5 is required for pyslim 1.1 initialization metadata;
the older predefined-focal runner uses a separate SLiM binary and remains intact.

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
bash scripts/bootstrap_origin_onset.sh
bash scripts/launch_segregating_introgression.sh run
```

For another host, configure `SLIM_BIN`, `GAMMA_NATIVE_DIR`, a compiled
`bin/gamma_smc`, and `SIM_OUTPUT_DIR`. The bootstrap uses uv. The launcher accepts
`plan`, `simulate`, `audit`, `decode`, `analyze`, or `run`; plan performs no
simulations. Read `run_status.json` for the phase and `simulate_status.json` for
the heartbeat, active tasks, completed count, and failures. Per-region logs
distinguish prehistory, focal choice, loss, sample absence, and successful runs.

The previous pilot is stopped and preserved at `sim/origin_onset_pilot_q1`.
Its `stop_receipt.json` preserves the status at supersession. Live launch and
validation receipts are stored with the new outputs.

## Validation and launch record

Sixteen focused tests passed. Four small developmental cases (neutral and s=0.03
at each onset) completed simulation, independent hash/genotype/ascertainment
auditing, and Gamma-SMC decoding. They use 110 kb, 20 sampled diploids, and 100
pairs and never enter production denominators. Three final full-size cases
completed and were audited: I50 s=0.03, I10 s=0, and I10 s=0.03. The I50 s=0
preflight was transferred into the full queue after 17 recorded rejected
attempts; those seeds are not replayed. No completion or FPR estimate is claimed
for that unfinished neutral preflight.

Twenty synthetic PNG/vector-PDF plot pairs were generated. The power and AF
layouts were visually inspected, and PDFs checked for embedded raster images.
Synthetic scores are only layout/calibration tests. The validation record is
`results/segregating_introgression_validation/validation_summary.json`.

The complete queue was launched in a hidden WSL process with four workers.
`launch_receipt.json` records process IDs and the configuration; live state is
in the status files. A successful queue proceeds to audit, decode, and analysis.
