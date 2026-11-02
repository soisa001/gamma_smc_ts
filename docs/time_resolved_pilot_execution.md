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
