# Focal origin/onset pilot: 10 targets against 1,000 nulls

Authorized 2026-09-24. Configuration: `configs/origin_onset_pilot.json`.
The design and biological assumptions are in `proposed_origin_onset_grid.md`.

The selection grid is 0.001 through 0.010, then 0.020 and 0.030. Four treatments
(I50, I10, D50, D10) yield 480 selected test regions. Each of three null families
has 1,000 calibration regions and 10 independent unselected test regions:
3,510 retained regions in total. I50 and I10 share the introgressed null family.
The two de novo families differ in focal mutation birth time (50 or 10 kya),
even though their background demographic history is the same. Their age-matched
focal null distributions must not be pooled.

The first ten existing selected I50 replicates at s=0.001--0.010 can supply
100 regions after checksum, model, runtime, seed, and genotype verification.
The old nearest-ordinary-variant neutral archive is not reused as a matched
focal-allele control. Missing work is about 3,410 retained regions, plus rejected
attempts needed for survival and sample observation.

All selected and unselected new regions go through the same SLiM focal-allele
engine. Introgressed focal alleles are fixed in the archaic source before the
50-kya pulse. De novo focal alleles are placed in exactly one EAS genome at
onset. Population survival and observation in the sampled panel are required;
fixation is retained. No terminal AF target or lower AF cutoff is applied.
No neutral test target is a calibration replicate. Role enters the seed identity.

Default resources are four single-threaded workers and a 26-GB total worker
address-space budget (6.5 GB per worker), with a 20-GB disk reserve and the
existing 3-TB volume-use ceiling. Coordinator/runtime overhead is additional.
Introgressed arms use the archived Q=5 scaling; de novo arms use Q=1 so one copy
means one copy at the original population size. Within each comparison the
null and target scaling match. Cross-origin numerical convergence is unproven.
The Q=1 single-copy survival-conditioned nulls are substantially more expensive
than the archived introgressed selected simulations; no completion ETA is assumed.

## Phases and saved outputs

- `plan`: write/check deterministic task, pair, position, parameter, and runtime
  manifests. No simulations or decoder calls.
- `simulate`: reuse audited compatible selected inputs and generate missing
  simulations. Save raw/cropped trees, focal carrier masks, accepted population
  AF trajectories, generated SLiM scripts, attempts, logs, and SHA-256 receipts.
- `audit`: independently read/hash saved trees and verify focal genotype/sample
  alignment, lengths, sample size, carrier masks, and observed allele frequency.
- `decode`: save Gamma-SMC posteriors and pair-class raw fraction profiles for
  truth and decoding. Check pair/position manifests and reproduce native all-pair
  summaries. No simulation generation is performed by this phase.
- `analyze`: require complete targets and nulls; calculate focal empirical p
  values separately for every null family, statistic, source, and cutoff; save
  tables and PNG/vector-PDF figures. No per-replicate best-cutoff tuning.
- `run`: execute simulate, audit, decode, analyze sequentially, resuming each
  phase from verified receipts. The study lock prevents concurrent runners.

The output root on this host is
`D:/phase2simselection/sim/origin_onset_pilot`, or
`/mnt/d/phase2simselection/sim/origin_onset_pilot` in WSL.
Read `run_status.json` for the overall active phase and any failure, and
`simulate_status.json`, `audit_status.json`, or `decode_status.json` for detailed
progress. Region-level `status.json` and logs identify slow or failed tasks.
Scientific cache identities depend on model inputs, sample/pair manifests,
runtime semantics, and relevant file hashes, not the repository commit.
Adding target replicates later does not reseed or invalidate completed work.

The empirical p is `(1 + number of null scores >= target score) / 1001`.
Use p<=0.05 and cutoffs 10,20,30,40,50 kya. ALT/ALT uses its own raw within-pair
fraction, never carrier-mass weighting. Missing ALT/ALT classes are NA and
produce no discovery; tables record availability. Report the number of calls
out of 10 and Wilson intervals. A 10-target pilot has 10-percentage-point
resolution and cannot precisely validate a 5% neutral false-positive rate.
The user's neutral-only error comparison is labeled neutral detection/FPR;
the 10:1,000 testing/calibration ratio is not selection prevalence for FDR.

Trajectories describe accepted survivors. No establishment/loss probability is
inferred from the rejected-attempt counts, which omit internal SLiM restores.
Legacy trees do not acquire historical trajectories retroactively. Full
spatial profiles are supplementary; all power/error calculations use only the
prespecified focal allele coordinate. Figures remain on disk, not in notebook
output. Independent simulation replicates, not pairs, are the sampling units.

## Clone/pull, set up, and run or resume (Linux/WSL)

This host already has compatible SLiM 4.2.2 and Gamma-SMC binaries. On another
host, provide `SLIM_BIN`, `GAMMA_NATIVE_DIR`, and a compiled `bin/gamma_smc` first;
do not use these paths as portable HPC binary locations. `SIM_OUTPUT_DIR` and
`LEGACY_SIM_ROOT` can be overridden; an absent legacy root simply disables reuse.
The bootstrap uses uv and an isolated `.venv-origin`, with simulation versions
matching the original archive. It does not replace another task's environment.

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
bash scripts/launch_origin_onset.sh run
```

Do not start a second runner while the existing study lock is held. A read-only
status cell is:

```bash
%%bash
set -euo pipefail
cat /mnt/d/phase2simselection/sim/origin_onset_pilot/run_status.json
cat /mnt/d/phase2simselection/sim/origin_onset_pilot/simulate_status.json
```

Focused development validation uses `tests/test_origin_onset.py`, small real
introgressed/de novo simulations in a separate validation output root, and
synthetic-data plotting checks. Validation-region dimensions and counts differ
from production and never enter its calibration or power denominators.

## Validation and launch record

Nine focused tests passed (eight unrelated cases deselected). Five small real
cases completed simulation, hash/genotype auditing, and Gamma-SMC decoding:
introgressed s=0, introgressed 10-kya s=0.03, de novo 50-kya s=0, de novo
50-kya s=0.03, and de novo 10-kya s=0.03. The latter verifies the predicted
focal truth-recency bound above its birth time. The small 20-haplotype D10
neutral check was stopped after one completed unobserved survivor and part of
a second attempt; that bounded validation cancellation is recorded separately
and is not a production failure. Production uses the authorized 400 haplotypes.

Thirty synthetic-data PNG/PDF plot pairs were generated; the full coefficient
grid was visually checked and a title/legend overlap corrected. PDFs were
checked for absence of embedded raster images. See
`results/origin_onset_pilot_validation/validation_summary.json` for the validation
parameters and checks. These validation scores do not enter production results.

A full-size preflight reused an audited I50 s=0.001 region and generated an
I50 s=0.02 region and matched introgressed null. The complete 3,510-region pilot
was then launched in a hidden WSL process. `launch_receipt.json` in the output
root records its process identities, resources, output path, and launch context.
Live progress is in the status files, not a frozen count in this document.
