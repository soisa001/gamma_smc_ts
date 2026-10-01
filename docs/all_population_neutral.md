# Neutral pointwise calibration across six populations

Authorized 2026-10-01. Run 1,000 retained neutral regions for each of AFR, AMR,
EAS, EUR, MID, and SAS (6,000 total). Each population uses its own pointwise
median PHLASH history. This is a set of six separate single-recipient models,
not a joint demographic model with migration between the six populations.
The same hypothetical 2% Neanderthal ancestry pulse at 50 kya is imposed on
every population, including AFR and AMR, as requested.

## Scientific specification

- No selection, no de novo treatment, no population/time rescaling.
- Neanderthal effective diploid size 3,600; split at 500 kya; no added bottleneck.
- Generation time 25 years; mutation rate 1.25e-8/bp/generation and uniform
  recombination rate 1e-8/bp/generation.
- Simulate 11 Mb, choose the nearest eligible archaic-derived polymorphism to
  the midpoint within the central 1 Mb, then crop to 10 Mb with the focal at 5 Mb.
- The allele must be segregating in the full recipient census immediately
  after introgression at 50 kya. There is no 10-kya ascertainment cohort.
- Preserve the previous observation conditioning: retry if the chosen allele
  is lost or absent from the present-day sample. Fixation after onset is kept.
  Do not retrospectively choose a different surviving allele in the same run.
- Sample 200 diploids and the same fixed 10,000 haplotype pairs. ALT/ALT is the
  raw within-carrier-pair fraction, without carrier-mass weighting.
- Score only the prespecified focal site at T=10,20,30,40,50 kya. Report true
  TMRCA and Gamma-SMC posterior-mean calls separately. Decoder modes are
  `exp10=accurate` and `backward_alignment=fixed`, with scaled mutation rate
  0.00075 and recombination/mutation ratio 0.8.
- PHLASH histories retain their stored generation/Ne scale. Extend the first
  supported Ne constantly from 100 generations to the present, and the last
  supported Ne into the deeper ancestral past. This conditions calibration on
  the median history; demographic uncertainty is not integrated out. The
  historical PHLASH fitting mutation rate is not asserted by these artifacts.

Mutation rate, sequences, ascertainment, and seeds match the previous EAS
campaign. Its compatible completed I50 calibration nulls are read and verified
in place, along with compatible focal truth/decoded scores. Unfinished EAS
regions use separate new directories; the older selected/I10 campaign stays
paused. EAS seeds are unchanged, while the population label distinguishes the
other populations' deterministic seed streams. Reuse depends on scientific
parameters, runtime, pair manifest, and output hashes, not Git revisions.

## What the cutoff means

For each population, pair statistic, and T separately, use
`p = (1 + number of null scores >= target score)/(1 + N_null)` and p<=0.05.
This is **site-level calibration**, not a 10-Mb maximum or a correction across
populations, pair statistics, or T cutoffs.

`pointwise_cutoffs.csv` reports the exact critical boundary. A finite target
score must strictly exceed that boundary. With 1,000 nulls, at most 49 null
scores may equal or exceed the target, so the boundary is the 50th-largest null
score. Ties can make rejection more conservative; if the boundary equals 1,
no bounded raw-fraction score can reject. Missing ALT/ALT classes are no-calls
and rank below available null scores, retaining the prior unconditional
denominator convention. Available-only interpolated q95 is also supplied as a
descriptive column, clearly separate from the testing boundary.

Partial tables explicitly show actual null counts, available pair-class counts,
and `provisional=True`. Early completed regions can favor faster trajectories;
do not treat a partial table as the final 1,000-replicate calibration. No
independent neutral test cohort is generated in this cutoff-only experiment,
so these outputs do not claim an independently validated FPR.

## Workflow and outputs

Configuration: `configs/all_population_neutral.json`. One pool shares 20 CPU
slots across all six populations. Its workers have an aggregate 180-GB
address-space allowance (9 GB each), reserving 20 GB within the requested
200-GB cap for the paused older campaign and orchestration. Numerical library
threads and decoder threads are one per worker. The launcher records runtime,
binary and resource hashes, all default parameters, seeds, attempts, and paths.

Phases are `plan`, `simulate`, `decode`, `report`, or `run`. `run` simulates,
audits, and scores each region before moving to another, updates provisional
tables every 20 completed jobs, then creates PNG/vector-PDF plots and a ZIP.
`report` can refresh plots from the completed calibration receipts. A study
lock prevents overlapping writers. The run halts new submissions after a
failure; receipts allow deterministic resumption after the cause is fixed.

Outputs include `run_status.json`, `cohort_counts.csv`, `pointwise_cutoffs.csv`,
`null_scores.csv`, `neutral_cutoffs.pdf`, individual PNG/PDF figures, and
`neutral_cutoffs_bundle.zip`. Per-region simulation and scoring logs preserve
provenance and rejection reasons. Development preflight outputs use a separate
directory and never enter production denominators.

On the configured WSL analysis host (set the output, Python, SLiM, and native
library paths for another host):

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
bash scripts/run_population_neutral.sh run
```
