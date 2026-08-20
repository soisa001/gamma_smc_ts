# CoSi2 EAS and Han framework

This additive bundle renders parameter files for CoSi2 2.4.0, pinned
to commit `5da717cb8bf7324fb961e772a0c9042385de381a`. It does **not** contain production simulations.

## Canonical cells

| Cell | s | mutation birth (generations ago) | selection onset (generations ago) | deterministic endpoint from gate center | present population AF gate |
|---|---:|---:|---:|---:|---:|
| `eas_s0p01` | 0.01 | 1,475 | 1,475 | 0.200 | 0.195-0.205 |
| `han_introgressed_s0p01` | 0.01 | 2,400 | 645 | 0.234 | 0.195-0.205 |

The EAS mutation-birth ages are one-copy deterministic fixed-point centers on
the checked-in PHLASH median, not inferred ages. EAS uses 25 years/generation.
The Han/CHB model uses 29 years/generation, and the mutation is born as one copy
in Neanderthal at the explicitly configurable donor age. Selection starts in
CHB at generation 645, the integer realization of
18.7 kya, when 5R19's continuous archaic migration stops.

## Critical Han distinction

`OutOfAfricaArchaicAdmixture_5R19` does **not** contain a 1.2% pulse. It models symmetric,
continuous Neanderthal-Eurasian migration at 8.25e-06
per generation from the 60.7-kya OOA split until 18.7 kya. The
Ragsdale-Gravel estimate of 1.2% +/-
0.6% is a genome-wide ancestry result,
not a per-locus migration parameter. The simple symmetric two-state contact
expectation for the catalog interval is
1.181%, consistent with that
estimate.

Following the requested ascertainment, this framework separately applies an
operational CHB allele-AF screen at generation 645 around
1.2%
(1.1%-1.3%). Native
`sweep_mult_standing` enforces only the present CHB population-AF range. Set
`COSI_SAVE_TRAJ`, then retain a Han output only when `selfreq_3` passes this
migration-end screen. This is a user-specified locus ascertainment rule, not a
claim that genome-wide ancestry and every introgressed allele have equal AF.
Use the CLI `validate-trajectory` phase before replaying an accepted trajectory.

The endpoint gate is a **population** AF interval. CoSi2 subsequently draws the
sampled carrier count binomially, so a 200-haplotype sample is not constrained
to exactly 20%.

CoSi2's `MSweep::makeSweepModel` constructs selected/unselected siblings for
every population in the base model and looks up a sampling request for each one.
The Han parameter file therefore emits explicit `sample_size 0` requests for
YRI, CEU/OOA, Neanderthal, and ArchaicAFR, while sampling
200 haplotypes from CHB. The zero requests are required CoSi2
bookkeeping; they do not add sampled chromosomes.

The default 2,400-generation donor birth is a transparent age-grid starting
point, not a calibrated age. It predates the full Neanderthal-contact interval.
Scan donor birth ages or supply a separately validated trajectory after the
smoke phase. With the fixed 645-generation selection duration and configured
gate center p=0.012, the deterministic endpoint
summary is s=0.01: p=0.234. Joint intermediate-plus-endpoint acceptance
can be extremely low. For the narrow 1.1%-1.3% default Han boundary gate, the
current feasibility estimate is roughly one joint acceptance per 28 million
trajectory draws; treat this as a planning estimate, not a calibrated rate.
`run_commands.sh` is a syntax smoke, does not implement an outer rejection loop,
and does not promise that its one Han draw will pass the boundary validator.

stdpopsim records `mutation_rate=None` for 5R19 because the model was calibrated
using recombination. The mutation rate and recombination map in these parameter
files are explicitly labeled study simulation overrides. The one-line constant
map starts at bp 1 because CoSi2 v2.4 rejects a zero start; its first rate covers
the region start and remains in force through the sequence end.

## Intended phases

1. `plan`: render and strictly validate this bundle; run no CoSi2 simulation.
2. Smoke: run one fixed-seed trajectory/coalescent replicate per cell from
   `run_commands.sh` and inspect stderr plus the saved trajectory.
3. Han gate: validate CHB AF at the migration-end/selection-start boundary.
4. Replay: strip `popsize_*` columns, validate every population-frequency
   column, and use `COSI_LOAD_TRAJ` for a later full-region replicate.
5. Scale only after smoke acceptance and runtime profiling. Keep raw ARG to
   tskit conversion as a separately tested downstream phase.

## Build CoSi2 and render this framework

```bash
git clone git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
git pull --ff-only origin AOU_run_opt
uv sync --frozen --extra selection
uv run python scripts/run_cosi2_framework.py plan

git clone --branch v2.4.0 --depth 1 https://github.com/broadinstitute/cosi2.git external/cosi2
test "$(git -C external/cosi2 rev-parse HEAD)" = "5da717cb8bf7324fb961e772a0c9042385de381a"
cd external/cosi2 && ./configure && make -j4 && cd ../..
```

Run commands are in `run_commands.sh`. They intentionally request one smoke
replicate per cell, use fixed seeds, and run from this directory because CoSi2
resolves `recomb_file` relative to the current working directory. All mutable
outputs are written outside this exact, immutable bundle under
`../cosi2_framework_runs/smoke/`; this lets `validate-trajectory` revalidate the
bundle inventory after a run. A smoke output is not accepted until the
trajectory and endpoint validators pass. Increase replicate counts only after
measuring the joint gate acceptance rate.

CoSi2 emits ms-format haplotypes and optional raw ARG edges, not native tskit
tree sequences. Any ARG-to-tskit bridge should be a separate validated phase.
