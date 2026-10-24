# Neutral SV mutation-age simulation pilot

This pilot estimates the age of a structural mutation from modern SNP
observations and phased carrier status. It compares inferred ages with known
simulation ages, and separately measures how much information an exact modern
genealogy would provide. The focal events are neutral, single-origin binary
SV surrogates, appropriate as an initial model of unique MEI-like events.
They are not simulated structural sequencing reads, inversions with suppressed
recombination, recurrent SVs, or selected alleles.

## Generative model

The model reuses the repository's EAS PHLASH median size history and its hybrid
introgression demography: 2% Neanderthal ancestry introduced at 50 kya, modern
and archaic separation at 500 kya, donor diploid Ne=3,600, generation time
25 years, mutation rate 1.25e-8 per base per generation, and mean recombination
rate 1e-8 per base per generation. The PHLASH input is checksum-validated.
These rounded demographic choices reproduce the project's specification;
they are not a newly fitted demographic model.

There are 160 independent 500-kb regions, each with 400 present-day EAS
haplotypes and two donor haplotypes sampled at 60 kya. The donor sequences
are used only for truth/ascertainment and never enter Gamma-SMC decoding.
Ancestry simulation uses 100 generations of discrete-time Wright-Fisher
followed by the standard coalescent, with recorded migration events.
Five-kb recombination hotspots recur every 100 kb with a fivefold multiplier;
rates are normalized to preserve the mean. Gamma-SMC uses its constant-rate
model, so the decoding comparison includes recombination-model mismatch.

An independent neutral candidate-event overlay is generated on the genealogy
at rate 5e-8, with continuous sites to avoid recurrent focal events. Candidate
order is randomly permuted. At most one event of each history is retained
per region, in its central 60%, conditional on at least two modern carriers
and two noncarriers:

- **Post-introgression:** age below 50 kya, origin in the modern population,
  and every carrier's locus descends through the introgression pulse.
- **Imported old:** age between 50 and 500 kya, mutation on a donor-population
  branch, and all modern carriers descend through that pulse. Population of
  mutation is traced using branch migration records, not merely node labels.
- **Ancestral retained:** age above 500 kya, observed in a sampled archaic
  chromosome and at least two modern carriers without introgression at that
  locus. Additional imported copies may coexist. This establishes ancestral
  retention; it does not impose a particular discordant species-tree topology.

Regions lacking a qualifying event are retained and logged, not resimulated
until they satisfy an artificial quota. This is a stratified, observed-event
benchmark and does not estimate natural prevalence or discovery FDR.

SNPs are generated with an independent JC69 overlay at the stated mutation
rate. Candidate events are not inserted as SNP emissions. Modern samples are
exported in a phased VCF using exact one-based conversion, while focal query
coordinates and pair manifests use zero-based indices. The idealized pilot
has correct phase, SV genotype and polarization; it does not include read
mapping or genotype-error simulation.

## Observation and decoding design

For each event, sample up to 96 carrier/carrier pairs and 96 carrier/noncarrier
pairs. For a separate oracle-ancestry sensitivity method, sample up to 96
opposite-SV-state pairs with both members on introgressed ancestry. Sampling
uses genotype and ancestry class, never true TMRCA or mutation age. Each method
uses its own event's pair set, even though pair unions are decoded efficiently
over the whole region.

The modern SNP VCF is decoded at each focal coordinate using the existing
Gamma-SMC fork. The shape/rate posterior parameters are retained. The shared
scaled mutation rate is 0.00075; conversion to years is
alpha/beta * theta/(2*mu) * generation_time.

Each region is decoded twice: a single-base focal mask, and a 2,001-bp mask
around each focal anchor. Both conditions use positive callable BED intervals.
Their difference measures sensitivity to losing nearby SNP observations,
not robustness to every kind of structural-variant error.

## Estimators and evaluation

The deterministic split is by independent region: replicate modulo four
0/1 is training, 2 is calibration, and 3 is testing. All events and both mask
conditions from one region remain together.

1. **Carrier mean:** average decoded carrier/carrier posterior mean, treated
   as a deliberately naive age estimator with no claimed age interval.
2. **Raw composite CDF:** multiply carrier CDFs and carrier/noncarrier survival
   functions, using a uniform age prior from 100 to 5,000,000 years. Evaluate
   on a 600-bin logarithmic grid with linear-age bin widths in the integral.
   Report median and 5th/95th percentiles. These are uncalibrated composite
   intervals, not asserted Bayesian credible intervals.
3. **Balanced composite CDF:** give each pair class total weight one by averaging
   log CDF/survival scores within class, then use the same prior and grid.
4. **Composite with oracle ancestry:** add introgressed carrier/noncarrier
   constraints to the raw composite, removing exact pair duplicates already
   present in the general carrier/noncarrier set. Distinct pairs still share
   genealogy; this remains an explicitly uncalibrated sensitivity score.
5. **Simulation-trained kNN:** 25 nearest training events after training-only
   standardization. Features are carrier/noncarrier counts, seven quantiles
   of decoded posterior means per pair class, and mean posterior probability
   of coalescing before 50 kya. Neither true ancestry labels nor coordinates,
   origin categories or true trees enter these features. Predict the median
   log age and 5th/95th percentiles among neighbors.
6. **Frequency-only kNN:** the same training procedure using only the two
   genotype counts. This tests whether TMRCA features add useful information.
7. **Trained composite CDF:** divide the raw log compatibility by a temperature
   from {1,3,10,30,100}, with either a uniform-age or log-uniform-age prior.
   Select the combination by mean log-density loss on training events only.
   Evaluate both its unexpanded interval and a region-calibrated expansion.
   This directly tests whether attenuation of correlated pair constraints is
   more useful than the fixed class-balanced weighting.

For each kNN version and the trained composite, a second interval is expanded in log age using separate
calibration regions. The calibration score is the maximum interval error over
events within a region, and the finite-sample 90% order statistic supplies the
expansion. This targets region-level simultaneous coverage under exchangeable
simulation regions, not conditional coverage for every origin category or
validity under a different demographic/ascertainment distribution.

Point errors use symmetric fold error: max(estimate/truth, truth/estimate).
Report medians, within-twofold accuracy, signed median bias, interval coverage,
and interval width. A pre/post call requires the entire interval to fall on
one side of 50 kya. Report its call rate, wrong-call fraction over all targets,
and error rate among calls. Results are stratified by history and include the
post-introgression/imported-old subset without the older ancestral class.
Uncertainty summaries resample independent regions, not pairwise comparisons
or individual events, in 500 bootstrap draws with seed 1729.

The exact-tree reference computes L=max carrier/carrier TMRCA and U=min
carrier/noncarrier TMRCA using all modern haplotypes. Each known mutation age
must fall in [L,U]. The reference reports how often [L,U] spans 50 kya and how
often a truly old imported allele has a young carrier crown. These are logical
tree bounds, not 90% posterior intervals.

## Reproducibility and files

Configuration: `configs/sv_mutation_age.json`. Entry point:
`python -m gamma_smc_aou.sv_mutation_age`. Launcher:
`scripts/run_sv_mutation_age.sh`. Phases are `simulate`, `audit`, `decode`,
`analyze`, and resumable `run`.

Default resources are four single-threaded workers and a total worker
address-space budget of 26 GiB. The existing seed base 20380101 and pair seed
1729 are retained; deterministic named stages generate independent ancestry,
candidate, observation and pair streams. Configuration, runtime versions,
source/decoder hashes, commands, masks, pair manifests, seeds and checksums
are recorded. Git commit is provenance, not a simulation-cache key.

Per-region outputs retain ancestry/SNP tree sequences, the observed VCF,
separate observation and truth manifests, raw Gamma posteriors, compact arrays,
decoder summaries and logs. Cached outputs are checksum-validated; an
incompatible or corrupt cache raises an error rather than silently overwriting
it. The audit rereads trees and verifies mutation branch bounds, carrier
identity, ancestry membership and pair classes.

Aggregate outputs include `metrics.tsv`, `metric_bootstrap.tsv`,
`predictions.tsv`, `features.tsv`, `oracle_bounds.tsv`, `oracle_summary.tsv`,
`truth_summary.tsv` and `calibration.json`. Figures are written to PNG and
editable vector PDF with the Agg backend. No plot is emitted to notebook or
chat output.

## Clone, install, run or resume on Linux

The native bootstrap is needed only on a machine without this repository's
compiled decoder. It installs rootless native dependencies; the Python setup
uses uv and the existing pinned origin-workflow constraints. On an HPC system
with supplied native dependencies, set `GAMMA_NATIVE_DIR` instead.

```bash
%%bash
set -euo pipefail
repo="$HOME/gamma_smc_ts"
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
if [ ! -x bin/gamma_smc ]; then
  bash scripts/bootstrap_uv.sh --skip-tests --python 3.12
fi
export UV_BIN="$(command -v uv || printf '%s' "$repo/.tools/uv-bin/uv")"
bash scripts/bootstrap_origin_onset.sh
export GAMMA_NATIVE_DIR="${GAMMA_NATIVE_DIR:-$repo/.native}"
export SIM_OUTPUT_DIR="$HOME/sim/sv_mutation_age_20261024"
bash scripts/run_sv_mutation_age.sh --phase run
```

Use `--replicates 4` for a smoke run; it is insufficient for training/calibration.
An initial small run can be extended with the same configuration and seeds.
No controlled or external participant data are accessed by this workflow.
