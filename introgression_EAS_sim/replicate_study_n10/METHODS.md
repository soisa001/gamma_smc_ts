# Replicated introgressed-sweep methods

## Scope

This campaign repeats the 27 AncientEurasia introgressed-sweep specifications
defined in [`../METHODS.md`](../METHODS.md). It requests ten independently
seeded slots per biological cell; it does not assume or guarantee that all ten
will pass the bounded search and acceptance checks. It changes replication,
restartability, and aggregation only; it does not change the demography,
selected-site model, sampling design, or Gamma-SMC estimands.

The completed bounded campaign had 270 requested slots: 155 accepted and
decoded across 24/27 represented cells, 104 `exhausted`, 11 `timed_out`, and 0
`failed`. Of the represented cells, 22 had at least two paired trajectories and
were eligible for the macro headline.

Matched neutral simulations are intentionally absent. Consequently, every
result is a descriptive within-selected-simulation comparison and cannot
estimate a p-value, false-positive rate, sensitivity, specificity, power, or
ROC AUC.

## Replicate unit and randomization

The replicate unit is one independently seeded, frequency-conditioned SLiM
trajectory and its 100-diploid Han sample. Diploid individuals, haplotype
pairs, genomic positions, and TMRCA cutoffs within a trajectory are not treated
as independent replicates.

Each base specification is expanded to stable identifiers
`<base_specification_id>__rep001` through `__rep010`. For external draw `d`, the
seed is deterministically derived from SHA-256 of
`base_seed:expanded_specification_id:d`, using fixed base seed 20240523. Adding
future replicate slots therefore does not change existing identifiers, seeds,
or caches.

## Prespecified design and acceptance

The full grid is:

- selection coefficients 0.01, 0.005, and 0.001;
- inclusive present-day Han/sample AF intervals `[floor, 0.95]`, with floors
  0.10, 0.20, ..., 0.90; and
- ten requested, independently seeded trajectory slots per selection-by-floor
  cell, or 270 requested slots; accepted counts may be smaller.

The AF contracts are nested floors, not disjoint bins. Results retain both the
requested floor and achieved sample AF. Acceptance additionally requires 100
sampled diploids, one focal mutation, at least two ref/ref individuals, one
heterozygote, and two alt/alt individuals. No threshold is relaxed on retry.

## Bounded search and missingness

For every incomplete slot, a local invocation sets the effective timeout to
`min(5 minutes, remaining time in the 30-minute cumulative budget)`. A slot
with no cumulative time remaining becomes `exhausted` without another launch.
Across resumptions, it also has ceilings of 40 launched internal trajectories
and 20 completed external sample-panel draws. The campaign uses uniform passes
over all 270 slots rather than stopping preferentially once an easy cell has
enough accepted trajectories.

Every launched seed is persisted before SLiM starts. Completed sample-panel
draws and interrupted internal conditioning trajectories are recorded
separately. A timeout advances the deterministic seed and is never counted as
a completed or accepted draw. Accepted tree checkpoints permit recovery after
a crash between tree return and final artifact materialization. Per-slot
cross-process locks prevent concurrent processes from consuming the same seed.

Conditioning failures and timeouts are nonrandom. All effect summaries are
therefore conditional on accepted trajectories, and cell coverage is always
reported alongside effects. Missing cells are retained explicitly rather than
silently omitted.

## Integrity and software contract

Simulation reuse requires the exact biological/runtime contract plus SHA-256
agreement for the tree sequence, sample and pair manifests, truth profiles,
and focal contrasts. Invocation histories are bound to the simulation contract
hash. Atomic writes use bounded retry for transient Windows/OneDrive access
denials and otherwise fail closed.

A checksum-current complete cache is revalidated and then treated as a no-op:
it does not consume another seed, launch SLiM, or rerun Gamma-SMC. Decoder
worker/thread settings are execution provenance rather than scientific cache
identity; the thread count is retained and validated in class completion
records under the declared output-invariance contract.

The production simulation engine is stdpopsim 0.3.0 with SLiM 4.2.2. The
runner rejects SLiM 5.2 because the unmodified stdpopsim-generated script is
incompatible with it. Gamma-SMC decoding uses the repository's Linux AVX2
binary and records its SHA-256 in every class decode contract and the campaign
manifest.

Some early tasks were created before per-task source hashes were recorded and
therefore retain `legacy_unknown` source-epoch evidence. Artifact validation
confirms current checksums, schemas, and parseability, but cannot prove semantic
equivalence across unbound implementation epochs. A frozen-source rerun is
required to remove that provenance limitation.

## Gamma-SMC decoding

Every accepted 10-Mb tree sequence is decoded for the overall sample and the
ref/ref, heterozygous, and alt/alt genotype classes. After a read-only
throughput audit, the completed local campaign used 24 concurrent decoder
workers with one thread per worker. These are execution/provenance settings;
the immutable plan may retain its original four-thread default, which is not a
claim about the completed execution and does not invalidate an otherwise
validated output cache. The decoder uses a 10-kb output stride and cache size,
mutation rate 1.25e-8, recombination-to-mutation ratio 0.8, 25 years per
generation, and thresholds 1, 4.5, 10, 20, 30, 40, and 50 kya.

Truth and Gamma-SMC are compared at the 5-Mb focal site. The two primary
internal effects are:

1. alt/alt minus ref/ref normalized area under the within-genotype TMRCA CDF
   integrated from `(0, 0)` through 50 kya, expected to be positive; and
2. alt/alt minus ref/ref mean TMRCA, expected to be negative.

The first quantity is an integrated CDF area, not classifier/ROC AUC.

## Replicate summaries

Aggregation starts from all 270 requested slots and retains every terminal
status, including accepted, exhausted, timed-out, and failed slots. It adds
available checksum-validated truth trajectories from accepted simulations and
Gamma-SMC/paired trajectories only from accepted simulations with valid
decodes. It preserves the exact base specification, selection coefficient, and
AF floor and never pools trajectories from different AF cells for cell-level
effects. Reported outputs include:

- requested, attempted, accepted, and decoded trajectories per cell;
- achieved AF and genotype-count distributions;
- per-cell mean, median, standard deviation, and seeded 5,000-draw trajectory
  bootstrap 95% intervals for each focal effect;
- Wilson 95% intervals for expected-sign rates;
- paired truth-versus-Gamma bias, MAE, RMSE, Pearson and Spearman association;
- truth expected-sign rate and Gamma expected-sign rate separately;
- raw truth/Gamma sign agreement; and
- Gamma expected-sign detection conditional only on trajectories whose tree
  truth has the expected sign.

Raw sign agreement may count jointly contrary signs as agreement. Conditional
detection never counts a jointly contrary pair as success. All uncertainty
uses trajectories as the sampling units. The headline cross-cell result is the
cell-equal macro mean over biological cells with at least two paired
accepted-and-decoded trajectories (22 cells in the completed campaign); no
confidence interval is computed for that macro mean. The trajectory-weighted
micro summary is secondary because cells with more accepted-and-decoded
trajectories contribute more weight. Outputs
are compact campaign reports and machine-readable tables rather than copies of
bulk intermediates. Figures use readable letter-size canvases, large text, and
both PNG and PDF formats.
