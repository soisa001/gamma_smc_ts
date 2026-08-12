# No-introgression EAS sweep study

## Status and scientific question

The implemented study asks how a single de novo beneficial allele changes
Gamma-SMC coalescence summaries under uncertainty in the East Asian
effective-population-size history, without archaic introgression. It exposes
separate plan, simulation, decode, and plot phases. The bounded local execution,
including incomplete cells and validation scope, is reported separately in
[`RUN_RESULTS.md`](RUN_RESULTS.md); this document defines the reusable run
contract rather than promoting a partial pilot to a full power study.

The pilot seeks one accepted simulation per design cell. It is descriptive,
not a power study: matched neutral simulations, calibrated null distributions,
and p-values are deliberately deferred. A cell that exhausts its acceptance
budget or declared timeout remains a failed cell; it is not retained as data.

## Design matrix

| Component | Pilot specification |
|---|---|
| Simulator | stdpopsim 0.3.0 with SLiM 4.2.2 and tree-sequence recording |
| Sequence | 10,000,000 bp; focal selected site at 5,000,000 bp |
| Sample | 100 present-day diploid EAS individuals (200 haplotypes) |
| Mutation rate | 1.25e-8 per bp per generation |
| Recombination rate | 1e-8 per bp per generation |
| Generation time | 25 years |
| Selection coefficient | 0.01, 0.005, 0.001 |
| Dominance coefficient | 0.5 |
| SLiM scaling | `Q=10` |
| Forward burn-in | `slim_burn_in=0.1` N generations, followed by stdpopsim recapitation |
| Origin | one de novo copy in the EAS population |
| Present-day selected-allele frequency | nominal 0.50; inclusive population and sampled-panel acceptance band 0.45-0.55 |
| Demographic sensitivity paths | empirical pointwise 2.5%, 50%, and 97.5% EAS Ne curves |
| Attempt budget | at most 20 completed external sample-panel draws per design cell; internal conditioning separately timeout-bounded |
| Accepted replicates | one per design cell if all acceptance checks pass |

The implemented frequency band remains an exploratory placeholder. It should
be replaced or expanded only after a separately documented survey of
frequencies for eligible archaic SNVs; the current 0.50 target is not claimed
to represent that empirical distribution.

## EAS demography and uncertainty

The exact checked-in [`resources/EAS.npz`](resources/EAS.npz) artifact has schema
`phlash.aou.log-ne-mvn/v1`. It stores 100 posterior-median PHLASH bootstrap Ne
curves on a common 10,000-point geometric grid from 100 through 40,000
generations ago. The implementation calculates the 0.025, 0.50, and 0.975
quantiles of `bootstrap_ne` independently at every time point and uses those
three curves as lower-envelope, median, and upper-envelope sensitivity
histories.

These are empirical **pointwise** intervals across 100 curves. They are not a
simultaneous confidence band, three coherent posterior draws, or known
biological truth. Outputs must retain the `q025`, `median`, or `q975` label rather
than describing all three paths as equally likely demographic histories. The
resource's arrays, checksum, and source commit are recorded in
[`resources/EAS.provenance.json`](resources/EAS.provenance.json) and
[`resources/README.md`](resources/README.md).

No migration or archaic-admixture event is added to these EAS Ne histories.
Calendar labels use 25 years per generation throughout.

## Focal-mutation origin and frequency conditioning

The selected allele is introduced as one copy at the central site, so this is a
hard sweep from a de novo mutation rather than selection on standing
variation. The origin generation is fixed by each tracked specification; it is
not sampled from a proposal window during SLiM execution.

The implemented origin ages come from the deterministic
`logit_fixed_frequency_grid` approximation under each time-varying Ne curve,
using a one-copy starting frequency and the nominal final frequency 0.50:

| Selection coefficient | `q025` | `median` | `q975` |
|---:|---:|---:|---:|
| 0.01 | 1,275 gen / 31.875 kya | 1,289 gen / 32.225 kya | 1,308 gen / 32.700 kya |
| 0.005 | 2,635 gen / 65.875 kya | 2,725 gen / 68.125 kya | 2,781 gen / 69.525 kya |
| 0.001 | 15,523 gen / 388.075 kya | 15,728 gen / 393.200 kya | 16,306 gen / 407.650 kya |

These values are fixed simulation inputs recorded in `specifications.tsv` and
`origin_age_contracts.json`, but they are approximate trajectory-derived ages,
not exact biological age estimates. This is the meaning of
`origin_age_is_exact=false` in the specification table.

At `Q=10`, a single introduced copy is in the rescaled population of roughly
`Ne/Q` diploids. The solver therefore uses the approximate initial frequency
`Q/(2Ne)`, before integral population-size rounding, rather than the unscaled
`1/(2Ne)`. stdpopsim also rounds each event to an integral scaled SLiM tick.
The specification records `Q`, the one-copy semantics, and the effective
rounded mutation and survival-check times. The unrounded ages remain requested
design inputs, not exact biological or forward-event ages.

For a draw to be retained, the focal allele must survive, the present-day
population frequency condition and observed 100-diploid sample frequency must
both lie in the inclusive 0.45-0.55 interval, and the sample must contain at
least two `ref/ref`, one heterozygous, and two `alt/alt` individuals. The runner
uses at most 20 completed external sample-panel draws with deterministic seeds.
A row is written only after stdpopsim's internally frequency-conditioned SLiM
call returns; internal losses and frequency misses retried inside that call are
not separate `external_draws.tsv` rows. A contract-bound ledger records each
launched seed and distinguishes an interrupted internal trajectory from a
completed sample-panel draw, so a retry advances rather than repeats the seed.
Completed draws record frequency, genotype counts, and acceptance; failure
never widens the interval.

This conditioning parallels the published fixed-frequency design while making
the mutation origin explicit. The published `cosi2` experiment used 100
replicates of 100 haploid 10-Mb samples at 1.25e-8 mutation and 1e-8
recombination per bp per generation, but noted that mutation age could not be
controlled in that mode. This pilot intentionally changes the demography,
selection grid, forward simulator, diploid sample, and replicate count.

The paper's separate SLiM experiment instead fixed mutation ages at 100, 500,
750, 1,000, 2,000, 3,000, 4,000, and 5,000 generations and did not control the
sample allele frequency. This study does not take the full Cartesian product of
that age grid with a 50% frequency condition. It uses the s- and
demography-specific trajectory-derived ages above, because many fixed-age /
fixed-frequency combinations would be effectively unreachable. A direct age
sensitivity replication remains a distinct follow-up design.

## Forward-simulation approximation and SLiM version

The runtime uses stdpopsim's `Q=10` rescaling and a short
`slim_burn_in=0.1` N-generation forward burn-in. stdpopsim then performs its
standard recapitation/rescaling step and overlays neutral mutations. This is an
exploratory approximation: SLiM rescaling is not equivalent to the unscaled
selected process, and a 0.1 N burn-in relies heavily on recapitation. Any
biological interpretation requires sensitivity checks against less aggressive
scaling and longer burn-in; none is claimed completed here.

The general repository prefix retains SLiM 5.2 for diagnostics, but unmodified
stdpopsim 0.3.0 scripts use the removed `Subpopulation.genomes` API under SLiM
5.2. The EAS study therefore defaults to a separate, version-checked SLiM 4.2.2
prefix. A compatibility-transformed SLiM 5.2 run is diagnostic only, not a
production study engine or accepted data source.

## Gamma-SMC decoding and summaries

Decode each accepted 10-Mb simulation at all of these calendar thresholds:

| Cutoff | Generations at 25 years/generation |
|---:|---:|
| 1 kya | 40 |
| 4.5 kya | 180 |
| 10 kya | 400 |
| 20 kya | 800 |
| 30 kya | 1,200 |
| 40 kya | 1,600 |
| 50 kya | 2,000 |

Pair manifests and focal-site genotypes must be fixed before decoding. Report
the posterior probability `P(TMRCA < x)` and posterior mean TMRCA for:

- all prespecified pairs;
- pairs drawn from homozygous-reference individuals (`ref/ref`); and
- pairs drawn from homozygous-selected individuals (`alt/alt`).

Heterozygotes may contribute to the overall summary but are excluded from the
strict `ref/ref` versus `alt/alt` contrast. Acceptance requires at least two
sampled diploids in each homozygous class and at least one heterozygote. Record
the number of individuals, haplotypes, and pairs in every class.

The expected sweep direction is **higher** `P(TMRCA < x)` and **lower** mean
TMRCA among selected-allele haplotypes near the focal site. A reversal is a
result to investigate, not the expected positive-selection signature. The
probability-versus-cutoff curve is a TMRCA CDF. If its area is summarized, call
it an integrated CDF over a stated cutoff range; it is not ROC AUC. ROC AUC
requires labeled positive and negative replicates and a varying decision
threshold.

`ref/ref` is an internal genotype-class comparator from the same selected
simulation, not a matched neutral replicate. Accordingly, this pilot will not
report neutral-calibrated p-values, false-positive rates, power, or ROC AUC.

## Produced figures and tables

Figures are written as both PNG and PDF with large, letter-page-readable text.
The implemented report set is:

1. the three pointwise EAS Ne sensitivity paths and proposed origin times;
2. a complete 3-by-3 simulation-status grid, including timeouts and last sample AF;
3. accepted sample frequency and completed-draw effort by design cell;
4. genome-wide `P(TMRCA < x)` and mean-TMRCA heatmaps by genotype class;
5. focal-site CDFs and `alt/alt - ref/ref` probability contrasts;
6. focal-site posterior mean TMRCA contrasts; and
7. tree-sequence truth versus Gamma-SMC estimates where the estimands match.

Machine-readable specifications, completions, manifests, and checksums carry
the full parameter and seed provenance. Each figure carries the dimensions
needed for its panel; the manifest connects it to the complete run contract.

## Reproducibility contract

The implementation separates planning, simulation/acceptance, integrity
checks, Gamma-SMC decode, and reporting into restartable phases. Cache validity
depends on the full design and sample/pair manifests plus checksums and
parseability of the relevant outputs, not merely on a Git commit. Logs preserve
resolved defaults, paths, software versions, resource hashes, commands, seeds,
attempt history, and validation status. Corrupt or incomplete artifacts must
never satisfy a cache check.

When a nonzero per-specification timeout is declared, expiry terminates that
worker and its SLiM descendant and writes a `timed_out` failure record. A
timeout, including a short diagnostic timeout, is failure metadata rather than
simulation data and must not enter acceptance tables, plots, or biological
interpretation. Exhausting the 20-draw acceptance budget is likewise failure,
not an accepted replicate.

## References

- Schweiger and Durbin. [Ultrafast genome-wide inference of pairwise coalescence times](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538485/).
- Shlyakhter et al. [cosi2: efficient simulator of exact and approximate coalescent with selection](https://doi.org/10.1093/bioinformatics/btu562).
- Haller and Messer. [SLiM](https://messerlab.org/slim/).
