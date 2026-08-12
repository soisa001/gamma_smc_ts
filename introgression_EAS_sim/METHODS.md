# Introgressed EAS sweep study

## Status and scientific question

The implemented study asks how Gamma-SMC coalescence summaries respond when a
genuinely Neanderthal-lineage allele enters the ancestor of Han through
introgression and is beneficial from the onset of that pulse. It exposes
separate plan, simulation, decode, and plot phases. The bounded local execution,
including accepted and timed-out cells, is reported in
[`RUN_RESULTS.md`](RUN_RESULTS.md); this document defines the reusable methods
contract.

The pilot seeks one accepted trajectory per requested specification. It is
descriptive: matched neutral simulations, null calibration, p-values,
false-positive rates, and power are deferred. A specification that exhausts
its acceptance budget or declared timeout remains failed; it is not retained
as data.

## Ancient Eurasian demographic model

The implementation uses the stdpopsim Homo sapiens demographic model
[`AncientEurasia_9K19`](https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_ancienteurasia_9k19)
and preserves its event ordering and population labels in the selection-aware
SLiM simulation.

The official stable catalog documents the following study-critical values:

| Quantity | Catalog value | Generations at 25 years/generation |
|---|---:|---:|
| Neanderthal-to-Eurasian admixture | 56.8 kya | 2,272 |
| Admixture proportion | 2.96% | not applicable |
| Han/Loschbour split requested by the catalog | 50.4 kya | 2,016 |
| Han/Loschbour split realized at `Q=10` | 50.25 kya | 2,010 |
| Present-day Han diploid Ne | 6,300 | not applicable |
| Human/Neanderthal split | 696 kya | 27,840 |
| Generation time | 25 years | 1 |
| Catalog mutation rate | 1.22e-8 per bp per generation | not applicable |

The catalog names the admixture recipient `Loschbour`. Because the pulse at
2,272 generations predates the Han/Loschbour split at 2,016 generations, that
population label denotes the common non-African ancestral branch at the pulse;
the introgressed lineage can therefore descend into Han.

For comparability with the Gamma-SMC sweep experiment, this study deliberately
configures a mutation rate of **1.25e-8**, rather than silently inheriting the
catalog's 1.22e-8. It uses 1e-8 recombination per bp per generation, a
10,000,000-bp sequence, a focal site at 5,000,000 bp, and 100 present-day Han
diploid samples (200 haplotypes). Resolved stdpopsim, SLiM, and
demographic-model versions and the serialized event table are part of the run
provenance contract.

## Archaic-specific focal allele and selection onset

The focal allele has one fixed origin on the Neanderthal lineage at 2,400
generations ago (60 kya), after the human/Neanderthal split and before the
56.8-kya introgression pulse. It must remain present in the Neanderthal source
from immediately after its origin through the pre-pulse check, and its source
frequency must be at least 0.10 at generation 2,273, exactly one generation
before the requested pulse at generation 2,272. Under the exploratory `Q=10`
scaling used here, stdpopsim rounds both events to the same 2,270-generation
SLiM tick while preserving callback ordering; this pilot therefore does not
retain a literal one-generation gap. The same site is ancestral in
modern-human lineages before gene flow, and recurrent mutation at the focal
site is disallowed. These conditions make the simulated allele
introgression-specific by construction and exclude incomplete lineage sorting
as its route into Han.

At the pulse, the allele becomes standing variation in the recipient branch.
Positive selection begins at that pulse, not at its older mutation-origin time.
Selection acts in the common recipient branch until the realized Han/Loschbour
split tick and continues in Han thereafter; the Han present-day frequency is
the conditioned outcome. stdpopsim 0.3 rounds extended events but truncates
positive catalog demographic times in its generated `time_to_tick()` helper.
The `Q=10` catalog split is therefore realized at 2,010 generations, and both
fitness callbacks are explicitly aligned to that tick so there is no neutral
gap across the population transition. Requested and realized times, the event
table, and the frequency-condition contracts are serialized separately. The
returned present-day tree sequence does not expose exact source frequency,
number of transferred copies, or the full allele-frequency trajectory, so none
of those quantities is claimed as an observed per-attempt result.

## Selection and frequency design

| Component | Pilot specification |
|---|---|
| Simulator | stdpopsim 0.3.0 with SLiM 4.2.2 and tree-sequence recording |
| Sample | 100 present-day Han diploid individuals (200 haplotypes) |
| Selection coefficient | 0.01, 0.005, 0.001 |
| Dominance coefficient | 0.5 |
| SLiM scaling | `Q=10` |
| Forward burn-in | `slim_burn_in=0.1` N generations, followed by stdpopsim recapitation |
| Mutation origin | Neanderthal, 2,400 generations / 60 kya |
| Pre-pulse source condition | Neanderthal AF >=0.10 at generation 2,273 |
| Selection onset | 56.8 kya (2,272 generations), at introgression |
| Han final-frequency intervals | inclusive `[target, 0.95]`, for targets 0.10, 0.20, ..., 0.90 |
| Attempt budget | at most 20 completed external sample-panel draws per design cell; internal conditioning separately timeout-bounded |
| Accepted replicates | one per design cell if all acceptance checks pass |

The nine final-frequency contracts are nested thresholds, not disjoint bins. A
frequency of 0.72 lies in the mathematical acceptance intervals for targets
0.10 through 0.70. Each target nevertheless has its own specification ID and
deterministically derived seed stream; results must not be interpreted as nine
independent frequency bins.

The 0.95 inclusive cap prevents fixation and supports a `ref/ref` comparator.
Both the present-day Han population condition and the observed 100-diploid
sample frequency must be within `[target, 0.95]`. Sample validation further
requires at least two `ref/ref`, one heterozygous, and two `alt/alt`
individuals. A frequency alone therefore does not guarantee acceptance. The
runner uses at most 20 completed external sample-panel draws. A row is written
only after stdpopsim's internally conditioned SLiM call returns, so internal
losses, fixation, and frequency misses are not separate rows. A contract-bound
ledger persists each launched seed and labels an interrupted internal
trajectory; a retry advances to a new deterministic seed. Completed draws
retain sampled frequency, genotype counts, and acceptance, and failure never
relaxes the frequency or genotype requirements.

## Forward-simulation approximation and SLiM version

The runtime uses stdpopsim's `Q=10` rescaling and a short
`slim_burn_in=0.1` N-generation forward burn-in. stdpopsim then performs its
standard recapitation/rescaling step and overlays neutral mutations. This is an
exploratory approximation: SLiM rescaling is not equivalent to the unscaled
selected process, and a 0.1 N burn-in relies heavily on recapitation. Any
biological interpretation requires sensitivity checks against less aggressive
scaling and longer burn-in; none is claimed completed here.

For this demography, `Q=10` reduces the declining Neanderthal population to as
few as eight simulated diploid individuals and can round a short growth phase
to no forward ticks. stdpopsim emits explicit warnings for both effects. They
are major exploratory-pilot limitations, not harmless logging noise.

Catalog demographic and selected-site extended events use different tick
conversion paths in stdpopsim 0.3. The study therefore computes the realized
positive-time demographic tick as `floor(time/Q)*Q` and uses that value for the
Loschbour-to-Han fitness transition. At `Q=1`, requested and realized split
times are both 2,016 generations.

The general repository prefix retains SLiM 5.2 for diagnostics, but unmodified
stdpopsim 0.3.0 scripts use the removed `Subpopulation.genomes` API under SLiM
5.2. The EAS study therefore defaults to a separate, version-checked SLiM 4.2.2
prefix. A compatibility-transformed SLiM 5.2 run is diagnostic only, not a
production study engine or accepted data source.

## Frequency-grid motivation and limits

The 0.10-0.90 grid is an exploratory sensitivity range, not a claimed empirical
distribution of archaic-specific Han SNV frequencies. Two literature examples
show why a broad range is useful but do not define that distribution:

- A Neanderthal-introgressed haplotype at
  [`DPEP1`](https://doi.org/10.1093/molbev/msv176) reached more than 50% in East
  Asians and was reported under positive selection. Its highlighted
  `rs460879-T` allele predates the human/Neanderthal divergence and was
  reintroduced after loss, so that allele does **not** satisfy this study's
  archaic-specific/no-ILS focal-site rule.
- A Han study of
  [adaptively selected cis-regulatory elements](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917166/)
  reported positively selected Neanderthal-derived enhancers near `HYAL1` and
  higher frequencies of two Neanderthal alleles in East Asians than in
  non-East-Asians. It supports the biological plausibility of high-frequency
  regulatory introgression but is not a genome-wide frequency calibration.

Replacing the exploratory grid requires a variant-level dataset with explicit
ancestral/derived orientation, archaic allele state, Han frequency, callability,
and an ILS/recurrent-mutation exclusion rule.

## Gamma-SMC decoding and summaries

Decode each accepted 10-Mb Han sample at these calendar thresholds:

| Cutoff | Generations at 25 years/generation |
|---:|---:|
| 1 kya | 40 |
| 4.5 kya | 180 |
| 10 kya | 400 |
| 20 kya | 800 |
| 30 kya | 1,200 |
| 40 kya | 1,600 |
| 50 kya | 2,000 |

Fix the sample and pair manifests before decoding. Report `P(TMRCA < x)` and
posterior mean TMRCA for all prespecified pairs, `ref/ref` pairs, and `alt/alt`
pairs. Heterozygotes may contribute to the overall summary but are excluded
from the strict homozygote contrast. Acceptance requires at least two sampled
diploids in each homozygous class and at least one heterozygote. Every summary
must carry the number of individuals, haplotypes, and pairs in its class.

Near a beneficial allele, selected-allele haplotypes are expected to have
**higher** `P(TMRCA < x)` and **lower** mean TMRCA than the internal reference
class. A probability-versus-cutoff curve is a TMRCA CDF. Any area under it must
be labeled an integrated CDF over the stated cutoff range, not ROC AUC. ROC AUC
would require labeled selected and neutral replicates and a varying decision
threshold.

`ref/ref` is an internal comparator exposed to the same demography,
introgression pulse, and selected simulation. It is not a matched neutral
replicate. No p-value, false-positive rate, power estimate, or neutral-calibrated
claim will be made in this pilot.

## Produced figures and tables

Figures are emitted in PNG and PDF with large, letter-page-readable text. The
implemented report set is:

1. the AncientEurasia event diagram with mutation origin, pulse, selection onset, and Han split;
2. a complete 3-by-9 status grid showing accepted and timed-out cells;
3. accepted sample frequency and completed-draw effort by selection coefficient and AF threshold;
4. genome-wide `P(TMRCA < x)` and mean-TMRCA heatmaps by genotype class;
5. focal-site CDFs and `alt/alt - ref/ref` probability contrasts;
6. focal-site posterior mean TMRCA contrasts; and
7. tree-sequence truth versus Gamma-SMC estimates where the estimands match.

The returned tree sequence does not expose exact donor AF, migrant-copy count,
or a full selected-allele trajectory, so the figures do not fabricate those
quantities. Machine-readable specifications, completions, manifests, and
checksums retain the full available parameter, sample, decoder, and seed
provenance.

## Reproducibility contract

The implementation separates planning, simulation/acceptance, integrity
checks, Gamma-SMC decode, and reporting into restartable phases. Cache validity
depends on design and sample/pair manifests plus relevant output checksums and
parseability, not merely on a Git commit. Logs preserve resolved defaults,
paths, versions, serialized demographic events, commands, seeds, rejection
attempts, and validation status. Corrupt or incomplete artifacts must never
satisfy a cache check.

When a nonzero per-specification timeout is declared, expiry terminates that
worker and its SLiM descendant and writes a `timed_out` failure record. A
timeout, including a short diagnostic timeout, is failure metadata rather than
simulation data and must not enter acceptance tables, plots, or biological
interpretation. Exhausting the 20-draw acceptance budget is likewise failure,
not an accepted replicate.

## References

- stdpopsim. [`AncientEurasia_9K19` stable catalog entry](https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_homsap_models_ancienteurasia_9k19).
- Kamm et al. [Ancient Eurasian demographic model](https://doi.org/10.1080/01621459.2019.1635482).
- Schweiger and Durbin. [Ultrafast genome-wide inference of pairwise coalescence times](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538485/).
- Ding et al. [Reintroduction of a homocysteine level-associated allele into East Asians by Neanderthal introgression](https://doi.org/10.1093/molbev/msv176).
- Zhang et al. [Adaptive selection of cis-regulatory elements in the Han Chinese](https://pmc.ncbi.nlm.nih.gov/articles/PMC10917166/).
- Haller and Messer. [SLiM](https://messerlab.org/slim/).
