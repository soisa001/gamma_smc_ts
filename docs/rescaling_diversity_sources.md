# Rescaling and introgressed diversity: source review

Reviewed 2026-09-24. This is a primary-source review and analytic check, not a
new simulation or a measurement of the running pilot. No running jobs or
simulation parameters were changed for this review.

**Subsequent decision:** the user then explicitly requested removal of rescaling
and a restart. `configs/origin_onset_pilot.json` now sets Q=1 for both origins.
The restart keeps the donor bottleneck and other biological parameters unchanged;
see `origin_onset_pilot_run.md`. This is not a completed Q-convergence experiment.

## Local model audit

The pilot's `fresh_power.model` uses archaic Ne=3,600 except for Ne=10 over
100 original generations immediately before the 50-kya, 2% pulse. The split
is at 700 kya. These assumptions match the recent EAS h400 array. The older
`run7_models.build_eas_with_archaic` instead uses constant archaic Ne=3,600,
with defaults of a 55-kya, 2.5% pulse in `run7_config.py`. Thus the immediate
baseline and older run7 model are different biological models, even though
both originally used Q=5. Fixing the focal allele in the archaic source does
not itself require making that source nearly clonal.

The imposed bottleneck is approximately 50--52.5 kya; by itself it does not
force coalescence below a 10-kya cutoff. Later recipient drift, ancestry
founder counts, ascertainment, and selection also affect carrier genealogies.
The research calculation below should not be interpreted as measured
present-day carrier nucleotide diversity or as attribution of the pilot's
thresholded scores to the donor bottleneck alone.

## What rescaling should preserve

Rescaling is intended to preserve population-scaled drift, mutation,
recombination, and selection. It does not simply discard four-fifths of the
polymorphisms at Q=5. The diffusion argument keeps products such as Ne*mu,
Ne*r, and Ne*s fixed while shortening elapsed simulation time. This is an
approximation; finite populations, initial allele counts, strong selection,
and long linked regions can violate its assumptions. The recent theoretical
review distinguishes approximately invariant statistics from quantities that
do not scale identically, and specifically discusses recombination over long
chromosomal regions. [Johri, Pouyet and Charlesworth, 2026](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012261)

## Verified stdpopsim implementation

The stable 0.3.0 engine source implements:

- Population sizes divided by Q, then rounded when creating/resizing populations.
- Event ages rounded to multiples of Q generations before scheduling.
- Selection coefficients, growth rates, and continuous migration rates multiplied by Q.
- Recombination transformed as `r_Q = (1 - (1 - 2*r)**Q) / 2`, approximately Q*r
  for small per-link rates. Admixture pulse proportions retain their original values.
- Node, migration, and mutation times multiplied by Q before recapitation.
- Recapitation using the demographic population sizes and original recombination map.
- Neutral mutations overlaid on these restored-time trees using the original
  contig mutation rate, retaining existing mutations; focal event sites are
  excluded from neutral overwriting.
- Warnings for populations below 50 individuals and advice to compare Q values.

Consequently, correctly postprocessed output does not retain branch lengths
that are five times too short. But mutation overlay cannot restore genealogical
branches that were already lost during the forward simulation.
[stdpopsim engine source](https://popsim-consortium.github.io/stdpopsim-docs/stable/_modules/stdpopsim/slim_engine.html)

The local package version and generated scripts should be checked against
these upstream rules when interpreting the pilot.

## Analytic check for a very small donor bottleneck

If a local donor bottleneck uses 10 diploids for 100 generations, Q=5 maps it
to 2 diploids for 20 simulated generations. This is an extreme finite-population
case despite the apparently modest overall Q. It also bounds the number of
physical donor gene copies available in any bottleneck generation to four.
The ancestry represented by those copies is not equivalent to saying the
biological model has only four original introgression founders.

As an **illustrative idealized diploid Wright-Fisher calculation**, ignoring
mutation and migration, the chance that two lineages remain distinct across
a constant-size bottleneck is `(1 - 1/(2*N))**duration`:

| Representation | N | Generations | Two-lineage noncoalescence probability |
|---|---:|---:|---:|
| Q=1 | 10 | 100 | 0.0059205292 (0.5921%) |
| Q=5 | 2 | 20 | 0.0031712119 (0.3171%) |

These numbers are derived analytically here, not empirical pilot results or an
exact calculation for SLiM's full pedigree, recombination, conditioning, or
migration scheme. They demonstrate both points: the intended bottleneck itself
already produces extremely strong coalescence, and additional finite-size
effects are plausible after scaling. Removing rescaling alone does not make
that donor bottleneck mild.

## Implications and proposed validation

Simulation studies find that the size and direction of rescaling bias depend
on the model and statistic; they do not justify assuming that Q=5 always
reduces diversity. A published systematic study recommends assessing the
specific model with unscaled replicates and reports that relatively small
comparison cohorts can reveal substantial biases.
[Dabi and Schrider, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11708920/)

The practical recommendation is to separate two checks:

1. **Model realism:** establish whether the donor bottleneck and focal marker
   construction intentionally make archaic carrier haplotypes nearly
   homogeneous; compare those assumptions with the established model.
2. **Numerical convergence:** compare fixed-seed independent replicate cohorts
   at Q=1 and Q=5 under unchanged biological parameters and ascertainment.
   Examine donor/carrier nucleotide diversity, haplotype diversity, true
   carrier TMRCA distributions, focal AF, tract/LD patterns, and both all-pair
   and ALT/ALT `frac_recent_T`. Include neutral and selected conditions, with
   attention to the largest requested coefficient: s=0.03 becomes 0.15 in
   the Q=5 forward model.

A useful first diagnostic on existing outputs is to compare mutation-based
carrier diversity with branch-based carrier diversity, and inspect the full
TMRCA distribution rather than only thresholded fractions. This separates
mutation-overlay concerns from genealogy concerns. A Q comparison must be
generated explicitly before declaring convergence; different demographic
or ascertainment arms are not substitutes. These checks are recommendations
only and were not executed by this source review.
