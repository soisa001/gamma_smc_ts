# Dating SV origin separately from its route into modern humans

Method sketch, 2026-10-24. Motivated by the user-provided ASHG Board 7021T
abstract and proposed comparisons of introgressed SV carriers and noncarriers.
This is a proposed analysis, not an implemented classifier or empirical result.

## Scientific target

Estimate two distinct properties for each precisely resolved SV allele:

1. When did this structural mutation arise relative to the relevant
   introgression event and the modern-human/archaic population separation?
2. Did the sampled human copies descend through archaic-to-modern gene flow,
   ancestral retention without that gene flow, or another route?

An ancestral polymorphism can both persist through lineage sorting and later
be imported through introgression. Consequently, mutation origin and the
transmission route of a particular copy are not mutually exclusive categories.
An old allele on an introgressed background is not necessarily a mutation that
originated in Neanderthals. Allow unresolved and mixed histories.

## Genealogical basis and the important asymmetry

Measure times backward from the present. Let A be the mutation age and I the
introgression time. For a single-origin, correctly polarized, biallelic SV with
no recurrence, reversion, or genotyping error, at its true local genealogy:

    L = max T(i,j), over distinct derived-carrier pairs
    U = min T(i,j), over derived-carrier/ancestral-noncarrier pairs
    L <= A <= U

L is the sampled carrier-clade root age; U is the age of its attachment to the
nearest sampled noncarrier lineage. Using only a subset of carriers and
noncarriers generally weakens these bounds. At least two carriers and one
noncarrier are needed for this complete bracket; a singleton lacks a
carrier-pair lower bound. Carrier/noncarrier mean differences are descriptive
and do not substitute for these bounds.

This is the concordant/discordant-pair principle underlying
[GEVA (Albers and McVean, 2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6992231/).
GEVA combines pairwise age constraints using a composite construction and
explicitly notes overconfidence in its approximate posteriors. The proposed
extension here uses Gamma-SMC to supply the pairwise distributions and adds
explicit ancestry-route discrimination for SVs.

The following are logical implications under the assumptions above, not
empirically validated calling rules:

- L > I implies the mutation predates introgression. Donor ancestry evidence
  is still required to conclude that these copies were imported.
- U < I implies the mutation is younger than that introgression event.
- L < I < U leaves the mutation on a branch spanning introgression and does
  not identify which side of the event it arose on.

For example, with illustrative ages I=50 kya, L=10 kya, U=150 kya, both a
20-kya mutation on introgressed ancestry and an 80-kya mutation subsequently
imported are compatible with the same sampled tree. These numbers illustrate
identifiability, not demographic estimates.

A neutral old SV can pass through one surviving introgressed founder, giving
a young modern carrier crown. Positive selection can also erase older carrier
lineages from the surviving sample. A long stem may remain, but it does not
locate the mutation along that stem. Multiple surviving donor carrier lineages
would improve the lower bound; their survival is not guaranteed. The branch
inequalities do not require neutrality, although estimating their uncertainty
and interpreting population-level distributions require demographic and
selection calibration.

## Phase 1: establish alleles, ancestry, and observation quality

For each assembly-resolved event, record exact breakpoints, insertion sequence
or orientation/copy-number state, genomic build, ancestral/derived polarization,
phased genotype and its confidence, sample population, and local-ancestry
posterior. Reference versus alternate is not automatically ancestral versus
derived. Treat distinct alleles at a multiallelic locus separately.

Partition present-day haplotypes into introgressed derived carriers (I+),
introgressed ancestral noncarriers (I-), nonintrogressed derived carriers (H+),
and nonintrogressed ancestral noncarriers (H-). Retain ancestry uncertainty
rather than treating tract calls as ground truth. Use local ancestry, not an
individual's overall ancestry or continent, to define these groups.

Prioritize identical, sequence-resolved MEI events for the pilot. Examine
inversions separately because suppressed recombination changes inference and
tract-length expectations. Duplication analyses need orthologous, uniquely
mapped sequence rather than paralog comparisons. Validate apparent archaic SV
absence against ancient coverage and callability; an uncallable genotype is
missing evidence, not an ancestral allele.

## Phase 2: decode SNP-based pairwise genealogies at SV anchors

Use phased SNP observations and callable invariant sequence. Use SV state to
label pairs, not as a synthetic SNP mismatch in the decoder. Apply comparable
callability masks across classes, including missing deletion interiors and
unreliable breakpoints; do not interpret absent calls as confident reference.
Keep a common mutation-rate/time calibration across carrier classes.

The local fork already exposes `--output_positions_file` with exact 0-based
coordinates, `--pairs_file` for explicit haplotype pairs, and `--threads`.
See `src/gamma_smc.cpp`, `src/io.h`, `src/data_processor.h`, and
`python/gamma_smc_aou/decoder.py`. These query coordinates are posterior-output
boundaries, not added mutation observations. Source inspection confirms this
interface; no new execution or decoder validation was performed for this note.

For insertions, query a homologous anchor and nearby callable flanks. For
deletions and inversions, query both breakpoints and suitable flanking anchors;
inspect interior estimates only where the comparison is homologous and valid.
Retain sufficient surrounding sequence for decoding and assess sensitivity to
flank extent and masks. A posterior at an exact coordinate can still obtain
most of its information from nearby sequence. Discordant flank results should
flag recombination or uncertainty, not be automatically averaged away.

The [Gamma-SMC documentation](https://github.com/regevs/gamma_smc) specifies a
Gamma posterior with shape alpha and rate beta, in units of 2Ne generations.
Thus the mean is alpha/beta. For a threshold g in generations, use
F_Gamma(g/(2Ne); alpha,beta), equivalently
F_Gamma(2*mu*g/theta; alpha,beta) when theta=4Ne*mu. Retain distribution
parameters and credible intervals, not just posterior means.

Decode these comparisons:

| Pairs | Purpose |
| --- | --- |
| I+ / I+ | Carrier ancestry within introgressed tracts; younger bound on SV age |
| I+ / I- | Separation of SV states within archaic ancestry; potentially informative older bound |
| I- / I- | Background coalescence and founder/selection context, not a direct mutation-age bound |
| Any derived carrier / any ancestral noncarrier | Stronger global older bound using all reliable groups |
| Carrier pairs including H+ | Older surviving carrier lineages; potential ancestral retention or transfer |

Balance sampling across populations and haplotype clusters. Include diverse
carrier lineages and close carrier/noncarrier neighbors. Retain exact pair
identities. An SV fixed among introgressed haplotypes has no I+ / I- comparison;
do not impute that missing evidence.

## Phase 3: estimate age compatibility with honest uncertainty

Gamma-SMC pairwise marginals do not constitute a joint local-tree posterior.
The [original Gamma-SMC paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538485/)
describes its pairwise approximation; joint tree construction is discussed as
a potential extension. Consequently, independently drawing all pairwise times
and taking their extrema does not produce a valid posterior on L and U.

A practical screening statistic is a GEVA-like compatibility curve:

    C(a) = product_CC F_ij(a)^w_ij
           * product_CN [1 - F_ij(a)]^v_ij

Here CC means derived-carrier pairs, CN means derived/ancestral pairs, and F is
the decoded TMRCA CDF. Use this as a composite compatibility score, not an
ordinary likelihood or an automatically calibrated posterior. Shared
haplotypes and shared sequence create substantial dependence; weights alone
do not guarantee calibrated uncertainty. Do not multiply SNP-derived ancestry
and TMRCA probabilities as independent evidence either.

Calibrate age-class scores and decision thresholds using simulations passed
through the complete ascertainment and decoding pipeline. For a smaller set
of high-priority SVs, use a jointly coherent, dated local-tree/ARG analysis
validated for the demographic setting, then place the structural mutation on
compatible branches while propagating tree and genotype uncertainty.

Treat the introgression time as uncertain and potentially lineage-specific
under multiple pulses. If using an externally supplied event interval,
evaluate age support on each side and retain the overlapping mass as ambiguous.
Any probabilities derived from a fitted classifier or age model need held-out
calibration and explicit prior assumptions.

## Phase 4: distinguish transmission through introgression from ILS

Age alone does not determine route. Combine the age evidence with local
modern/archaic topology, tract extent in genetic distance, linked archaic SNPs,
and validated archaic SV genotype/sequence when available. Joint modern and
archaic genealogies have been used to study introgression and ILS, for example
in [Schaefer et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8284891/).

Under a simple isolation model without later exchange, an ILS modern/archaic
coalescence must predate their population separation. A younger coalescence
therefore argues against that ILS-only model, subject to dating uncertainty
and demographic misspecification. An older coalescence does not exclude
introgression: the introgressing donor and the sequenced archaic can themselves
have a deep common ancestor. Tract length is supporting evidence conditional
on local recombination and selection; a long inversion is not a generic
introgression certificate.

The inspected Gamma-SMC interface does not expose an ancient-sample age model.
Use a method validated for serially sampled genomes for modern/archaic dating,
or implement and validate that extension separately. For ancient age s and
MRCA age t, total mutation-bearing length is 2t-s, not 2t, and t must be at
least s. Merely correcting a posterior mean does not repair the full HMM.

Report separate evidence for these histories:

| History | Interpretation |
| --- | --- |
| SV predates gene flow; carrier lineage imported | Supports an introgressed SV, whether donor-derived or older ancestral variation |
| SV younger than gene flow; carrier background imported | Supports a new mutation on introgressed ancestry |
| Ancestral SV retained without focal archaic-to-modern transmission | Consistent with ILS/ancestral retention |
| Ancestral SV with both retained and imported copies | Mixed history; allele-wide binary label is inadequate |
| Modern standing SV transferred onto archaic ancestry, recurrent SV, or uncertain genealogy | Alternative or unresolved history |

## Phase 5: establish what is identifiable before large-scale application

Simulate the same demographic background while tracking the true mutation
branch and focal lineage's migration route. Include: old imported SVs with one
and multiple surviving founders; new post-introgression mutations; ancestral
polymorphisms retained with and without subsequent importation; standing modern
SVs transferred onto an introgressed background; and recurrent/error cases.
Use selection-aware simulation for selection scenarios rather than labeling
neutral simulations as selected after the fact.

Vary mutation age around the admixture interval, donor frequency, present-day
carrier frequency, donor structure, pulse number, recombination, and selection
onset/strength. Model inversion recombination separately. Apply realistic
phasing, missingness, SV genotyping, and local-ancestry errors. Match observed
ascertainment: surviving SVs selected for overlap with called introgressed
tracts. Include an oracle-truth stage to separate intrinsic nonidentifiability
from decoder and genotype error.

Use independent simulation replicates for calibration and evaluation. Do not
treat pairs from the same locus as independent replicates. Report confusion
rates, calibration, and the fraction left unresolved, stratified by allele
count, SV class, and number of surviving donor lineages. Compare against
carrier-mean TMRCA, tract/LD evidence alone, and age evidence alone to quantify
the gain from each component.

Implementation should remain phased and resumable: manifests/QC, targeted
decoding, scoring, calibration, then candidate review. Default to four threads
and a 26-GB memory budget, preserve existing RNG seeds (or fixed iteration
seeds 0,1,2,... for new replicates), and log every parameter including defaults.
Reuse outputs based on sample/pair/query manifests, relevant inputs, masks,
parameters and verified output correctness, rather than invalidating solely
because a repository commit changed. Use uv for Python tooling. Future plots
should be written as PNG and editable-vector PDF with large text, inspected
locally, and not displayed in notebook or chat output.

Per-SV output should include allele identity and polarization, group counts,
callability, pair coverage, age support/intervals, introgression-time assumptions,
ancestry-route evidence, archaic genotype callability, flank consistency,
recurrence flags, calibration stratum, and the unresolved probability/category.

## Initial deliverable

Start with the sequence-resolved MEI candidates from the supplied abstract,
including both the 79 candidate introgressed events and the 405 older candidates
as unvalidated candidates, not truth labels. Run a targeted simulation pilot
first. The decisive evaluation is whether I+ / I- information distinguishes
new mutations from old SVs whose surviving carrier crowns are young. An
honest unresolved category for stems spanning admixture is a required outcome.
