# Proposed origin, onset, and TMRCA grid

The active design is now [the segregating introgressed-allele experiment](segregating_introgression_run.md):
100 targets per cell, natural focal choice at 50 or 10 kya, and de novo deferred.
The earlier grid below is retained as history.

Pilot authorized 2026-09-24: start with 10 test targets per cell and 1,000
calibration nulls per matching specification. The initial design below is now
implemented in `configs/origin_onset_pilot.json`; simulation, audit, decoding,
and focal analysis are separate resumable phases. Existing trajectories are
not replayed. The user confirmed that de novo means one new EAS ALT copy at
selection onset.
The user confirmed the final positive coefficient is 0.03 and specified
10 independent targets versus 1,000 calibration nulls. For introgressed
targets, both calibration nulls and unselected test targets use the same
unselected archaic-allele specification; only their replicate identities differ.

## Review of the saved array

- `configs/eas_q02_50k_h400_array.json`: EAS PHLASH demography; nominal 2%
  introgression at 50 kya; selection starts at that pulse; s=0.001 through
  0.010 in steps of 0.001; 100 observed surviving selected replicates per s;
  1,000 unconditional neutral regions; 400 haplotypes and 10,000 sampled pairs.
- `sim/eas_q02_h400/simulation_status.json` outside the repository reports
  2,000/2,000 complete, audit phase, no failures, updated 2026-09-18. The older
  array overview's 1,947/2,000 paused statement is historical, not current.
  This review read that receipt; it did not rehash the full simulation archive.
- Saved figure packs document decoded cohorts at s=0.002 and s=0.005, each
  with 100 selected regions and 1,000 neutrals. At T=50 kya and p<=0.05 their
  all-pair detection is respectively 5/100 and 61/100. The original lab-meeting
  overview predates the s=0.002 addition. Full-grid decoded availability must
  be inventoried from receipts before future work is scheduled.
- Existing cutoffs are 5, 10, 20, 30, 40, 50 kya. The array config already has
  alpha=0.05; figure scripts also produced 0.01 and 0.001 sensitivity panels.
- Existing selected simulations condition on survival and sample observation;
  fixation is retained. Neutral simulations anchor to an ordinary segregating
  variant, and subsequent archaic-marker analyses use their own ascertainment.
  They are not automatically matched neutral focal alleles of a specified age.
- Current primary error-rate panels report positional FPR, not discovery FDR.
- The score called decoded `frac_recent_T` thresholds each pair's posterior
  mean TMRCA. It is not the mean posterior probability P(TMRCA<T | data).
- Spatial all-pair results exist for s=0.002, but that arm's saved archaic
  labels cover only the selected allele. This suffices to define pairs by the
  focal allele; it does not supply local archaic marker labels everywhere.

Source files: `python/gamma_smc_aou/fresh_power.py`,
`python/gamma_smc_aou/pair_class_profiles.py`,
`scripts/decode_saved_s002.py`, `scripts/lab_meeting_figures.py`, and
`docs/results/lab_meeting_s00{2,5}_all_pairs_20260921/`.

## Biological grid

| Treatment | Focal allele origin | Selection starts | Frequency at selection onset |
|---|---|---:|---|
| I50 | Archaic allele enters in the nominal 2% pulse at 50 kya | 50 kya | Approximately 2%, stochastic pulse realization |
| I10 | Same introgressed allele and pulse at 50 kya | 10 kya | Outcome of 40 kyr of neutral drift; do not reset to 2% |
| D50 | One new EAS copy at 50 kya | 50 kya | 1/(2 N_EAS(50 kya)) in the simulated diploid population |
| D10 | One new EAS copy at 10 kya | 10 kya | 1/(2 N_EAS(10 kya)) in the simulated diploid population |

Confirmed common positive s grid:

`0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.020, 0.030`.

Use the same 12 coefficients in all four treatments to make comparisons direct.
The previously proposed 0.015, 0.040, and 0.050 are removed; 0.300 was a typo.
Include s=0 matched controls. Keep h=0.5 and the same fitness convention.
With 25-year generations, 50 kya is 2,000 generations and 10 kya is 400.
Fivefold higher s is a rough duration-based comparison, not a promise of equal
power or equal final AF, especially when initial copy numbers differ.

Working demographic assumption: keep the existing EAS demographic background,
including background archaic admixture, in all arms; D50/D10 change the focal
allele's origin to an EAS mutation instead of an archaic allele. The focal copy
is placed on a randomly sampled EAS haplotype without an ancestry filter.
Removing the background pulse entirely would be a separate demographic
comparison, not an implicit change bundled with the focal-origin factor.

For this initial pilot, use 10 observed surviving selected replicates per cell:
4 x 12 x 10 = 480. Replicate IDs and seeds are stable if expanded later.
Retain fixation and do not impose a terminal AF target. This is explicitly a
conditional detection experiment, not power per mutation introduced.

Three matched neutral focal-allele families are retained: introgressed at 50 kya,
de novo at 50 kya, and de novo at 10 kya. For EACH family, use 1,000 observed
calibration nulls plus 10 additional independent observed s=0 test targets.
The introgressed specification is explicitly confirmed by the user; the same
matched-origin design carries forward to the previously requested de novo arms.
I50 and I10 share the same s=0 family because there is no selection onset under
neutrality. Share its calibration cohort and independent test cohort across s
and the two introgressed onset arms; do not regenerate identical controls for
each heatmap cell. Repeated use is not independent replication.

| Test | Targets | Reference for empirical p values |
|---|---|---|
| Selection detection | 10 selected focal alleles per positive-s cell | 1,000 matching unselected focal alleles |
| Neutral detection/error | 10 additional unselected focal alleles per null family | The same 1,000 matching unselected focal alleles |

For the introgressed comparison, every unselected focal allele in both rows is
archaic and introduced by the same pulse. There is no ordinary non-archaic or
no-introgression comparator in this updated design. Preserve matching simulation
specifications, ascertainment, demographic history, pulse, and allele definition;
use disjoint seeds and replicate identities for neutral targets and nulls.
Match final observation ascertainment across selected and unselected cohorts;
record allele age, origin, population, and initial copies rather than substituting
a nearby unrelated SNP. Survival/observation conditioning remains the draft
default, not a claim of detection per original introduction.

Nominal retained total: 3,510 (480 selected + 3,000 calibration nulls + 30
independent neutral test targets). Following the authorized Q=1 restart on
2026-09-24, existing Q=5 I50 trees cannot supply this cohort. Only verified
compatible Q=1 regions may be retained. This is not an attempted-simulation budget:
rare neutral de novo survivors can require many attempts. Reuse of individual
neutral artifacts requires proof of matching focal identity and ascertainment;
the existing 1,000 unconditional neutral regions do not satisfy that by default.
Keep those existing neutrals only as a separately labeled historical reference
unless their compatibility with the new focal-allele specification is verified.

To test the specific hypothesis that a neutral de novo introduction rarely
reaches high AF, also plan an unbiased fixed-attempt trajectory cohort with
losses retained. Its size and cost should be set after an authorized feasibility
phase, before examining comparative power. Record AF=0 for loss and distinguish
population survival from observation in the 400-haplotype sample. Do not infer
unconditional survival probabilities from accepted trees or from outer retry
counts: the current SLiM survival-conditioning machinery also rejects internally.
The 3,510 retained-tree count excludes this separately budgeted trajectory work.
This optional attempted-origin study is not needed for the requested conditional
focal comparison and is not included in its launch scope.

## Analysis grid (reuses each biological replicate)

| Dimension | Values |
|---|---|
| TMRCA threshold T | 10, 20, 30, 40, 50 kya |
| Test cutoff | p<=0.05 only |
| Pair population | All 10,000 manifest pairs; ALT/ALT subset defined by the focal allele |
| Score sources, reported separately | True pairwise TMRCA; Gamma-SMC decoded posterior-mean hard calls |
| Source-specific calibration | Truth against truth nulls; decoded against decoded nulls |
| Primary endpoint | Prespecified focal position, matched-position null |
| Spatial diagnostic | Same pair identities across positions; display +/-500 kb, save full 10 Mb |
| AF baseline | Focal ALT frequency, with its own matched neutral calibration |

Thus there are 4 x 12 x 5 x 2 = 480 selected power cells per source, or 960
across the two separately reported sources. There are also 3 x 5 x 2 = 30
neutral-detection cells per source. These are repeated analyses of the same
simulations, not independent experiments or simulation jobs. AF has 48 selected
cells and is not replicated artificially over T.

For pair set C at coordinate z:

`R_C(T,z) = count{pairs in C with TMRCA(z)<T} / count{pairs in C}`.

Use C=all pairs and C=ALT/ALT pairs. The latter is the within-carrier fraction,
without multiplying by carrier-pair frequency, AF, or AF squared; no REF/REF
subtraction, ALT/REF penalty, carrier-mass gate, AF matching, or minimum positive
pair-count filter is part of the primary score. Mathematically it conditions on
both haplotypes carrying ALT, but it does not weight by how common ALT is.

At zero ALT/ALT pairs, raw R_ALT_ALT is undefined: export NA, never invent a
zero raw fraction. Mark the test unevaluable and count it as no discovery in
the enrolled-cohort detection/FPR calculation (equivalently p=1 for that rule).
Also report evaluable fraction and detection among evaluable replicates. A
single sampled carrier yields no pair; even two carriers can have their one
pair absent from the fixed 10,000-pair manifest. Keep pair counts in all tables.
This availability effect is especially relevant to de novo arms. Do not pool
carrier pairs across replicates to estimate power or representative distributions.

Keep the same existing pair manifest (seed 1729) as the primary comparison.
At low AF, targeted extra carrier pairs would change the sampling design;
such a sensitivity analysis must be declared separately and applied to nulls
and selected arms alike, rather than silently repairing missing pairs.

For decoded R, substitute each pair's decoded posterior-mean hard call for
true TMRCA. Preserve exact numerical conventions in reused results. Existing
native extraction uses a >= boundary comparison while true TMRCA uses strict
<; record and check boundary cases rather than silently changing cached values.
If posterior CDF averages are wanted later, define a separate score and null.

## Scientific predictions and limits

De novo starting frequency is much smaller than 2%. Under a matched neutral
process without survival ascertainment this should reduce the chance of high
terminal AF. It does not guarantee more power: selected copies also start
rarer, may be lost, and have less time to rise in D10. Among surviving observed
neutral alleles the AF distribution is different from that among all origins.
Show both populations of trajectories when making this claim.

There is a particularly sharp prediction for raw carrier recency. With one
mutation origin and no recurrent/back mutation at the focal site, all present
carriers descend from that copy. Their focal-site TMRCA is no older than its
origin. Therefore for T strictly greater than the birth time, true ALT/ALT
R(T) is 1 whenever a carrier pair is available, under both s=0 and s>0.
D10 must therefore saturate at focal T=20,30,40,50 kya; the 10-kya boundary
depends on exact event timing/ties. D50 has the corresponding boundary near
50 kya. This is a deduction from the proposed model, not a simulated result.
Decoded scores may depart from this due to inference error. The raw fraction
can still differ spatially away from the mutation, where recombination matters.
Do not add sub-10-kya cutoffs silently; flag them as an optional future extension
if focal D10 ALT/ALT discrimination becomes an objective.

Initial frequency 1/(2N) is also the convention in the
[msprime sweep API](https://tskit.dev/msprime/docs/stable/api.html#msprime.SweepGenicSelection).
The implementation must use the actual simulated population size at the event,
not the sample size of 200 diploids. The user authorized removal of rescaling
on 2026-09-24: all arms now use Q=1. One copy therefore means one original-size
EAS genome. Q=5 archived introgressed trees cannot be reused in this cohort.
The archaic bottleneck remains 10 diploids for 100 generations; its biological
effect on diversity is separate from numerical rescaling.

## Calibration, detection, and false discoveries

For each treatment/null family, score, source, and T, rank against matching
neutral calibration scores, using the conservative upper-tail empirical p:

`p = (1 + number of calibration scores >= test score) / (n_calibration + 1)`.

Use all 1,000 calibration nulls to rank each of the 10 selected targets and
each of the 10 independent neutral test targets. There is no 800/200 split
and no cross-validation fold partition in this revised design. Never include
a test target in its own reference cohort. Keep biological/ascertainment
parameters identical between neutral targets and calibration nulls, and check
that input identities and seeds are disjoint. Reuse the same cohort assignments
for truth and decoded scores, all pairs and ALT/ALT, and every T.

The existing figure workflow used 400 calibration neutrals per fold; regenerate
legacy comparisons under the new 1,000-null rule before interpreting changes as
biological. Minimum attainable p is 1/1001. Raw-score differences determine
separate thresholds; p<=0.05 stays fixed. Retain ties conservatively. Missing
ALT/ALT scores remain NA in exports; for ranks only, an unavailable reference
score is a no-call outcome ordered below every finite score, retaining the
1,000-null denominator. An unavailable target always has p=1. Report availability
for nulls and targets and evaluable-only summaries alongside enrolled-cohort rates.

- Detection/power: TP/N_selected, with NA pair scores treated as no discovery
  in the enrolled conditional cohort. Also give the evaluable-only denominator.
- Neutral detection/FPR: FP/10 among the independent unselected test targets,
  using the same ascertainment and missing-pair rule as the selected targets.
  Report the count as well as the percentage and uncertainty interval.

Terminology: the user calls this neutral-only error experiment FDR. Its requested
fraction of unselected targets detected is conventionally FPR, so plot it as
"Neutral detection rate (FPR)" and explain the correspondence in captions.
The 10 targets : 1,000 nulls ratio specifies testing and calibration cohort
sizes, not disease/selection prevalence among tested loci. The earlier assumed
pi1=1/11 mixed-discovery FDR calculation is removed. It must not be inferred
from this calibration ratio. With only null test targets, any discoveries are
all false: FDP=1 if there is a discovery (and conventionally 0 if there is none);
the formal batch FDR is then the probability of at least one discovery, not FP/10.
No separate mixed-prevalence FDR estimate is requested by this clarification.

Replicate-level bootstrap of the full calibration/evaluation procedure can
include null-fit uncertainty; pair-level resampling cannot replace independent
simulations. With 10 targets, both power and FPR change in 10-percentage-point steps. This pilot cannot precisely validate a nominal 5% FPR; report counts and wide intervals.

Treat all five T values as prespecified separate comparisons. Do not take the
best T per replicate or call a replicate positive if any cell has p<=0.05
without calibrating that combined rule. Any future regional/gene maximum must
likewise be calibrated against the same window summary in neutral regions.

## Figure plan

1. Design table showing origin, mutation/pulse time, selection onset, initial
   copies, s, and cohort counts, including eligible and ineligible ALT/ALT tests.
2. Two separate decoded detection heatmap figures: all pairs and raw ALT/ALT.
   Each has four treatment panels, rows=s and columns=T. Fixed scale 0-100%.
   Repeat as separate true-TMRCA results, without mixing score sources.
3. Held-out FPR panels with one common 0-10% scale across all treatments and
   methods. If any value exceeds 10%, expand the common scale globally; never
   clip it. A shared neutral arm should be labeled as shared, not independent
   evidence repeated across s. Show evaluable-only FPR as a diagnostic.
4. The user's neutral-only error comparison is shown by figure 3, not a second
   mixed-prevalence FDR heatmap. Include raw detected counts out of 10 and
   confidence intervals alongside that shared-scale neutral-detection figure.
5. AF distributions versus s, with matched s=0 controls, fixation fractions,
   common 0-100% AF axes, and AF-only detection curves at p<=0.05. Separately
   plot AF vs raw R_ALT_ALT, annotated with pair counts, to inspect redundancy
   and saturation without weighting the score by AF.
6. AF trajectories for new simulations: median and central intervals, loss
   and fixation fractions, origin/onset annotations. Keep attempted-origin
   and surviving-observed cohorts separate. Save trajectories during future
   generation. Historical trajectories were not saved for the existing archive;
   do not imply they exist or replay simulations without authorization.
7. Raw R(T) distributions at the focal position for every treatment and s,
   overlaying matched neutrals, with fixed 0-1 y scale and all five cutoffs.
   Weight replicates equally, not by their numbers of carrier pairs.
8. Representative spatial raw profiles for every treatment and s: select
   neutral and selected replicates closest to their cohort's median final AF,
   break ties by replicate ID, never choose on significance. Plot all pairs and
   focal-ALT/ALT pairs separately, with all five T curves and fixed 0-1 axes.
   Export IDs, AF, pair counts, truth and decoded profiles, and selection rule.
   If the chosen example has no carrier pair, show that absence explicitly;
   an additional evaluable example may be labeled as such, never substituted
   silently. Display +/-500 kb and retain full-region profile tables.

Save all figures to disk as PNG and vector PDF, with large letter-page-readable
text, horizontal legends above panels, editable PDF text, and no rasterized
data/colorbars. Set MPLBACKEND=Agg before importing pyplot. Visually inspect
saved renders without displaying figures inline or emitting notebook plots.

## Future implementation phases and checks

1. The positive-s grid and independent 10-target/1,000-null pilot design are now
   confirmed. Retain the documented demographic-background and survival/observation
   defaults. The accompanying CSV expands biological cells and separates neutral
   target and calibration roles; it is not a runnable simulation configuration.
   The user has now explicitly authorized starting this pilot.
2. Implement/audit focal origin, matched nulls,
   event timing, single-copy semantics, survival bookkeeping, retained
   trajectories, and pair-availability behavior. Use focused synthetic tests
   and a small real end-to-end validation before scale. The restarted pilot uses Q=1 throughout.
3. Generate only missing compatible inputs, audit saved artifacts, then decode
   only missing compatible profiles. Do not combine simulation and decoding
   implicitly. Keep originals immutable and use separate output roots per
   scientific treatment.
4. Extract raw scores, calibrate, produce plots, and verify labels, scales,
   counts, vector outputs, and manifests.

Keep study seed 20380101 and pair seed 1729. Extend new task seed identities
with origin, introduction time, onset, s, cohort role, replicate and attempt; preserve all
legacy seeds exactly. Include actual model/engine/scaling settings in cache
identity, not just task labels. Cache validity depends on sample/pair manifests,
scientific inputs, relevant output hashes and semantic correctness, not raw
Git commit version. Log paths, all defaults, seeds, counts, losses, retries,
versions, parameter hashes and reuse decisions. Verify readable trees, sample
and position alignment, focal genotype/carrier equality, checksums, finite
monotone scores, pair denominators, and separation of calibration and test cohorts.

New resource defaults: 4 total worker threads, 26 GB RAM. Do not inherit the
old launcher's 20-worker/300-GB settings. The authorized launcher is `scripts/launch_origin_onset.sh`; see `docs/origin_onset_pilot_run.md` for setup, resume, and status commands.
