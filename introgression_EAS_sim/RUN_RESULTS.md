# Introgressed EAS sweep: bounded local run

## Execution outcome

This report freezes the exploratory local campaign completed on 2026-08-11.
It is one accepted trajectory per reachable design cell, not a replicate-based
power study. Of the 27 requested selection-by-frequency specifications, 13
produced an accepted 100-diploid Han panel and all 13 were decoded with
Gamma-SMC. The other 14 cells remain explicit `timed_out` / `not_simulated`
rows rather than being omitted or relaxed.

| Selection coefficient | Accepted AF-floor cells | Timed-out AF-floor cells |
|---:|---|---|
| 0.01 | 10%, 40% | 20%, 30%, 50%, 60%, 70%, 80%, 90% |
| 0.005 | 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80% | 90% |
| 0.001 | 10%, 20%, 40% | 30%, 50%, 60%, 70%, 80%, 90% |

Accepted sampled selected-allele frequencies ranged from 0.115 to 0.91. The
intervals are nested `[floor, 0.95]` contracts, not disjoint AF bins. The exact
completed sample-panel draws, interrupted internal trajectories, final
invocation time, last completed sample AF, and failure reason are in
`results/simulation_status.tsv`; decoder coverage is in
`results/decode_status.tsv`.

The corrected campaign used an initial bounded pass followed by one
deterministic five-minute retry for each of the 14 unfilled cells. That retry
added no accepted panels, so execution stopped at the declared pilot boundary.
A deterministic trajectory screen flagged all `s=0.001` cells and
`s=0.005, AF>=0.90` as difficult under the fixed archaic-source contract, but
it was diagnostic rather than an exclusion gate: three `s=0.001` cells did
accept.

## Corrected population-transition timing

The AncientEurasia catalog requests the Han/Loschbour split at 2,016
generations. stdpopsim 0.3's positive-time demographic conversion realizes
that split at 2,010 generations when `Q=10`. Both the demographic transition
and the selected-site fitness callback now use the realized 2,010-generation
tick. The final data therefore contain no unintended ten-generation neutral
gap between selection in the ancestral recipient and selection in Han.

The archaic-only focal origin remains at 2,400 generations, the requested
source check is at 2,273 generations, and introgression plus selection onset
are requested at 2,272 generations. The latter two events collide at the same
2,270-generation scaled tick while retaining callback order, so this pilot
does not preserve their literal one-generation separation.

## Descriptive focal-site results

The expected sweep direction is positive `alt/alt - ref/ref` normalized
integrated TMRCA-CDF area and negative mean-TMRCA difference. At the focal
site, both signs occurred together in 10/13 accepted simulations using
tree-sequence truth and in 10/13 using Gamma-SMC. Considered separately, the
CDF-area sign agreed in 10/13 and the mean-TMRCA sign in 12/13 for both
sources.

The same three cells had non-positive CDF-area contrasts in truth and
Gamma-SMC: `s=0.005` with AF floors 10%, 30%, and 70%. Only the 10%-floor cell
also had a positive, rather than sweep-direction, mean-TMRCA contrast.

For the prespecified representative cell `intro_s0p005_af50`:

| Source | Normalized integrated CDF area, alt-ref | Mean TMRCA, alt-ref |
|---|---:|---:|
| Tree-sequence truth | +0.129537 | -7,746.41 generations (-193,660.19 years) |
| Gamma-SMC | +0.247507 | -23,691.10 generations (-592,277.59 years) |

These are internal genotype-class contrasts from a selected simulation. The
area is the normalized integral of the TMRCA CDF over the requested 1-50-kya
cutoffs; it is **not** ROC AUC. There is no matched neutral null in this pilot,
so these values do not define a p-value, false-positive rate, or detection
power.

## Runtime and integrity scope

The accepted simulations were generated on Windows with Python 3.11.15,
stdpopsim 0.3.0, tskit 1.0.3, and SLiM 4.2.2. Gamma-SMC decoding ran under WSL
with the existing Linux AVX2 binary whose SHA-256 is
`e57182c105e216c5a3c4afad991d7f51343dbf21d4c1889b38af97a3c52e6d40`.
The bootstrap rebuilds and verifies a decoder for a fresh checkout; this local
result does not claim that the checked binary was freshly compiled during this
campaign.

The final result manifest contains 13 plotted specifications, 13 complete
decoder rows, 14 timeout rows, 23 artifact-relative paths, and SHA-256
checksums. All declared artifacts were re-read and checksum-verified after
aggregation; compressed result tables passed gzip CRC validation. Bulk tree
sequences and decoder intermediates remain restartable under the ignored
`work/` directory.

## Pilot limitations

- `Q=10` and `slim_burn_in=0.1` are speed-oriented approximations. In this
  demography, rescaling reduces the declining Neanderthal population to roughly
  8-20 simulated diploids and collapses a short growth phase to zero forward
  ticks.
- The serialized forward-event contract establishes an archaic-only focal
  origin and no recurrence. Mutation-node population metadata is not used as
  proof of origin, and the returned tree sequence cannot recover exact donor
  AF, migrant-copy count, or the full selected trajectory.
- The 10%-90% grid is exploratory. No qualifying genome-wide distribution of
  archaic-specific, no-ILS Han SNV frequencies was identified for replacing it
  in this run.
- One trajectory per accepted cell in this pilot is suitable for method
  inspection only. A completed exploratory extension requested ten seeded slots
  per cell: 270 requested, 155 accepted and decoded across 24/27 represented
  cells (22 with at least two paired trajectories for the macro headline), 104
  exhausted, 11 timed out, and 0 failed. See
  [`replicate_study_n10/RUN_RESULTS.md`](replicate_study_n10/RUN_RESULTS.md).
  Those replicate summaries remain conditional on accepted trajectories;
  less-aggressive scaling, longer burn-in, matched neutral simulations, and
  null-calibrated p-values remain follow-up work.

The machine-readable source of truth is
`results/results_manifest.json`; methods and the reusable run contract are in
[`METHODS.md`](METHODS.md).
