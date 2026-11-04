# Three EAS selection episodes

## Completed experiment and interpretation

The new simulation/inference phase used 24 workers for 48.109 minutes. It completed 5,879 new proposals: 5,541 selected proposals (131–132 per treatment) and 338 neutral regions. Combined with 3,505 verified prior proposals, all 9,384 submitted units finished without failures or incomplete work. The analysis contains 493 surviving, present-day segregating selected regions and 676 independent neutral regions across four matched pools. Each pool has 20 development, 99 calibration, and 50 held-out regions. No treatment reached the target of 100 accepted survivors; all power estimates below are conditional on survival and use the actual denominator.

Posterior mass improves detection over posterior medians, particularly for introgression. For 10-kya onset, mass detects 36/56 at s=0.02 and 56/59 at s=0.03, versus 13/56 and 28/59 using medians. Effects at s≤0.01 remain near background detection rates. For selection ongoing since 50 kya, mass detects 9/33 at s=0.005 and 28/40 at s=0.01. That episode spans the entire scan and therefore cannot test outside-time false calls.

For the 40–20 kya episode, narrow-bin mass detects 0/26 at s=0.005, 0/25 at s=0.01, and 10/15 at s=0.02 in introgression simulations. At s=0.02, five of 15 regions have near-site calls entirely outside the true period, and eight of 15 have at least one outside-time component. Mean temporal F1 at the exact selected coordinate is only 0.0646, despite region F1 of 0.7407. Wider 1/2/5/10-kyr windows improve joint detection to 15/15 and temporal F1 to 0.1724; five of 15 still contain an outside-time component. Thus overlap-based detection is appreciably easier than accurate temporal localization. The wider scan detects only 1/25 of the s=0.01 introgression episodes.

Held-out posterior-mass FPR is 2/50 (4%; Wilson 95% interval 1.1–13.5%) for the shared introgression pool and the recent/50-kya de novo pools, and 0/50 (upper 95% bound 7.1%) for the 40–20 kya de novo pool. The primary and wider scans have the same observed FPR here, but their union is not calibrated. These are per-1-Mb-region error rates, not genome-wide family-wise error control.

The exact all-19,900-pair true-tree audit identifies a subsampling limitation. At recent de novo s=0.03, narrow-bin true-tree detection rises from 13/26 with 1,000 pairs to 18/26 with all pairs. At recent introgression s=0.02 it rises from 3/56 to 22/56; wider windows with all true pairs detect 46/56. This is a diagnostic benchmark, not an optimal-detector upper bound. Posterior inference in this pilot still uses 1,000 unfiltered pairs. Full audit results are in `all_pairs_truth_performance.tsv` and `all_pairs_truth_replicates.tsv`.

The absence of accepted older strong-selection simulations is a conditioning result, not zero estimated detection power. For example, among 132 introgression proposals with selection ongoing since 50 kya at s=0.05, 129 fixed and three were lost. None remained segregating. The rejected proposals and exact accepted counts are retained below and in `proposal_summary.tsv`.

Figures are saved as PNG and vector PDF: `power_de_novo`, `power_introgression`, `time_profiles_de_novo`, `time_profiles_introgression`, and `older_episode_example`. The example was chosen after observing outcomes and is illustrative. All tests and execution-source checks are recorded in `validation_checks.json`; the exact executed source and current CLI source are preserved separately.

Selection 10 kya to present, 50 kya to present, and 40–20 kya; s = 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05; target 100 surviving segregating replicates per treatment (42 treatments).
EAS PHLASH history, unscaled WF, h=0.5, 25 years/generation, 1 Mb, 100 diploids. De novo starts at one copy at onset. Introgression injects 1% donor haplotypes at 50 kya, with a fixed donor allele, human–Neanderthal split 500 kya and donor Ne 3,600. Selected loss/fixation rejected; neutral regions unconditioned with nearest segregating focal label.
Detection scans all 49 adjacent one-kyr intervals from 1 to 50 kya, every 10 kb. Development nulls standardize interval coalescence fractions. Independent calibration maxima control the whole-region scan at nominal 5% per estimator; held-out nulls estimate FPR. All three estimators use the same 1,000 unfiltered pairs (of 19,900).
Joint TPR requires a call within 100 kb of the selected site and overlapping the true selection interval. Spatial TPR ignores time. Temporal F1 is the Dice overlap at the exact focal coordinate. A wholly outside temporal call is reported separately from held-out neutral FPR. Region F1 depends on the selected/null mixture.

Completed proposals: 9384; selected regions analyzed: 493.

| Window | Origin | s | Accepted / target | Prefix proposals | Lost | Fixed |
|---|---|---:|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 2/100 | 132 | 130 | 0 |
| 10 kya → present | de_novo | 0.005 | 13/100 | 528 | 515 | 0 |
| 10 kya → present | de_novo | 0.01 | 15/100 | 528 | 513 | 0 |
| 10 kya → present | de_novo | 0.02 | 20/100 | 528 | 508 | 0 |
| 10 kya → present | de_novo | 0.03 | 26/100 | 528 | 502 | 0 |
| 10 kya → present | de_novo | 0.04 | 10/100 | 132 | 122 | 0 |
| 10 kya → present | de_novo | 0.05 | 7/100 | 132 | 125 | 0 |
| 10 kya → present | introgression | 0.001 | 8/100 | 132 | 124 | 0 |
| 10 kya → present | introgression | 0.005 | 38/100 | 528 | 490 | 0 |
| 10 kya → present | introgression | 0.01 | 43/100 | 528 | 485 | 0 |
| 10 kya → present | introgression | 0.02 | 56/100 | 528 | 472 | 0 |
| 10 kya → present | introgression | 0.03 | 59/100 | 527 | 468 | 0 |
| 10 kya → present | introgression | 0.04 | 14/100 | 132 | 117 | 1 |
| 10 kya → present | introgression | 0.05 | 13/100 | 132 | 117 | 2 |
| 50 kya → present | de_novo | 0.001 | 0/100 | 132 | 132 | 0 |
| 50 kya → present | de_novo | 0.005 | 0/100 | 132 | 132 | 0 |
| 50 kya → present | de_novo | 0.01 | 1/100 | 132 | 131 | 0 |
| 50 kya → present | de_novo | 0.02 | 0/100 | 132 | 131 | 1 |
| 50 kya → present | de_novo | 0.03 | 0/100 | 132 | 128 | 4 |
| 50 kya → present | de_novo | 0.04 | 0/100 | 132 | 130 | 2 |
| 50 kya → present | de_novo | 0.05 | 0/100 | 132 | 125 | 7 |
| 50 kya → present | introgression | 0.001 | 13/100 | 132 | 119 | 0 |
| 50 kya → present | introgression | 0.005 | 33/100 | 132 | 99 | 0 |
| 50 kya → present | introgression | 0.01 | 40/100 | 132 | 69 | 23 |
| 50 kya → present | introgression | 0.02 | 0/100 | 132 | 36 | 96 |
| 50 kya → present | introgression | 0.03 | 0/100 | 132 | 14 | 118 |
| 50 kya → present | introgression | 0.04 | 0/100 | 132 | 10 | 122 |
| 50 kya → present | introgression | 0.05 | 0/100 | 132 | 3 | 129 |
| 40–20 kya | de_novo | 0.001 | 0/100 | 132 | 132 | 0 |
| 40–20 kya | de_novo | 0.005 | 0/100 | 132 | 132 | 0 |
| 40–20 kya | de_novo | 0.01 | 4/100 | 132 | 128 | 0 |
| 40–20 kya | de_novo | 0.02 | 3/100 | 132 | 128 | 1 |
| 40–20 kya | de_novo | 0.03 | 0/100 | 132 | 129 | 3 |
| 40–20 kya | de_novo | 0.04 | 0/100 | 132 | 126 | 6 |
| 40–20 kya | de_novo | 0.05 | 0/100 | 132 | 127 | 5 |
| 40–20 kya | introgression | 0.001 | 9/100 | 132 | 123 | 0 |
| 40–20 kya | introgression | 0.005 | 26/100 | 132 | 106 | 0 |
| 40–20 kya | introgression | 0.01 | 25/100 | 132 | 105 | 2 |
| 40–20 kya | introgression | 0.02 | 15/100 | 132 | 103 | 14 |
| 40–20 kya | introgression | 0.03 | 0/100 | 131 | 102 | 29 |
| 40–20 kya | introgression | 0.04 | 0/100 | 131 | 92 | 39 |
| 40–20 kya | introgression | 0.05 | 0/100 | 131 | 96 | 35 |

| Window | Origin | s | Truth joint TPR | Median joint TPR | Mass joint TPR | Mass region F1 | Mass temporal F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 0/2 | 0/2 | 0/2 | 0.000 | 0.000 |
| 10 kya → present | de_novo | 0.005 | 0/13 | 0/13 | 0/13 | 0.000 | 0.000 |
| 10 kya → present | de_novo | 0.01 | 0/15 | 0/15 | 0/15 | 0.000 | 0.000 |
| 10 kya → present | de_novo | 0.02 | 0/20 | 0/20 | 1/20 | 0.087 | 0.000 |
| 10 kya → present | de_novo | 0.03 | 13/26 | 7/26 | 7/26 | 0.400 | 0.116 |
| 10 kya → present | de_novo | 0.04 | 9/10 | 7/10 | 7/10 | 0.737 | 0.452 |
| 10 kya → present | de_novo | 0.05 | 7/7 | 7/7 | 7/7 | 0.875 | 0.783 |
| 10 kya → present | introgression | 0.001 | 0/8 | 0/8 | 0/8 | 0.000 | 0.000 |
| 10 kya → present | introgression | 0.005 | 0/38 | 0/38 | 1/38 | 0.049 | 0.000 |
| 10 kya → present | introgression | 0.01 | 0/43 | 0/43 | 2/43 | 0.085 | 0.000 |
| 10 kya → present | introgression | 0.02 | 3/56 | 13/56 | 36/56 | 0.766 | 0.402 |
| 10 kya → present | introgression | 0.03 | 20/59 | 28/59 | 56/59 | 0.957 | 0.657 |
| 10 kya → present | introgression | 0.04 | 8/14 | 14/14 | 14/14 | 0.933 | 0.709 |
| 10 kya → present | introgression | 0.05 | 7/13 | 9/13 | 13/13 | 0.929 | 0.680 |
| 50 kya → present | de_novo | 0.001 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.005 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.01 | 1/1 | 0/1 | 1/1 | 0.500 | 0.000 |
| 50 kya → present | de_novo | 0.02 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.03 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.04 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.05 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.001 | 0/13 | 0/13 | 0/13 | 0.000 | 0.000 |
| 50 kya → present | introgression | 0.005 | 21/33 | 0/33 | 9/33 | 0.409 | 0.044 |
| 50 kya → present | introgression | 0.01 | 38/40 | 4/40 | 28/40 | 0.800 | 0.288 |
| 50 kya → present | introgression | 0.02 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.03 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.04 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.05 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.001 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.005 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.01 | 0/4 | 0/4 | 0/4 | 0.000 | 0.000 |
| 40–20 kya | de_novo | 0.02 | 2/3 | 0/3 | 1/3 | 0.500 | 0.017 |
| 40–20 kya | de_novo | 0.03 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.04 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.05 | NA | NA | NA | NA | NA |
| 40–20 kya | introgression | 0.001 | 0/9 | 0/9 | 0/9 | 0.000 | 0.000 |
| 40–20 kya | introgression | 0.005 | 0/26 | 0/26 | 0/26 | 0.000 | 0.000 |
| 40–20 kya | introgression | 0.01 | 3/25 | 0/25 | 0/25 | 0.000 | 0.000 |
| 40–20 kya | introgression | 0.02 | 10/15 | 1/15 | 10/15 | 0.741 | 0.065 |
| 40–20 kya | introgression | 0.03 | NA | NA | NA | NA | NA |
| 40–20 kya | introgression | 0.04 | NA | NA | NA | NA | NA |
| 40–20 kya | introgression | 0.05 | NA | NA | NA | NA | NA |

| Window | Origin | Method | Neutral false-positive regions | 95% Wilson interval |
|---|---|---|---:|---:|
| 10 kya → present | de_novo | truth | 3/50 | 0.021–0.162 |
| 10 kya → present | de_novo | median | 2/50 | 0.011–0.135 |
| 10 kya → present | de_novo | mass | 2/50 | 0.011–0.135 |
| 10 kya → present | introgression | truth | 0/50 | 0.000–0.071 |
| 10 kya → present | introgression | median | 2/50 | 0.011–0.135 |
| 10 kya → present | introgression | mass | 2/50 | 0.011–0.135 |
| 50 kya → present | de_novo | truth | 5/50 | 0.043–0.214 |
| 50 kya → present | de_novo | median | 2/50 | 0.011–0.135 |
| 50 kya → present | de_novo | mass | 2/50 | 0.011–0.135 |
| 50 kya → present | introgression | truth | 0/50 | 0.000–0.071 |
| 50 kya → present | introgression | median | 2/50 | 0.011–0.135 |
| 50 kya → present | introgression | mass | 2/50 | 0.011–0.135 |
| 40–20 kya | de_novo | truth | 1/50 | 0.004–0.105 |
| 40–20 kya | de_novo | median | 0/50 | 0.000–0.071 |
| 40–20 kya | de_novo | mass | 0/50 | 0.000–0.071 |
| 40–20 kya | introgression | truth | 0/50 | 0.000–0.071 |
| 40–20 kya | introgression | median | 2/50 | 0.011–0.135 |
| 40–20 kya | introgression | mass | 2/50 | 0.011–0.135 |

## Posterior mass: temporal errors near the selected site

Denominators are surviving selected regions. “Only outside” means a call within 100 kb but no temporal overlap. “Any outside component” also counts an outside component alongside a correct call. Outside-time error is unobservable for the 50-kya ongoing episode because it spans the entire scan.

| Window | Origin | s | Only outside time | Any outside component |
|---|---|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 0/2 | 0/2 |
| 10 kya → present | de_novo | 0.005 | 0/13 | 0/13 |
| 10 kya → present | de_novo | 0.01 | 0/15 | 0/15 |
| 10 kya → present | de_novo | 0.02 | 0/20 | 0/20 |
| 10 kya → present | de_novo | 0.03 | 0/26 | 0/26 |
| 10 kya → present | de_novo | 0.04 | 0/10 | 0/10 |
| 10 kya → present | de_novo | 0.05 | 0/7 | 0/7 |
| 10 kya → present | introgression | 0.001 | 0/8 | 0/8 |
| 10 kya → present | introgression | 0.005 | 0/38 | 0/38 |
| 10 kya → present | introgression | 0.01 | 0/43 | 0/43 |
| 10 kya → present | introgression | 0.02 | 0/56 | 0/56 |
| 10 kya → present | introgression | 0.03 | 0/59 | 0/59 |
| 10 kya → present | introgression | 0.04 | 0/14 | 0/14 |
| 10 kya → present | introgression | 0.05 | 0/13 | 0/13 |
| 50 kya → present | de_novo | 0.01 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.001 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.005 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.01 | Unobservable | Unobservable |
| 40–20 kya | de_novo | 0.01 | 0/4 | 0/4 |
| 40–20 kya | de_novo | 0.02 | 2/3 | 3/3 |
| 40–20 kya | introgression | 0.001 | 0/9 | 0/9 |
| 40–20 kya | introgression | 0.005 | 3/26 | 3/26 |
| 40–20 kya | introgression | 0.01 | 4/25 | 4/25 |
| 40–20 kya | introgression | 0.02 | 5/15 | 8/15 |

## Secondary scan: widths 1, 2, 5, 10 kyr

All windows start on the same 1-kyr grid. Each estimator is calibrated against maxima over all widths and positions; the primary/secondary union is not calibrated.

| Window | Origin | s | Mass joint TPR | Mass temporal F1 | Neutral false-positive regions |
|---|---|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 0/2 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.005 | 0/13 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.01 | 0/15 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.02 | 1/20 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.03 | 7/26 | 0.149 | 2/50 |
| 10 kya → present | de_novo | 0.04 | 7/10 | 0.482 | 2/50 |
| 10 kya → present | de_novo | 0.05 | 7/7 | 0.683 | 2/50 |
| 10 kya → present | introgression | 0.001 | 0/8 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.005 | 1/38 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.01 | 2/43 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.02 | 36/56 | 0.403 | 2/50 |
| 10 kya → present | introgression | 0.03 | 56/59 | 0.574 | 2/50 |
| 10 kya → present | introgression | 0.04 | 14/14 | 0.642 | 2/50 |
| 10 kya → present | introgression | 0.05 | 13/13 | 0.598 | 2/50 |
| 50 kya → present | de_novo | 0.001 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.005 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.01 | 1/1 | 0.000 | 2/50 |
| 50 kya → present | de_novo | 0.02 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.03 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.04 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.05 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.001 | 0/13 | 0.000 | 2/50 |
| 50 kya → present | introgression | 0.005 | 9/33 | 0.059 | 2/50 |
| 50 kya → present | introgression | 0.01 | 28/40 | 0.344 | 2/50 |
| 50 kya → present | introgression | 0.02 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.03 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.04 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.05 | NA | NA | 2/50 |
| 40–20 kya | de_novo | 0.001 | NA | NA | 0/50 |
| 40–20 kya | de_novo | 0.005 | NA | NA | 0/50 |
| 40–20 kya | de_novo | 0.01 | 0/4 | 0.000 | 0/50 |
| 40–20 kya | de_novo | 0.02 | 2/3 | 0.165 | 0/50 |
| 40–20 kya | de_novo | 0.03 | NA | NA | 0/50 |
| 40–20 kya | de_novo | 0.04 | NA | NA | 0/50 |
| 40–20 kya | de_novo | 0.05 | NA | NA | 0/50 |
| 40–20 kya | introgression | 0.001 | 0/9 | 0.000 | 2/50 |
| 40–20 kya | introgression | 0.005 | 0/26 | 0.000 | 2/50 |
| 40–20 kya | introgression | 0.01 | 1/25 | 0.000 | 2/50 |
| 40–20 kya | introgression | 0.02 | 15/15 | 0.172 | 2/50 |
| 40–20 kya | introgression | 0.03 | NA | NA | 2/50 |
| 40–20 kya | introgression | 0.04 | NA | NA | 2/50 |
| 40–20 kya | introgression | 0.05 | NA | NA | 2/50 |

## Limits

- First 100 acceptances in uninterrupted proposal order per treatment; report actual deadline-limited counts.
- Null introgression pool shared across time windows because s=0 throughout; neutral de novo pools match forward duration.
- Primary: 49 adjacent intervals [1,2),...,[49,50) kya. Secondary: widths 1/2/5/10 kyr with 1-kyr starts, jointly calibrated over widths and positions. No data for [0,1) kya. Truth used only for evaluation.
- Conservative +1 empirical p-values from independent calibration-region maxima over all 100 positions and all candidate intervals; alpha .05 separately per estimator and scan. The union of the two scan variants is not calibrated.
- Fixed 1,000 unfiltered pairs out of 19,900; neither an exhaustive all-pairs posterior analysis nor independent pair tests.
- True-tree benchmarks apply the same interval statistic to true TMRCAs; they are not an optimal selection detector or a theoretical upper bound on power.
- Near-site mistiming is a call within 100 kb with zero temporal overlap; near-site off-target flags any wholly outside connected component, even alongside a correct component. Center versions restrict evaluation to the exact selected coordinate.
- Ongoing 50-kya selection occupies the full observed range, so outside-time false calls cannot be assessed there.
- center_peak_kya is the midpoint of the highest standardized window even when it is not significant; it is not by itself an inferred selection date.
- Region F1 depends on selected/null cohort sizes; balanced_region_f1 standardizes to equal selected and neutral prevalence. Temporal F1 compares center support to the observable selection interval.
- Reuse is restricted to prior 10-kya proposals and their nulls; previous simulation compute is outside the current budget.

## Illustrative dating error

`40_20k/de_novo/s0.02/selected_00039` is the first analyzed 40–20 kya region with true-tree joint detection and a posterior-mass center call wholly outside the selected period. This example was selected after inspecting outcomes and is illustrative; aggregate performance above supplies the denominators. The three panels in `older_episode_example.pdf` use the same 1,000 pairs, with independently calibrated thresholds for each estimator.
