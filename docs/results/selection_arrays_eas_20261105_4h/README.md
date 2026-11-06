# Three EAS selection episodes

The resumed continuation finished on November 5, 2026 (Central time) after **3 h 57 min 26 s of active simulation time**, excluding the overnight pause. It added **40,015 completed proposals and 1,323 analyzed selected survivors** to the earlier cohort. The final analysis includes **1,816 selected regions and 676 neutral regions**; **12/42 treatments reached 100 survivors**. All 49,399 submitted proposals completed, with zero failures or pending jobs. Nine accepted in-flight overshoots are excluded by the first-100 rule.

Earlier completed regions were reused across all three episodes within this run directory. Their seeds, tree/posterior/score/receipt hashes, calibration arrays, and every earlier replicate metric were verified unchanged. The original experiment also imported prior 10-kya simulations from a separate run. The four-hour continuation budget excludes that earlier compute. See `CONTINUATION.md` and `continuation_provenance.json` for accounting and the nested before/after comparison.

The larger cohort supports three findings:

- Posterior mass substantially improves detection over the median for recent introgressed selection: 66/100 versus 19/100 at s=0.02, and 94/100 versus 49/100 at s=0.03. Recent weak selection remains difficult: mass detects 1/100 at s=0.005 and 4/100 at s=0.01. Recent de novo power rises from 2/49 at s=0.02 to 30/80, 46/60, and 80/82 at s=0.03, 0.04, and 0.05.
- Older ongoing introgressed selection is detectable at weaker coefficients: mass detects 25/100 at s=0.005 and 69/100 at s=0.01. However, the ongoing 50-kya episode covers the entire observed time axis, so this experiment cannot test outside-period temporal error for that episode.
- The terminated 40–20 kya episode remains difficult to date. At s=0.02 with introgression, mass detects overlap in 51/100 regions, but 36/100 have near-site calls only outside the true period, and 47/100 have at least one wholly outside component. Mean exact-center temporal F1 is 0.049. The secondary 1/2/5/10-kyr scan raises overlap detection to 68/100 and temporal F1 to 0.129, with 19/100 wholly mistimed regions. Broader overlap does not establish accurate interval recovery.

Held-out posterior-mass neutral FPR is 2/50 (4%; Wilson 95% CI 1.1–13.5%) for each applicable pool, except the 40-kya de novo null, which has 0/50 (upper 95% bound 7.1%). The introgression null pool is shared across episodes. These are per-1-Mb-region error rates, separately calibrated for each estimator and scan, not a genome-wide error guarantee. Selected-region mistiming and neutral-region false positives answer different questions and are reported separately.

| Primary posterior-mass example | Joint TPR | Balanced region F1 | Event F1 | Mean center temporal F1 |
|---|---:|---:|---:|---:|
| Recent introgression, s=0.02 | 66/100 | 0.776 | 0.600 | 0.396 |
| Recent introgression, s=0.03 | 94/100 | 0.949 | 0.662 | 0.643 |
| 50-kya ongoing introgression, s=0.005 | 25/100 | 0.388 | 0.350 | 0.051 |
| 50-kya ongoing introgression, s=0.01 | 69/100 | 0.798 | 0.636 | 0.256 |
| 40–20 kya introgression, s=0.02 | 51/100 | 0.658 | 0.338 | 0.049 |

Balanced region F1 uses equal selected/neutral prevalence; event F1 penalizes unmatched connected call components; temporal F1 measures the called time support at the exact focal coordinate. A low temporal F1 can reflect absent or spatially displaced calls, narrow coverage, or mistiming. The separate mistiming counts are needed to identify actual outside-period errors. Rates are conditional on present-day whole-population segregation, and zero surviving replicates make power unestimable.

The exact **19,900-pair true-TMRCA audit** covers all 2,492 analyzed regions, reusing 1,169 verified cached regions and computing 1,323 new ones. For introgression, its primary TPR is 58/100 and 92/100 for 50-kya ongoing s=0.005 and 0.01, and 21/100 and 73/100 for 40–20 kya s=0.01 and 0.02. This leaves scope to improve inferred-TMRCA performance for these episodes. The true-tree scan is the same statistic applied to true times, not an optimal detector or a theoretical upper bound. The posterior scan remains a fixed 1,000-pair analysis; exhaustive all-pairs posterior decoding was not performed.

All five figures were visually checked and saved as PNG plus vector PDF with embedded editable TrueType fonts. Twelve relevant tests passed. The fixed-seed audit found three individual phase-seed collisions among rejected proposals only; no analyzed region was involved and all complete simulation seed vectors were unique. Seeds were preserved. `validation_checks.json` records the checks, and `output_hashes.json` covers the complete saved report.

Selection 10 kya to present, 50 kya to present, and 40–20 kya; s = 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05; target 100 surviving segregating replicates per treatment (42 treatments).
EAS PHLASH history, unscaled WF, h=0.5, 25 years/generation, 1 Mb, 100 diploids. De novo starts at one copy at onset. Introgression injects 1% donor haplotypes at 50 kya, with a fixed donor allele, human–Neanderthal split 500 kya and donor Ne 3,600. Selected loss/fixation rejected; neutral regions unconditioned with nearest segregating focal label.
Detection scans all 49 adjacent one-kyr intervals from 1 to 50 kya, every 10 kb. Development nulls standardize interval coalescence fractions. Independent calibration maxima control the whole-region scan at nominal 5% per estimator; held-out nulls estimate FPR. All three estimators use the same 1,000 unfiltered pairs (of 19,900).
Joint TPR requires a call within 100 kb of the selected site and overlapping the true selection interval. Spatial TPR ignores time. Temporal F1 is the Dice overlap at the exact focal coordinate. A wholly outside temporal call is reported separately from held-out neutral FPR. Region F1 depends on the selected/null mixture.

Completed proposals: 49399; selected regions analyzed: 1816.

| Window | Origin | s | Accepted / target | Prefix proposals | Lost | Fixed |
|---|---|---:|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 20/100 | 1226 | 1206 | 0 |
| 10 kya → present | de_novo | 0.005 | 35/100 | 1622 | 1587 | 0 |
| 10 kya → present | de_novo | 0.01 | 37/100 | 1622 | 1585 | 0 |
| 10 kya → present | de_novo | 0.02 | 49/100 | 1622 | 1573 | 0 |
| 10 kya → present | de_novo | 0.03 | 80/100 | 1622 | 1542 | 0 |
| 10 kya → present | de_novo | 0.04 | 60/100 | 1225 | 1165 | 0 |
| 10 kya → present | de_novo | 0.05 | 82/100 | 1225 | 1143 | 0 |
| 10 kya → present | introgression | 0.001 | 100/100 | 1169 | 1069 | 0 |
| 10 kya → present | introgression | 0.005 | 100/100 | 1257 | 1157 | 0 |
| 10 kya → present | introgression | 0.01 | 100/100 | 1222 | 1122 | 0 |
| 10 kya → present | introgression | 0.02 | 100/100 | 1014 | 914 | 0 |
| 10 kya → present | introgression | 0.03 | 100/100 | 985 | 885 | 0 |
| 10 kya → present | introgression | 0.04 | 98/100 | 1225 | 1122 | 5 |
| 10 kya → present | introgression | 0.05 | 91/100 | 1225 | 1111 | 23 |
| 50 kya → present | de_novo | 0.001 | 2/100 | 1225 | 1223 | 0 |
| 50 kya → present | de_novo | 0.005 | 1/100 | 1225 | 1223 | 1 |
| 50 kya → present | de_novo | 0.01 | 8/100 | 1225 | 1216 | 1 |
| 50 kya → present | de_novo | 0.02 | 0/100 | 1225 | 1206 | 19 |
| 50 kya → present | de_novo | 0.03 | 0/100 | 1225 | 1193 | 32 |
| 50 kya → present | de_novo | 0.04 | 0/100 | 1225 | 1181 | 44 |
| 50 kya → present | de_novo | 0.05 | 0/100 | 1225 | 1164 | 61 |
| 50 kya → present | introgression | 0.001 | 100/100 | 765 | 665 | 0 |
| 50 kya → present | introgression | 0.005 | 100/100 | 351 | 251 | 0 |
| 50 kya → present | introgression | 0.01 | 100/100 | 325 | 171 | 54 |
| 50 kya → present | introgression | 0.02 | 1/100 | 1225 | 286 | 938 |
| 50 kya → present | introgression | 0.03 | 0/100 | 1225 | 139 | 1086 |
| 50 kya → present | introgression | 0.04 | 0/100 | 1225 | 86 | 1139 |
| 50 kya → present | introgression | 0.05 | 0/100 | 1225 | 35 | 1190 |
| 40–20 kya | de_novo | 0.001 | 5/100 | 1225 | 1220 | 0 |
| 40–20 kya | de_novo | 0.005 | 6/100 | 1225 | 1219 | 0 |
| 40–20 kya | de_novo | 0.01 | 12/100 | 1225 | 1213 | 0 |
| 40–20 kya | de_novo | 0.02 | 21/100 | 1225 | 1200 | 4 |
| 40–20 kya | de_novo | 0.03 | 4/100 | 1225 | 1192 | 29 |
| 40–20 kya | de_novo | 0.04 | 0/100 | 1225 | 1175 | 50 |
| 40–20 kya | de_novo | 0.05 | 0/100 | 1225 | 1170 | 55 |
| 40–20 kya | introgression | 0.001 | 100/100 | 1219 | 1119 | 0 |
| 40–20 kya | introgression | 0.005 | 100/100 | 582 | 482 | 0 |
| 40–20 kya | introgression | 0.01 | 100/100 | 471 | 369 | 2 |
| 40–20 kya | introgression | 0.02 | 100/100 | 967 | 743 | 124 |
| 40–20 kya | introgression | 0.03 | 4/100 | 1225 | 940 | 281 |
| 40–20 kya | introgression | 0.04 | 0/100 | 1225 | 900 | 325 |
| 40–20 kya | introgression | 0.05 | 0/100 | 1225 | 896 | 329 |

| Window | Origin | s | Truth joint TPR | Median joint TPR | Mass joint TPR | Mass region F1 | Mass temporal F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 0/20 | 0/20 | 0/20 | 0.000 | 0.000 |
| 10 kya → present | de_novo | 0.005 | 0/35 | 0/35 | 1/35 | 0.053 | 0.000 |
| 10 kya → present | de_novo | 0.01 | 0/37 | 0/37 | 1/37 | 0.050 | 0.000 |
| 10 kya → present | de_novo | 0.02 | 0/49 | 0/49 | 2/49 | 0.075 | 0.000 |
| 10 kya → present | de_novo | 0.03 | 48/80 | 26/80 | 30/80 | 0.536 | 0.172 |
| 10 kya → present | de_novo | 0.04 | 53/60 | 48/60 | 46/60 | 0.852 | 0.534 |
| 10 kya → present | de_novo | 0.05 | 82/82 | 80/82 | 80/82 | 0.976 | 0.791 |
| 10 kya → present | introgression | 0.001 | 0/100 | 0/100 | 0/100 | 0.000 | 0.000 |
| 10 kya → present | introgression | 0.005 | 0/100 | 0/100 | 1/100 | 0.019 | 0.000 |
| 10 kya → present | introgression | 0.01 | 0/100 | 0/100 | 4/100 | 0.075 | 0.000 |
| 10 kya → present | introgression | 0.02 | 3/100 | 19/100 | 66/100 | 0.786 | 0.396 |
| 10 kya → present | introgression | 0.03 | 31/100 | 49/100 | 94/100 | 0.959 | 0.643 |
| 10 kya → present | introgression | 0.04 | 30/98 | 64/98 | 93/98 | 0.964 | 0.641 |
| 10 kya → present | introgression | 0.05 | 37/91 | 68/91 | 88/91 | 0.972 | 0.643 |
| 50 kya → present | de_novo | 0.001 | 0/2 | 0/2 | 0/2 | 0.000 | 0.000 |
| 50 kya → present | de_novo | 0.005 | 1/1 | 0/1 | 0/1 | 0.000 | 0.000 |
| 50 kya → present | de_novo | 0.01 | 8/8 | 0/8 | 8/8 | 0.889 | 0.371 |
| 50 kya → present | de_novo | 0.02 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.03 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.04 | NA | NA | NA | NA | NA |
| 50 kya → present | de_novo | 0.05 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.001 | 3/100 | 1/100 | 0/100 | 0.000 | 0.000 |
| 50 kya → present | introgression | 0.005 | 58/100 | 1/100 | 25/100 | 0.394 | 0.051 |
| 50 kya → present | introgression | 0.01 | 92/100 | 8/100 | 69/100 | 0.807 | 0.256 |
| 50 kya → present | introgression | 0.02 | 1/1 | 0/1 | 1/1 | 0.500 | 0.658 |
| 50 kya → present | introgression | 0.03 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.04 | NA | NA | NA | NA | NA |
| 50 kya → present | introgression | 0.05 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.001 | 0/5 | 0/5 | 0/5 | 0.000 | 0.000 |
| 40–20 kya | de_novo | 0.005 | 0/6 | 0/6 | 0/6 | 0.000 | 0.000 |
| 40–20 kya | de_novo | 0.01 | 0/12 | 0/12 | 0/12 | 0.000 | 0.000 |
| 40–20 kya | de_novo | 0.02 | 15/21 | 0/21 | 14/21 | 0.800 | 0.085 |
| 40–20 kya | de_novo | 0.03 | 3/4 | 0/4 | 4/4 | 1.000 | 0.197 |
| 40–20 kya | de_novo | 0.04 | NA | NA | NA | NA | NA |
| 40–20 kya | de_novo | 0.05 | NA | NA | NA | NA | NA |
| 40–20 kya | introgression | 0.001 | 0/100 | 0/100 | 0/100 | 0.000 | 0.000 |
| 40–20 kya | introgression | 0.005 | 2/100 | 0/100 | 0/100 | 0.000 | 0.000 |
| 40–20 kya | introgression | 0.01 | 22/100 | 0/100 | 3/100 | 0.057 | 0.003 |
| 40–20 kya | introgression | 0.02 | 73/100 | 2/100 | 51/100 | 0.667 | 0.049 |
| 40–20 kya | introgression | 0.03 | 2/4 | 0/4 | 2/4 | 0.500 | 0.164 |
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
| 10 kya → present | de_novo | 0.001 | 0/20 | 0/20 |
| 10 kya → present | de_novo | 0.005 | 0/35 | 0/35 |
| 10 kya → present | de_novo | 0.01 | 0/37 | 0/37 |
| 10 kya → present | de_novo | 0.02 | 0/49 | 0/49 |
| 10 kya → present | de_novo | 0.03 | 0/80 | 0/80 |
| 10 kya → present | de_novo | 0.04 | 0/60 | 0/60 |
| 10 kya → present | de_novo | 0.05 | 0/82 | 0/82 |
| 10 kya → present | introgression | 0.001 | 0/100 | 0/100 |
| 10 kya → present | introgression | 0.005 | 0/100 | 0/100 |
| 10 kya → present | introgression | 0.01 | 0/100 | 0/100 |
| 10 kya → present | introgression | 0.02 | 0/100 | 0/100 |
| 10 kya → present | introgression | 0.03 | 0/100 | 0/100 |
| 10 kya → present | introgression | 0.04 | 0/98 | 0/98 |
| 10 kya → present | introgression | 0.05 | 0/91 | 0/91 |
| 50 kya → present | de_novo | 0.001 | Unobservable | Unobservable |
| 50 kya → present | de_novo | 0.005 | Unobservable | Unobservable |
| 50 kya → present | de_novo | 0.01 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.001 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.005 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.01 | Unobservable | Unobservable |
| 50 kya → present | introgression | 0.02 | Unobservable | Unobservable |
| 40–20 kya | de_novo | 0.001 | 0/5 | 0/5 |
| 40–20 kya | de_novo | 0.005 | 0/6 | 0/6 |
| 40–20 kya | de_novo | 0.01 | 0/12 | 0/12 |
| 40–20 kya | de_novo | 0.02 | 6/21 | 9/21 |
| 40–20 kya | de_novo | 0.03 | 0/4 | 2/4 |
| 40–20 kya | introgression | 0.001 | 1/100 | 1/100 |
| 40–20 kya | introgression | 0.005 | 4/100 | 4/100 |
| 40–20 kya | introgression | 0.01 | 16/100 | 17/100 |
| 40–20 kya | introgression | 0.02 | 36/100 | 47/100 |
| 40–20 kya | introgression | 0.03 | 1/4 | 1/4 |

## Secondary scan: widths 1, 2, 5, 10 kyr

All windows start on the same 1-kyr grid. Each estimator is calibrated against maxima over all widths and positions; the primary/secondary union is not calibrated.

| Window | Origin | s | Mass joint TPR | Mass temporal F1 | Neutral false-positive regions |
|---|---|---:|---:|---:|---:|
| 10 kya → present | de_novo | 0.001 | 0/20 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.005 | 1/35 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.01 | 1/37 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.02 | 2/49 | 0.000 | 2/50 |
| 10 kya → present | de_novo | 0.03 | 30/80 | 0.202 | 2/50 |
| 10 kya → present | de_novo | 0.04 | 46/60 | 0.504 | 2/50 |
| 10 kya → present | de_novo | 0.05 | 80/82 | 0.657 | 2/50 |
| 10 kya → present | introgression | 0.001 | 0/100 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.005 | 1/100 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.01 | 4/100 | 0.000 | 2/50 |
| 10 kya → present | introgression | 0.02 | 66/100 | 0.392 | 2/50 |
| 10 kya → present | introgression | 0.03 | 94/100 | 0.556 | 2/50 |
| 10 kya → present | introgression | 0.04 | 93/98 | 0.563 | 2/50 |
| 10 kya → present | introgression | 0.05 | 88/91 | 0.558 | 2/50 |
| 50 kya → present | de_novo | 0.001 | 0/2 | 0.000 | 2/50 |
| 50 kya → present | de_novo | 0.005 | 0/1 | 0.000 | 2/50 |
| 50 kya → present | de_novo | 0.01 | 8/8 | 0.434 | 2/50 |
| 50 kya → present | de_novo | 0.02 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.03 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.04 | NA | NA | 2/50 |
| 50 kya → present | de_novo | 0.05 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.001 | 0/100 | 0.000 | 2/50 |
| 50 kya → present | introgression | 0.005 | 25/100 | 0.068 | 2/50 |
| 50 kya → present | introgression | 0.01 | 69/100 | 0.307 | 2/50 |
| 50 kya → present | introgression | 0.02 | 1/1 | 0.744 | 2/50 |
| 50 kya → present | introgression | 0.03 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.04 | NA | NA | 2/50 |
| 50 kya → present | introgression | 0.05 | NA | NA | 2/50 |
| 40–20 kya | de_novo | 0.001 | 0/5 | 0.000 | 0/50 |
| 40–20 kya | de_novo | 0.005 | 0/6 | 0.000 | 0/50 |
| 40–20 kya | de_novo | 0.01 | 0/12 | 0.000 | 0/50 |
| 40–20 kya | de_novo | 0.02 | 18/21 | 0.213 | 0/50 |
| 40–20 kya | de_novo | 0.03 | 4/4 | 0.318 | 0/50 |
| 40–20 kya | de_novo | 0.04 | NA | NA | 0/50 |
| 40–20 kya | de_novo | 0.05 | NA | NA | 0/50 |
| 40–20 kya | introgression | 0.001 | 0/100 | 0.000 | 2/50 |
| 40–20 kya | introgression | 0.005 | 0/100 | 0.000 | 2/50 |
| 40–20 kya | introgression | 0.01 | 5/100 | 0.008 | 2/50 |
| 40–20 kya | introgression | 0.02 | 68/100 | 0.129 | 2/50 |
| 40–20 kya | introgression | 0.03 | 2/4 | 0.232 | 2/50 |
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
- Same-manifest continuations reuse verified completed proposals from every regime. Cross-root historical reuse is restricted to compatible prior 10-kya proposals and their nulls; previous simulation compute is outside the current budget.

## Illustrative dating error

`40_20k/de_novo/s0.02/selected_00039` is the first analyzed 40–20 kya region with true-tree joint detection and a posterior-mass center call wholly outside the selected period. This example was selected after inspecting outcomes and is illustrative; aggregate performance above supplies the denominators. The three panels in `older_episode_example.pdf` use the same 1,000 pairs, with independently calibrated thresholds for each estimator.

## Reproduction on Linux / WSL / HPC

The following notebook command cell starts a fresh deterministic experiment with a four-hour simulation cap, then validates and analyzes it. Rerunning with the same output path resumes verified completed proposals; it does not overwrite them based on the Git revision. Existing native binaries are retained to avoid changing the manifest contract. A fresh run has no historical simulation cache, so its time-limited sample counts can differ from this continuation. Analysis runs after the simulation cap. No plots are displayed inline and no data are pushed by this cell.

```bash
%%bash
set -euo pipefail
REPO=gamma_smc_ts_selection
if [[ -d "$REPO/.git" ]]; then
  git -C "$REPO" pull --ff-only
else
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
fi
cd "$REPO"

# Phase 1: install the pinned environment only when it is absent.
if [[ ! -x .venv/bin/python || ! -x .native/bin/slim || ! -x bin/gamma_smc ]]; then
  bash scripts/bootstrap_uv.sh --skip-tests
fi
if [[ -x .tools/uv-bin/uv ]]; then UV="$PWD/.tools/uv-bin/uv"; else UV="$(command -v uv)"; fi
export LD_LIBRARY_PATH="$PWD/.native/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MPLBACKEND=Agg
OUT="$PWD/../sim/selection_arrays_eas_fresh"
REPORT="$PWD/../results/selection_arrays_eas_fresh"
mkdir -p "$OUT" "$REPORT"

# Phase 2: fixed-seed simulations, with existing valid receipts reused.
"$UV" run --frozen --all-extras python scripts/run_selection_arrays.py \
  --out "$OUT" --slim "$PWD/.native/bin/slim" --decoder "$PWD/bin/gamma_smc" \
  --workers 24 --seconds 14400 --drain-seconds 150 --accepted 100 --attempts 10000 \
  2>&1 | tee -a "$OUT/console.log"

# Phase 3: integrity checks, inference metrics, and exact true-pair benchmark.
"$UV" run --frozen --all-extras python scripts/report_selection_arrays.py \
  --root "$OUT" --out "$REPORT" 2>&1 | tee -a "$REPORT/validation.log"
"$UV" run --frozen --all-extras python scripts/audit_selection_arrays.py \
  --root "$OUT" --out "$REPORT" --workers 24 2>&1 | tee -a "$REPORT/audit.log"
```

To rebuild this particular report from the retained local raw outputs, use `scripts/report_selection_arrays.py --root /mnt/d/phase2simselection/sim/selection_arrays_20261104 --out <report-directory>`, followed by the all-pairs audit. This requires the raw outputs at their recorded paths; the repository contains summary arrays and provenance, not the full simulated trees/posteriors. The prior truth cache can be copied into the report directory together with its proof before the audit; the audit verifies the algorithm, tree, and cache hashes before reuse. The historical continuation arguments and absolute deadline are recorded in `schedule.json`.

The pause-aware comparison command is:

```bash
uv run --frozen --all-extras python scripts/summarize_selection_continuation.py \
  --before docs/results/selection_arrays_eas_20261104 \
  --after docs/results/selection_arrays_eas_20261105_4h \
  --started-epoch 1793914147 --previous-active-seconds 1375.5174491405487
```

This last command rewrites `CONTINUATION.md`, its two comparison tables, and `continuation_provenance.json`; it does not run simulations. The generated report and comparison files can be rebuilt independently of the manual interpretation paragraphs above.
