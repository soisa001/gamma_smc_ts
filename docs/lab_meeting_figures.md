# EAS lab meeting figures - completed cohort, 21 September 2026

The [15-figure collection](results/lab_meeting_20260921/README.md) uses all
1,000 selected regions (100 at each s from 0.001 to 0.010) and 1,000 neutral
regions. It supersedes the 947-selected snapshot in the earlier positional
and iHS summaries. Gamma-SMC comparisons remain limited to the already
decoded s=0.005 cohort and neutrals. This figure task started no simulations
and no new Gamma-SMC decoding; it extracted the final 53 focal records from
the completed archive and refreshed iHS from those saved regions.

All-pair and ALT/ALT carrier scores use `frac_recent_T`. Carrier mass is
`n_recent_ALT_ALT / n_all_pairs`, with no ALT/REF penalty or REF/REF subtraction.
The primary endpoint compares a fixed position with that same position in
neutral regions; error rates are positional FPR, not FDR among discoveries.
The fixed 100-kb iHS window is separately calibrated.

At T=50 kya and nominal p<=0.05, completed-cohort results are:

| s | AF power | Carrier-mass power | All-pair power |
|---|---:|---:|---:|
| 0.001 | 43% | 47% | 5% |
| 0.002 | 70% | 70% | 19% |
| 0.003 | 79% | 79% | 41% |
| 0.004 | 92% | 92% | 62% |
| 0.005 | 100% | 100% | 86% |
| 0.006 | 98% | 98% | 89% |
| 0.007 | 99% | 99% | 94% |
| 0.008 | 98% | 98% | 96% |
| 0.009 | 100% | 100% | 99% |
| 0.010 | 100% | 100% | 100% |

The respective neutral FPRs are 5.1%, 5.0%, and 5.4%. At s=0.001,
T=20-kya carrier mass detects 54/100, compared with 43/100 using AF alone
(12 carrier-only detections, one AF-only detection). This comparison across
cutoffs is exploratory. At s=0.005, decoded carrier-mass power at T=50 kya
is 100%, versus 61% for decoded all-pair recency, with FPRs 5.1% and 5.0%.

Suggested main sequence: 01, 02, 03, 04, 05, 07, 08, 09, 10. Figure 06 checks
neutral calibration; figures 11-15 cover example regions, p=0.001 sensitivity,
carrier genealogies, demographic assumptions, and pair-class diagnostics.
Each figure has a high-resolution PNG, a vector PDF with editable embedded
text, a caption, and supporting CSVs. The combined PDF is 15 pages in 16:9
layout. Source hashes and checks accompany the figures. All plotted PDFs
were rendered and visually inspected.

The presentation ZIP is on D: at
`D:/phase2simselection/sim/eas_lab_meeting_20260921/EAS_lab_meeting_figure_pack.zip`.
The same presentation figures and tables are committed here without the
duplicate ZIP or raw tree sequences. See the package's
[reproduction instructions](results/lab_meeting_20260921/REPRODUCE.md) for
the build, optional summary refresh, and verification/package phases.

Interpretation limits are in the caption index: selected survival conditioning,
an unconditional positional null, oracle archaic markers, a strong shared
archaic bottleneck, uniform rates, and no genome-wide discovery-FDR claim.
The genealogy diagnostic does not establish a count of introgressing founders
or demonstrate a gain from multiple independent introgressing haplotypes.
