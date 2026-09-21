# EAS lab meeting figures - completed cohort, 21 September 2026

The [15-figure collection](results/lab_meeting_20260921/README.md) presents
Gamma-SMC decoded TMRCA statistics for its primary results. Decoding is
available for 100 selected regions at s=0.005 and all 1,000 neutral regions.
Other selected arms have not been decoded, so the presentation does not show
decoded TMRCA power over the full selection grid. AF and iHS use simulated
genotypes and retain all 1,000 selected regions (100 per s=0.001 to 0.010).
Figure 07 alone includes true TMRCA as an explicitly labeled validation
comparison. No new simulations or Gamma-SMC decoding were started.

All-pair and ALT/ALT carrier scores use `frac_recent_T`. Carrier mass is
`n_recent_ALT_ALT / n_all_pairs`, with no ALT/REF penalty or REF/REF subtraction.
The primary endpoint compares a fixed position with that same position in
neutral regions; error rates are positional FPR, not FDR among discoveries.
The fixed 100-kb iHS window is separately calibrated.

At s=0.005 and T=50 kya, results using the decoded TMRCA scores are:

| Method | Power at p<=0.05 | Positional FPR at p<=0.05 | Power at p<=0.01 |
|---|---:|---:|---:|
| Archaic AF (genotypes) | 100% | 5.1% | 99% |
| Gamma-SMC carrier mass | 100% | 5.1% | 98% |
| Gamma-SMC all-pair recency | 61% | 5.0% | 41% |
| Positional iHS (genotypes) | 67% | 4.8% | 50% |

Carrier power remains 100% at nominal 5% for all tested T cutoffs. At nominal
1%, carrier power is 97-99% versus AF power of 99%. This cohort therefore
does not establish a decoded coalescence gain at weaker selection. The separate
pooled p<=0.001 sensitivity gives 92% carrier power, 90% AF power, and 23%
all-pair power; its extreme neutral tail is too small for precise validation.
The explicit validation figure compares 61% decoded all-pair power with 86%
using true TMRCA, with carrier power at 100% for both sources.

Suggested main sequence: 01, 02, 03, 04, 05, 06, 08, 09, 10. Figure 07 is
truth-versus-decoding validation; figures 11-15 cover example regions,
p=0.001 sensitivity, decoded carrier-score components, demographic
assumptions, and decoded pair-class diagnostics.
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
The score-component diagnostic does not identify introgressing founders or
demonstrate a gain from multiple independent introgressing haplotypes.
Figures 03_a and 03_b remain deferred: no historical AF trajectories were
saved, and the user declined replaying the simulations. The s=0.006 AF tail
passed its integrity audit; all selected replicates remain included.
