# Pairwise TMRCA by the focal archaic allele

The alt/alt versus ref/ref profiles show a strong within-carrier recent-coalescence
pattern that is obscured by the all-pair average, especially for selection
beginning at 10 kya. This is a descriptive analysis of the selected simulations;
it does not establish power or a neutral false-positive rate for a new detector.

All 200 selected regions were profiled: 100 with onset at 50 kya and 100 with
onset at 10 kya, using s=0.005, EAS, 400 haplotypes, 10,000 existing sampled
pairs, and the unchanged 2% introgression pulse at 50 kya. No simulation or HMM
decoding was rerun. Original seeds, sequences, pair manifests, and saved gamma
posteriors were reused.

## Exactly what alt/alt means

The focal selected SNP is treated as a known archaic-origin allele, as specified
by the simulation design. Each member of an arbitrary haplotype pair contributes
its allele at that SNP:

- **alt/alt:** both haplotypes carry the alternate archaic allele.
- **ref/ref:** both carry the reference allele at that SNP.
- **alt/ref:** one carries each allele.
- **all pairs:** the complete original 10,000-pair panel.

These classes are not restricted to within-individual pairs or to homozygous
individuals. Heterozygous individuals can contribute haplotypes to alt/alt and
ref/ref pairs with other individuals. A focused test explicitly checks this
case. Reference status at this marker is not a claim that a chromosome has no
archaic ancestry elsewhere.

Pair membership is fixed at the focal SNP, and those same pairs are followed
across the entire 10 Mb region. Thus the distance profiles show how the
association with that marker fades with genomic distance. They do not redefine
the pair classes at each stride.

## Figures

Each figure has true TMRCA in the top row and decoded TMRCA in the bottom row;
columns compare selection onset at 50 and 10 kya. Green is alt/alt, blue is
ref/ref, orange is alt/ref, and the dashed gray curve is all pairs.

The primary files are:

- [Pairwise ages, SNP-centered zoom](results/pair_classes/analysis/mean_tmrca_years_zoom.pdf)
  and [full 10 Mb region](results/pair_classes/analysis/mean_tmrca_years_full.pdf).
- [Fraction recent at T=50 kya, zoom](results/pair_classes/analysis/frac_recent_50000_zoom.pdf)
  and [full region](results/pair_classes/analysis/frac_recent_50000_full.pdf).
- [Fraction recent at T=10 kya, zoom](results/pair_classes/analysis/frac_recent_10000_zoom.pdf)
  and [full region](results/pair_classes/analysis/frac_recent_10000_full.pdf).

PNG copies are saved alongside the PDFs. Equivalent fraction plots are also
provided for T=5, 20, 30 and 40 kya. There are 14 figures, each saved in both
formats. The zoom extends 500 kb on each side; the full profiles use every
existing 10 kb stride. Age plots use a logarithmic age axis and report actual
ages, including older ancestry. All recent-coalescence thresholds remain at
or below 50 kya.

Within each region, mean pairwise age and `frac_recent_T` are calculated
separately in each pair class. Curves then average equally across available
regions. The decoded age diagnostic is the mean of the individual pairs'
posterior mean ages. Decoded `frac_recent_T` uses the existing native
posterior-mean hard-call criterion. No averaged posterior-mass statistic is used.

Shading is a pointwise 95% Student-t interval for the mean across independent
simulated regions. It does not treat overlapping haplotype pairs as independent
observations, and it is not a simultaneous interval for the entire profile.

## At the archaic SNP

The values below are mean class-specific fractions of pairs with TMRCA below
50 kya, expressed as percentages. All-pair values average all 100 regions per
onset. Class-specific availability is described below.

| Selection onset | TMRCA source | Alt/alt | Ref/ref | Alt/ref | All pairs |
|---|---|---:|---:|---:|---:|
| 50 kya | True | 91.2% | 30.2% | 0.0% | 59.5% |
| 50 kya | Decoded | 68.0% | 27.8% | 0.046% | 45.2% |
| 10 kya | True | 97.2% | 22.3% | 0.0% | 20.7% |
| 10 kya | Decoded | 78.9% | 22.2% | 0.075% | 19.3% |

The corresponding mean pairwise ages at the SNP are:

| Selection onset | TMRCA source | Alt/alt mean age | Ref/ref mean age | Alt/ref mean age |
|---|---|---:|---:|---:|
| 50 kya | True | 35.8 kya | 409.4 kya | 1,301.4 kya |
| 50 kya | Decoded | 56.3 kya | 405.3 kya | 1,217.3 kya |
| 10 kya | True | 22.2 kya | 464.1 kya | 1,403.8 kya |
| 10 kya | Decoded | 46.6 kya | 421.9 kya | 1,210.0 kya |

For the 10-kya onset arm, the contrast also appears at T=10 kya: true alt/alt
`frac_recent_10000` is 12.74%, versus 0.90% for ref/ref. Decoding gives 25.11%
and 2.44%, respectively. A selection onset at 10 kya does not imply that all
surviving carrier lineages coalesce after that date; the selected allele was
already present on introgressed variation.

The average true-TMRCA contrast fades with distance. For the 10-kya onset arm,
alt/alt `frac_recent_50000` falls from 97.2% at the SNP to 56.6% at 50 kb,
38.7% at 100 kb, 24.6% at 250 kb and 20.8% at 500 kb. These off-center values
average the left and right positions. Ref/ref is 20.2% at 500 kb. The contrast
in the 50-kya onset arm is also largely gone by 250–500 kb. This is the decay
of an average association profile, not a tract boundary assigned to each pair.

## Why the all-pair average is weak for the later onset

For each region and stride, the following identity holds exactly:

\[
F_{\mathrm{all}}(T)=
w_{AA}F_{AA}(T)+w_{RR}F_{RR}(T)+w_{AR}F_{AR}(T),
\]

where each weight is the fraction of sampled pairs in that class. In the
10-kya onset cohort, only 7.01% of sampled pairs are alt/alt on average;
64.50% are ref/ref and 28.49% are alt/ref. The mixed pairs have very old TMRCA
at the marker and contribute essentially no recent pairs at these cutoffs.

At T=50 kya, the exact mean contributions to the later-onset all-pair fraction
are:

| Source | Alt/alt contribution | Ref/ref contribution | Alt/ref contribution | Total |
|---|---:|---:|---:|---:|
| True | 6.4769 percentage points | 14.1822 pp | 0 pp | 20.6591% |
| Decoded | 4.9609 pp | 14.2814 pp | 0.0240 pp | 19.2663% |

Contributions are computed within each region before averaging. Multiplying
the cohort-average weights by the cohort-average conditional fractions would
not generally reproduce this table, because weights and fractions covary.

Thus a requirement for an unusually high all-pair fraction can discard a
region that has a strong recent-coalescence contrast among its archaic-allele
carriers. These figures support evaluating pair-class information directly.
Neutral introgressed carrier pairs can also be related; matched neutral
archaic-pair calibration is still required before interpreting this contrast
as selection detection or claiming improved specificity.

## Empty classes and verification

All 200 regions remain in the output. The fixed 10,000-pair manifest contains
no ref/ref pair in three onset-50k regions and no alt/alt pair in two onset-10k
regions. Those class summaries are missing, not zero. Accordingly, the ref/ref
mean for onset 50 kya uses 97 regions, and the alt/alt mean for onset 10 kya uses
98. Other displayed class means use 100 regions. No fixed or nearly fixed
region is discarded; all-pair statistics always retain the full cohort.

Three focused tests passed: arbitrary cross-individual class assignment,
explicit true pair ages with empty-class retention, and padded raw-posterior
replay against native mean counts, including rejection of a truncated stream.
All 200 regions reproduce the original all-pair true and decoded counts exactly.

The independent audit checks 1,600,000 profile rows, verifies 400,000 complete
pair-class partitions, and recomputes 112,000 summary means. It also verifies
input-profile and output-artifact hashes. Headless figure checks confirm text
within the canvas, legends outside panels, and vector PDFs with editable text.
Plots were not displayed in chat.

Full per-region profiles are stored under
`D:/phase2simselection/sim/eas_pair_classes/profiles`. Compact figure/data copies
are in [the artifact directory](results/pair_classes/artifact_manifest.json),
including [focal summaries](results/pair_classes/analysis/focal_summary.csv),
[pair counts](results/pair_classes/analysis/pair_counts.csv), and
[the independent audit](results/pair_classes/audit.json).

## Rerun in the configured WSL environment

The phases are `profiles`, `analyse`, and `run`; completed profiles are reused
when their scientific inputs, computation, and output hashes match. `analyse`
regenerates the summaries and figures and runs the independent audit. It does
not rerun simulation or HMM decoding.

```bash
%%bash
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
cd "$repo"
git pull --ff-only
bash scripts/launch_pair_class_profiles.sh run
```

This assumes the existing simulation/posterior inputs and native dependencies.
The launcher supports `SIM_BASELINE_DIR`, `SIM_LATER_DIR`,
`SIM_COMPARISON_DIR`, and `PAIR_CLASS_OUTPUT_DIR` path overrides. Twenty workers
are used, and existing RNG seeds are retained; this analysis uses no new RNG.
