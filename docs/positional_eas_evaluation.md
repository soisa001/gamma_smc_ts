**Primary endpoint updated to a matched position or gene window.**

**Completed-cohort update (21 September 2026):** the tables below are the
earlier 947-selected snapshot. See the [lab meeting figure collection](lab_meeting_figures.md)
for the completed 1,000-selected cohort, refreshed iHS, and presentation figures.

The user's clarified question is whether a particular position or gene has a
signal unusual under neutrality at that same position or within that same
gene-sized window. A peak elsewhere in the 10-Mb neutral region is irrelevant
to this endpoint. The earlier regional-maxima results remain valid for their
different, whole-region question; they are not the primary power estimates
for this positional analysis.

The saved EAS simulations use 400 haplotypes, 10,000 fixed sampled pairs,
a 2% introgression pulse and immediate selection onset at 50 kya, generation
time 25 years, mutation rate 1.25e-8 and recombination rate 1e-8 per base per
generation. The carrier score is `n_recent_ALT_ALT / n_all_pairs`, equivalently
ALT/ALT `frac_recent_T` times the ALT/ALT pair fraction. Decoded counts use
the existing posterior-mean TMRCA hard calls; no mean posterior-mass statistic
is used.

At the focal position, only **43/1,000 neutrals (4.3%)** have an assigned
eligible archaic marker at AF >=18%, versus 909/1,000 with such a marker
somewhere in 10 Mb. The mean assigned-marker AF at the focal position is
2.066%, counting absent markers as zero. Thus the earlier 90.9% figure does
not answer the user's local question.

For the carrier/AF test, evaluate TMRCA at 5 Mb on the original 10-kb grid.
Assign the nearest eligible archaic marker within 5 kb, as in the existing
scan. All 1,000 neutral positions remain in the null; a position without a
nearby archaic marker has carrier and AF score zero. There are 169/1,000
neutral positions with an assigned archaic marker. This is an unconditional
FPR across pre-specified genomic positions, not an FPR conditional on already
having selected an observed archaic SNP. All-pair TMRCA remains defined at
every position. The selected allele lies exactly at 5 Mb in every selected
region. Selected simulations continue to be conditioned on allele survival
and observation, including fixation.

In each of the five existing folds, compare the test-position score with the
same position in 400 independent neutral calibration regions. Test 200 other
neutral positions per fold; each of the 1,000 neutrals is tested once.
Use `(1 + count(calibration_score >= test_score)) / 401`. No maximum over
positions, chromosomes or T cutoffs enters these p-values. The minimum
resolvable rank is 1/401; p=0.001 cannot be evaluated with this split.

For a pre-specified gene interval, any chosen summary (including a maximum)
should instead be calibrated against that same summary in a matching neutral
interval. A peak outside the interval never enters that calibration. Testing
each position separately and calling a gene if any position passes would be
a different rule and needs its own gene-level calibration. Neither endpoint
by itself estimates FDR among discoveries; the reported quantity is FPR under
neutrality.

**True-TMRCA results, T=50 kya, nominal p<=0.05.**

| s | Selected positions | AF power | ALT/ALT carrier-mass power | All-pair power |
|---|---:|---:|---:|---:|
| 0.001 | 92 | 42.4% | 46.7% | 5.4% |
| 0.002 | 93 | 71.0% | 71.0% | 19.4% |
| 0.003 | 95 | 78.9% | 78.9% | 40.0% |
| 0.004 | 94 | 92.6% | 92.6% | 60.6% |
| 0.005 | 100 | 100.0% | 100.0% | 86.0% |
| 0.006 | 94 | 97.9% | 97.9% | 90.4% |
| 0.007 | 94 | 98.9% | 98.9% | 93.6% |
| 0.008 | 95 | 97.9% | 97.9% | 95.8% |
| 0.009 | 95 | 100.0% | 100.0% | 98.9% |
| 0.010 | 95 | 100.0% | 100.0% | 100.0% |

Held-out positional FPR is 5.1% for AF, 5.0% for carrier mass, and 5.4% for
all-pair recency. At nominal p<=0.01, carrier-mass power for s=0.001 through
0.005 is 17.4%, 45.2%, 61.1%, 79.8%, and 99.0%, with 0.9% positional FPR.
AF power is 15.2%, 43.0%, 60.0%, 78.7%, and 99.0%, also with 0.9% FPR.
This change in calibration recovers substantial weak-selection power.
Carrier mass and AF remain similar under the current model at T=50 kya.

Shorter recency cutoffs provide a larger exploratory gain for the weakest
selection arm. At T=20 kya and p<=0.05, carrier-mass power is 54.3%, 73.1%,
85.3%, 95.7%, and 100% for s=0.001 through 0.005, at 5.1% positional FPR.
For s=0.001, this is 50/92 detections versus 39/92 for AF: twelve additional
detections and one loss (paired exact p=0.00342, unadjusted across the time
cutoffs and selection arms examined). At T=50 kya the paired difference is
five gains and one loss (p=0.219). All six cutoffs are reported separately;
there is no per-replicate search over T in the p-value calculation.

For the existing decoded s=0.005 arm only, carrier mass has 100% power and
5.1% positional FPR at p<=0.05. The all-pair decoded statistic has 61% power
and 5.0% positional FPR, compared with 86% true all-pair power. No additional
Gamma-SMC decoding was run. The other selection arms have true-TMRCA results
only in this comparison. All six saved T cutoffs are retained in the tables.

As a separate p<=0.001 sensitivity analysis, pool all 1,000 neutrals for
selected-position ranks and compare each neutral with the other 999. These
scores require no fitted normalization. True carrier-mass power at T=50 kya
is 2.2%, 16.1%, 38.9%, 59.6%, and 86.0% for s=0.001 through 0.005; AF power
is 2.2%, 18.3%, 41.1%, 62.8%, and 90.0%. The decoded s=0.005 carrier result
is 92%. The leave-one-out neutral call rate is 1/1,000 for these scores, but
only one extreme rank is available: this is not precise validation of a
0.1% tail. These pooled results are separate from the five-fold estimates.

This remains a preliminary analysis of the 947 saved selected regions;
53 planned selected simulations are unfinished. Simulation generation is
paused. The implementation reproduces the old s=0.005 focal truth values
and independently recomputes all saved p-values from the calibration vectors.

**Signal as a function of distance.** The already decoded s=0.005 cohort also
supports pointwise evaluation at every 10-kb stride across the entire 10-Mb
crop. At each stride, the null uses that same coordinate in the neutral
replicates; the computation never searches for an external maximum. At
T=50 kya and p<=0.05, decoded carrier-mass detection is:

| Distance from selected allele | Detection | Neutral positional FPR |
|---|---:|---:|
| -200 kb | 31% | 4.9% |
| -100 kb | 60% | 4.4% |
| -50 kb | 77% | 4.7% |
| 0 | 100% | 5.1% |
| +50 kb | 82% | 4.6% |
| +100 kb | 64% | 4.7% |
| +200 kb | 46% | 5.3% |

These linked-position detections describe the spatial extent of the signal,
not separate causal variants or localization false discoveries. They also
do not give the probability of at least one call in a gene. The focal counts
exactly reproduce the independent focal evaluation for both sources and all
thirteen AF/TMRCA statistics at both significance cutoffs. The full table is
saved under `eas_positional_distance_h400/metrics.csv`; selected distances
are exported to [`results/positional_distance_h400`](results/positional_distance_h400).

**What counts as an archaic marker.** Empirically, the intended ascertainment
is ALT present in Neanderthal or Denisovan and absent from an African
outgroup. The simulation uses an oracle origin rule, rather than sampling
these reference genomes: a single biallelic mutation, observed ALT, mutation
age at least the pulse age and younger than the archaic-modern split, and
all ALT carriers inheriting archaic pulse ancestry at that site. Between the
split and pulse the modeled populations are isolated. Therefore a retained
post-split mutation on the archaic ancestral route was not segregating in
the modern population before introgression. Merely sitting on an introgressed
segment would not establish that fact. This conservative rule excludes all
pre-split mutations, including some older variants that could have been lost
from the modern population before introgression.

The empirical Neanderthal/Denisovan-present, African-absent rule is an
approximation to this origin criterion: ancestral polymorphism and finite
reference sampling can alter the retained set. A literal empirical
ascertainment benchmark would require reference/outgroup sampling assumptions
that are not part of the current EAS-only array.

For a truly introgression-specific marker, `marker_AF <= local_archaic_AF`.
The genome-wide or region-average pulse proportion is not an upper bound on
local ancestry. The previously reported high neutral AF peaks were already
restricted to the oracle marker set; they were not caused by including
pre-split modern/archaic shared mutations. Those peaks outside the focal
position are excluded from the present calibration.

**Completed iHS comparison.** Raw extraction finished for all 1,947 saved
regions in 36.9 minutes with 20 workers. Every cached raw output passed its
hash check, and the analysis independently reconstructed the saved p-values.
At p<=0.05:

| s | Positional iHS power | Fixed 100-kb window iHS power |
|---|---:|---:|
| 0.001 | 10.9% | 8.7% |
| 0.002 | 12.9% | 19.4% |
| 0.003 | 32.6% | 36.8% |
| 0.004 | 53.2% | 58.5% |
| 0.005 | 67.0% | 81.0% |
| 0.006 | 64.9% | 71.3% |
| 0.007 | 43.6% | 75.5% |
| 0.008 | 32.6% | 63.2% |
| 0.009 | 27.4% | 53.7% |
| 0.010 | 14.7% | 47.4% |

Held-out local FPR is 4.8% for positional iHS and 5.0% for the fixed window.
At p<=0.01, positional power for s=0.001 through 0.005 is 1.1%, 3.2%, 18.9%,
30.9%, and 50.0%, at 0.7% FPR. Window power is 3.3%, 6.5%, 23.2%, 37.2%, and
58.0%, at 1.3% FPR. Thus conventional iHS does not outperform local archaic
AF or carrier mass under this model. The window and pointwise columns are
different endpoints and should be described separately.

Strong-arm iHS power declines as fixation and very high target AF become
common. For s=0.010, 46/95 selected alleles are fixed and only 2/95 focal
alleles themselves meet the finite-score/core-MAF rule; the positional method
can still use a nearby scorable SNP. These diagnostics are consistent with
iHS being an incomplete-sweep statistic, but are not a decomposition of all
causes of its power loss. Fixed/undefined cases were retained in every power
denominator. The checked iHS tables and audit are in
[`results/ihs_eas_h400`](results/ihs_eas_h400).

**iHS methods.** Compute conventional iHS from all biallelic
segregating variants across each saved 11-Mb tree sequence, with the 10-Mb
crop defining scored positions and the flanks supporting EHH integration.
Use the known uniform recombination rate, minimum core MAF 5%, EHH cutoff
0.05, and no missing-data gap rescaling because the sequences are fully
observed. Scores whose EHH does not decay before the simulated edge are
undefined. Raw scores are `ln(iHH_ALT/iHH_REF)`. Retain rare flanking variants
for haplotype construction; the MAF filter applies only to core scoring.

Standardize raw iHS in twenty derived-frequency bins using 400 neutral fit
regions, separate from the 400 calibration and 200 held-out neutral regions
in each fold. Selected regions never fit normalization. The positional
statistic is absolute standardized iHS at the nearest scorable core within
5 kb of 5 Mb, with ties resolved to the left; an undefined position scores
zero. The gene-window statistic is the fraction of scorable cores with
absolute iHS > 2 in the fixed 100-kb window centered at 5 Mb, requiring at
least twenty scorable cores. Each statistic is calibrated against its
identical local counterpart in neutral regions. No external peak enters
the p-value. These are separate endpoints, not a combined search over
window definitions. iHS does not require archaic marker annotation for this
conventional benchmark.

The implementation follows the frequency-standardized EHH contrast described
in [Voight et al. (2006)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072),
using the raw-score sign and integration options documented by
[scikit-allel](https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/selection.html).

The iHS environment is isolated from the simulation/decoder runtime and
pinned in `configs/requirements-ihs.txt`. Three focused tests verify direct
EHH integration and score sign, genetic-map scaling invariance, frequency
bin arithmetic, and local/window score behavior. A five-region extraction
smoke run passed, including a selected replicate with a fixed focal allele.

**Reproduce without resuming simulations.**

The following WSL cell uses the saved studies and existing simulation
environment on this machine. Raw tree sequences are not stored in Git.

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
cd "$repo"
git pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
"${HOME}/.local/bin/uv" --no-config run --no-project --python .venv/bin/python \
  python scripts/evaluate_positional_array.py
"${HOME}/.local/bin/uv" --no-config run --no-project --python .venv/bin/python \
  python scripts/summarize_positional_evaluation.py
bash scripts/launch_ihs_benchmark.sh run
"${HOME}/.local/bin/uv" --no-config run --no-project --python .venv-ihs/bin/python \
  python scripts/evaluate_positional_distance.py
"${HOME}/.local/bin/uv" --no-config run --no-project --python .venv-ihs/bin/python \
  python scripts/export_positional_summary.py --analysis all
```

Positional carrier/AF outputs are in
`/mnt/d/phase2simselection/sim/eas_positional_h400`; iHS intermediates and
results are in `/mnt/d/phase2simselection/sim/eas_ihs_h400`. Raw iHS inputs
and cached outputs are hash-verified. Cache identities depend on the saved
tree, relevant parameters, implementation and library versions, independently
of unrelated Git changes. Neither launcher invokes simulation generation or
Gamma-SMC decoding.

Checked compact tables are exported to
[`docs/results/positional_eas_evaluation`](results/positional_eas_evaluation).
`metrics.csv` holds the five-fold local results; `paired_af_comparison.csv`
holds paired detection gains/losses; `pooled_rank_metrics.csv` holds the
separate p=0.001 sensitivity analysis. Large raw trees and per-site scores
remain in the simulation directory.
