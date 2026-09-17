**Results.** Weighting recent coalescence by archaic-allele carrier support
substantially improved detection of selected regions relative to the all-pair
TMRCA scan in the current EAS benchmark. The final comparison comprised 100
selected and 1,000 neutral isolated 10 Mb regions, with a nominal 2% archaic
introgression pulse at 50 kya, selection beginning immediately at the pulse,
and selection coefficient s=0.005. Each region contained 400 sampled haplotypes
and was evaluated using the same panel of 10,000 haplotype pairs. Selected
alleles were required to survive and be observed in the sample; every accepted
selected replicate was retained. None of these 100 alleles was fixed in the
sample. Their frequencies ranged from 25.75% to 99.25%, with mean 78.39% and
median 82.13%. The benchmark therefore represents incomplete, predominantly
high-frequency sweeps. These results do not establish performance across the
earlier s=0.001-0.010 grid with the final carrier-based detector.

The archaic ALT/ALT score measures the fraction of all sampled pairs that both
carry a candidate archaic allele and have recent TMRCA. At the primary time
cutoff of 50 kya and a calibrated regional p<=0.05, the decoded score detected
83/100 selected regions on the 10 kb grid, while calling 48/1,000 neutral regions.
Descriptive 95% Wilson intervals were 74.5-89.1% for power and 3.64-6.31% for
the neutral call fraction. Scanning every eligible archaic site yielded similar
performance: 84/100 selected and 49/1,000 neutral regions called. Thus scanning
every archaic site did not provide a clear improvement over the 10 kb grid in
this cohort.

| Statistic at T=50 kya | Candidate positions | Selected regions called | Neutral regions called |
|---|---|---:|---:|
| All-pair recent fraction | 10 kb grid | 14/100 (14%) | 48/1,000 (4.8%) |
| Archaic ALT/ALT score | 10 kb grid | 83/100 (83%) | 48/1,000 (4.8%) |
| Archaic allele frequency alone | Grid-assigned archaic markers | 84/100 (84%) | 49/1,000 (4.9%) |
| All-pair recent fraction | Every archaic site | 33/100 (33%) | 49/1,000 (4.9%) |
| Archaic ALT/ALT score | Every archaic site | 84/100 (84%) | 49/1,000 (4.9%) |
| Archaic allele frequency alone | Every archaic site | 85/100 (85%) | 52/1,000 (5.2%) |

All rows use separately calibrated regional maxima. The allele-frequency
benchmark has no TMRCA cutoff. Restricting the all-pair scan to archaic sites
changes both the candidate set and its placement; it is not simply a finer
version of the uniform grid.

The difference between true and decoded TMRCA depended on the statistic. On
the grid, the all-pair score at T=50 kya had 41% power with true TMRCA but only
14% with decoded TMRCA, at neutral call fractions of 5.2% and 4.8%, respectively.
For the archaic ALT/ALT score, true and decoded power were both 83%, with
neutral call fractions of 5.0% and 4.8%. At exact archaic sites, true and decoded
ALT/ALT-score power were 81% and 84%, respectively. This agreement concerns
calibrated detection performance; it does not establish accurate individual
TMRCA estimates. Across the six tested time cutoffs from 5 to 50 kya, decoded
ALT/ALT-score power remained 79-84% on the grid and 79-87% at exact sites, with
neutral call fractions of 4.1-5.1%. Decoded all-pair grid power was 12-18% across
the same cutoffs. These cutoff comparisons were exploratory, and the primary
reported cutoff remained 50 kya.

At thresholds trained to target approximately 70% power, the archaic ALT/ALT
score called 71/100 selected and 13/1,000 neutral regions on the grid, or 69/100
and 14/1,000 at exact sites. The grid estimates had descriptive 95% Wilson
intervals of 61.5-79.0% for power and 0.76-2.21% for the neutral call fraction.
Thresholds were determined from selected training regions and were not adjusted
to force 70% power in the held-out regions. All-pair thresholds trained toward
the same power target called many more neutral regions.

| Decoded statistic, T=50 kya where applicable | 10 kb grid: power / neutral call fraction | Every archaic site: power / neutral call fraction |
|---|---:|---:|
| All-pair recent fraction | 71% / 58.8% | 70% / 28.0% |
| Archaic ALT/ALT score | 71% / 1.3% | 69% / 1.4% |
| Frequency-weighted ALT/ALT minus REF/REF | 72% / 2.0% | 72% / 2.4% |
| Archaic allele frequency alone | 70% / 1.4% | 71% / 1.4% |

Carrier support was essential to the successful ALT/ALT statistic. At exact
sites and T=50 kya, the unweighted within-carrier recent fraction reached one
in 970/1,000 decoded neutral regions and 99/100 selected regions. Across
regions, the median minimum number of ALT/ALT pairs supporting a regional
maximum was one. A rare neutral marker could therefore attain the same
unweighted score as a selected allele carried by many haplotypes. With
conservative treatment of these ties, the unweighted regional ALT/ALT statistic
made no calls at p<=0.05. Weighting by n_ALT/ALT divided by the total pair count
distinguished sparse support from evidence involving a substantial fraction of
the sampled panel.

Neither REF/REF subtraction nor an ALT/REF penalty provided a demonstrated
advantage over the simpler carrier score. At regional p<=0.05, removing REF/REF
subtraction changed decoded grid power from 81% to 83% without changing the
4.8% neutral call fraction. Near the 70% power target, the simpler score removed
seven neutral grid calls and ten exact-site calls, without adding new neutral
calls, while power decreased by one and three percentage points. The exact-site
neutral comparison had a nominal paired p=0.00195 and Holm-adjusted p=0.0313
across the 16 comparisons in that analysis; the grid comparison had nominal
p=0.0156 and adjusted p=0.234. These conditional comparisons are exploratory
because thresholds are shared across cross-validation folds and the same cohort
informed development of the score. They do not establish superiority on new
simulation draws or equivalence of power.

The allele-frequency benchmark was similarly effective to the ALT/ALT score
at both operating points. Consequently, the large improvement over all-pair
TMRCA detection cannot yet be attributed to genealogical recency beyond the
information supplied by archaic carrier frequency. The current evidence supports
the simpler ALT/ALT score as a practical detector under the specified model,
while leaving its incremental benefit over allele frequency unresolved.

The neutral call fraction is a regional false-positive rate, not the proportion
of discoveries that are false. For example, the grid operating point with 71%
power and a 1.3% neutral call fraction implies a discovery FDR of 14.1% and F1
of 0.777 if 10% of tested regions are selected. At 1% prevalence, the same rates
imply FDR=64.4% and F1=0.474. These are illustrative calculations under assumed
prevalence, not estimates for an empirical genome-wide scan. The present results
are conditional on allele survival and observation, oracle archaic markers,
fixed rates, and isolated 10 Mb regions. They require confirmation on new seeds
after freezing the final method.

**Methods.** We simulated an EAS recipient population using a piecewise-constant
effective-population-size history derived from the pointwise median of the
tracked PHLASH EAS bootstrap fits. The first supported estimate, at 100
generations, was extended to the present; the last supported estimate was
extended into the older tail. Both selected and neutral models included an
archaic lineage splitting 700 kya, an archaic effective population size of 3,600,
and an imposed bottleneck to size 10 for the 100 generations immediately before
the 50-kya pulse. This bottleneck is a model assumption of the benchmark. The
nominal pulse proportion was 0.02; realized ancestry was not conditioned to a
frequency band. Mutation and recombination rates were fixed at 1.25e-8 and
1e-8 per base per generation, respectively, with 25 years per generation.

Neutral regions were generated with msprime 1.4.2. Selected regions were
generated with SLiM 4.2.2 through stdpopsim 0.3.0, using rescaling factor five
and a burn-in parameter of 0.1. A focal allele was placed at fixation in the
archaic source after the split. Selection began in EAS at introgression and
continued to the present, with unscaled selection coefficient s=0.005 and
dominance h=0.5. Selected trajectories were conditioned on population survival
of the allele, and sampled panels lacking the allele were rejected and retried
using deterministic attempt seeds. Accepted replicates were not filtered by
allele frequency or by detection outcome; fixation was allowed. This estimates
detection conditional on an observable surviving allele, not the unconditional
probability that a new introgressed allele survives and is detected. The use of
different simulation engines for neutral and selected regions remains an
assumption requiring a matched-engine validation.

Each simulation covered 11 Mb and supplied a 10 Mb crop. In selected replicates,
the crop placed the selected allele at 5 Mb. In neutral replicates, the crop
placed the observed segregating biallelic variant nearest the original center
at 5 Mb, without a minor-allele-frequency threshold. This neutral crop anchor
was not required to be archaic. We sampled 200 diploid individuals, yielding
400 haplotypes, and used 10,000 unique unordered pairs drawn without replacement
from all possible haplotype-index pairs. The fixed pair-selection seed was
1729. Simulation seeds were deterministically derived from study seed 20380101,
replicate identity, simulation component, and retry attempt. Overlapping pairs
were not treated as independent observations for significance testing.

Archaic markers were identified using simulation ancestry truth. An eligible
marker had one mutation, was biallelic with observed ALT, had mutation age from
the pulse to before the archaic split, and had archaic pulse ancestry in every
ALT-carrying haplotype at that position. The same criteria were applied to both
cohorts. Recurrent or multiallelic histories were excluded, and fixed ALT
markers were retained. Ancestry was recovered by deterministic replays that
were checked against the original genealogies and sample identities; the
original mutations and decoder inputs were retained. The blind regional scan
did not use the identity of the selected allele. Regions without eligible
markers remained in the denominator and received score zero; there were 16
such neutral regions.

True pairwise TMRCA was extracted from the local trees. Decoded TMRCA was
obtained from the local Gamma-SMC implementation using all original input
variants, with output at the union of grid and eligible-marker positions.
The decoder retained its embedded constant-size flow field, scaled mutation
rate theta=0.00075, and supplied rho/theta=0.8. These nominal units correspond
to a diploid reference N0=15,000 and supplied rho=0.0006. The simulation's
variable demographic history was not incorporated into this flow field.
The existing accurate exponential calculation and corrected backward-output
alignment were retained. True and decoded statistics were calibrated against
their own corresponding neutral distributions. This calibration evaluates
the decoder's statistic under the simulated null; it does not establish
absolute age accuracy or recover information lost during decoding.

We used only hard pairwise recent-coalescence calls. For true TMRCA, a pair was
recent when its age was strictly below T. For decoded TMRCA, we applied the
native posterior-mean threshold convention, including its float32 shape
clipping and inclusive equality rule. Specifically, in scaled units the test
was beta*T >= clip(alpha, 2^-16, 2^16). The statistic was `frac_recent_T`, not
the average posterior probability of TMRCA below T. The evaluated cutoffs were
5, 10, 20, 30, 40 and 50 kya; 50 kya is the current primary reporting cutoff.
Decoded class counts were checked against the native output, and true grid
counts were checked against the original tree-based profiles.

For marker j, let r_(k,T) indicate a recent-coalescence call for pair k, let
A_j be the set of sampled pairs whose two haplotypes both carry ALT at j,
and let N=10,000 be the total number of sampled pairs. We defined

\[
f_{AA,T}(j)=\frac{\sum_{k\in A_j}r_{k,T}}{|A_j|},\qquad
S_T(j)=\frac{|A_j|}{N}f_{AA,T}(j)
=\frac{1}{N}\sum_{k\in A_j}r_{k,T}.
\]

Thus ALT/ALT denotes two arbitrary carrier haplotypes, not a diploid homozygote.
The score was zero when no sampled ALT/ALT pairs existed. At fixation, the
score reduced to the all-pair recent fraction. The final score contained no
REF/REF subtraction, ALT/REF penalty, or additional all-pair gate. The all-pair
benchmark used the recent fraction across all N pairs. The allele-frequency
benchmark used the maximum candidate archaic ALT frequency in the region.

The primary scan evaluated TMRCA at 0, 10,000, ..., 9,990,000 bp. At each grid
position, carrier membership was supplied by the nearest eligible archaic
marker within 5 kb, resolving equal distances toward the left. Positions
without an assigned marker supplied no carrier score, while the all-pair
baseline could use every grid position. A secondary scan evaluated each
eligible archaic marker at its exact coordinate. For each method and cutoff,
the regional statistic was the maximum over its candidate positions. One
position sufficed; the final score did not require consecutive significant
strides. The full positional search was therefore included before calibration.

We retained five fixed whole-region folds. In each fold, 200 neutral and 20
selected regions were held out, 400 separate neutral regions supplied
calibration scores, and 80 selected regions were available for training.
The remaining 400 neutral regions supported gate-fitting comparisons; the
final ungated score had no gate to fit. Each region was tested once. For a
test-region score R and the 400 neutral calibration scores R_b, we calculated

\[
p=\frac{1+\sum_{b=1}^{400}\mathbf{1}(R_b\ge R)}{401}.
\]

Ties were included conservatively. We evaluated p<=0.01, 0.05 and 0.10,
with p<=0.05 as the primary nominal threshold. Cutoffs were evaluated separately
for the fixed methods; the reported p-values did not represent an uncalibrated
maximum across time cutoffs. The minimum obtainable held-out rank p was 1/401,
so this design cannot validate a p=0.001 operating point. The final formula was
chosen after exploring this cohort and was not prospectively validated merely
by retaining these folds.

For the approximately 70%-power operating point, each fold's threshold was
the lower 30th percentile of its 80 selected training-region scores, computed
separately for each method, position scheme and TMRCA source at T=50 kya.
Held-out regions were called when their score was at least this threshold and
strictly positive. This thresholding did not impose a nominal p-value or a
specified false-positive rate. Held-out power and neutral calls were then
measured without further adjustment.

Power was the proportion of selected 10 Mb regions called anywhere in the
region. Regional false-positive rate Q was the proportion of neutral regions
called anywhere. Detection was not required to lie within a distance window
around the selected SNP. False-positive labels were not assigned according to
whether a called peak included the selected haplotype. A separate oracle
focal-SNP diagnostic was not pooled with these blind-scan results. If P denotes
power and pi the assumed prevalence of selected regions, illustrative
discovery FDR and F1 were calculated as

\[
\mathrm{FDR}(\pi)=\frac{(1-\pi)Q}{\pi P+(1-\pi)Q},\qquad
F_1(\pi)=\frac{2\pi P}{\pi(1+P)+(1-\pi)Q}.
\]

Wilson intervals described count uncertainty conditional on the saved
decisions; they did not propagate shared threshold-fitting uncertainty.
Carrier mass and the REF/REF contrast were compared by two-sided exact
McNemar calculations using discordant region calls, with Holm adjustment
across the 16 displayed comparisons. These are conditional exploratory
diagnostics, since fitted thresholds are shared and earlier method searches
are outside that adjustment.

The workflow used 20 workers, a 300 GB RAM allowance, deterministic seeds,
and manifest- and hash-verified caching. No simulations or decoding were rerun
for this synthesis. Input parameters, pair identities, fold assignments,
calibration scores, predictions, and integrity checks are retained on D and in
the repository's compact result exports. The variable-Ne decoder discussed
in the research note has not been implemented or evaluated for these results.
The separate audit's transition-rate convention discrepancy also remains to
be resolved in an independent decoder benchmark; no rate correction was
silently applied to this cohort.

**Supporting records and reproduction.** The principal data sources are the
[regional metrics](results/allele_class_ablation/metrics.csv),
[70%-target metrics](results/allele_class_ablation/target70_metrics.csv),
[saturation diagnostics](results/allele_class_ablation/saturation.csv),
[focal allele-frequency inventory](results/allele_focal_comparison/regions.csv),
and [paired REF/REF comparison](results/alt_alt_ref_comparison/report.md).
The corresponding five source data files were checked against their artifact
manifests when preparing this text. The current
[score preset](../configs/archaic_alt_alt_signal.json),
[simulation configuration](../configs/eas_q02_50k_h400.json),
[simulation implementation](../python/gamma_smc_aou/fresh_power.py),
[demographic construction](../python/gamma_smc_aou/eas_sweep_models.py),
[ancestry-label implementation](../python/gamma_smc_aou/joint_scan_labels.py),
[profile extraction](../python/gamma_smc_aou/joint_scan_profiles.py),
[regional calibration](../python/gamma_smc_aou/allele_class_ablation.py),
and [flow-field audit](variable_ne_flow_field.md) specify implementation details.

The following WSL cell updates a dedicated checkout and regenerates the compact
current-score summary and paired comparison from the existing D-drive cache.
It does not simulate or decode new regions.

```bash
%%bash
set -euo pipefail
repo=gamma_smc_results_writeup
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
cd "$repo"
git pull --ff-only
uv --no-config run --no-project python scripts/summarize_archaic_alt_alt_signal.py \
  --source /mnt/d/phase2simselection/sim/eas_allele_class_ablation \
  --out /mnt/d/phase2simselection/sim/eas_results_writeup/current_score
uv --no-config run --no-project --with scipy python scripts/compare_alt_alt_ref_subtraction.py \
  --source /mnt/d/phase2simselection/sim/eas_allele_class_ablation \
  --out /mnt/d/phase2simselection/sim/eas_results_writeup/ref_comparison
```
