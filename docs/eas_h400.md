# Increasing the EAS panel to 400 haplotypes

This experiment changes the sample from 100 to 200 diploid individuals, giving
400 distinct haplotypes. It retains 10,000 distinct sampled pairs so that the
comparison measures increasing the number of haplotypes at a fixed pair budget.
The focused arm is s=0.005, with 100 selected and 1,000 neutral simulations.
All sample-fixed selected alleles remain in the experiment.

The configuration is `configs/eas_q02_50k_h400.json`. All other scientific
conditions match the fresh 200-haplotype study: EAS, 2% introgression and immediate
selection at 50 kya, fixed mutation and recombination rates, simulation of 11 Mb
and site-aligned cropping to an isolated 10 Mb region. The primary statistic
remains `frac_recent_T` based on posterior mean TMRCA. No earlier handoff
simulation data are used.

The same fixed seed convention is retained. The larger sample is simulated
anew; these are not nested samples of the same saved genealogy. Original
200-haplotype data remain available for comparison. The new output root is
`D:/phase2simselection/sim/eas_q02_h400`. Resource allowances remain 20 workers,
300 GB RAM subject to availability, and 3 TB total simulation storage on D.

## Analysis

The existing paired truth/decoded grid evaluates the same 540 threshold, stride,
gap, and run rules per statistic source, using a matched null for each panel.
It retains inclusive p<=alpha and removes each neutral replicate from its own
reference. The p=0.001 tail remains resolution-limited with 1,000 neutrals.

`panel_size_comparison.py` compares focal neutral standard deviation, selected
minus neutral mean, average profile roughness, and decoded-versus-true error.
Roughness includes real genealogical variation and must not alone be interpreted
as decoder noise. It also evaluates whole-region maximum `frac_recent_T`, with
and without normalization at each position, separately for the six time cutoffs.

Five folds split complete independent replicates within each panel and arm.
Each fold reserves 200 neutral and 20 selected regions for evaluation. Of the
800 remaining neutrals, 400 fit the normalization and 400 calibrate the null
distribution of the region maximum. These groups do not overlap. This avoids
using the same neutral observations to fit normalization and calibrate its
maximum. Each cutoff and score is evaluated at regional nominal alpha 0.01,
0.05, and 0.10. A separate threshold targets at least 70% sensitivity in the 80
training selected regions. Fixed samples are retained in training and testing.

Every primary prediction is whether there is at least one qualifying call
anywhere in the region. No distance-to-allele criterion is used. No method is
selected using held-out outcomes; comparing or choosing the best reported rule
is exploratory. These folds are recomputed on the focused s=0.005 cohort, so
individual results need not exactly equal the previous all-arm cross-validation.

## Invocation

Use the existing pinned Linux/WSL uv environment and native dependencies described
in [the fresh-study setup](fresh_eas_power.md). The baseline raw simulations must
already be available locally. The scripts log parameters and phase receipts and
validate saved simulation and profile hashes. Interrupted runs resume completed
work rather than regenerate it. Analysis never starts before simulation success
when using the full `run` phase.

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
bash scripts/launch_eas_h400.sh smoke
bash scripts/launch_eas_h400.sh run
```

The `analyse` phase recomputes only comparisons from completed profiles; `status`
validates completed simulation receipts. Detailed outputs stay on D; compact
tables, provenance and PNG/vector PDF figures can be staged in the repository.
