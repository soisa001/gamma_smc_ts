# Provisional neutral-model CDFs, 2026-10-05

Verified completed receipts contain 161 regions: AFR 5, AMR 11, EAS 119, EUR 10,
MID 9, SAS 7. Simulations remain paused. These tables cover the focal site only,
with archaic-allele segregation required immediately after the 50-kya pulse
and observation required today. Each population is below 1,000 planned nulls.

`four_cutoffs_decoded_all_pairs.csv` supplies the four descriptive percentages
per population: 95th/99th percentiles of the raw recent-pair fraction at
T=10,000/50,000 years. Values use the inverse empirical CDF (nearest rank),
not smoothed or extrapolated tails. For every population except EAS, both
percentiles equal the observed maximum because of the small sample size.

`cutoffs_all_methods.csv` also covers true TMRCA and raw ALT/ALT fractions.
It distinguishes empirical percentiles from exact conservative empirical-p
boundaries: p=(1+#null>=target)/(n+1), reject when the finite target strictly
exceeds the boundary. Sample sizes below 19/99 cannot attain p<=0.05/0.01.
The file flags these unsupported thresholds and bounded-score saturation.
Missing ALT/ALT scores are omitted from CDFs, but remain in the testing
denominator under the existing no-call convention. Counts are supplied.

The completed subset can favor faster trajectories. Neither the descriptive
tails nor the fitted boundaries establish the final null calibration.

Local figures and bundle:
`sim/all_population_neutral_20261001/cdf_20261005/`. There are eight overview
plots (two time thresholds, two pair statistics, two TMRCA sources), 24
population-specific two-panel pages, PNG/vector PDFs, and one combined PDF.
Source receipts and both score-file hashes were checked before plotting.

Reproduce without restarting simulations, on the configured analysis host:

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
~/.local/bin/uv --no-config run --no-project --python .venv-origin/bin/python \
  python scripts/plot_population_null_cdfs.py \
  --root /mnt/d/phase2simselection/sim/all_population_neutral_20261001 \
  --out /mnt/d/phase2simselection/sim/all_population_neutral_20261001/cdf_20261005
```
