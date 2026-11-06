# One selected simulation versus 100 neutrals

Open [the two-page PDF](selected_vs_100_neutrals.pdf): posterior-median and posterior-mass bin curves, then pointwise versus scan-wide p-values.

Example: `10k/introgression/s0.01/selected_00020`, position 420000 bp (selected allele at 500000 bp), highlighted mass bin [1, 2] kya. Present AF 57.14%. Selection s=0.01 from 10 kya to present; EAS demography, 1% Neanderthal pulse at 50 kya. Both methods use the same selected replicate, genomic position, and 1000 pairs. Each panel highlights its own maximum standardized score within 1–10 kya at that position; highlighted bins may differ.

The example is the first previously detected posterior-mass replicate in frozen order. The position/bin maximizes its standardized mass score within ±100 kb and 1–10 kya. This is an illustration chosen to show a call, not an unbiased performance estimate. The original scan correction covers positions and times within each region, not selection among 100 replicates for illustration.

The 100 displayed neutrals are all 50 held-out regions plus the first 50 calibration regions by replicate. They are distinct simulations, with input hashes verified; this combined set is not a new independent validation cohort. Both selected and neutral curves use the same genomic coordinate, not each neutral's nearest-segregating label.

For each time bin k is the number of neutral values **equal to or above** the selected value. Raw tail fraction k/100; conservative empirical p=(k+1)/101. Thus at most 4/100 gives p≤0.05; 5/100 gives p=0.0594. The dashed curve is the 96th ascending order statistic, which the selected curve must strictly exceed. It is a pointwise rank threshold, not a confidence band. Posterior-median curves count pairs with medians in the bin; posterior-mass curves sum pairwise posterior probabilities in the bin. Neither plot is a pairwise TMRCA credible interval.

The second page preserves the original scan-wide p-values from 99 independent calibration-region maxima, standardized using 20 development regions and searching 100 genomic positions × 49 time bins. No threshold was refit for the selected example. No 0–1 kya bin is evaluated. Detailed values are in `timewise_statistics.tsv` and `summary.json`.

## Reproduce from the verified summary cache

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
uv run --frozen python scripts/plot_selected_vs_100_neutrals.py --workers 4
```
