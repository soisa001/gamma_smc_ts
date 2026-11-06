# Selection at 40–20 kya, s=0.02: complete introgression cohort

Open [the four-page figure series](selection_episode_series.pdf). All pages are also saved individually as PNG and vector PDF; no inline plotting is used.

These are the first 100 present-day segregating survivors in the frozen proposal order, not a subset chosen for attractive trajectories or significant results. The corresponding de novo treatment has 21 survivors and is not included here. The 100 introgression survivors required 967 proposals: 743 lost the focal allele and 124 fixed it. The null comparison uses the 50 held-out neutral introgression regions, without focal-allele survival or AF matching. Both cohorts have the same 1% Neanderthal pulse at 50 kya and EAS demographic history.

## Figure guide

1. [Allele-frequency trajectories](01_allele_frequency.pdf): all 100 whole-population trajectories, the pointwise median, and AF distributions at introgression, selection onset, selection cessation, and today. Selection is shaded; trajectory time runs from past to present. The selected cohort is conditioned on present-day segregation.
2. [Selected versus neutral TMRCAs](02_selected_vs_neutral_tmrca.pdf): cumulative coalescence fractions and adjacent 1-kyr bin masses at the same fixed 500-kb coordinate in both cohorts. Lines are means; bands are pointwise 10th–90th percentiles across regions, not confidence or posterior credible intervals. Right-column y-scales vary by estimator. Mass is not conditioned on coalescence before 50 kya.
3. [Spatial and temporal excess](03_spatial_coalescence_excess.pdf): mean selected-minus-held-out-neutral coalescence mass across the 1-Mb region, sampled every 10 kb. All panels share the same color scale. Dashed lines mark 20 and 40 kya; dotted lines mark the ±100-kb detection tolerance. These are descriptive excesses, not significance maps.
4. [Called intervals for every survivor](04_called_intervals.pdf): each row is one replicate, ordered by present-day AF from low to high and aligned across methods. Blue marks significant support at any scanned position within ±100 kb. Grey marks the true selected period. This union can combine support from different positions and does not preserve spatial connected-component identity.

## Allele-frequency milestones

| Time | Median AF | 10th–90th percentiles |
|---|---:|---:|
| 50 kya: introgression | 1.00% | 1.00%–1.00% |
| 40 kya: selection starts | 2.17% | 0.56%–4.23% |
| 20 kya: selection ends | 96.79% | 87.83%–98.83% |
| Present | 95.75% | 84.24%–99.22% |

## Timing and detection

| Method | Overlapping near-site call / 100 | Only mistimed / 100 | Any wholly outside component / 100 | Center temporal F1 | Mean focal excess peak |
|---|---:|---:|---:|---:|---:|
| True TMRCAs · 19,900 pairs | 73 | 23 | 46 | 0.070 | 38.5 kya |
| Posterior median · 1,000 pairs | 2 | 12 | 12 | 0.000 | 11.5 kya |
| Posterior mass · 1,000 pairs | 51 | 36 | 47 | 0.049 | 2.5 kya |

The true mean focal excess peaks at 38.5 kya, inside the simulated episode. The inferred profiles put substantial mass at younger times; their mean excess peaks are 11.5 kya for the median and 2.5 kya for posterior mass. These are descriptive maxima of the mean selected-minus-neutral curve, not fitted onset/end dates or per-replicate significance calls. True TMRCAs use all 19,900 pairs; inference retains the original fixed 1,000 pairs, so the panels do not isolate inference error from pair subsampling.

The primary scan covers [1,2), …, [49,50) kya; [0,1) is absent. Each method uses 20 independent development regions and 99 calibration regions, with conservative +1 empirical p-values against maxima over all 100 genomic positions and all 49 intervals. The 50 plotted neutral regions are held out from calibration. No thresholds were refitted to the selected cohort or the plotted null cohort. All 450 region/method metric records reproduce the frozen reports.

A region can have both a correctly overlapping component and a wholly outside component. “Only mistimed” counts regions with near-site calls but no temporal overlap; “any wholly outside component” also includes errors alongside a correct call. Center temporal F1 refers to the exact focal coordinate, whereas the plotted call union uses ±100 kb.

## Reproduce the figures without rerunning simulations

The checked cache contains the plotted CDFs and all 100 AF trajectories. Cache validity depends on the source report, cohort identities, and file hashes, not Git commit identity. On a fresh clone, the saved cache supports rerendering without the large raw trees or posterior files. Initial cache construction requires the raw paths recorded in the cohort manifest. There are no RNG draws in this plotting workflow.

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
uv run --frozen python scripts/plot_selection_episode_series.py --workers 4
```

The script separates verified data extraction/cache reuse from plotting and metric validation. `series_data.npz` contains the full plotting arrays; `af_milestones.tsv` and `tmrca_profiles.tsv` contain tabular summaries; `plot_row_order.tsv` maps each displayed row to its simulation; `verified_calls.tsv` contains the recomputed metrics. `input_manifest.json`, `provenance.json`, `pdf_validation.json`, and `output_hashes.json` record provenance and integrity checks.
