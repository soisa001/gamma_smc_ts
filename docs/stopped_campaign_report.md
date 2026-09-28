# Stopped introgression campaign report

The user stopped simulation generation on 2026-09-28 and requested a comprehensive
report. The frozen snapshot contains 2,740 regions: 2,351 selected regions,
194 independent neutral tests, and 195 calibration nulls. I50 has 97 nulls and
96 independent neutral targets; I10 has 98 of each. Selected cells have 98
replicates except I50 s=0.001, which has 97. No replacement simulations are run.

The planned grid remains 100 targets per treatment and 1,000 calibration nulls
per onset. This stopped snapshot is **provisional**: both null calibration and
the selection of completed versus unfinished simulations limit interpretation.
Completion speed depends on survival/retry behavior, so this is not a random
subsample of the planned grid.

## Computation

`stopped_campaign_report.py` freezes the task manifest and complete simulation
receipts. It checks scientific identity, accepted seeds, and receipt hashes.
Analysis inputs (cropped tree sequence, carrier mask, and AF trajectory) are
hash-checked and focal genotypes independently compared with the carrier mask.
This does not claim a new full ascertainment re-audit.

The same 10,000 pairs and original seeds are retained. True focal TMRCA uses
strict `< T`. Gamma-SMC uses the established native posterior-mean rule, including
its float32 arithmetic and `<= T` threshold convention. No carrier-mass weighting
is applied. All-pairs and ALT/ALT statistics are separate. Carrier-pair membership
is fixed by the focal allele across the spatial profile. Missing ALT/ALT pairs
produce an unavailable score and a no-call; they remain in target denominators.
Unavailable null scores rank below all available scores, matching the existing
pipeline's unconditional no-call convention.

The decoder sees the complete 10-Mb region. Its output is requested only at the
focal site, plus a +/-500-kb profile on the existing 10-kb grid for one
representative of each treatment/control cohort. Representatives are chosen
closest to median final AF, with task-ID tie breaks; rounding distances to 12
decimal places removes floating-point noise in exact AF-count ties. They are
individual examples, not cohort means. Temporary VCFs use identical one-based
coordinates and individual order to the standard streaming decoder input.

Calibration uses actual available independent nulls, separately by onset,
source, pair class, and cutoff: `p=(1 + # null scores >= target)/(1 + n_null)`.
Alpha remains 0.05; ties are conservative. With 97/98 nulls, the smallest p-values
are 1/98 and 1/99, and the largest attainable significant p-values are 4/98 and
4/99. Wilson intervals condition on the fitted null and exclude uncertainty in
null estimation. Each cutoff is a separate test, without a best-cutoff inference.

Independent unselected introgressed test regions estimate **FPR**, the
neutral-only endpoint the user called FDR. Discovery FDR for a mixed population
is not inferred from the artificial target:null ratio.

Truth and decoding are explicit, restartable phases with per-region input/output
hash receipts. Reports include decoded comparisons only after decoding the full
frozen snapshot; early decoder completions never define calibration cohorts.
`calibration_diagnostics.csv` records score saturation, missing pairs, p-value
resolution, and whether any bounded target score could achieve p<=0.05.

## Artifacts

The report directory includes a browsable HTML index, a combined vector PDF,
separate PNG/vector PDF figures, and a ZIP of report deliverables. CSV tables
contain exact cohort counts, focal raw scores, per-target p-values, calls,
power/FPR with intervals and coverage, AF and retry summaries, representative
IDs, calibration diagnostics, and unfinished task IDs. Provenance and output
checksums accompany the bundle. Plots are saved to disk, not notebook output.
Every FPR heatmap uses the same fixed 0–100% scale; power also uses 0–100%.

## Reproduce or resume analysis only

Use the configured uv environment and existing simulation output on the analysis
host. This launcher does not invoke SLiM or resume simulation generation.

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
export SIM_OUTPUT_DIR=/mnt/d/phase2simselection/sim/segregating_introgression_h400
export REPORT_OUTPUT_DIR="$SIM_OUTPUT_DIR/report_20260928"
export REPORT_WORKERS=20
bash scripts/report_stopped_campaign.sh all
```

Use `report` instead of `all` to rebuild figures and tables from completed
analysis caches. On another host, configure `ORIGIN_PYTHON`, `GAMMA_NATIVE_DIR`,
and the input/output paths as for the original campaign.
