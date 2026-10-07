# Provisional neutral CDFs: completed simulations only

Each observation is one prespecified focal site from an independent retained
neutral region, not every position within the 10-Mb region. Focal alleles are
archaic-derived, segregating immediately after the 2% pulse at 50 kya, and
observed today; subsequent fixation is retained. No new simulations or decoding
were run by this reporting script. The simulation runner is not modified.

The four main numbers are percentages of pairs with TMRCA below 10,000 or
50,000 years, not TMRCA ages. The blue step is the focal population ECDF; gray
steps are the other populations. Orange markers give descriptive 95th and 99th
percentiles using the inverse empirical CDF (nearest-rank order statistic).
They are placed on the nominal 95%/99% guides; a finite ECDF can jump past these
levels. No smoothing or fitted tail model is used. Missing ALT/ALT pair classes
are excluded from the displayed CDF and their counts are reported.

The descriptive percentiles are NOT the exact significance thresholds.
cutoffs_all_methods.csv separately reports conservative p=(1+#null>=score)/(n+1)
boundaries for p<=0.05 and p<=0.01; a finite target score must strictly exceed
the boundary. Missing null pair classes remain in that test's denominator,
matching the campaign convention. Fewer than 19 nulls cannot attain p<=0.05;
fewer than 99 cannot attain p<=0.01. Saturation at score=1 can also prevent
rejection. Unsupported thresholds are flagged, never replaced with percentiles.

All populations are below the planned 1,000 retained nulls. Small-sample tail
percentiles often equal the observed maximum. Completed trajectories may be
biased toward faster runs. These plots are provisional descriptions, not a
finished null calibration. The all-pairs decoded tables/figures are the primary
view; truth and ALT/ALT are supplied separately.

All source calibration receipts and truth/decoded score-file hashes were checked.
PNG and editable vector PDF outputs share the same figure code.

This snapshot uses the mapped-autosomal-region MVN-demography campaign. The
CDFs pool completed regions. Demographic-history CDF bands are not inferred
from singleton histories; history_coverage.csv reports replication. Completion
order can favor faster histories, so incomplete-cohort tails remain provisional.

## Archived snapshot: 2026-10-07

535 completed regions: AFR 88, AMR 90, EAS 89, EUR 90, MID 89, SAS 89.
The combined PDF has 36 pages, including allele-frequency distributions with
95th-percentile markers. Simulation state was unchanged by report generation.
The full local bundle also includes individual PNG/vector PDF figures, focal
scores, source receipts, and SHA-256 provenance manifests.

ALT/ALT at 50 kya reaches a 100% upper-tail boundary for every population in
both true and decoded TMRCA. Even a target score of 100% cannot reject at
p <= 0.05 with this snapshot. These calibration nulls do not provide an
independent neutral-target false-positive estimate or selected-treatment power.

Validation: all 96 cutoff rows checked against sorted raw scores and direct
empirical p-values; source and output hashes checked; both ZIP archives passed
CRC validation; all 36 PDF pages contain vector content without embedded raster
images. Representative summary and plot pages were rendered and inspected.

## Reproduce from local campaign outputs

Run from a Linux/WSL checkout. The campaign path must contain the completed
calibration receipts and score files; those large intermediate files are not
included in Git. Use a new output directory for a new snapshot. Add
`--reuse-snapshot` to redraw exactly an existing report's verified receipt list.
Neither command starts simulations.

```bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only origin AOU_run_opt
: "${CAMPAIGN:?Set CAMPAIGN to the local campaign directory}"
: "${REPORT:?Set REPORT to the report output directory}"
uv run --extra selection --with pypdf --with reportlab \
  python scripts/plot_population_null_cdfs.py --root "$CAMPAIGN" --out "$REPORT"
uv run --extra selection --with pypdf --with reportlab \
  python scripts/summarize_genomic_null_snapshot.py --campaign "$CAMPAIGN" \
  --snapshot "$REPORT" --report-date 2026-10-07
```
