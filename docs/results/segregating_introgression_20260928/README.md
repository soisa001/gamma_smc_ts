# Available results at the 2026-09-28 stop

2,740 completed regions: 2,351 selected, 194 independent neutral tests, and
195 calibration nulls. The results table contains true-TMRCA, Gamma-SMC decoded
TMRCA, and allele-frequency results. Both TMRCA sources cover all 2,740 retained
regions, with 546 metric rows in total. The comparison table joins truth and
decoded rates for the same treatment endpoints.

Calibration uses 97 I50 and 98 I10 nulls, rather than 1,000 each. All estimates
are provisional. Completion-time selection and calibration uncertainty limit
interpretation. Neutral detection estimates FPR, not mixed-prevalence FDR.

True all-pairs examples: at s=0.01, I50/T=50 kya detects 92/98 (93.9%), and
I10/T=10 kya detects 51/98 (52.0%). True ALT/ALT detection is zero at T=20–50 kya:
matched nulls often score 1, so even the maximum bounded score cannot achieve
p<=0.05 under conservative empirical ties. See the main report-methods document
for exact definitions, missing-pair handling, and reproducible commands.

Local deliverables are in
`sim/segregating_introgression_h400/report_20260928/`: `index.html`,
`comprehensive_report.pdf`, `report_bundle.zip`, individual PNG/vector PDFs,
raw scores, predictions, and provenance receipts. The allele-frequency plots
show the exact empirical upper 5% calibration boundary: AF must exceed the line
to achieve p<=0.05, with conservative treatment of ties. Calibration-null IQR
shading remains omitted.

After completing this snapshot report, simulation generation was authorized
to resume at 20 workers with a 200-GB memory cap. These tables preserve the
stopped snapshot and do not incorporate subsequent simulation completions.
