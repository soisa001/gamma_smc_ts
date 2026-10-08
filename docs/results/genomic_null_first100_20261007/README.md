# Empirical neutral CDFs: frozen completed cohort

Each observation is one prespecified focal site from an independent retained
neutral region, not every position within the 10-Mb region. Focal alleles are
archaic-derived, segregating immediately after the 2% pulse at 50 kya, and
observed today; subsequent fixation is retained. No new simulations or decoding
were run by this reporting script. The simulation runner is not modified.

The four main numbers are percentages of pairs with TMRCA below 10,000 or
50,000 years, not TMRCA ages. The blue step is the focal population ECDF; gray
steps are the other populations. Orange markers give empirical 95th and 99th
percentiles using the inverse empirical CDF (nearest-rank order statistic).
They are placed on the nominal 95%/99% guides; a finite ECDF can jump past these
levels. No smoothing or fitted tail model is used. Missing ALT/ALT pair classes
remain in the denominator below finite score support; their CDF mass appears
at the left edge. These are no-calls, not biological zero scores.

The user's specified p-value is count(null >= observed) / n, without any +1
adjustment. Ties count. Undefined targets are no-calls. With n=100, p<=0.05
requires strictly exceeding sorted score 95; p<=0.01 requires strictly
exceeding sorted score 99. Equality at those boundaries does not reject.
The primary percentiles use all n nulls, with null no-calls below finite
support. Available-only percentiles are separately labeled in the CSV.
The CSV also exports tied counts/fractions, p at the cutoff and p at score 1.
When all null scores equal the observation, p=1. At score 1, p is the fraction
of null scores equal to 1. Saturation can prevent rejection at either level.

All populations are below the planned 1,000 retained nulls. With 100 nulls,
the p grid has 0.01 steps and the 99th percentile is the second-largest score;
the 1% tail is imprecisely estimated. A zero empirical tail is possible and
does not establish zero population probability. These plots are provisional
calibration results. The all-pairs decoded tables/figures are the primary
view; truth and ALT/ALT are supplied separately.

All source calibration receipts and truth/decoded score-file hashes were checked.
PNG and editable vector PDF outputs share the same figure code.

This snapshot uses the mapped-autosomal-region MVN-demography campaign. The
CDFs pool completed regions. Demographic-history CDF bands are not inferred
from singleton histories; history_coverage.csv reports replication. Completion
order can favor faster histories, so incomplete-cohort tails remain provisional.

## Verified first-100 milestone

This frozen report contains exactly 100 simulations in each population,
rep0000 through rep0099; all 600 receipts and source score hashes were checked.
Later completed simulations are excluded. The campaign continues running.
The PDF has 37 vector pages: four summary pages, allele-frequency distributions,
and 32 population CDF pages. Individual PNG/PDF plots, focal scores, source
receipts and complete output hashes are available in the local first_100_each
bundle beside the campaign. No calibration-null IQRs were added.

Validation: ten focused empirical-tail tests passed, including equality,
saturation, tied cutoffs and missing carrier pairs. All 96 cutoff rows were
independently checked against sorted scores and inclusive tail counts at and
around every observed score. Source/output hashes and ZIP CRC checks passed.
Summary pages and representative CDF/allele-frequency pages were rendered and
visually inspected. The PDF contains no embedded raster images.

For reproduction commands and the no-call convention, see
../../genomic_population_neutral.md, section 'First 100 per population'.
