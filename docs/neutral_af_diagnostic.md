**Why an 18% archaic allele can fail the regional scan.**

This diagnostic reads the saved 1,000 neutral ancestry replays and existing
regional scores. No simulations or decoding were started. The array remains
paused. Each replay and sample/ancestry mapping was checked against its saved
SHA-256 receipt, and its seed against the validated simulation inventory.

At the pre-specified 5.5-Mb coordinate in the original 11-Mb simulation, mean
sampled archaic ancestry is 2.06625%, consistent with the nominal 2% pulse.
However, 836/1,000 regions have no sampled archaic ancestry at that position.
Among the 164 regions where ancestry is present, its mean frequency is
12.5991%. These observations obey the empirical identity
`mean(F) = P(F>0) * mean(F | F>0)`: 0.0206625 = 0.164 * 0.1259909.

| Neutral summary | Number evaluated | Mean frequency | Fraction at least 18% |
|---|---:|---:|---:|
| Local ancestry at a fixed position, including zero ancestry | 1,000 | 2.07% | 4.4% |
| Same position, conditional on ancestry being present | 164 | 12.60% | 26.8% |
| Nearest observed eligible archaic SNP to the focal coordinate | 984 | 9.30% | 15.9% |
| Maximum eligible archaic SNP frequency on the 10-kb grid across 10 Mb | 1,000 | 38.00% | 90.9% |

The first two rows measure total local pulse ancestry, not the frequency of
an ascertained SNP. The last two rows use the existing eligible archaic SNP
definition. Sixteen neutral regions have no eligible archaic markers; these
remain in the regional scan with score zero. The listed fractions are
descriptive empirical tail fractions, not newly cross-validated p-values.

Thus, 18% is in the upper 4.4% at a position specified independently of the
data under the unconditional ancestry null. It is common after restricting
to surviving ancestry, and very common when searching for the strongest
marker across 10 Mb. These are different inferential questions. The previous
selection array was conditioned on survival and observation; its s=0.001
mean of 17.8913% cannot be compared with 2% as if the latter were the null
mean among observed surviving alleles. Its median selected AF is 13.875%.
A mean shift across many replicates is also different from reliably
identifying selection in one region.

The neutral regional maximum AF has median 38.25%, pooled 95th percentile
64.75%, and maximum 88%. The existing five-fold p<=0.05 rule requires a
score strictly above 62.75%, 66.75%, 67.00%, 65.25%, or 62.75%, depending on
the calibration fold. These are the exact 20th-largest values among 400
calibration regions; equality is rejected because ties count in the tail.
The paired carrier-mass boundaries at T=50 kya are 0.3448, 0.3891, 0.4188,
0.3824 and 0.3448.

An allele at 18% can contribute at most approximately 0.18 squared, or
3.24%, to carrier mass, even if every carrier pair has TMRCA below T. For
400 haplotypes and exactly 72 carriers, the all-possible-pairs weight is
72*71/(400*399)=0.03203; the actual score uses the existing 10,000 sampled
pairs. This weighting makes a young, moderate-frequency allele compete
against much higher-frequency neutral tracts. It is a limitation of this
particular regional score, not proof that the underlying partial sweep
cannot be detected by another statistic.

The saved neutral genealogy supplies a concrete explanation for the broad
local distribution. At this coordinate there are, on average, 8.96 ancestral
lineages remaining just before 50 kya, despite sampling 400 present-day
haplotypes. Pulse migration is sampled on those ancestral branches; one
migrant lineage can have many present-day descendants. Among the 164
ancestry-present regions, 90.85% have one sampled-carrier lineage at that
boundary. Increasing the present-day panel reduces sampling noise but does
not eliminate variation caused by these shared genealogies. These numbers
are specific to the current EAS demographic model.

**Comparison with iHS.** The original iHS method integrates the spatial decay
of extended haplotype homozygosity on the two allelic backgrounds and
standardizes the ratio by allele frequency. It therefore evaluates the
unusual persistence of a haplotype relative to its frequency. Its authors
also reported limited power for low-frequency sweeps in their HapMap-based
benchmark, so an 18% allele does not guarantee easy detection. The study
also used clusters of extreme scores in windows rather than only the most
extreme isolated SNP. See [Voight et al. (2006)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072).

An appropriate next comparison would evaluate iHS and a carrier-coalescence
score calibrated conditional on AF and spatial persistence, using the same
introgressed neutral simulations and the same regional error definition.
No iHS power has been measured in this array. Such a comparison should
retain the 10-Mb regional calibration if that remains the inferential goal.

Machine-readable summary: [summary.json](results/neutral_af_diagnostic/summary.json).
The full fixed-position per-region diagnostics remain at
`/mnt/d/phase2simselection/sim/eas_neutral_af_diagnostic_20260918`.

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
  python scripts/diagnose_neutral_af.py
```
