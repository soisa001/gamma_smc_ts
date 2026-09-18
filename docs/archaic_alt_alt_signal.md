# Archaic ALT/ALT signal

The score definition is unchanged. The current primary calibration is
[matched-position or gene-window evaluation](positional_eas_evaluation.md);
the older whole-10-Mb maximum results below are historical comparisons.

For the immediate-onset introgression experiment, the primary statistic is the contribution of recently coalescing
archaic-carrier pairs to the complete haplotype-pair panel. It was previously
reported as the **carrier mass** score; the calculation is already implemented
and evaluated. Its matched-position method is `mass_50000`; the earlier
whole-region method was `mass_g0_r1_T50000`.

## Definition

At an archaic allele j and time cutoff T, define

\[
f_{AA,T}(j)=\widehat P(T_{\mathrm{MRCA}}<T\mid\mathrm{ALT/ALT},j)
=\frac{n_{\mathrm{recent\ ALT/ALT},T}(j)}{n_{\mathrm{ALT/ALT}}(j)}.
\]

The **archaic ALT/ALT signal** is

\[
S_T(j)=f_{AA,T}(j)\frac{n_{\mathrm{ALT/ALT}}(j)}{N_{\mathrm{pairs}}}
=\frac{n_{\mathrm{recent\ ALT/ALT},T}(j)}{N_{\mathrm{pairs}}}.
\]

Thus S is the fraction of all sampled haplotype pairs that both carry the
archaic alternate allele and have recent TMRCA. The probability notation is
conditional on ALT/ALT membership. It denotes the empirical `frac_recent_T`
within that class, not a second probability weight applied to all pairs.

For decoded results, threshold each pair's posterior-mean TMRCA using the
existing native hard-call convention, then count recent pairs. Do not average
posterior mass. For truth results, use the true local TMRCA. ALT/ALT refers to
two arbitrary haplotypes carrying the archaic allele, not diploid homozygotes.

For example, if 2,500 of 10,000 sampled pairs are ALT/ALT and 80% of those
pairs coalesce within T, S=0.25*0.80=0.20: 20% of the complete pair panel
supports recent coalescence among archaic-allele carriers.

The score uses no REF/REF subtraction, ALT/REF penalty, or additional all-pair
gate. Reference-class summaries may be retained as secondary diagnostics.
If there are no ALT/ALT pairs, S=0. At fixation, all pairs are ALT/ALT and the
score becomes the all-pair recent fraction. No AF floor is imposed, and all
selected regions remain in the denominator.

## Manuscript wording

> For each candidate archaic allele, we quantified recent coalescence among
> pairs of haplotypes carrying the archaic alternate allele. We defined the
> archaic ALT/ALT signal as the fraction of carrier pairs with estimated TMRCA
> below a specified time cutoff, weighted by the proportion of all sampled
> haplotype pairs that carried the allele on both haplotypes. Equivalently,
> the statistic measures the fraction of all sampled pairs that both share
> the archaic allele and coalesce within the specified time interval. This
> combines the recency of carrier genealogies with their representation in
> the sample. We assessed significance against neutral simulations subjected
> to the same archaic-marker ascertainment, pair sampling, and positional
> evaluation. At a pre-specified position, the null uses that same position
> in each neutral replicate. For a pre-specified gene interval, the null
> uses the same statistic within a matching neutral interval.

## Scan and interpretation

Use the existing fixed pair panel. The primary scan evaluates every 10 kb,
assigning the nearest eligible archaic marker within 5 kb to define carrier
membership. The comparison scan evaluates every eligible archaic marker at
its exact coordinate. Marker-free positions cannot contribute signal.
The current primary endpoint is the focal position. A pre-specified gene
window can be scored separately and calibrated against an identical neutral
window. The earlier whole-region comparison used the maximum S over all
evaluated candidates in 10 Mb; it answers a different question.

The primary cutoff is T=50 kya. The all-pair benchmark remains
`frac_recent_T` across the complete pair panel. The ALT/ALT signal is its
contribution from pairs carrying the candidate archaic allele: it requires
the same pairs to satisfy both the carrier and recent-coalescence conditions.
S is an evidence score, not a selection p-value. Obtain local p-values by
ranking the local score against the same position/window in independent
neutral calibration regions, including ties. Peaks outside the target window
are excluded.

Current scope is EAS, a 2% pulse at 50 kya with selection starting immediately,
s=0.001 through 0.010, fixed mutation/recombination rates, 400 haplotypes and 10,000 sampled
pairs, and isolated 10 Mb regions. The later-onset arm is excluded. Primary
whole-region defaults are recorded in [the historical analysis preset](../configs/archaic_alt_alt_signal.json).

## Historical whole-region performance

These are the already audited results for this exact score, without a separate
all-pair gate. No new simulations or decoding were needed for the reframing.
At T=50 kya and regional p<=0.05, decoded power is **83/100** selected regions
with **48/1,000** neutral calls on the 10 kb grid; every-site scanning gives
**84/100** and **49/1,000**, respectively.

With thresholds selected from training regions to target approximately 70%
power, held-out decoded results are **71/100 selected and 13/1,000 neutral**
calls on the grid, or **69/100 selected and 14/1,000 neutral** calls at every
site. The [generated tables](results/archaic_alt_alt_signal/tables.md) also
include true-TMRCA results. The [full cutoff metrics](results/archaic_alt_alt_signal/metrics.csv)
and [70%-target metrics](results/archaic_alt_alt_signal/target70_metrics.csv)
are available as machine-readable summaries.

The statistic intentionally depends on carrier frequency as well as recency.
The AF-only benchmark remains relevant to testing whether TMRCA adds detection
value. This score was adopted after examining the existing-cohort ablation;
the reported cross-validated estimates are exploratory, not prospective
confirmation on new seeds. The oracle archaic labels and existing simulation
model remain assumptions of the benchmark.

The earlier [joint-score report](joint_scan_results.md) and
[class-ablation report](allele_class_ablation_results.md) document development
and comparisons; their historical primary labels are not the current score.

## Reproduce the current summary

On the configured WSL machine with cached analysis on D:

```bash
%%bash
cd /mnt/d/phase2simselection/code
git pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
"${HOME}/.local/bin/uv" --no-config run --no-project --python .venv/bin/python \
  python scripts/summarize_archaic_alt_alt_signal.py
```

The summary verifies source hashes and the completed analysis audit, then
selects the existing score's rows. The underlying simulation and decoding
artifacts retain their original identities and seeds.
