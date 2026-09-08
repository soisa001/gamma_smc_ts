# Regional localization metrics for the fresh EAS study

This analysis reads only the completed fresh study in `sim/eas_q02`. It checks
every decoded profile against its completion-receipt SHA-256 and byte count,
checks the sample manifest, exact coordinates, pair-count fractions, and cutoff
ordering, and reproduces all 6,000 published selected focal p-values. It does not
simulate, decode, or use old simulation outputs.

## What is scored

Each replicate has 1,000 profile positions at 10 kb spacing over the cropped
10 Mb. The statistic is exclusively `frac_recent_T`: the fraction of pairs whose
decoded posterior mean TMRCA is below T. T is separately 5, 10, 20, 30, 40, or
50 kya; these are statistic cutoffs, not different simulated selection onsets.
The selection onset is 50 kya in every replicate.

The primary operating point is an unadjusted, strict p < 0.05 at every stride.
p < 0.01 is also reported. Each position is compared with the same position in
the 1,000 independent neutral replicates. This retains decoder boundary effects
and the allele-centered design instead of assuming every position is identically
distributed. Ties count against significance:

`p_selected = (1 + number of neutral values >= selected value) / 1001`.

For evaluating neutral regions, the tested replicate is removed from its own
reference: `p_neutral = (1 + number of other neutral values >= value) / 1000`.
The small denominator difference is retained explicitly. Strict p < 0.001
cannot be evaluated with these leave-one-out neutral ranks and is not reported.
The 1,000 linked strides are not treated as 1,000 independent null replicates.

## Localization truth and metrics

The default target is the positions within inclusive ±100 kb of the selected
allele at 5 Mb (21 profile positions). This carries forward the old handoff's
localization tolerance as an explicit analysis assumption. It is not a measured
causal sweep footprint. Selection can affect flanking genealogies outside that
interval, so a call labeled false by this rule need not be biologically unrelated
to selection. Sensitivity results use ±50, ±250, and ±500 kb as well.

- **Stride scoring:** TP is a significant stride inside the target; FP is a
  significant stride outside it; FN is a nonsignificant stride inside it.
- **Peak scoring:** Merge only consecutive significant strides. Score either
  all peaks or peaks of at least five strides. Each replicate has one true event.
  A peak overlaps the target if at least one of its significant strides is inside
  it. At most one called peak matches that event; extra peaks are FP, including
  duplicate calls overlapping the target. If no peak matches, FN is one.
- **FDR:** Average the replicate-level `FP / max(TP + FP, 1)` across all 100
  selected simulations in an arm. No-call simulations contribute zero. This
  differs from the old handoff's average restricted to called simulations.
- **F1:** Average replicate-level `2 TP / (2 TP + FP + FN)` across all 100 selected
  simulations. No-call simulations contribute zero. Both metrics are also saved
  with counts pooled across replicates, under distinct column names.
- **Neutral regions:** Every call is false. Report mean false-call count and the
  fraction of neutral regions with any call. The latter is the all-null regional
  FDR / family-wise error rate. Neutral F1 is not an informative endpoint.

Bootstrap intervals use 2,000 resamples of whole selected replicates, with the
original study's fixed seed. They are conditional on the estimated neutral
reference and omit its Monte Carlo uncertainty.

These are localization metrics in isolated 10 Mb regions. Per-stride alpha is
not an FDR-adjusted threshold. No genome-wide error-control or prevalence-mixed
FDR estimate is implied, and the selected and neutral replicate counts are not
used as a biological prevalence assumption.

## Run the analysis again

The raw study must already exist locally; it is not distributed through Git.
Use the study's pinned uv environment as documented in `fresh_eas_power.md`.
This command is analysis only and never reruns simulations:

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
uv --no-config run --no-project --python .venv/bin/python \
  python -m gamma_smc_aou.fresh_region_metrics \
  --study /mnt/d/phase2simselection/sim/eas_q02 --workers 20
```

Outputs are under `analysis/region`: per-replicate CSVs, selected and neutral
summary CSVs, annotated FDR/F1 figures in PNG and vector PDF, and provenance with
all input-profile and output hashes. The primary figures use p < 0.05 and ±100 kb.
The summarized CSVs and figures are small enough to publish in Git; raw trees
and full simulation profiles remain on D.
