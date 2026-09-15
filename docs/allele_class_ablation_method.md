# ALT/ALT and REF/REF ablation without an ALT/REF penalty

This comparison was specified after the original joint-scan results, in response
to the request to remove the mixed-pair penalty and compare carrier and
noncarrier recency. It is exploratory in the existing fresh simulation cohort.
Only the 1,000 neutral and 100 selected EAS regions with selection starting at
the 50-kya introgression pulse are used. Existing seeds, folds, original
400-haplotype/10,000-pair profiles, and all regions are retained.

At each eligible archaic SNP, fAA and fRR are the ALT/ALT and REF/REF hard-call
`frac_recent_T` values. They are not averaged posterior probabilities. REF
means noncarrier of that marker, not an independently simulated neutral allele.
Power for either score means detecting a selected 10 Mb region using that
class; false-positive rate means calling a neutral 10 Mb region.

The requested primary comparisons have no all-pair gate or frequency weight:

- ALT/ALT alone: fAA.
- REF/REF alone: fRR.
- Carrier contrast: max(fAA - fRR, 0).

To isolate the ALT/REF penalty from other changes, also compare
wAA*max(fAA-fRR,0) with the original
wAA*max(fAA-fRR,0)*(1-fAR), where wAA=nAA/Npairs. Include the simpler
carrier mass wAA*fAA, all-pair recency, and archaic AF as benchmarks.
Each carrier/class score is evaluated ungated and with the original median
all-pair gate. The latter holds the other parts of the original detector fixed.

Evaluate every eligible archaic site (primary placement) and the existing
10 kb grid with its nearest eligible marker within 5 kb. Mask unassigned grid
positions for every class/contrast score; their stored placeholder REF class
is not real noncarrier evidence. The all-pair grid baseline uses all grid
positions. Score each region by its maximum over candidate positions. Only
single-position support is tested here. T=5,10,20,30,40,50 kya; T=50 kya is the
primary cutoff. No cutoff/run selection is performed.

An empty class contributes no stand-alone evidence. With no ALT/ALT pairs the
contrast is zero. With no REF/REF pairs use zero as its contrast baseline,
retaining complete sweeps. No AF or minimum pair-count filter is imposed.
Report how often the maximum unweighted class fraction saturates at one, and
the supporting pair counts; rare-allele saturation is an anticipated concern,
not grounds for excluding regions afterward.

Keep the previous whole-region fold assignments: 400 fit neutrals, 400 separate
calibration neutrals, 200 test neutrals, and 80 selected training/20 selected
test regions per fold. Recompute and verify the median gate from fit neutrals.
Apply conservative upper-tail rank p=(1+#calibration scores>=score)/401.
Report p<=.01,.05,.1. For fixed T=50-kya rules and AF, also set thresholds from
the lower 30th percentile of selected training scores and measure held-out
power and neutral call fraction. Zero evidence never qualifies. Preserve ties;
do not randomize them to force nominal type-I error or a 70% power target.

Compare the saved original scores for all overlapping unchanged rules. Check
that exact-site true-TMRCA scores are identical with and without the ALT/REF
factor. Independently audit every new calibration rank, call count, and
70%-target training threshold. Existing simulations and decoding are not rerun.
