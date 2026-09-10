# Archaic pair-class confirmation after an all-pair scan

This specification was written before examining the new joint detector's power
or false-positive results. The primary comparison is
`joint_g50_r1_T50000`, separately for true and posterior-mean TMRCA calls.
Other gates, times, runs, and the training-selected detector are exploratory.

## Region, labels, and pair classes

Reuse only the fresh EAS 400-haplotype experiment: 1,000 neutral 10 Mb regions,
100 selected regions with s=0.005 and onset 50 kya, and 100 with onset 10 kya.
The original 2% pulse, fixed mutation/recombination rates, 10,000 pair manifest,
seeds, sampled-allele survival condition, and every selected replicate remain.
No genotype, tree, or original decoding input is replaced.

Recover neutral ancestry by replaying the original seed with migration recording.
Simplified node, edge, individual and population tables must equal the original.
Every pulse migration must have represented sample descendants throughout its
recorded interval; an orphan record fails verification rather than becoming an
absence call. Recording extra unary migration nodes was tested and rejected:
that setting changed the genealogies. It is not used in the analysis.

For selected simulations, replay the original SLiM seed with one callback that
remembers every archaic individual at the pulse. Remembering consumes no random
draws. Original sampled haplotypes are mapped by their SLiM genome identifiers.
After removing the additional remembered samples, every pedigree edge over
every genomic interval below the archaic split, and its node time and population,
must match the original. The recovered focal ALT mask must match the saved mask.
The replay is used only to recover ancestry; its mutations are never used.

An eligible marker is an observed, biallelic, single-mutation ALT allele in the
original simulation. Its mutation age lies between the pulse and the archaic
split, and every ALT haplotype inherits archaic pulse ancestry at that position.
The same rule applies to neutral and selected regions. This is an oracle panel
of archaic-origin variants, not an inferred archaic SNP panel. Variants with
ambiguous recurrent/multiallelic histories are excluded. Fixed ALT markers are
retained. A region with no usable markers or no sampled ALT/ALT pairs stays in
the denominator and gets no joint call.

The classes are ALT/ALT, REF/REF and ALT/REF among the fixed arbitrary haplotype
pairs. They do not refer to homozygous individuals. A marker defines its own pair
classes. The true selected coordinate and mutation type are never used to select
scan candidates or rank a region.

## Two ways to place the scan

1. **10 kb grid:** evaluate TMRCA at 0, 10,000, ..., 9,990,000 bp. Attach the nearest
   eligible archaic marker within 5 kb, breaking equal distances toward the left.
   An empty grid location cannot contribute a joint score. The all-pair baseline
   can still use every grid location.
2. **Every observed archaic site:** evaluate TMRCA and pair classes at the exact
   position of every eligible marker. This includes fixed ALT markers when
   present, preserving the requirement to retain complete sweeps.

Decode the union of grid and marker positions with the unchanged Gamma-SMC
binary and all original VCF variants. Use the native posterior-mean hard call,
with its exact float32 clipping and tie rule. Do not average posterior mass.
Verify extracted counts against the native summary, and true grid counts against
the original truth profile. Record any decoded-grid change induced by adding
output positions, rather than silently assuming output-density invariance.

## Score and first-stage gate

For cutoff T, let fAA, fRR, fAR, and fAll be `frac_recent_T` in the three classes
and all pairs. Let wAA = nAA / Npairs. The confirmation score is

    J_T = wAA × max(fAA_T − fRR_T, 0) × (1 − fAR_T).

It rewards a substantial contribution of recent ALT/ALT pairs, greater recent
coalescence among carriers than noncarriers, and little recent coalescence in
mixed pairs. Its scale reflects carrier frequency as well as TMRCA; that is
intentional. Archaic allele frequency alone is an explicit ablation, so better
performance than the all-pair baseline cannot be attributed to TMRCA without
also checking this ablation.

The primary first-stage gate requires fAll at the same position and T to exceed
the median of pooled neutral grid values from the **fit** regions. Thus it is a
permissive local screen, not a 5% whole-region significance requirement. An AND
rule after an already strict all-pair regional call cannot recover regions that
the first stage missed. Gates at the 80th and 95th neutral percentiles and no
gate are included as comparisons. At fixation, missing REF/REF and ALT/REF
fractions use zero comparator baselines, making J_T=fAll_T; the region is retained.
With no ALT/ALT pairs, J_T=0. No pair-count independence or binomial pair-level
standard error is assumed.

Scores are calculated at T=5,10,20,30,40,50 kya. The primary T is 50 kya. An
additional score takes the maximum over these six cutoffs before calibration.

For spatial evidence, each physical 10 kb bin receives its maximum candidate
score. The region score is the maximum score for one bin, or the maximum of the
minimum scores in two or three consecutive bins. Multiple adjacent SNPs inside
one bin cannot satisfy a multi-bin requirement. The primary uses one bin.
The full site/time/run search belonging to a detector is completed before its
region score is calibrated.

## Held-out power and neutral-region call fractions

Use five fixed whole-replicate folds seeded with 20380101. Each outer fold has
400 neutral fit regions, 400 neutral calibration regions, and 200 neutral test
regions. Per onset, 80 selected regions are available for training and 20 are
held out. The local gate uses fit neutrals only. Calibration uses the upper-tail
rank p=(1+number of calibration scores >= candidate score)/(1+400), including ties.
Report calls at p<=0.01, 0.05, and 0.10. Every replicate is tested once.

The trained joint detector chooses gate, cutoff and run length using only fit
neutrals and the 160 selected training regions, maximizing their pooled training
power at 5%. It then uses the independent calibration and test regions. Training
selection is repeated inside each outer fold. No empirical onset label is needed
to apply that detector to a test region.

An additional operating point sets a score threshold from the lower 30th
percentile of selected training scores, separately for each onset, and evaluates
the resulting power and neutral call fraction on the outer test fold. Zero
evidence is never called. This measures the cost of aiming at 70% power; it does
not assert that the target will be reached out of sample.

Power is the fraction of selected regions with **any call anywhere in the 10 Mb
region**. The user's operational “FDR” is the fraction of neutral regions with
any call, statistically a regional false-positive rate. No distance-to-causal
window defines success or failure. Conventional discovery FDR depends on the
prevalence pi of selected regions:

    FDR(pi) = (1 − pi) × FPR / [pi × power + (1 − pi) × FPR].

Report illustrative FDR and F1 for pi=0.5, 0.1 and 0.01. These are conditional on
the survival/observation ascertainment and the isolated-region simulation model,
not a guarantee for a whole-genome empirical scan. Wilson intervals summarize
held-out counts; they do not include dependence from shared fitted thresholds.

These are exploratory, cross-validated estimates within the current simulation
cohort. The selected pair-class summaries from this cohort informed development
of the score, so the cohort is not an untouched prospective validation set.
The primary formula was frozen before inspecting its new whole-region joint
power/null results. A final chosen method should be confirmed with new seeds.
The benchmark also inherits the existing msprime-neutral / scaled-SLiM-selected
simulation design and perfect archaic-marker labels. The selected allele is
aligned to the central grid coordinate by the original simulation protocol;
the every-site comparison is less dependent on that grid alignment.

The new score explicitly uses archaic carrier frequency. Earlier findings about
the all-pair statistic being insensitive to realized ancestry do not establish
the same property for this score. Its null is the specified neutral 2% pulse
model, without purifying selection against archaic ancestry.

With 400 calibration regions the smallest valid rank p is 1/401, so this held-out
design cannot resolve p=0.001. A separately marked exploratory table uses all
1,000 neutrals to calibrate fixed ungated scores at p<=0.001. A selected score
must exceed every neutral score. Neutral leave-one-out ranks in that table are
constrained by construction and are **not independent false-positive validation**.

## Reproducibility

Run `scripts/launch_joint_scan.sh labels`, then `profiles`, then `analyse`, or
`run` for all phases. Work uses 20 workers, one native thread per worker, and D:
storage. Input manifests, artifact hashes, replay equality proofs, full command
lines, parameters, fold assignments, calibration references, scores and held-out
predictions are saved. Existing phase outputs are reused only after identity and
integrity checks. Large intermediates stay under `sim/eas_joint_scan`; compact
results and PNG/vector PDF figures are copied to the repository for review.

Figures are saved without inline output. The plotting code checks canvas bounds,
legend separation from data panels, and vector PDF content.
