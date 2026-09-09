# Mean, median, posterior-mass calls and delayed selection

This experiment compares hard-call definitions of `frac_recent_T` on identical
simulated sequences and pair manifests. It also adds selection beginning at
10 kya, after an unchanged 2% introgression pulse at 50 kya. The later-onset arm
is an illustrative timing experiment, not a fitted model of EPAS1.

## Fixed scientific design

- EAS; 400 haplotypes (200 diploid individuals); 10,000 distinct haplotype pairs.
- Fixed mutation rate 1.25e-8 and recombination rate 1e-8 per base per generation.
- The decoder keeps theta=0.00075, rho/theta=0.8 and 25 years per generation.
  No alternative Ne, mutation rate, recombination rate, or flow field is tested.
- TMRCA cutoffs: 5, 10, 20, 30, 40 and 50 kya. No cutoff exceeds 50 kya.
- Selection coefficient s=0.005 and dominance=0.5 for both onset arms.
- 100 selected replicates for onset 50 kya, 100 new selected replicates for onset
  10 kya, and one shared panel of 1,000 neutral replicates.
- Simulate 11 Mb and crop/site-align an isolated 10 Mb region. Selected alleles
  must be observed in the sample; fixation is retained. Neutral alignment uses
  the nearest segregating site without an allele-frequency floor.
- The later arm has neutral drift between introgression and selection. The
  survival condition spans the pulse to the present, so these are power estimates
  conditional on observing the selected allele, not establishment probabilities.
- RNG seed base 20380101 and pair seed 1729 are unchanged. The existing seed
  formula includes selection onset, so later-onset replicates have distinct,
  reproducible genealogies. They are not paired genealogical counterfactuals.

The new selected simulations use `configs/eas_q02_10k_h400.json`. It specifies
zero additional neutrals because onset has no effect on the neutral demographic
model. Baseline files under `sim/eas_q02_h400` remain immutable and are checked
against their completion receipts. They are the fresh study from this task,
not the old handoff simulation data.

## Four decoded fractions and true TMRCA

| Label | A pair counts as recent when |
|---|---|
| mean | Posterior mean TMRCA is below T; the established baseline |
| median | Posterior median TMRCA is at most T, equivalently P(TMRCA < T) >= 0.5 |
| prob80 | P(TMRCA < T) >= 0.8 |
| prob90 | P(TMRCA < T) >= 0.9 |
| truth | True pair TMRCA is below T |

Each result is the number of passing pairs divided by 10,000. None uses
`mean_p_lt_T`. Mean and truth profiles are reused by reference with checksum
validation. One pass of the same decoder binary saves the per-pair gamma
posteriors. A small native helper reuses `src/recent_stats.h` to compute all four
hard-call summaries from those posteriors. This avoids three repeated HMM passes.
For every replicate, both the native pass and replayed mean must reproduce the
original mean's exact integer counts. The helper checks raw-stream length,
compression integrity and posterior validity; the wrapper checks pair and
position manifests. `--no_recent_probability` skips computation of unused
average probabilities; it does not alter the HMM or hard-call criteria.

Small native tests compare all four rules with direct incomplete-gamma
evaluation and check posterior replay, padded pair lanes and truncated streams.
The smoke phase checks all pair counts, valid fractions, age-cutoff monotonicity,
mean <= median, and prob90 <= prob80 <= median on one replicate per arm. Full
analysis checks posterior-mass nesting over all 1,200 input regions. Different
hard-call rules always receive separately matched neutral calibration.

## Regional evaluation

The primary endpoint is a qualifying call anywhere in the entire 10 Mb region.
All 100 selected replicates per onset are retained, including any fixed samples.
No selected-carrier labels or distance-to-focal tolerance determine a call.

Two scores are prespecified at each age cutoff: the maximum raw fraction and
the maximum fraction standardized at each position. Five deterministic folds
split whole regions. In each fold, 400 neutrals fit normalization, another 400
calibrate the regional score, and 200 are held out. Selected regions split
80 for training and 20 held out. The same fold roles are used for each hard-call
rule and onset comparison. The shared neutral sample is not counted as a new
independent sample for each onset.

Each score/cutoff is evaluated at nominal regional alpha 0.01, 0.05 and 0.10.
A separate policy sets the threshold to attain at least 70% sensitivity among
training selected regions and measures the held-out selected and neutral calls.
The simple raw maximum does not need normalization, but uses the same calibration
split for comparability. All ranks include ties and the finite-sample +1 term.

A secondary grid compares p<=0.05, 0.01 and 0.001; reporting strides 10, 20 and
50 kb; zero or one allowed missing stride; and runs of 1, 3, 5, 10 or 20 significant
strides. Neutral ranks leave the replicate itself out. Each grid entry reports
regional any-call rates. Its p=0.001 tail is resolution-limited by 1,000 neutrals;
the grid is exploratory and is not itself a calibrated regional-level 0.001 test.

These are internal per-rule comparisons. Choosing a winner across all displayed
cutoffs, sources and evidence rules requires separate selection and confirmation.
Matched null calibration addresses the decoder's fixed-Ne mismatch under the
simulation model. It does not guarantee recovery of information lost in decoding.

## Phases, output and validation

`fresh_power --phase simulate` writes simulations, mean decoding and true profiles
without the original focal-only analysis. It supports onset at or after the pulse
in forward time; population survival and observation conditioning are unchanged.
The existing onset-50k configuration and seeds retain their original behavior.

`scripts/launch_recent_calls.sh` supports `simulate`, `smoke`, `decode`, `analyse`
and `run`. `run` finishes and validates all selected simulations, then decodes and
analyses the full comparison. `smoke` runs new hard-call definitions only for
replicate zero; the later-onset smoke simulation must first be available.

Simulation output: `D:/phase2simselection/sim/eas_q02_onset10k_h400`.
Comparison output: `D:/phase2simselection/sim/eas_recent_calls`.
The comparison's `analysis` folder contains regional call summaries, the complete
held-out predictions, fitted thresholds, secondary grid, allele frequencies,
diagnostics, provenance and figures saved as PNG and editable-text vector PDF.
Figures are not displayed inline or emitted as notebook output.
The independent audit recomputes the 960 regional summaries from 528,000 unique
held-out predictions, checks output hashes, and verifies that the 192 mean/truth
baseline rows reproduce the previously completed 400-haplotype comparison.

Twenty workers run one decoder thread each; simulation and decoding phases do
not overlap their worker pools. WSL currently exposes about 188 GiB RAM, below the
user's 300 GB ceiling. The storage guard conservatively reserves 20 GB within the
3 TB allowance using whole-D-volume usage. No per-commit cache invalidation is
used: scientific inputs, pair/position manifests and output hashes determine
whether completed work can be reused. Each phase uses a filesystem lock.

Using the existing pinned WSL/Linux uv environment and native dependencies
described in [the fresh-study setup](fresh_eas_power.md):

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
bash scripts/launch_recent_calls.sh run
```

Paths can be overridden with `SIM_BASELINE_DIR`, `SIM_LATER_DIR` and
`SIM_COMPARISON_DIR`. A new clone still requires the existing simulation inputs;
they are not stored in Git. The separate [flow-field research note](variable_ne_flow_field.md)
describes a possible variable-Ne decoder without changing this experiment.
