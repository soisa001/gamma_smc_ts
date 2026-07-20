# High-frequency s=0.05 selected validation

## Design

- One 10 Mb diploid simulation under the saved two-epoch demography: ancestral
  `Ne=10,000`, recent `Ne=20,000`, with growth 100 generations ago.
- A single-copy selected allele was introduced at 5 Mb 180 generations
  (4,500 years) ago with additive `s=0.05` (`h=0.5`).
- Mutation and recombination rates were `1.25e-8` and `1e-8` per bp per
  generation. The final pseudo-data contain 2,000 diploid samples (4,000
  haplotypes); the statistic averages the 2,000 within-individual haplotype
  pairs.
- The selected trajectory was rejection-sampled with 20 worker processes until
  population AF was at least 0.30. Rejected trajectories reuse one independently
  generated neutral ancestry because the neutral ancestry and forward allele
  trajectory are independent. Only the accepted trajectory was materialized.
- The existing 100 neutral truth simulations and 100 official Gamma-SMC v0.2
  decodes were reused unchanged. All Gamma-SMC summaries use a 1 kb output grid.

## Selected trajectory and tree-sequence truth

- Seed-order acceptance occurred at zero-based attempt 2,514 (2,515 attempts to
  acceptance; 2,520 screens completed in 20-way batches).
- Final population AF was **0.379525** and sampled AF was **0.38275**. The sample
  contained 274 hom-alt, 983 heterozygous, and 743 hom-ref diploids.
- At the selected base, true `P(TMRCA < 4,500 years)` across all 2,000
  within-individual pairs was **0.138**. No neutral truth replicate was as large,
  giving `(0+1)/(100+1) = 0.00990099`.
- The two-epoch theoretical neutral expectation was 0.006479. The 100 neutral
  simulations had mean 0.006555, bias 0.0000761, RMSE 0.001669, and a 0.04
  leave-one-out fraction at `p<=0.05`.
- At 5 Mb, hom-alt mean TMRCA was 161.4 generations (95% CI 159.0-163.9), versus
  19,008.0 generations (95% CI 18,380.8-19,656.6) for hom-ref diploids.

## Official Gamma-SMC decode

- Container: `docker.io/regevsch/gamma_smc:v0.2` with `--only_within`, scaled
  mutation rate 0.0005, recombination/mutation ratio 0.8, and
  `--output_at_stride 1000`.
- The container read 53,188 segregating sites, decoded 2,000 pairs at 10,000
  positions in **12.43 seconds**, and used **0.796 GB** peak memory. Posterior
  reduction took 1.40 seconds.
- At 5 Mb, Gamma-SMC estimated `P(TMRCA < 4,500 years)=0.091691`; the decoded
  neutral mean was 0.023462 (95% interval 0.015309-0.037047). No decoded neutral
  replicate was as large, so the Gamma-SMC Monte Carlo p-value was also
  **0.00990099**.
- Across the full profile, decoded-versus-truth bias was +0.01427, MAE 0.01594,
  RMSE 0.01752, and Pearson correlation 0.8522. The absolute probability is
  biased upward in neutral regions and downward at the selected center, but the
  matched decode-to-decode calibration still ranks the selected center beyond
  all 100 nulls.

## Spatial scan

- There were 2,440/9,999 complete grid windows with raw pointwise `p<0.05`,
  merged into 37 tracts. The tract spanning the selected site was
  **4.7125-5.4165 Mb** (704 windows; minimum p=0.00990099).
- In the +/-500 kb interval, 980/1,001 windows had raw pointwise `p<0.05`.
- Simulation-calibrated window-density p-values were **0.00990099** both for the
  full 10 Mb scan (2,440 selected versus null median 369) and the +/-500 kb scan
  (980 selected versus null median 29).
- No window had BH q<0.05. With only 100 null simulations, the smallest possible
  p-value is 1/101, so ordinary BH correction over 9,999 windows cannot reject.
  Use the simulation-calibrated regional-density statistic (or substantially
  more null simulations) for region-level inference; treat the 37 pointwise
  tracts as descriptive breakpoints, not 37 independent discoveries.

## Main outputs

- `truth_null_and_selected_pvalue.png`: theoretical/empirical truth null and
  selected truth p-value.
- `gamma_smc_center_null_and_selected_pvalue.png`: decoded neutral distribution
  and selected Gamma-SMC center p-value.
- `selected_decoded_vs_neutral_pvalues.png`: full and +/-500 kb decoded profile,
  neutral interval, pointwise p-values, and tract breakpoints; legend is outside.
- `decoded_scan_density_calibration.png`: global and local regional-density nulls.
- `selected_allele_frequency_trajectory.png`: realized AF from origin to present
  with present (`0`) at the left.
- `selected_hom_alt_vs_hom_ref_tmrca.png`: hom-alt/hom-ref truth means and 95% CIs.
- `decoded_study_metrics.json`: machine-readable decode, error, and p-value summary.
- `selected_rejection_attempts.tsv`: complete AF screening audit.
