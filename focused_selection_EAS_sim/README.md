# Focused EAS selected-versus-neutral campaign

This is the isolated, nonsynced campaign root for the focused Gamma-SMC study.
The immutable default plan contains twelve selected biological cells and six
matched neutral banks:

- EAS PHLASH median demography without introgression, at final AF 10%, 20%,
  and 30%.
- `stdpopsim` `AncientEurasia_9K19`, sampled in Han after the modeled
  Neanderthal pulse, at final AF 10%, 20%, and 30%.
- Each demography-by-AF cell requests 10 accepted trajectories for each of
  `s=0.01` and `s=0.005`: 120 selected units. Each demography-by-AF cell has
  one bank of 100 independently seeded matched neutral trajectories: 600
  neutral units. The same neutral bank is used for both prespecified s
  comparisons because the null has no selection coefficient; it is not copied
  and counted twice. The campaign therefore contains 720 simulations.

The plan uses a 10-Mb sequence centered at 5 Mb, 25 years per generation, and
thresholds 1, 4.5, 10, 20, 30, 40, and 50 kya. `execution_units.tsv` fixes the
unit identities and seeds. `study_design.json` fixes the analysis-facing
contract and its checksum. Re-running `plan` validates both without rewriting
them.

```powershell
uv sync --frozen --all-extras
uv run --no-sync python scripts/run_focused_selection_campaign.py plan `
  --repo-root C:\Users\Lenovo\codex\gamma_smc_ts
```

All high-churn files belong below `focused_selection_EAS_sim/work/`, which is
ignored by Git. The CLI rejects campaign roots outside the checkout or inside
any path component containing `OneDrive`.

The `simulate`, `decode`, and `analyze` phases are implemented as idempotent,
checksummed stages. They preserve these contracts:

1. Neutral units use independently seeded `msprime` ancestries under the cell's
   exact demography: DTWF for the recent 200 generations, then
   `SmcPrimeApproxCoalescent`. The 1-kya and 4.5-kya endpoints (40 and 180
   generations) therefore lie wholly in the recent DTWF interval; ancestry
   older than 200 generations remains an explicit coalescent-with-recombination
   approximation rather than the exact Hudson coalescent. Both neutral and selected draws
   now begin with a 500-diploid candidate pool that must fall in the target AF
   band (target +/- 2.5 percentage points). Only then is a 100-diploid panel
   sampled uniformly from all subsets with exactly 20, 40, or 60 alternate
   copies and at least two hom-ref, one heterozygote, and two hom-alt diploids.
   For EAS neutral units, each ancestry is the requested 10 Mb and ascertainment
   is restricted to branches in the fixed 5-Mb marginal tree, weighted by
   mutation-time length. For Han, the executor scans a 20-Mb ancestry by its
   Loschbour-to-Neanderthal pulse-migration records, ascends those lineages to
   the archaic branch crossing the fixed 2,400-generation mutation age, weights
   nonoverlapping candidates by genomic span, and crops a 10-Mb window. This
   migration-indexed construction permits one older archaic origin to feed one
   or several pulse lineages while retaining an auditable no-ILS route. Both
   use fixed batches of four ancestries and the same exact-k panel sampler.
   Neutral units cannot be duplicated across AF cells and counted as independent.
2. Selected units use the established stdpopsim/SLiM model contracts. EAS is a
   de novo mutation; Ancient Eurasia uses an archaic-specific introgressed
   allele with selection beginning at introgression. Both demographies remove
   the two terminal population-AF conditions from the SLiM event stream and use
   the same external candidate-pool band plus exact-k panel ascertainment. Han
   additionally conditions the selected copy on survival along the post-pulse
   Loschbour-to-Han recipient path. Its target- and s-specific fixed selection
   endpoint must come from the disjoint terminal-AF calibration: production
   refuses to start unless `calibration/han_selection_end_frozen.json`, its
   Q=5 production contract, production-seed digest, complete six-cell mapping,
   and every child-artifact checksum validate. Selection is neutral after the
   frozen endpoint.
   The external 500-diploid pool gate is enforced after strict selected-mutation
   identity validation and before exact-panel sampling; a high-frequency pool
   cannot be forced down to a target-frequency panel and accepted.
3. Gamma-SMC runs once per unit, retaining the raw within-individual pair
   posterior. Overall, hom-ref, heterozygous, and hom-alt summaries are derived
   from that same output instead of launching four decoders.
4. Analysis retains one focal row per unit, source, genotype class, and
   threshold. It reports internal selected-simulation contrasts, matched
   selected-versus-neutral contrasts, ROC AUC, descriptive neutral
   leave-one-out ranks, conditional power, tree-truth accuracy, and one-sided
   upper-tail Monte Carlo p-values with the +1 correction.
5. An empirical input table must retain population/demography, AF, genotype
   class, threshold, statistic, value, pair count, and source provenance so the
   simulation-null p-value has an auditable comparison target.

Each planned selected unit is one accepted replicate, not one draw. It has a
deterministic disjoint retry stream capped at 50 draws, with a five-minute
wall-clock timeout per draw and a 45-minute cumulative timeout. Thus the plan
contains exactly ten accepted immutable unit IDs per biological cell when
complete; there are no replacement/backfill IDs. At the prespecified minimum
20% calibration hit rate, exhausting all 50 draws has probability below
1.5e-5 per unit. Biological draw exhaustion, cumulative timeout exhaustion,
and execution failure are reported separately. The selected stage uses the
stdpopsim-compatible SLiM 4.2.2 executable. SLiM
5.2 is installed for benchmarking but is not used because stdpopsim 0.3.0 does
not generate a compatible script. The selected forward model uses `Q=5` and
`slim_burn_in=0.1`; results are therefore exploratory until scaling and burn-in
sensitivity runs are added.

From WSL in the repository root, run the long phases separately so each has its
own durable log:

```bash
mkdir -p focused_selection_EAS_sim/logs

.venv-wsl/bin/python scripts/run_focused_selection_campaign.py simulate \
  --repo-root . --simulation-class neutral --workers 20 \
  2>&1 | tee focused_selection_EAS_sim/logs/simulate_neutral.log

.venv-wsl/bin/python scripts/run_focused_selection_campaign.py simulate \
  --repo-root . --simulation-class selected --workers 4 \
  --slim-bin .native-stdpopsim/bin/slim --slim-scaling-factor 5 \
  --han-selection-calibration \
    focused_selection_EAS_sim/calibration/han_selection_end_frozen.json \
  --selected-max-draws 50 --selected-draw-timeout-minutes 5 \
  --selected-cumulative-timeout-minutes 45 \
  2>&1 | tee focused_selection_EAS_sim/logs/simulate_selected.log

.venv-wsl/bin/python scripts/run_focused_selection_campaign.py decode \
  --repo-root . --workers 20 --threads 1 --decoder-bin bin/gamma_smc \
  2>&1 | tee focused_selection_EAS_sim/logs/decode.log

.venv-wsl/bin/python scripts/run_focused_selection_campaign.py analyze \
  --repo-root . \
  2>&1 | tee focused_selection_EAS_sim/logs/analyze.log

.venv-wsl/bin/python scripts/run_focused_selection_campaign.py report \
  --repo-root . \
  2>&1 | tee focused_selection_EAS_sim/logs/report.log
```

`results/simulation_status.tsv` and `results/decode_status.tsv` expose unit-level
coverage. Every simulation completion binds the scientific contract, source
files, dependency versions, input resources, SLiM binary where applicable, and
output checksums. Every decode completion binds the simulation inputs,
Gamma-SMC binary/settings, raw posterior, and compact output checksums.

The report command can also be run during a partial campaign. It writes the
deterministic `RUN_RESULTS.md` from validated completion records, reports exact
AF/genotype QC, attempts/runtime, and coverage, and labels unavailable analyses
as pending. Once an analysis completion manifest exists, the generator validates
all manifested checksums and table schemas and rejects an incomplete primary
result grid instead of summarizing inconsistent finalized tables. Empirical
p-values remain explicitly unavailable until a provenance-complete empirical
input table is supplied.
