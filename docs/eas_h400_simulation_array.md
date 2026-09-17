**Simulation archive for later decoding.** This extends the existing fresh EAS
400-haplotype cohort to s=0.001 through 0.010 in increments of 0.001, with 100
accepted selected replicates per coefficient and the same 1,000 neutral
replicates. The completed s=0.005 arm and neutrals are verified and reused;
900 selected replicates are newly generated. No old handoff simulations or
200-haplotype replicates are used.

The array was subsequently paused at the user's request with 1,947/2,000
regions validated. See the [paused-array evaluation](eas_h400_paused_evaluation.md)
for the saved counts, preliminary true-TMRCA results and evaluation-only
commands. Generation and decoding remain paused.

The configuration is [eas_q02_50k_h400_array.json](../configs/eas_q02_50k_h400_array.json).
Selection starts immediately at the 50-kya, nominal 2% introgression pulse.
All other model inputs are unchanged: PHLASH EAS history, archaic split and
bottleneck, dominance 0.5, fixed mutation/recombination rates 1.25e-8/1e-8,
25-year generations, SLiM rescaling factor five, 11 Mb simulation and a
site-aligned 10 Mb crop. All surviving observed selected alleles are retained,
including fixation. There is no terminal-frequency target or AF filter.

Each replicate saves the complete mutation-bearing `simulation.trees`, the
cropped `decoded_input.trees`, `focal_carriers.npy`, accepted-attempt metadata,
checksums, and execution logs. The root pair manifest retains the original
10,000 pairs from 400 haplotypes (seed 1729). The study seed remains 20380101;
adding coefficients does not change any existing replicate seed.
These inputs preserve true genealogies, genotypes and focal carrier labels for
future decoding. Genome-wide oracle archaic-marker labels can subsequently be
recovered using the existing ancestry-replay workflow; they are not newly
computed by this archive phase.

The output directory is `/mnt/d/phase2simselection/sim/eas_q02_h400` in WSL,
or `D:\phase2simselection\sim\eas_q02_h400` on Windows. The archive phase uses
20 single-threaded workers, a maximum RAM allowance of 300 GB (subject to WSL's
available memory), and the existing conservative 3 TB volume-use limit with a
20 GB reserve. Expanding this cohort updates the root study/sample manifests
and retains exact earlier versions under `manifest_history/`. The scientific
fingerprint, pair manifest and completed simulation files remain compatible.
Earlier analysis outputs describe their original cohort and are not regenerated.

**Phases.** The dedicated launcher accepts:

- `archive` (default): generate only missing simulations, then independently
  reload and audit all saved simulations.
- `simulate-only`: generate missing simulations and validate their saved inputs.
- `simulation-status`: verify and inventory existing simulations, without
  generating missing work.
- `audit-simulations`: verify all planned simulations; fail if any are missing.

The existing general `fresh_power --phase simulate` retains its historical
meaning (simulation plus decoding, without analysis) for older launchers.
The new dedicated launcher deliberately does not accept that phase or `run`.
No decoder or all-pair truth-profile calculation is called for new simulations.
The existing decoder binary is hashed only for compatibility with the saved
study fingerprint; it is not executed during simulation-only work.

`simulation_status.json` records progress, per-arm completion, errors, and
`decoding_requested: false`. The historical `run_status.json` continues to
describe the earlier decoded run. `sample_manifest.json` includes all planned
2,000 region identities. `simulation_inventory.csv` indexes verified saved
regions with absolute input paths, selected frequencies, accepted seeds,
checksums, and whether a prior decoded completion exists. The inventory is
written after startup validation and at completion or a handled failure; the
status JSON is the live progress record during the run.

Validation checks the scientific fingerprint and task/seed identity, all three
simulation-file hashes, readable raw/cropped trees, exact sample counts and
sequence lengths, diploid ordering, original and shifted focal coordinates,
and equality of focal genotypes to the saved carrier mask and frequency. A
checksum-consistent but biologically inconsistent carrier mask fails validation.
Missing or corrupt completed inputs stop the workflow. Interrupted simulations
resume from validated phase receipts; incomplete simulation attempts are
restarted with the original deterministic seed. A study lock prevents competing
runners. Restarting WSL interrupts computation, but completed receipts remain
reusable.

**Run or resume in WSL.** This uses the existing local runtime and data disk;
override `SLIM_BIN`, `GAMMA_NATIVE_DIR`, and `UV_BIN` when using another host.
The general environment setup is documented in
[fresh_eas_power.md](fresh_eas_power.md).

```bash
%%bash
set -euo pipefail
repo=/mnt/d/phase2simselection/code
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
cd "$repo"
git pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
bash scripts/launch_eas_h400_array.sh archive
```

While a runner holds the study lock, read `simulation_status.json` and
`logs/simulation_array.log` directly. After it exits, the following verifies
the full archive without decoding:

```bash
%%bash
cd /mnt/d/phase2simselection/code
bash scripts/launch_eas_h400_array.sh audit-simulations
```

Ten focused `test_fresh_power.py` tests passed before launch, including a real
small neutral simulation through the new CLI phase with decoder/truth-profile
functions replaced by failure sentinels, an idempotent restart, an independent
audit, preserved historical status/manifests, fixed-allele retention, and
rejection of corrupt or inconsistent saved inputs. No full test suite was run.

**Hypothesis for later decoding: founder diversity and complementary evidence.**
When selection expands one introgressed haplotype, linked archaic allele
frequencies can already provide substantial power. Carrier coalescence and
allele frequency may then describe much of the same expansion. With several
introgressed backgrounds, a joint signal of archaic carriage and recent
coalescence may add information beyond allele frequency alone. This is a
hypothesis to test, not an established gain in the present results.

Several backgrounds can mean different biological situations: multiple
introgressing founders bearing the same selected allele, different archaic
alleles carried on different backgrounds, or recombined descendants of one
founder. These should not be treated as interchangeable. At a single marker,
the current carrier-mass score is
`M_T = (n_alt_alt / n_total) * frac_recent_T_alt_alt`.
At fixed allele frequency and pair sampling, any gain comes from differences
in the fraction of carrier pairs with TMRCA below T. Distinct successful
founders can also produce older between-founder carrier pairs and weaken this
recent-coalescence fraction. Thus, more backgrounds do not by themselves
guarantee greater power for the current score.

The current regional detector takes the maximum score across positions; it
does not combine evidence from multiple markers or require spatially sustained
carrier coalescence. A later aggregation method could test whether several
moderate, spatially coherent signals help, with neutral simulations calibrating
their linkage and dependence. An aggregated carrier score must be compared
against an allele-frequency baseline using the same marker set, physical
windows, and aggregation opportunities, in addition to the existing maximum-AF
baseline.

The 2% pulse in this array does not specify exactly one contributing haplotype.
The model fixes the selected allele in the archaic source before introgression,
and the source bottleneck can reduce diversity without enforcing one founder.
The later analysis should measure ancestry contribution and relatedness rather
than assume the number of founders from the pulse fraction. Where retained
ancestry records or deterministic replays permit, record the number of
contributing pulse lineages, their descendant weights (and effective diversity),
and the distribution of carrier coalescence times.

The planned comparisons are AF-only versus carrier mass, followed by any
proposed regional aggregation, using truth and decoded `frac_recent_T` at
T <= 50 kya. Evaluate power at matched neutral 10 Mb region call rates,
stratify by measured founder diversity, and use AF-matched comparisons to
isolate information beyond frequency. Keep any score tuning separate from
evaluation. Linked archaic SNPs are not independent replications. Explicit
founder-diversity simulation arms would be a separate follow-up with matched
neutral controls; this note does not change the running array's model or start
decoding.
