# Reproduce the s=0.002 figure packs

Repository: git@github.com:soisa001/gamma_smc_ts.git, branch AOU_run_opt.
Run in WSL with the existing D: archive at `/mnt/d/phase2simselection/sim`.
The repository includes summarized results, not tree sequences or raw posteriors.

```bash
git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
bash scripts/launch_s002_figures.sh build
bash scripts/launch_s002_figures.sh package
```

The launcher uses uv and the existing `.venv`. Phases are idempotent:
`decode` decodes the 100 saved s=0.002 trees with 20 single-thread workers,
hash-checks cached results, and never starts simulations or ancestry replays.
`analyze` calculates focal calibration and all-pair spatial power/FPR.
`build` renders both packs from those analyses; `package` validates PDF/PNG
structure and verifies ZIP contents. Run decode/analyze only if their saved
outputs need to be built. Original simulation and pair seeds are retained.

Inputs: eas_q02_h400, eas_s002_saved_decoding, eas_joint_scan,
eas_allele_class_ablation, eas_ihs_h400, eas_lab_meeting_20260921/analysis,
and eas_lab_meeting_s002_20260921/analysis. Exact paths and hashes are in
figure_provenance.json; decoder inputs, defaults and outputs have manifests
and per-replicate receipts in eas_s002_saved_decoding.

Saved s=0.002 archaic carrier labels are known only at the selected allele.
No spatial carrier score is produced. The all-pairs-only pack excludes
carrier-mass plots. Spatial scores use all 10,000 sampled pairs at every stride.
After changing plots, render the PDFs and visually inspect every page; manual
review is valid only for the hashes recorded in quality_checks.json.
