# Reproduce the s=0.005 all-pairs-only pack

Repository: git@github.com:soisa001/gamma_smc_ts.git, branch AOU_run_opt.
Run in WSL on the machine with the saved D: archive at
`/mnt/d/phase2simselection/sim`.

```bash
git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
bash scripts/launch_s005_all_pair_figures.sh build
bash scripts/launch_s005_all_pair_figures.sh package
```

The launcher uses uv and the existing `.venv`. Both phases are idempotent.
`build` verifies source hashes, checks the selected coefficient and onset,
independently recomputes focal power/FPR from saved integer pair counts,
checks agreement with the spatial center, and renders 10 figures.
`package` checks PDF/PNG structure, writes SHA-256 inventories and verifies
the ZIP. There is no simulation, replay, or decoder phase.

Required saved inputs: eas_joint_scan, eas_allele_class_ablation,
eas_positional_distance_h400, and eas_lab_meeting_20260921/analysis.
Exact inputs, hashes and parameters are in figure_provenance.json. Plotting
uses four profile readers and one thread per numerical-library process.
The simulation and pair RNG seeds are retained; plotting introduces no
randomness. Raw trees and per-pair posteriors remain on D:; the repository
contains figures, summarized data, and provenance.

Render and visually inspect every PDF page after a plot edit. Manual review
in quality_checks.json applies only to the exact recorded PDF hashes.
