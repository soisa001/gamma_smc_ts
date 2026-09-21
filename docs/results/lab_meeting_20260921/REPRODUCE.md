# Reproduce this figure collection

Repository: git@github.com:soisa001/gamma_smc_ts.git, branch AOU_run_opt.
The scripts run in WSL Ubuntu and use the saved study on D: at
`/mnt/d/phase2simselection/sim`. Source data are not included in Git.

```bash
git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
bash scripts/launch_lab_meeting_figures.sh build
bash scripts/launch_lab_meeting_figures.sh package
```

`build` renders the saved analyses; `package` checks all output PDFs/PNGs,
writes SHA-256 inventories, and verifies the ZIP. Both are idempotent.
Use `refresh` instead of `build` only when focal summaries need to be rebuilt:
it checks the completed archive, reuses verified focal/iHS caches, refreshes
positional calibration, and builds figures. It never starts SLiM or Gamma-SMC.
Focal extraction and iHS use 20 workers; plotting uses four profile readers
and one thread per numerical-library process. Seeds are preserved from the
simulation manifest and no figure uses random jitter or resampling.

Required saved inputs include `eas_q02_h400`, `eas_h400_pause_evaluation_20260917`,
`eas_allele_class_ablation`, `eas_joint_scan`, `eas_ihs_h400`,
`eas_positional_distance_h400`, and this collection's `analysis` directory.
The exact files and hashes used to construct figures are recorded in
`figure_provenance.json`. Original parameter and analysis audits accompany
the package. ZIP contents are only figures, summarized simulation data, captions,
and provenance; tree sequences and per-site raw decoding files stay on D:.

After changing figure code, render the combined PDF with Poppler and inspect
every page. The PDF/PNG checks cannot replace visual inspection. A manual
review recorded in `quality_checks.json` applies only to its exact PDF hashes.
