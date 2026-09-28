# Introgression and selection manuscript draft

- [Methods and Results](introgression_selection_draft.md): editable working manuscript with mean posterior P(TMRCA < 50,000 years), detailed introgression methods, and locally supported results. Selection simulations are excluded.
- [Evidence and open items](evidence_and_open_items.md): sources, historical/current discrepancies, and exact production outputs needed to replace remaining fields.
- [Verified aggregate summary](verified_summary.json): numerical checks and source fingerprints. No participant-level records are included.
- `package_manifest.json` and `verify_draft.py`: portable integrity check for this document package. This check does not reproduce controlled-data analyses.

The biological results were not manufactured from test fixtures or transferred from the legacy 4,500-year scan. Manuscript placeholders deliberately remain where the final production evidence is unavailable.

To retrieve the draft and verify its files, run this code-only cell. It does not access controlled buckets or execute a cohort analysis.

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts_manuscript"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --single-branch --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
  git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
cd "$REPO"
uv run --no-project python docs/manuscript/verify_draft.py
```

The validation checks UTF-8 text integrity, SHA-256 hashes, numerical count identities, and the requested time-scale conversion. The original local-source audit remains under the workspace's `manuscript_review_20260928/` directory; source hashes in the aggregate summary permit later reconciliation with those local files.
