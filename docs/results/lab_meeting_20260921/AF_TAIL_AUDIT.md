# s=0.006 lower-AF tail: integrity audit

All 1,000 selected replicates passed checks of source identity, seed/attempt,
arm, file hashes, raw/cropped focal genotypes, carrier masks, and census counts.
No replicate was removed or replaced.

The s=0.006 arm has 10/100 replicates below 50% sample AF. Its lowest run,
rep0020, has 5.75% sample AF (23/400) and 5.22785% population AF (787/15,054).
Its low endpoint therefore exists in the simulated population itself. The
s=0.005 arm also has 10/100 below 50%, but its minimum is higher (25.75%).
Median AF is 82.125% at s=0.005 and 89.5% at s=0.006. Different arms use
independent seeded paths, so a higher median need not imply a higher minimum.

These checks found no evidence that the tail is broken. They do not reconstruct
its history or exclude every possible modelling issue. Removing valid low-AF
outcomes would alter the simulated distribution, so they remain in all figures.

No historical AF time series was saved in the 1,000 archived selected runs.
The user chose to skip replaying them, so 03_a and 03_b are deferred. The
existing 15 figures are unchanged. No historical simulation replay or
replacement was started; the original archive remains unchanged.

Audit files:
- analysis_audit/af_tail_audit.json
- analysis_audit/af_tail_selected_endpoint_audit.csv
- analysis_audit/af_tail_s006_lower_tail.csv
- analysis_audit/af_tail_af_tail_by_s.csv

Audit command (WSL, existing checkout, no simulations):

```bash
~/.local/bin/uv --no-config run --no-project --python .venv/bin/python python scripts/audit_selected_af_tail.py
```
