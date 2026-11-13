"""Author-only support diagnostics from the frozen public-data fit."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
with np.load(ROOT/"agent_submission"/"fitted_model.npz") as fit:
    exposure = fit["context_counts"]
    alt = fit["destination_counts"]
    if alt.shape == (4096, 4):
        mask = np.ones_like(alt, dtype=bool)
        mask[np.arange(4096), np.arange(4096)//16%4] = False
        alt = alt[mask].reshape(4096, 3)
    assert alt.shape == (4096, 3)
    result = {
        "source": "agent_submission/fitted_model.npz; counts learned solely from public calibration",
        "calibration_sites": int(exposure.sum()),
        "calibration_mutations": int(alt.sum()),
        "six_base_contexts": 4096,
        "substitution_channels": int(alt.size),
        "mean_sites_per_context": float(exposure.mean()),
        "mean_mutations_per_context": float(alt.sum(axis=1).mean()),
        "mean_events_per_channel": float(alt.mean()),
        "zero_observation_channels": int((alt == 0).sum()),
        "fraction_zero_observation_channels": float((alt == 0).mean()),
    }
(ROOT/"support_diagnostics.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
