# EAS selected calibration quarantine

This directory preserves the immutable EAS calibration plan and the 120
`screen20` failures produced on 2026-08-13 before any SLiM simulation began.
All 120 ledger rows failed while rebuilding the extended-event contract because
the executor derived the AF half-width from TSV-rounded bounds. For example,
the canonical lower bound `0.07500000000000001` was read as `0.075`, producing
a numerically equivalent but byte-distinct event contract.

- Original plan manifest SHA-256: `56ad08f5e23e5aa5ea0a214fd06a9c5f29b7d2f5d7404bc061f7c963622185ac`
- Original manifest contract SHA-256: `4f44c0836ce14aaad3949a04e6f70a66050b515d2807c6ddd86e53c25baaab5e`
- Original `screen20_plan.tsv` SHA-256: `8644b88256f6adb53d3a0c2965b2ef334f05428c5c35f0b03b70b3cb855c323b`
- Original `screen20_ledger.tsv` SHA-256: `bda657e6df14c9bd575e60a28c62fe32b6161f203713d94a0bb7fe4278050249`
- Planned calibration-module SHA-256: `45d40d3023d6739fe8efc3de0a62f87aef04e764c4354942cf4dacda89978873`
- Preserved per-draw `last_failure.json` records: 120
- Preserved per-draw lock/provenance records: 120

The plan must be regenerated because its source binding intentionally changes
when the executor is fixed. The old artifacts remain here as evidence only and
must not be resumed or mixed with the regenerated plan.
