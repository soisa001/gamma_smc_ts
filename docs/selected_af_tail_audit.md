# Audit of the s=0.006 allele-frequency tail

The audit found no corrupt, mislabelled, or inconsistent endpoint records.
All 1,000 selected replicates passed checks of their original seed and attempt,
arm identity, production source hashes, simulation fingerprint, file hashes,
raw/cropped focal genotypes, sampled-carrier masks, and recorded census counts.
No replicates were removed or replaced.

The s=0.006 tail contains 10/100 replicates below 50% sample AF. The lowest,
rep0020 (seed 1089285226), has 23/400 sampled carriers (5.75%) and 787/15,054
population copies (5.22785%). Thus its low sample frequency agrees with a low
population frequency; it is not merely a bad draw from a high-AF population.
The remaining low-tail records also agree with their population endpoints.

| s | Median sample AF | 5th percentile | Replicates below 50% |
|---|---:|---:|---:|
| 0.005 | 82.125% | 47.7375% | 10/100 |
| 0.006 | 89.5% | 24.4125% | 10/100 |
| 0.007 | 96.25% | 63.625% | 3/100 |

These arms use different deterministic seeds, rather than paired stochastic
paths. The s=0.006 median rises as expected while its lower tail extends
farther than that of s=0.005. Endpoint integrity supports retaining these
outcomes. It does not establish the precise historical cause of the long tail
or rule out every possible modelling issue.

The archive stores present-day sample frequencies, present-day population
counts, and archaic mutation-placement metadata. None of the 1,000 final tree
sequences contains a historical frequency trajectory; their samples are all
at the present. The generator records placement and final counts, not
generation-by-generation population frequencies. Simplified final genealogies
do not provide the full population census through time.

Recovering actual paths would require replaying the original seeds with
frequency logging and checking the reproduced endpoints and tree tables.
The user declined those replays, so figures 03_a and 03_b were deferred and
the existing 15 figures remain unchanged. A logger passed a SLiM initialization
dry run, but no historical simulation replay or replacement was started.

The [audit report](results/lab_meeting_20260921/AF_TAIL_AUDIT.md) and its source
tables are included in the updated figure ZIP. The original archive remains
unchanged. To rerun only the read-only endpoint audit, from the existing WSL
checkout:

```bash
~/.local/bin/uv --no-config run --no-project --python .venv/bin/python \
  python scripts/audit_selected_af_tail.py
```
