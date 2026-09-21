# EAS s=0.002: decoded presentation figures

Two separate collections use 100 saved selected regions at s=0.002 and the
existing 1,000 decoded neutral regions:

- `docs/results/lab_meeting_s002_20260921`: 15 figures comparing the focal
  carrier score, all-pair recency, AF and iHS. Full-grid AF/iHS plots are
  genotype-based context, with s=0.002 highlighted.
- `docs/results/lab_meeting_s002_all_pairs_20260921`: 10 figures using only
  all-pair TMRCA scores, including matched-position FPR and spatial profiles.

Every primary TMRCA plot uses Gamma-SMC `frac_recent_T`, with posterior-mean
TMRCA hard calls. Figure 07 alone includes true TMRCA as a validation reference.
All cutoffs are at or below 50 kya. No simulations or ancestry/trajectory
replays were run. The previous s=0.005 collection is preserved.

## Focal results at T=50 kya

| Method | Power, p<=0.05 | Neutral FPR | Power, p<=0.01 | Neutral FPR |
|---|---:|---:|---:|---:|
| Archaic AF | 70% | 5.1% | 41% | 0.9% |
| Decoded ALT/ALT carrier mass | 69% | 5.1% | 38% | 1.0% |
| Decoded all-pair recency | 5% | 5.0% | 3% | 1.4% |

Positional iHS power is 12% and 3% at these two thresholds. Across the six
TMRCA cutoffs, decoded all-pair power is 3-11% at nominal 5%; its highest
observed value is at 20 kya, an exploratory comparison. At 50 kya, true
all-pair power is 19% versus 5% after decoding. True focal carrier-mass power
is 70% versus 69% after decoding. These counts do not establish a gain over AF.

The pooled p<=0.001 sensitivity gives 17% AF, 16% decoded carrier-mass and 0%
decoded all-pair power. It uses a different calibration size and has only one
extreme neutral rank; it is not precise validation of a 0.1% tail.

## Scope and verification

Archaic carrier labels exist at the saved selected allele, but not across the
rest of each s=0.002 region. Spatial carrier plots are omitted at the user's
request. The all-pair spatial analysis is available at all 1,000 10-kb strides;
unavailable carrier labels are never interpreted as absence of ancestry.

The original 400-haplotype panel, 10,000-pair manifest, mutation/recombination
rates, simulation seeds and decoder defaults are retained. Decoder outputs
and input trees have SHA-256 receipts. Saved allele genotypes match the carrier
labels exactly; extracted all-pair recent counts match native Gamma-SMC output;
pair classes partition the same panel. The independent focal truth summary and
the center of the spatial evaluation must reproduce the focal counts exactly.

Calibration is at the fixed 5-Mb coordinate, or separately at each spatial
coordinate. Five-fold evaluation uses 400 calibration neutrals and 200 held-out
neutrals per fold. Carrier calibration retains marker-free neutral positions
with score zero. FPR is the fraction of neutral positions called, not discovery
FDR. All selected replicates remain in the denominator, including fixation;
power is conditional on the original survival/observation rule.

## Reproduction on the machine with the saved D: archive

```bash
git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
cd gamma_smc_ts
bash scripts/launch_s002_figures.sh decode
bash scripts/launch_s002_figures.sh analyze
bash scripts/launch_s002_figures.sh build
bash scripts/launch_s002_figures.sh package
```

The launcher uses uv, 20 single-thread decoder workers, and separate phases.
`decode` checks and reuses completed caches. With analyses already present,
only `build` and `package` are needed. Raw trees and posteriors remain on D:;
Git contains figures, summarized results, provenance, and reproducible code.
