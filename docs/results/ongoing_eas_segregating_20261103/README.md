# Present-day-segregating EAS selection pilot

Selected simulations require strictly 0 < focal allele frequency at present < 1. Neutral genealogies are unconditioned; their focal label is the nearest present-day segregating variant. This supersedes the exploratory rejection of neutral focal alleles.
Calibration uses 10 independent development, 20 calibration and 10 held-out neutral regions per origin, from complete blocks of the original main and supplemental schedules. The 10 kb grid is unchanged; whole-region maxima and FPR are unaffected by the neutral focal label. The variant is polymorphic in the sampled 200 haplotypes, which establishes present-day segregation.
Selection starts at 50 kya and remains enabled to present. The corrected coefficient grid is 0.005, 0.01, 0.02, 0.03. Stronger recent-onset comparisons (around 5 kya) are reserved for a separate experiment; no recent-onset result is claimed here.
The introgression pulse is 1% at 50 kya, with 500 kya human–Neanderthal divergence, Neanderthal Ne=3,600, and an allele fixed in donors. De novo starts from one copy. EAS PHLASH median demography; 25 years/generation; h=0.5; Q=1; 1 Mb; 100 diploids; 1,000 unfiltered pairs; 10 kb stride; cutoffs 1–50 kya.
Rejected trajectories stop at absorbing loss/fixation. This preserves exact rejection conditioning because no recurrent focal mutation is introduced. Previously completed proposals are reused with receipt and artifact hashes. Cohorts keep their originally assigned seeds; accepted simulations are never moved between cohorts.
Analysis uses only the fully completed proposal prefix within each origin/cohort/coefficient stratum; later completions are recorded but excluded. This avoids choosing individual replicates by completion order, but deadline-dependent stopping and small cohorts remain limitations. Wilson intervals are exploratory.
Calibration uses a maximum across the genomic/time grid, independent null development, and conservative +1 upper-tail p-values. At least 19 independent calibration nulls are needed to resolve p<=0.05. Missing selected cohorts are not reported as zero power.
All scanned times lie inside the selection interval, so temporal overlap does not demonstrate dating accuracy. Temporal F1 measures [1,50) kya coverage. Whole-population segregation does not require both alleles to appear in the sampled 100 diploids.

Completed proposals: 865; usable prefix proposals: 817; accepted selected proposals: 30.

| Origin | Cohort | s | Completed prefix | Accepted | Lost | Fixed |
|---|---|---:|---:|---:|---:|---:|
| de_novo | selected | 0.005 | 63 | 0 | 63 | 0 |
| de_novo | selected | 0.01 | 51 | 0 | 51 | 0 |
| de_novo | selected | 0.02 | 63 | 0 | 62 | 1 |
| de_novo | selected | 0.03 | 63 | 0 | 60 | 3 |
| introgression | selected | 0.005 | 54 | 17 | 37 | 0 |
| introgression | selected | 0.01 | 46 | 13 | 22 | 11 |
| introgression | selected | 0.02 | 62 | 0 | 16 | 46 |
| introgression | selected | 0.03 | 62 | 0 | 11 | 51 |

## Performance availability

Conditional TPR, neutral FPR, region F1 and Wilson 95% intervals are in performance.tsv.
The exploratory neutral-rejection outcomes remain in proposal_summary.tsv for provenance; they are not used for calibration or FPR in this final analysis.

| s | Posterior-mass 1-kyr TPR | Neutral FPR | Region F1 | Temporal coverage F1 |
|---:|---:|---:|---:|---:|
| 0.005 | 4/17 | 0/10 | 0.381 | 0.060 |
| 0.01 | 4/13 | 0/10 | 0.471 | 0.017 |
| 0.02 | Not estimable (0 accepted) | 0/10 | NA | NA |
| 0.03 | Not estimable (0 accepted) | 0/10 | NA | NA |
