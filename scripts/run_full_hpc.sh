#!/usr/bin/env bash
set -euo pipefail

# Override these through the scheduler environment for each population/region.
REPO="${REPO:-$PWD}"
OUT="${OUT:-$REPO/sim_results/full}"
N_SIMS="${N_SIMS:-1000}"
N_DIPLOIDS="${N_DIPLOIDS:-2000}"
LENGTH="${LENGTH:-1000000}"
MU="${MU:-1.25e-8}"
R="${R:-1e-8}"
NE="${NE:-10000}"
SEED="${SEED:-1729}"

extra_sim_args=()
[[ -n "${HISTORIES:-}" ]] && extra_sim_args+=(--histories "$HISTORIES")
[[ -n "${MUTATION_MAP:-}" ]] && extra_sim_args+=(--mutation-map "$MUTATION_MAP")
[[ -n "${RECOMBINATION_MAP:-}" ]] && extra_sim_args+=(--recombination-map "$RECOMBINATION_MAP")

python -m gamma_smc_aou.cli simulate \
  --output-dir "$OUT" --replicates "$N_SIMS" --diploids "$N_DIPLOIDS" \
  --length "$LENGTH" --mutation-rate "$MU" --recombination-rate "$R" \
  --ne "$NE" --seed "$SEED" --save-trees "${extra_sim_args[@]}"

# Decode empirical data (one within-individual pair per diploid).
python -m gamma_smc_aou.cli decode \
  --executable "${GAMMA_SMC:-gamma_smc}" --input "$EMPIRICAL_VCF" --input-format vcf \
  --output "$OUT/empirical.within.tsv" --theta "$THETA" \
  --rho-over-theta "$RHO_OVER_THETA" --mutation-rate "$MU"

# Decode each saved tree replicate as a separate scheduler array job, for example:
# gamma-smc-aou decode --input "$OUT/trees/replicate_${TASK_ID}.trees" --input-format trees \
#   --output "$OUT/decoded/replicate_${TASK_ID}.tsv" --theta "$THETA" \
#   --rho-over-theta "$RHO_OVER_THETA" --mutation-rate "$MU"

# After the array completes:
# gamma-smc-aou calibrate --observed "$OUT/empirical.within.tsv" \
#   --sim-glob "$OUT/decoded/*.tsv" --observed-length "$LENGTH" \
#   --simulation-length "$LENGTH" --match nearest --output "$OUT/site_scan.tsv"
