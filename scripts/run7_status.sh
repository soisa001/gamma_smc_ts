#!/bin/bash
cd /home/mew/AllOfUs_Phase2/gamma_smc_ts || exit 1
echo "=== per-arm completed replicates ==="
for d in neutral eas_introgressed_s0p002 eas_introgressed_s0p003 \
         eas_introgressed_s0p005 eas_introgressed_s0p01 \
         denovo_s0p005 denovo_s0p01 denovo_s0p02 denovo_s0p05; do
  n=$(find "sim_results_run7/$d" -name endpoint.json 2>/dev/null | wc -l)
  printf "  %-26s %3d\n" "$d" "$n"
done
echo
echo "=== processes ==="
printf "  main run7 : %d\n" "$(pgrep -cf 'run_run7\.py' || true)"
printf "  denovo    : %d\n" "$(pgrep -cf 'run_run7_denovo\.py' || true)"
echo
echo "=== failures ==="
find sim_results_run7 -name endpoint.json -exec grep -l '"status": "failed"' {} + 2>/dev/null | wc -l
echo
echo "=== recent AF by arm ==="
.venv/bin/python - <<'PY'
import json, glob, collections
import numpy as np
by = collections.defaultdict(list)
for p in glob.glob("sim_results_run7/*/selected/replicates/*/endpoint.json") + \
         glob.glob("sim_results_run7/neutral/replicates/*/endpoint.json"):
    r = json.load(open(p))
    if r.get("status") != "completed":
        continue
    by[r["arm_id"]].append(r["final_allele_frequency"]["sample_af"])
for arm in sorted(by):
    v = np.array(by[arm])
    print("  %-26s n=%3d  AF mean %.3f median %.3f" % (arm, v.size, v.mean(), np.median(v)))
PY
