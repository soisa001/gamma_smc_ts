#!/bin/bash
cd /home/mew/AllOfUs_Phase2/gamma_smc_ts || exit 1
for d in neutral recent_denovo postintro_denovo introgressed_pulse introgressed_pulse_weak; do
  printf "  %-24s %3d\n" "$d" "$(find "sim_results_run8/$d" -name endpoint.json 2>/dev/null | wc -l)"
done
printf "  procs: %s\n" "$(pgrep -cf 'run_run8\.py')"
printf "  disk : %s\n" "$(df -h /home/mew | tail -1 | awk '{print $4" free"}')"
