#!/bin/bash
cd /home/mew/AllOfUs_Phase2/gamma_smc_ts || exit 1
total=0
for d in neutral pulse_s0p003 pulse_s0p004 pulse_s0p005 pulse_delayed_s0p005 \
         denovo_s0p003 denovo_s0p004 denovo_s0p005 denovo_delayed_s0p01 denovo_delayed_s0p05; do
  n=$(find "sim_results_run9/$d" -name endpoint.json 2>/dev/null | wc -l)
  total=$((total+n))
  printf "  %-24s %3d\n" "$d" "$n"
done
echo "  TOTAL $total / 1200"
printf "  procs %s  load %s  mem %s\n" "$(pgrep -cf 'run_run9\.py')" \
  "$(cut -d' ' -f1 /proc/loadavg)" "$(free -g | awk 'NR==2{print $3"G/"$2"G"}')"
