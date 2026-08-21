#!/bin/bash
cd /home/mew/AllOfUs_Phase2/gamma_smc_ts || exit 1
cp -f /mnt/c/Users/Lenovo/AllOfUs_Phase2/gamma_smc_ts/python/gamma_smc_aou/run7_analysis.py python/gamma_smc_aou/
cp -f /mnt/c/Users/Lenovo/AllOfUs_Phase2/gamma_smc_ts/scripts/run7_figures.py scripts/
sed -i 's/\r$//' python/gamma_smc_aou/run7_analysis.py scripts/run7_figures.py
export MPLBACKEND=Agg
.venv/bin/python scripts/run7_figures.py
