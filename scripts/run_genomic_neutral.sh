#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
root="${GENOMIC_OUTPUT_DIR:-${repo_dir}/../sim/genome_neutral_20261005}"
native_dir="${GAMMA_NATIVE_DIR:-/home/mew/AllOfUs_Phase2/gamma_smc_ts/.native}"
export LD_LIBRARY_PATH="${native_dir}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export MPLBACKEND=Agg PYTHONUNBUFFERED=1 UV_NO_CONFIG=1 UV_LINK_MODE=copy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
phase="${1:-plan}"
if [ "$#" -gt 0 ]; then shift; fi
cd "$repo_dir"
python_cmd=("${UV_BIN:-${HOME}/.local/bin/uv}" --no-config run --no-project --python "${ORIGIN_PYTHON:-${repo_dir}/.venv-origin/bin/python}" python)
if [[ "$phase" == prepare ]]; then
  mkdir -p "$root/maps"
  urls=(
    'https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw'
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz'
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
  )
  for url in "${urls[@]}"; do
    target="$root/maps/${url##*/}"
    if [[ ! -s "$target" ]]; then
      curl --fail --location --retry 3 "$url" --output "$target.tmp"
      mv -- "$target.tmp" "$target"
    fi
  done
  "${python_cmd[@]}" scripts/draw_population_histories.py --out "$root/histories"
  exec "${python_cmd[@]}" scripts/prepare_genomic_null_regions.py --maps "$root/maps" --out "$root/regions"
fi
exec "${python_cmd[@]}" -m gamma_smc_aou.genomic_neutral --out "$root/campaign" \
  --regions "$root/regions" --histories "$root/histories" \
  --slim "${SLIM_BIN:-${native_dir}/bin/slim}" --phase "$phase" "$@"
