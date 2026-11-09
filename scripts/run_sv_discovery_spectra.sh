#!/usr/bin/env bash
# Run in the controlled environment containing the BCF. No upload or git push.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 INPUT.bcf [OUTPUT_DIRECTORY]" >&2
    exit 2
fi
command -v bcftools >/dev/null
command -v uv >/dev/null
INPUT=$(realpath "$1")
OUT=${2:-sv_discovery_spectra}
mkdir -p "$OUT"
OUT=$(realpath "$OUT")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export MPLBACKEND=Agg
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
exec > >(tee "$OUT/run.log") 2>&1

date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
printf 'input=%s\noutput=%s\nthreads=4\nlength_bins=30\n' "$INPUT" "$OUT"
bcftools --version
uv --version

echo 'Phase 1: export sample manifest, header, and site annotations locally.'
# Sequential full-file query checks decoding through EOF; no index is required.
# Refresh these small exports on every run, so no stale input cache is reused.
bcftools view -h "$INPUT" > "$OUT/header.vcf.tmp"
mv "$OUT/header.vcf.tmp" "$OUT/header.vcf"
bcftools query -l "$INPUT" > "$OUT/samples.txt.tmp"
mv "$OUT/samples.txt.tmp" "$OUT/samples.txt"
bcftools query -u \
    -f '%CHROM\t%POS\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%INFO/AC\n' \
    "$INPUT" > "$OUT/sites.tsv.tmp"
mv "$OUT/sites.tsv.tmp" "$OUT/sites.tsv"

echo 'Phase 2: validate annotations and save aggregate plots/tables.'
uv run --script "$SCRIPT_DIR/plot_sv_discovery_spectra.py" \
    --sites-tsv "$OUT/sites.tsv" --samples-file "$OUT/samples.txt" \
    --source-bcf "$INPUT" --out "$OUT" --length-bins 30
date -u +"completed_utc=%Y-%m-%dT%H:%M:%SZ"
