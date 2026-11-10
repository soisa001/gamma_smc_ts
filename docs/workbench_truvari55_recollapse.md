# Stringent recollapse with Truvari 5.5.0

Run on the Workbench VM containing the staged BCF. This cell builds and installs
Truvari 5.5.0 with uv and then runs the independent 90%/70% experiments. It does not
access the controlled bucket, upload data, or display plots.

The builder uses Python 3.12, builds a local static zlib to supply the previously
missing `zlib.h`, and uses `uv build --wheel` for Truvari, pywfa, and pytabix.
No sudo is needed. A C compiler, make, ar, uv, and bcftools must already be on PATH.
The pywfa wheel is built locally because the distributed wheel did not contain a
usable Python 3.12 extension in our installation. Wheels are specific to the build
machine despite pywfa's upstream `py3-none-any` filename; do not share this cache
between machines. Source archives and cached output files are SHA256 checked.

```bash
%%bash
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"

BASE=/home/jupyter/ashg_revised
REPO="$BASE/gamma_smc_ts_svplots"
BCF="$BASE/v3_callset_comparison/inputs/v3_ultralong_stringent.bcf"
OUT="$BASE/v3_stringent_recollapse_90_70"
TOOLS="$BASE/tools/truvari55_uv"
mkdir -p "$OUT"
exec > >(tee -a "$OUT/setup_and_run_55.log") 2>&1
date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'

test -s "$BCF"
command -v uv
command -v bcftools
if [[ ! -d "$REPO/.git" ]]; then
    git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
    test "$(git -C "$REPO" branch --show-current)" = AOU_run_opt
    git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
git -C "$REPO" rev-parse HEAD

uv --no-config run --script "$REPO/scripts/build_truvari55.py" \
    --out "$TOOLS" --workers 4

"$TOOLS/venv/bin/python" "$REPO/scripts/compare_stringent_recollapse.py" \
    --bcf "$BCF" --out "$OUT" --workers 4 --refdist 500
date -u '+finished_utc=%Y-%m-%dT%H:%M:%SZ'
```

Both thresholds start from identical partitions, use `--gt all --keep common`,
and set `--sizemax -1` to include SVs above 50 kb. Literal unambiguous DEL/INS use
both sequence and size similarity. Symbolic calls are partitioned by full ALT
type and use size/position matching with `--pctseq 0`; these results do not imply
sequence similarity. Ambiguous or unsupported alleles pass through unchanged.
Types and representations never collapse together. BND is excluded.

`comparison.tsv` reports input, retained, and collapsed records by type and
representation, including counts above 50 kb. Per-sample ALT dosage and
carrier-record totals must be preserved. `run_manifest.json`, per-job
`complete.json` and `truvari.log`, and the builder's `build_manifest.json` and
`installed.txt` record versions, parameters, commands, checksums, and paths.

Existing source BCFs are reused. Partitions from the old implementation are
rebuilt once because pysam 0.24 can infer an interval end from SVLEN even when END
is invalid; the updated classifier checks the explicit span. Old collapse outputs
are not accepted as 5.5 results. Subsequent unchanged runs reuse validated results.

Local validation: a complete source build, a repeated cached build, and two
integration tests using actual Truvari 5.5.0 and bcftools passed. Synthetic cases
cover >50 kb literal and symbolic alleles, distinct 90%/70% outcomes, shared-carrier
protection, type separation, invalid spans, ambiguous bases, genotype totals, and
cache reuse. Controlled cohort data have not been run locally.
