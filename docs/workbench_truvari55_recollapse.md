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
and set `--sizemax -1` to include SVs above 50 kb.

| Type/representation | 90% experiment | 70% experiment |
| --- | --- | --- |
| DEL with valid intervals | Size AND reciprocal overlap >=90%; sequence off | Size AND reciprocal overlap >=70%; sequence off |
| Sequence-resolved INS | Size AND sequence similarity >=90% | Size AND sequence similarity >=70% |
| Symbolic INS | Pass through without collapse | Pass through without collapse |
| Ambiguous/unsupported alleles | Pass through without collapse | Pass through without collapse |
| Symbolic DUP/INV | Existing size/position criteria, size >=90% | Existing size/position criteria, size >=70% |

Symbolic calls stay partitioned by full ALT type. Types and representations never
collapse together. BND is excluded. Pass-through preserves each partition record's
alleles and genotypes; the workflow still adds audit IDs and recalculates count
annotations as before. No reference FASTA is needed.

`comparison.tsv` reports input, retained, and collapsed records by type and
representation, including counts above 50 kb, `pctovl`, and `action` (collapse or
passthrough). Progress messages print the active criteria. Per-sample ALT dosage and
carrier-record totals must be preserved. `run_manifest.json`, per-job
`complete.json` and `truvari.log`, and the builder's `build_manifest.json` and
`installed.txt` record versions, parameters, commands, checksums, and paths.

Wait for an existing run to finish before running this cell again. Existing source
BCFs and validated Truvari 5.5 partitions are reused. DEL results now go to
`pct90/DEL.symbolic_size_overlap` and `pct70/DEL.symbolic_size_overlap`; symbolic
INS go to `INS.symbolic_passthrough`. Previous size/position-only directories are
left on disk as a baseline but are not included in the updated `comparison.tsv`.
Unchanged jobs, including completed sequence-resolved INS jobs, remain reusable.

Local validation: a complete source build, a repeated cached build, and two
integration tests using actual Truvari 5.5.0 and bcftools passed. Synthetic cases
cover >50 kb literal and symbolic alleles, distinct 90%/70% outcomes, shared-carrier
protection, type separation, invalid spans, ambiguous bases, genotype totals, and
cache reuse. Equal-length DELs with 80% overlap merge only at 70%; those with 25%
overlap merge in neither experiment. Symbolic INS pass through even when their
lengths/positions would previously have allowed merging. Controlled cohort data
have not been run locally.
