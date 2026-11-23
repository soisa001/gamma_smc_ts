#!/usr/bin/env bash
set -euo pipefail
umask 077
cd "$(dirname "${BASH_SOURCE[0]}")/.."
base="$HOME/ashg_revised"
original="${POPPED_FULL_ROOT:-$HOME/popped_full_same_original_pos_p70}"
fixed="${REPAIRED_ROOT:-$HOME/popped_full_same_original_pos_p70_gt_union}"
long="${LENIENT_RUN:-$base/v3_lenient_recollapse_same_start_90_70}"
out="${UNIFIED_ROOT:-$HOME/unified_popped_lenient_p70}"
jobs="${JOBS:-4}"
mkdir -p "$out"
exec > >(tee -a "$out/workbench.log") 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH="$HOME/.local/bin:$PATH"
export POPPED_BILLING_PROJECT="${POPPED_BILLING_PROJECT:-${GOOGLE_PROJECT:-}}"
uv_bin="$(command -v uv || true)"
for tool in .tools/uv-bin/uv "$HOME/tlr_gt_audit_code/gamma_smc_ts/.tools/uv-bin/uv"; do
    if [[ -x "$tool" ]]; then uv_bin="$tool"; break; fi
done
if [[ -z "$uv_bin" ]]; then
    bash scripts/bootstrap_singleton_tools.sh
    uv_bin=.tools/uv-bin/uv
fi
chroms=("$@")
if (( ! ${#chroms[@]} )); then chroms=(chr{1..21}); fi
merge=(); repair=()
for c in "${chroms[@]}"; do
    [[ "$c" =~ ^chr([1-9]|1[0-9]|2[01])$ ]] || exit 2
    if [[ ! -s "$original/$c/final/done.json" || ! -s "$original/$c/final/old_to_new.jsonl.gz" ]]; then merge+=("$c"); fi
    if [[ ! -s "$fixed/$c/done.json" || ! -s "$fixed/$c/$c.popped.samepos.p70.gt_union.vcf.gz" ]]; then repair+=("$c"); fi
done
if (( ${#merge[@]} )); then
    POPPED_FULL_JOBS="$jobs" POPPED_FULL_ROOT="$original" bash scripts/run_popped_full_batch_workbench.sh --chromosomes "${merge[@]}"
    # Recheck repair receipts against regenerated memberships.
    repair=("${chroms[@]}")
fi
if (( ${#repair[@]} )); then
    "$uv_bin" run --no-project --python 3.12 scripts/repair_popped_cluster_gts.py \
        --batch-root "$original" --out "$fixed" --jobs "$jobs" --chromosomes "${repair[@]}"
fi
if ! PYTHONPATH=scripts "$uv_bin" run --no-project --python 3.12 --with pysam==0.24.1 --with google-crc32c==1.8.0 \
    python -c 'from build_unified_sv_panel import long_inputs; from pathlib import Path; import sys; long_inputs(Path(sys.argv[1]))' "$long"; then
    bcf="$base/v3_callset_comparison/inputs/v3_ultralong_lenient.bcf"
    if [[ ! -s "$bcf" ]]; then
        mkdir -p "$(dirname "$bcf")"
        # The nested prefix is also a directory-marker object. List from its parent.
        staging="$out/runtime/stage"
        mkdir -p "$staging"
        gcloud storage rsync gs://rw-migration-aou-rw-8cbcf26c/resources "$staging" \
            --recursive --exclude '^(?!ultralong_bnd_update/v3_ultralong_lenient\.bcf(?:\.csi)?$).*' \
            --checksums-only --do-not-decompress ${GOOGLE_PROJECT:+--billing-project="$GOOGLE_PROJECT"}
        for suffix in '' .csi; do
            test -s "$staging/ultralong_bnd_update/v3_ultralong_lenient.bcf$suffix"
            cp --reflink=auto "$staging/ultralong_bnd_update/v3_ultralong_lenient.bcf$suffix" "$bcf$suffix.partial"
            mv "$bcf$suffix.partial" "$bcf$suffix"
        done
    fi
    tools="$base/tools/truvari55_uv"
    if ! command -v bcftools >/dev/null; then
        bash scripts/bootstrap_singleton_tools.sh
        export PATH="$PWD/.native/bin:$PWD/.tools/singleton-native/bin:$PATH"
    fi
    "$uv_bin" run --script scripts/build_truvari55.py --out "$tools" --workers "$jobs"
    "$tools/venv/bin/python" scripts/compare_stringent_recollapse.py --bcf "$bcf" --out "$long" --same-start --workers "$jobs" --refdist 500
fi
"$uv_bin" run --no-project --python 3.12 scripts/build_unified_sv_panel.py --original "$original" \
    --repaired "$fixed" --long "$long" --out "$out" --jobs "$jobs" --recover-short --chromosomes "${chroms[@]}"
