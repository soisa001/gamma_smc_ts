#!/usr/bin/env bash
set -Eeuo pipefail

# Restart-safe All of Us Researcher Workbench runner for the empirical,
# random within-population haplotype-pair Gamma-SMC scan. Inputs remain in controlled GCS buckets;
# decoded summaries, controlled sample manifests, and plots go to the workspace bucket.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AOU="$REPO/scripts/aou.sh"

ALL_POPS=(AFR AMR EAS EUR MID SAS)
ALL_CHROMOSOMES=({1..22})

CHR_SPEC=""
POPS_SPEC=""
LOCAL_ROOT="${AOU_GAMMA_LOCAL_ROOT:-/home/jupyter/gamma_smc_workbench}"
INPUT_PREFIX="${AOU_GAMMA_INPUT_PREFIX:-}"
OUTPUT_PREFIX="${AOU_GAMMA_OUTPUT_PREFIX:-}"
BILLING_PROJECT="${AOU_GAMMA_BILLING_PROJECT:-${GOOGLE_PROJECT:-}}"
BCF_TEMPLATE="${AOU_GAMMA_BCF_TEMPLATE:-}"
INDEX_TEMPLATE="${AOU_GAMMA_INDEX_TEMPLATE:-}"
MASK_TEMPLATE="${AOU_GAMMA_MASK_TEMPLATE:-}"
ANCESTRY_URI="${AOU_GAMMA_ANCESTRY_URI:-gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv}"
QC_EXCLUSIONS_URI="${AOU_GAMMA_QC_EXCLUSIONS_URI:-gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv}"
RELATEDNESS_EXCLUSIONS_URI="${AOU_GAMMA_RELATEDNESS_EXCLUSIONS_URI:-gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv}"
GENE_ANNOTATION_URI="${AOU_GAMMA_GENE_ANNOTATION_URI:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz}"
THREADS="${AOU_GAMMA_THREADS:-12}"
THETA="${AOU_GAMMA_THETA:-0.00075}"
RHO_OVER_THETA="${AOU_GAMMA_RHO_OVER_THETA:-0.8}"
MUTATION_RATE="${AOU_GAMMA_MUTATION_RATE:-1.29e-8}"
GENERATION_TIME="${AOU_GAMMA_GENERATION_TIME:-25}"
THRESHOLD_YEARS="${AOU_GAMMA_THRESHOLD_YEARS:-4500}"
RECENT_CALL="${AOU_GAMMA_RECENT_CALL:-mean}"
OUTPUT_STRIDE="${AOU_GAMMA_OUTPUT_STRIDE:-10000}"
CACHE_SIZE="${AOU_GAMMA_CACHE_SIZE:-1000}"
PAIR_BLOCK="${AOU_GAMMA_PAIR_BLOCK:-256}"
N_RANDOM_PAIRS="${AOU_GAMMA_N_RANDOM_PAIRS:-100000}"
PAIRS_SEED="${AOU_GAMMA_PAIRS_SEED:-1729}"
EXCLUDE_WITHIN=0
SIGNAL_FRACTION="${AOU_GAMMA_SIGNAL_FRACTION:-0.02}"
MERGE_GAP="${AOU_GAMMA_MERGE_GAP:-20000}"
PROFILE_HALF_WIDTH="${AOU_GAMMA_PROFILE_HALF_WIDTH:-500000}"
VARIANT_HALF_WIDTH="${AOU_GAMMA_VARIANT_HALF_WIDTH:-100000}"
MIN_GENOTYPE_PAIRS="${AOU_GAMMA_MIN_GENOTYPE_PAIRS:-20}"
EXP10="${AOU_GAMMA_EXP10:-accurate}"
BACKWARD_ALIGNMENT="${AOU_GAMMA_BACKWARD_ALIGNMENT:-fixed}"
TOP_N="${AOU_GAMMA_TOP_N:-100}"
HIT_BIN_SIZE="${AOU_GAMMA_HIT_BIN_SIZE:-1000000}"
GENE_CONTEXT_FLANK="${AOU_GAMMA_GENE_CONTEXT_FLANK:-500000}"
ZOOM_YMAX="${AOU_GAMMA_ZOOM_YMAX:-0.04}"
HIT_LABEL_MIN_FRACTION="${AOU_GAMMA_HIT_LABEL_MIN_FRACTION:-0.02}"
MASK_ENABLED=1
UPLOAD=1
KEEP_INPUTS=0
FORCE=0
DRY_RUN=0
ALLOW_DIRTY=0

usage() {
    cat <<'EOF'
Usage: scripts/run_aou_workbench.sh -chr CHR|all -pops POP|all [options]

Required selection (case-insensitive; comma lists are also accepted):
  -chr, --chr 1|all         Run one chromosome, a comma list, or autosomes 1-22.
  -pops, --pops AFR|all     Run one population, a comma list, or all six
                            (AFR, AMR, EAS, EUR, MID, SAS).

Cloud and local paths:
  --input-prefix URI        Optional value for {input_prefix} in custom templates
  --output-prefix URI       Default: gs://rw-migration-aou-rw-fa99430f/
                            gamma_smc/results (WORKSPACE_BUCKET overrides)
  --billing-project ID      Requester-pays billing project. Default:
                            AOU_GAMMA_BILLING_PROJECT or GOOGLE_PROJECT
  --local-root PATH         Default: /home/jupyter/gamma_smc_workbench
  --bcf-template TEMPLATE   Default: AoU lrWGS phase-2 bubble-split chr BCF
  --index-template TEMPLATE Default: {bcf}.csi
  --mask-template TEMPLATE  Default: gs://rw-migration-aou-rw-fa99430f/
                            hardmask.hg38.v4.over99.bed
  --ancestry-uri URI        Default: v9 ancestry_preds.tsv
  --qc-exclusions-uri URI   Default: v9 QC flagged_samples.tsv
  --relatedness-exclusions-uri URI
                            Default: v9 relatedness_flagged_samples.tsv
  --gene-annotation-uri URI Default: GENCODE v50 basic GRCh38 GTF
  --no-mask                 Decode without a callable-region BED.

Decoder parameters:
  --threads N               Default: 12
  --theta X                 Default: 0.00075
  --rho-over-theta X        Default: 0.8
  --mutation-rate X         Default: 1.29e-8
  --generation-time X       Default: 25
  --threshold-years X       Default: 4500
  --recent-call RULE        Default: mean (mean or median)
  --output-at-stride N      Default: 10000 bp
  --cache-size N            Default: 1000 bp
  --pair-block N            Default: 256
  --n-random-pairs N        Default: 100000 haplotype pairs within each pop
  --pairs-seed N            Default: 1729
  --exclude-within          Exclude the same person's two haplotypes from draw
  --top-n N                 Default: 100 whole-genome windows per statistic
  --hit-bin-size N          Merge top-window hits through consecutive 1 Mb bins
  --gene-context-flank N    List protein-coding genes within +/-500000 bp
  --zoom-ymax X             Separate genome plot y ceiling; default: 0.04
  --hit-label-min-fraction X
                            Label ranked peaks strictly above X; default: 0.02

Candidate analysis:
  --signal-fraction X       Strict screen threshold; default: 0.02
  --merge-gap N             Merge signal-window gaps up to 20000 bp
  --profile-half-width N    Pair-TMRCA profile +/-500000 bp around each peak
  --variant-half-width N    Rank variants +/-100000 bp around each peak
  --min-genotype-pairs N    Require >=20 ref/ref and >=20 alt/alt decoded pairs

Run control:
  --force                   Ignore a valid completion record and decode again.
  --keep-inputs             Retain staged chromosome BCFs after upload.
  --no-upload               Do not upload aggregate outputs or plots.
  --allow-dirty             Permit tracked local source changes.
  --dry-run                 Print resolved work without accessing GCS.
  -h, --help                Show this help.

Default controlled inputs:
  gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/
    panel_bubble_split_vcf/aou_lr_phase2_v1.chr{chr}.bubble.split.bcf[.csi]
  ancestry_pred_other from the v9 ancestry table; QC and related samples removed
  gs://rw-migration-aou-rw-fa99430f/hardmask.hg38.v4.over99.bed

Default workspace outputs:
  gamma_smc/results/AFR/chromosomes/...
  gamma_smc/results/AFR/plots/...
  gamma_smc/results/summary/{scope}/...

Templates recognize {input_prefix}, {chr}, and (for the index) {bcf}. The
full-panel BCF is phased and may contain SNVs, indels, and SVs;
Gamma-SMC retains only biallelic segregating SNPs after population/QC subsetting.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

need_value() {
    [[ $# -ge 2 ]] || die "$1 requires a value"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -chr|--chr) need_value "$@"; CHR_SPEC="$2"; shift 2 ;;
        -pops|--pops) need_value "$@"; POPS_SPEC="$2"; shift 2 ;;
        --input-prefix) need_value "$@"; INPUT_PREFIX="$2"; shift 2 ;;
        --output-prefix) need_value "$@"; OUTPUT_PREFIX="$2"; shift 2 ;;
        --billing-project) need_value "$@"; BILLING_PROJECT="$2"; shift 2 ;;
        --local-root) need_value "$@"; LOCAL_ROOT="$2"; shift 2 ;;
        --bcf-template) need_value "$@"; BCF_TEMPLATE="$2"; shift 2 ;;
        --index-template) need_value "$@"; INDEX_TEMPLATE="$2"; shift 2 ;;
        --mask-template) need_value "$@"; MASK_TEMPLATE="$2"; shift 2 ;;
        --ancestry-uri) need_value "$@"; ANCESTRY_URI="$2"; shift 2 ;;
        --qc-exclusions-uri) need_value "$@"; QC_EXCLUSIONS_URI="$2"; shift 2 ;;
        --relatedness-exclusions-uri) need_value "$@"; RELATEDNESS_EXCLUSIONS_URI="$2"; shift 2 ;;
        --gene-annotation-uri) need_value "$@"; GENE_ANNOTATION_URI="$2"; shift 2 ;;
        --threads) need_value "$@"; THREADS="$2"; shift 2 ;;
        --theta) need_value "$@"; THETA="$2"; shift 2 ;;
        --rho-over-theta) need_value "$@"; RHO_OVER_THETA="$2"; shift 2 ;;
        --mutation-rate) need_value "$@"; MUTATION_RATE="$2"; shift 2 ;;
        --generation-time) need_value "$@"; GENERATION_TIME="$2"; shift 2 ;;
        --threshold-years) need_value "$@"; THRESHOLD_YEARS="$2"; shift 2 ;;
        --recent-call) need_value "$@"; RECENT_CALL="${2,,}"; shift 2 ;;
        --output-at-stride) need_value "$@"; OUTPUT_STRIDE="$2"; shift 2 ;;
        --cache-size) need_value "$@"; CACHE_SIZE="$2"; shift 2 ;;
        --pair-block) need_value "$@"; PAIR_BLOCK="$2"; shift 2 ;;
        --n-random-pairs) need_value "$@"; N_RANDOM_PAIRS="$2"; shift 2 ;;
        --pairs-seed) need_value "$@"; PAIRS_SEED="$2"; shift 2 ;;
        --exclude-within) EXCLUDE_WITHIN=1; shift ;;
        --signal-fraction) need_value "$@"; SIGNAL_FRACTION="$2"; shift 2 ;;
        --merge-gap) need_value "$@"; MERGE_GAP="$2"; shift 2 ;;
        --profile-half-width) need_value "$@"; PROFILE_HALF_WIDTH="$2"; shift 2 ;;
        --variant-half-width) need_value "$@"; VARIANT_HALF_WIDTH="$2"; shift 2 ;;
        --min-genotype-pairs) need_value "$@"; MIN_GENOTYPE_PAIRS="$2"; shift 2 ;;
        --top-n) need_value "$@"; TOP_N="$2"; shift 2 ;;
        --hit-bin-size) need_value "$@"; HIT_BIN_SIZE="$2"; shift 2 ;;
        --gene-context-flank) need_value "$@"; GENE_CONTEXT_FLANK="$2"; shift 2 ;;
        --zoom-ymax) need_value "$@"; ZOOM_YMAX="$2"; shift 2 ;;
        --hit-label-min-fraction) need_value "$@"; HIT_LABEL_MIN_FRACTION="$2"; shift 2 ;;
        --no-mask) MASK_ENABLED=0; shift ;;
        --force) FORCE=1; shift ;;
        --keep-inputs) KEEP_INPUTS=1; shift ;;
        --no-upload) UPLOAD=0; shift ;;
        --allow-dirty) ALLOW_DIRTY=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

[[ -n "$CHR_SPEC" ]] || die "-chr is required (use an autosome or all)"
[[ -n "$POPS_SPEC" ]] || die "-pops is required (use a population or all)"

WORKSPACE_BUCKET_VALUE="${WORKSPACE_BUCKET:-gs://rw-migration-aou-rw-fa99430f}"
if [[ -z "$OUTPUT_PREFIX" && -n "$WORKSPACE_BUCKET_VALUE" ]]; then
    OUTPUT_PREFIX="${WORKSPACE_BUCKET_VALUE%/}/gamma_smc/results"
fi
if [[ -n "$INPUT_PREFIX" ]]; then
    [[ "$INPUT_PREFIX" == gs://* ]] || die "input prefix must be a gs:// URI"
fi
if [[ "$UPLOAD" -eq 1 ]]; then
    [[ -n "$OUTPUT_PREFIX" ]] || die \
        "set WORKSPACE_BUCKET or pass --output-prefix gs://.../gamma_smc/results"
    [[ "$OUTPUT_PREFIX" == gs://* ]] || die "output prefix must be a gs:// URI"
fi
INPUT_PREFIX="${INPUT_PREFIX%/}"
OUTPUT_PREFIX="${OUTPUT_PREFIX%/}"
[[ -n "$BCF_TEMPLATE" ]] || \
    BCF_TEMPLATE='gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr{chr}.bubble.split.bcf'
[[ -n "$INDEX_TEMPLATE" ]] || INDEX_TEMPLATE='{bcf}.csi'
[[ -n "$MASK_TEMPLATE" ]] || \
    MASK_TEMPLATE='gs://rw-migration-aou-rw-fa99430f/hardmask.hg38.v4.over99.bed'
for source_uri in "$ANCESTRY_URI" "$QC_EXCLUSIONS_URI" \
    "$RELATEDNESS_EXCLUSIONS_URI"; do
    [[ "$source_uri" == gs://* ]] || die "controlled input must be a gs:// URI: $source_uri"
done

for integer_setting in THREADS OUTPUT_STRIDE CACHE_SIZE PAIR_BLOCK TOP_N \
    N_RANDOM_PAIRS PROFILE_HALF_WIDTH VARIANT_HALF_WIDTH MIN_GENOTYPE_PAIRS \
    HIT_BIN_SIZE; do
    value="${!integer_setting}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "$integer_setting must be a positive integer"
done
[[ "$PAIRS_SEED" =~ ^[0-9]+$ ]] || die "PAIRS_SEED must be a nonnegative integer"
[[ "$MERGE_GAP" =~ ^[0-9]+$ ]] || die "MERGE_GAP must be a nonnegative integer"
[[ "$GENE_CONTEXT_FLANK" =~ ^[0-9]+$ ]] || die \
    "GENE_CONTEXT_FLANK must be a nonnegative integer"
awk -v value="$SIGNAL_FRACTION" 'BEGIN { exit !(value >= 0 && value <= 1) }' || \
    die "SIGNAL_FRACTION must be in [0,1]"
awk -v value="$ZOOM_YMAX" 'BEGIN { exit !(value > 0 && value <= 1) }' || \
    die "ZOOM_YMAX must be in (0,1]"
[[ "$RECENT_CALL" =~ ^(mean|median)$ ]] || die \
    "--recent-call must be mean or median"

declare -a CHROMOSOMES=()
if [[ "${CHR_SPEC,,}" == "all" ]]; then
    CHROMOSOMES=("${ALL_CHROMOSOMES[@]}")
else
    declare -A SEEN_CHROMOSOMES=()
    IFS=',' read -r -a requested_chromosomes <<< "$CHR_SPEC"
    for chromosome in "${requested_chromosomes[@]}"; do
        chromosome="${chromosome//[[:space:]]/}"
        chromosome="${chromosome#chr}"
        chromosome="${chromosome#CHR}"
        [[ "$chromosome" =~ ^[0-9]+$ ]] || die "invalid chromosome: $chromosome"
        chromosome_number="$((10#$chromosome))"
        (( chromosome_number >= 1 && chromosome_number <= 22 )) || die \
            "chromosome must be between 1 and 22: $chromosome"
        if [[ -z "${SEEN_CHROMOSOMES[$chromosome_number]:-}" ]]; then
            CHROMOSOMES+=("$chromosome_number")
            SEEN_CHROMOSOMES[$chromosome_number]=1
        fi
    done
fi

declare -a POPULATIONS=()
if [[ "${POPS_SPEC,,}" == "all" ]]; then
    POPULATIONS=("${ALL_POPS[@]}")
else
    declare -A SEEN_POPS=()
    IFS=',' read -r -a requested_pops <<< "$POPS_SPEC"
    for population in "${requested_pops[@]}"; do
        population="${population//[[:space:]]/}"
        population="${population^^}"
        supported=0
        for candidate in "${ALL_POPS[@]}"; do
            [[ "$population" == "$candidate" ]] && supported=1
        done
        [[ "$supported" -eq 1 ]] || die "unsupported population: $population"
        if [[ -z "${SEEN_POPS[$population]:-}" ]]; then
            POPULATIONS+=("$population")
            SEEN_POPS[$population]=1
        fi
    done
fi

expand_template() {
    local template="$1" population="$2" chromosome="$3" bcf="${4:-}"
    local value="$template"
    value="${value//\{input_prefix\}/$INPUT_PREFIX}"
    value="${value//\{POP\}/$population}"
    value="${value//\{pop\}/${population,,}}"
    value="${value//\{chr\}/$chromosome}"
    value="${value//\{bcf\}/$bcf}"
    printf '%s\n' "$value"
}

print_plan() {
    echo "Gamma-SMC Workbench plan"
    echo "  populations: ${POPULATIONS[*]}"
    echo "  chromosomes: ${CHROMOSOMES[*]}"
    echo "  local root: $LOCAL_ROOT"
    echo "  requester-pays billing project: ${BILLING_PROJECT:-<unset>}"
    echo "  decoder: threads=$THREADS, theta=$THETA, rho/theta=$RHO_OVER_THETA"
    echo "  statistic: recent_call=$RECENT_CALL, threshold=$THRESHOLD_YEARS years"
    echo "  grid/cache: stride=$OUTPUT_STRIDE bp, cache=$CACHE_SIZE bp"
    echo "  pair draw: $N_RANDOM_PAIRS random haplotype pairs/pop, seed=$PAIRS_SEED, exclude_within=$EXCLUDE_WITHIN"
    echo "  candidates: fraction>$SIGNAL_FRACTION, merge_gap=$MERGE_GAP bp, profile=+/-$PROFILE_HALF_WIDTH bp, variants=+/-$VARIANT_HALF_WIDTH bp"
    echo "  ranked-hit report: top_n=$TOP_N, bins=$HIT_BIN_SIZE bp, gene_flank=+/-$GENE_CONTEXT_FLANK bp, zoom_ymax=$ZOOM_YMAX, label_min=$HIT_LABEL_MIN_FRACTION"
    echo "  gene annotation: $GENE_ANNOTATION_URI"
    echo "  ancestry: $ANCESTRY_URI (column ancestry_pred_other)"
    echo "  QC exclusions: $QC_EXCLUSIONS_URI"
    echo "  relatedness exclusions: $RELATEDNESS_EXCLUSIONS_URI"
    for chromosome in "${CHROMOSOMES[@]}"; do
        input_uri="$(expand_template "$BCF_TEMPLATE" "AFR" "$chromosome")"
        mask_uri="$(expand_template "$MASK_TEMPLATE" "AFR" "$chromosome")"
        echo "  chr$chromosome full-panel input: $input_uri"
        [[ "$MASK_ENABLED" -eq 1 ]] && echo "  chr$chromosome mask: $mask_uri"
    done
    if [[ "$UPLOAD" -eq 1 ]]; then
        for population in "${POPULATIONS[@]}"; do
            echo "  $population output: $OUTPUT_PREFIX/$population/{chromosomes,plots}/"
        done
        echo "  combined report: $OUTPUT_PREFIX/summary/{scope}/"
    fi
}

print_plan
if [[ "$DRY_RUN" -eq 1 ]]; then
    exit 0
fi

command -v gcloud >/dev/null 2>&1 || die "gcloud is required to stage workspace objects"
[[ -n "$BILLING_PROJECT" ]] || die \
    "set GOOGLE_PROJECT, AOU_GAMMA_BILLING_PROJECT, or --billing-project for requester-pays GCS access"
[[ -x "$AOU" ]] || die "runner is missing or not executable: $AOU"
[[ -x "$REPO/bin/gamma_smc" ]] || die \
    "Gamma-SMC is not built; run scripts/bootstrap_uv.sh first"
BCFTOOLS="${BCFTOOLS_BIN:-$REPO/.native/bin/bcftools}"
[[ -x "$BCFTOOLS" ]] || die "bcftools is not installed; rerun scripts/bootstrap_uv.sh"
grep -qw avx2 /proc/cpuinfo || die "Gamma-SMC requires an AVX2-capable CPU"

CODE_COMMIT="$(git -C "$REPO" rev-parse HEAD)"
if [[ "$ALLOW_DIRTY" -eq 0 ]] && {
    ! git -C "$REPO" diff --quiet || ! git -C "$REPO" diff --cached --quiet
}; then
    die "tracked repository files are modified; commit them or use --allow-dirty"
fi
if [[ "$ALLOW_DIRTY" -eq 1 ]] && {
    ! git -C "$REPO" diff --quiet || ! git -C "$REPO" diff --cached --quiet
}; then
    dirty_digest="$(git -C "$REPO" diff --binary HEAD | sha256sum | cut -d' ' -f1)"
    CODE_COMMIT="$CODE_COMMIT-dirty-$dirty_digest"
fi

mkdir -p "$LOCAL_ROOT"
LOCAL_ROOT="$(cd "$LOCAL_ROOT" && pwd)"
command -v flock >/dev/null 2>&1 || die "flock is required for restart-safe local locking"
exec 9>"$LOCAL_ROOT/.run_aou_workbench.lock"
flock -n 9 || die "another Workbench runner is already using $LOCAL_ROOT"

cloud_exists() {
    gcloud storage objects describe "$1" \
        --billing-project "$BILLING_PROJECT" >/dev/null 2>&1
}

object_fingerprint() {
    local uri="$1" description
    if ! description="$(
        gcloud storage objects describe "$uri" \
            --format='value(generation,size,crc32c_hash)' \
            --billing-project "$BILLING_PROJECT"
    )"; then
        die "required cloud object is absent or inaccessible: $uri"
    fi
    [[ -n "$description" ]] || die "could not fingerprint cloud object: $uri"
    printf '%s\n' "$description"
}

staged_object_is_current() {
    local uri="$1" destination="$2" fingerprint="$3"
    local source_record="${destination}.cloud-source"
    local expected observed generation expected_size crc32c
    expected="$(printf '%s\n%s' "$uri" "$fingerprint")"
    [[ -s "$destination" && -f "$source_record" ]] || return 1
    observed="$(<"$source_record")"
    IFS=$'\t ' read -r generation expected_size crc32c <<< "$fingerprint"
    [[ "$observed" == "$expected" && "$expected_size" =~ ^[0-9]+$ ]] && \
        [[ "$(stat -c '%s' "$destination")" == "$expected_size" ]]
}

stage_object() {
    local uri="$1" destination="$2" fingerprint="$3"
    local source_record="${destination}.cloud-source" temporary="${destination}.partial"
    local partial_record="${destination}.partial.cloud-source"
    local expected partial_observed="" generation expected_size crc32c
    expected="$(printf '%s\n%s' "$uri" "$fingerprint")"
    IFS=$'\t ' read -r generation expected_size crc32c <<< "$fingerprint"
    if staged_object_is_current "$uri" "$destination" "$fingerprint"; then
        echo "Reusing staged object: $destination"
        return 0
    fi
    mkdir -p "$(dirname "$destination")"
    # This is a runner-owned cache entry and cannot satisfy the current source
    # contract. Removing it first avoids needing space for two full panel BCFs.
    rm -f -- "$destination" "$source_record"
    if [[ -f "$partial_record" ]]; then
        partial_observed="$(<"$partial_record")"
    fi
    if [[ -e "$temporary" && "$partial_observed" != "$expected" ]]; then
        rm -f -- "$temporary" "$partial_record"
    fi
    printf '%s\n%s\n' "$uri" "$fingerprint" > "$partial_record"
    echo "Staging $uri"
    if ! gcloud storage cp "$uri" "$temporary" \
        --billing-project "$BILLING_PROJECT"; then
        return 1
    fi
    mv -f -- "$temporary" "$destination"
    printf '%s\n%s\n' "$uri" "$fingerprint" > "${source_record}.partial.$$"
    mv -f -- "${source_record}.partial.$$" "$source_record"
    rm -f -- "$partial_record"
}

stage_gene_annotation() {
    local uri="$1" destination="$2"
    local source_record="${destination}.source-uri"
    local temporary="${destination}.partial.$$"
    if [[ "$uri" == gs://* ]]; then
        local fingerprint
        fingerprint="$(object_fingerprint "$uri")"
        stage_object "$uri" "$destination" "$fingerprint"
        return 0
    fi
    if [[ "$uri" == https://* ]]; then
        if [[ -s "$destination" && -f "$source_record" ]] && \
            [[ "$(<"$source_record")" == "$uri" ]]; then
            echo "Reusing staged gene annotation: $destination"
            return 0
        fi
        command -v curl >/dev/null 2>&1 || die \
            "curl is required to stage the GENCODE gene annotation"
        mkdir -p "$(dirname "$destination")"
        rm -f -- "$temporary"
        echo "Staging $uri"
        if ! curl --fail --location --retry 3 --output "$temporary" "$uri"; then
            rm -f -- "$temporary"
            return 1
        fi
        mv -f -- "$temporary" "$destination"
        printf '%s\n' "$uri" > "${source_record}.partial.$$"
        mv -f -- "${source_record}.partial.$$" "$source_record"
        return 0
    fi
    [[ -s "$uri" ]] || die "gene annotation is absent or empty: $uri"
}

restore_if_present() {
    local uri="$1" destination="$2"
    local temporary="${destination}.partial.$$"
    [[ -s "$destination" ]] && return 0
    cloud_exists "$uri" || return 1
    mkdir -p "$(dirname "$destination")"
    rm -f -- "$temporary"
    if ! gcloud storage cp "$uri" "$temporary" \
        --billing-project "$BILLING_PROJECT"; then
        rm -f -- "$temporary"
        return 1
    fi
    mv -f -- "$temporary" "$destination"
}

remote_matches_file() {
    local uri="$1" local_file="$2"
    local temporary="${local_file}.remote.$$"
    [[ -s "$local_file" ]] || return 1
    rm -f -- "$temporary"
    if ! gcloud storage cp "$uri" "$temporary" \
        --billing-project "$BILLING_PROJECT" >/dev/null 2>&1; then
        rm -f -- "$temporary"
        return 1
    fi
    if cmp -s -- "$temporary" "$local_file"; then
        rm -f -- "$temporary"
        return 0
    fi
    rm -f -- "$temporary"
    return 1
}

make_upload_stage() {
    local stage_root="$LOCAL_ROOT/tmp/upload-staging"
    mkdir -p "$stage_root"
    mktemp -d "$stage_root/stage.XXXXXX"
}

safe_remove_upload_stage() {
    local stage="$1" stage_root resolved_stage resolved_root
    stage_root="$LOCAL_ROOT/tmp/upload-staging"
    resolved_stage="$(realpath -m -- "$stage")"
    resolved_root="$(realpath -m -- "$stage_root")"
    case "${resolved_stage}/" in
        "${resolved_root}/"stage.*/) rm -rf -- "$resolved_stage" ;;
        *) die "refusing to remove upload staging path outside $resolved_root: $stage" ;;
    esac
}

rsync_directory() {
    local source_directory="$1" remote_directory="$2"
    shift 2
    echo "Synchronizing $source_directory -> $remote_directory"
    # Checksums, rather than local/cloud mtimes, make an unchanged rerun a
    # metadata comparison instead of another upload. Extra remote objects are
    # deliberately retained; a scoped run must never delete another scope.
    gcloud storage rsync "$source_directory" "$remote_directory" \
        --recursive --checksums-only "$@" \
        --billing-project "$BILLING_PROJECT"
}

rsync_files() {
    local remote_directory="$1" stage source
    shift
    [[ "$#" -gt 0 ]] || die "rsync_files requires at least one source"
    stage="$(make_upload_stage)"
    for source in "$@"; do
        if ! ln -- "$source" "$stage/$(basename "$source")"; then
            safe_remove_upload_stage "$stage"
            return 1
        fi
    done
    if ! rsync_directory "$stage" "$remote_directory"; then
        safe_remove_upload_stage "$stage"
        return 1
    fi
    safe_remove_upload_stage "$stage"
}

rsync_single_file() {
    local source="$1" remote_directory="$2"
    rsync_files "$remote_directory" "$source"
}

upload_completed_chromosome() {
    local remote_directory="$1" summary="$2" run_json="$3" pairs_manifest="$4"
    local sample_list="$5" sample_audit="$6" decode_log="$7" bitmatrix="$8"
    local regions="$9" positions="${10}" candidate_directory="${11}"
    local completion="${12}" stage source
    local -a outputs=(
        "$summary" "$run_json" "$pairs_manifest" "$sample_list" "$sample_audit"
        "$bitmatrix" "${bitmatrix}.meta" "$regions" "$positions"
    )
    stage="$(make_upload_stage)"
    for source in "${outputs[@]}"; do
        if ! ln -- "$source" "$stage/$(basename "$source")"; then
            safe_remove_upload_stage "$stage"
            return 1
        fi
    done
    if [[ -s "$decode_log" ]] && \
        ! ln -- "$decode_log" "$stage/$(basename "$decode_log")"; then
        safe_remove_upload_stage "$stage"
        return 1
    fi
    # Hard links stage metadata only: no result data are duplicated locally.
    # Keeping the basename preserves chrN.candidates/ in Cloud Storage and
    # prevents different chromosomes from overwriting one shared directory.
    if ! cp -al -- "$candidate_directory" "$stage/"; then
        safe_remove_upload_stage "$stage"
        return 1
    fi
    if ! rsync_directory "$stage" "$remote_directory"; then
        safe_remove_upload_stage "$stage"
        return 1
    fi
    safe_remove_upload_stage "$stage"
    # Completion is the commit marker and is deliberately synchronized last.
    rsync_single_file "$completion" "$remote_directory"
}

safe_clear_run() {
    local result_directory="$1"
    shift
    case "${result_directory}/" in
        "${LOCAL_ROOT}/results/"*) ;;
        *) die "refusing to clear outputs outside local result root: $result_directory" ;;
    esac
    rm -f -- "$@"
}

safe_clear_candidate_directory() {
    local candidate_directory="$1"
    case "${candidate_directory}/" in
        "${LOCAL_ROOT}/results/"*) ;;
        *) die "refusing to clear candidate outputs outside local result root: $candidate_directory" ;;
    esac
    rm -rf -- "$candidate_directory"
}

check_staging_capacity() {
    local input_uri="$1" input_fingerprint="$2" local_input="$3"
    local index_uri="$4" index_fingerprint="$5" local_index="$6"
    local input_generation input_size input_crc index_generation index_size index_crc
    local available reclaimable=0 reserve=$((20 * 1024 * 1024 * 1024))
    local required capacity
    IFS=$'\t ' read -r input_generation input_size input_crc <<< "$input_fingerprint"
    IFS=$'\t ' read -r index_generation index_size index_crc <<< "$index_fingerprint"
    [[ "$input_size" =~ ^[0-9]+$ && "$index_size" =~ ^[0-9]+$ ]] || die \
        "could not read BCF/index sizes from their GCS metadata"
    available="$(df -PB1 "$LOCAL_ROOT" | awk 'NR==2 {print $4}')"
    [[ "$available" =~ ^[0-9]+$ ]] || die "could not determine free staging space"
    required=$reserve
    if ! staged_object_is_current "$input_uri" "$local_input" "$input_fingerprint"; then
        required=$((required + input_size))
        [[ -f "$local_input" ]] && \
            reclaimable=$((reclaimable + $(stat -c '%s' "$local_input")))
    fi
    if ! staged_object_is_current "$index_uri" "$local_index" "$index_fingerprint"; then
        required=$((required + index_size))
        [[ -f "$local_index" ]] && \
            reclaimable=$((reclaimable + $(stat -c '%s' "$local_index")))
    fi
    capacity=$((available + reclaimable))
    if (( capacity < required )); then
        die "insufficient local disk for chromosome staging: need $((required / 1073741824)) GiB including a 20-GiB reserve, have $((capacity / 1073741824)) GiB"
    fi
    echo "Disk preflight: $((required / 1073741824)) GiB required including reserve; $((capacity / 1073741824)) GiB available."
}

check_candidate_capacity() {
    local n_positions="$1" n_pairs="$2"
    local available reserve=$((5 * 1024 * 1024 * 1024)) required
    # Worst case while converting: raw alpha/beta (8 bytes/cell), an
    # uncompressed float32 mean array, and its compressed output (4+4 bytes).
    required=$((n_positions * n_pairs * 16 + reserve))
    available="$(df -PB1 "$LOCAL_ROOT" | awk 'NR==2 {print $4}')"
    [[ "$available" =~ ^[0-9]+$ ]] || die "could not determine candidate-analysis disk space"
    if (( available < required )); then
        die "insufficient disk for targeted pair-TMRCA replay: need $((required / 1073741824)) GiB including reserve, have $((available / 1073741824)) GiB"
    fi
    echo "Candidate disk preflight: $((required / 1073741824)) GiB worst-case; $((available / 1073741824)) GiB available."
}

metadata_directory="$LOCAL_ROOT/inputs/metadata"
local_ancestry="$metadata_directory/ancestry_preds.tsv"
local_qc_exclusions="$metadata_directory/flagged_samples.tsv"
local_relatedness_exclusions="$metadata_directory/relatedness_flagged_samples.tsv"
ancestry_fingerprint="$(object_fingerprint "$ANCESTRY_URI")"
qc_exclusions_fingerprint="$(object_fingerprint "$QC_EXCLUSIONS_URI")"
relatedness_exclusions_fingerprint="$(
    object_fingerprint "$RELATEDNESS_EXCLUSIONS_URI"
)"
stage_object "$ANCESTRY_URI" "$local_ancestry" "$ancestry_fingerprint"
stage_object "$QC_EXCLUSIONS_URI" "$local_qc_exclusions" \
    "$qc_exclusions_fingerprint"
stage_object "$RELATEDNESS_EXCLUSIONS_URI" "$local_relatedness_exclusions" \
    "$relatedness_exclusions_fingerprint"
if [[ "$GENE_ANNOTATION_URI" == gs://* || "$GENE_ANNOTATION_URI" == https://* ]]; then
    gene_annotation_name="$(basename "${GENE_ANNOTATION_URI%%\?*}")"
    local_gene_annotation="$LOCAL_ROOT/inputs/annotations/$gene_annotation_name"
    stage_gene_annotation "$GENE_ANNOTATION_URI" "$local_gene_annotation"
else
    local_gene_annotation="$(realpath -- "$GENE_ANNOTATION_URI")"
    stage_gene_annotation "$local_gene_annotation" "$local_gene_annotation"
fi

for chromosome in "${CHROMOSOMES[@]}"; do
    echo
    echo "=== Full panel chromosome $chromosome ==="
    input_uri="$(expand_template "$BCF_TEMPLATE" "PANEL" "$chromosome")"
    index_uri="$(expand_template "$INDEX_TEMPLATE" "PANEL" "$chromosome" "$input_uri")"
    mask_uri="$(expand_template "$MASK_TEMPLATE" "PANEL" "$chromosome")"
    [[ "$input_uri" == gs://* ]] || die "BCF template did not resolve to gs://: $input_uri"
    [[ "$index_uri" == gs://* ]] || die "index template did not resolve to gs://: $index_uri"
    if [[ "$MASK_ENABLED" -eq 1 ]]; then
        [[ "$mask_uri" == gs://* ]] || die \
            "mask template did not resolve to gs://: $mask_uri"
    fi
    input_fingerprint="$(object_fingerprint "$input_uri")"
    index_fingerprint="$(object_fingerprint "$index_uri")"
    mask_fingerprint=""
    if [[ "$MASK_ENABLED" -eq 1 ]]; then
        mask_fingerprint="$(object_fingerprint "$mask_uri")"
    fi

    input_directory="$LOCAL_ROOT/inputs/panel/chr$chromosome"
    local_input="$input_directory/$(basename "$input_uri")"
    local_index="$input_directory/$(basename "$index_uri")"
    bcf_samples="$input_directory/chr$chromosome.bcf_samples.txt"
    local_mask_source="$LOCAL_ROOT/inputs/masks/source/$(basename "$mask_uri")"
    local_mask="$LOCAL_ROOT/inputs/masks/callable/chr$chromosome.callable.bed"
    mask_audit="$LOCAL_ROOT/inputs/masks/callable/chr$chromosome.callable.audit.json"
    input_staged=0

    for population in "${POPULATIONS[@]}"; do
        echo
        echo "--- $population chromosome $chromosome ---"
        population_summary_dir="$LOCAL_ROOT/results/$population/chromosomes"
        mkdir -p "$population_summary_dir"
        summary="$population_summary_dir/chr$chromosome.gamma_smc.tsv"
        run_json="$summary.run.json"
        pairs_manifest="$population_summary_dir/chr$chromosome.pairs.tsv"
        bitmatrix="$population_summary_dir/chr$chromosome.recent.bits"
        sample_list="$population_summary_dir/chr$chromosome.samples.txt"
        sample_audit="$population_summary_dir/chr$chromosome.samples.audit.json"
        candidate_regions="$population_summary_dir/chr$chromosome.candidate_regions.tsv"
        candidate_positions="$population_summary_dir/chr$chromosome.candidate_positions.txt"
        candidate_directory="$population_summary_dir/chr$chromosome.candidates"
        candidate_manifest="$candidate_directory/candidate_analysis.json"
        candidate_summary="$candidate_directory/candidate_decode.gamma_smc.tsv"
        candidate_raw_directory="$LOCAL_ROOT/tmp/$population/chr$chromosome"
        candidate_raw="$candidate_raw_directory/candidate_decode.posteriors.zst"
        completion="$population_summary_dir/chr$chromosome.complete.json"
        decode_log="$population_summary_dir/chr$chromosome.decode.log"
        remote_directory="$OUTPUT_PREFIX/$population/chromosomes"

        validation_args=(
            workbench-validate
            --population "$population"
            --chromosome "$chromosome"
            --input-uri "$input_uri"
            --input-fingerprint "$input_fingerprint"
            --local-input "$local_input"
            --index-uri "$index_uri"
            --index-fingerprint "$index_fingerprint"
            --local-index "$local_index"
            --summary "$summary"
            --run-json "$run_json"
            --pairs-manifest "$pairs_manifest"
            --bitmatrix "$bitmatrix"
            --candidate-regions "$candidate_regions"
            --candidate-manifest "$candidate_manifest"
            --sample-list "$sample_list"
            --sample-audit "$sample_audit"
            --ancestry-uri "$ANCESTRY_URI"
            --ancestry-fingerprint "$ancestry_fingerprint"
            --local-ancestry "$local_ancestry"
            --qc-exclusions-uri "$QC_EXCLUSIONS_URI"
            --qc-exclusions-fingerprint "$qc_exclusions_fingerprint"
            --local-qc-exclusions "$local_qc_exclusions"
            --relatedness-exclusions-uri "$RELATEDNESS_EXCLUSIONS_URI"
            --relatedness-exclusions-fingerprint "$relatedness_exclusions_fingerprint"
            --local-relatedness-exclusions "$local_relatedness_exclusions"
            --completion "$completion"
            --theta "$THETA"
            --rho-over-theta "$RHO_OVER_THETA"
            --mutation-rate "$MUTATION_RATE"
            --generation-time "$GENERATION_TIME"
            --threshold-years "$THRESHOLD_YEARS"
            --recent-call "$RECENT_CALL"
            --output-at-stride "$OUTPUT_STRIDE"
            --cache-size "$CACHE_SIZE"
            --threads "$THREADS"
            --pair-block "$PAIR_BLOCK"
            --n-random-pairs "$N_RANDOM_PAIRS"
            --pairs-seed "$PAIRS_SEED"
            --signal-fraction "$SIGNAL_FRACTION"
            --merge-gap "$MERGE_GAP"
            --profile-half-width "$PROFILE_HALF_WIDTH"
            --variant-half-width "$VARIANT_HALF_WIDTH"
            --minimum-genotype-pairs "$MIN_GENOTYPE_PAIRS"
            --exp10 "$EXP10"
            --backward-alignment "$BACKWARD_ALIGNMENT"
            --code-commit "$CODE_COMMIT"
        )
        if [[ "$EXCLUDE_WITHIN" -eq 1 ]]; then
            validation_args+=(--exclude-within)
        fi
        if [[ "$MASK_ENABLED" -eq 1 ]]; then
            validation_args+=(
                --mask-uri "$mask_uri"
                --mask-fingerprint "$mask_fingerprint"
                --local-mask-source "$local_mask_source"
                --local-mask "$local_mask"
                --mask-audit "$mask_audit"
            )
        fi

        if [[ "$FORCE" -eq 0 && "$UPLOAD" -eq 1 ]]; then
            if [[ -s "$completion" ]] || \
                cloud_exists "$remote_directory/$(basename "$completion")"; then
                restore_if_present "$remote_directory/$(basename "$completion")" \
                    "$completion" || true
                restore_if_present "$remote_directory/$(basename "$summary")" \
                    "$summary" || true
                restore_if_present "$remote_directory/$(basename "$run_json")" \
                    "$run_json" || true
                restore_if_present "$remote_directory/$(basename "$pairs_manifest")" \
                    "$pairs_manifest" || true
                restore_if_present "$remote_directory/$(basename "$sample_list")" \
                    "$sample_list" || true
                restore_if_present "$remote_directory/$(basename "$sample_audit")" \
                    "$sample_audit" || true
                restore_if_present "$remote_directory/$(basename "$candidate_regions")" \
                    "$candidate_regions" || true
                restore_if_present \
                    "$remote_directory/$(basename "$candidate_directory")/$(basename "$candidate_manifest")" \
                    "$candidate_manifest" || true
            fi
        fi

        remote_outputs_ok=1
        if [[ "$UPLOAD" -eq 1 ]]; then
            remote_required=(
                "$remote_directory/$(basename "$summary")"
                "$remote_directory/$(basename "$run_json")"
                "$remote_directory/$(basename "$pairs_manifest")"
                "$remote_directory/$(basename "$sample_list")"
                "$remote_directory/$(basename "$sample_audit")"
                "$remote_directory/$(basename "$bitmatrix")"
                "$remote_directory/$(basename "${bitmatrix}.meta")"
                "$remote_directory/$(basename "$candidate_regions")"
                "$remote_directory/$(basename "$candidate_positions")"
                "$remote_directory/$(basename "$candidate_directory")/$(basename "$candidate_manifest")"
            )
            for remote_artifact in "${remote_required[@]}"; do
                if ! cloud_exists "$remote_artifact"; then
                    remote_outputs_ok=0
                    break
                fi
            done
        fi
        if [[ "$FORCE" -eq 0 && "$remote_outputs_ok" -eq 1 ]] && \
            "$AOU" "${validation_args[@]}" --check-only >/dev/null 2>&1; then
            echo "Validated completion; skipping decode."
            if [[ "$UPLOAD" -eq 1 ]] && ! remote_matches_file \
                "$remote_directory/$(basename "$completion")" "$completion"; then
                echo "Completion marker differs; checksum-syncing the chromosome outputs."
                upload_completed_chromosome "$remote_directory" "$summary" \
                    "$run_json" "$pairs_manifest" "$sample_list" "$sample_audit" \
                    "$decode_log" "$bitmatrix" "$candidate_regions" \
                    "$candidate_positions" "$candidate_directory" "$completion"
            fi
            continue
        fi

        # A decode can finish all outputs and then be interrupted during final
        # validation or upload. Rebuild the completion marker from those local
        # artifacts before clearing anything, so a retry does not recompute a
        # valid chromosome.
        if [[ "$FORCE" -eq 0 ]] && \
            "$AOU" "${validation_args[@]}" >/dev/null 2>&1; then
            echo "Recovered complete local outputs; skipping decode."
            if [[ "$UPLOAD" -eq 1 ]]; then
                upload_completed_chromosome "$remote_directory" "$summary" \
                    "$run_json" "$pairs_manifest" "$sample_list" "$sample_audit" \
                    "$decode_log" "$bitmatrix" "$candidate_regions" \
                    "$candidate_positions" "$candidate_directory" \
                    "$completion"
            fi
            continue
        fi

        if [[ "$input_staged" -eq 0 ]]; then
            check_staging_capacity "$input_uri" "$input_fingerprint" "$local_input" \
                "$index_uri" "$index_fingerprint" "$local_index"
            stage_object "$input_uri" "$local_input" "$input_fingerprint"
            stage_object "$index_uri" "$local_index" "$index_fingerprint"
            if [[ "$MASK_ENABLED" -eq 1 ]]; then
                stage_object "$mask_uri" "$local_mask_source" "$mask_fingerprint"
            fi
            mkdir -p "$input_directory"
            bcf_samples_temporary="${bcf_samples}.tmp.$$"
            rm -f -- "$bcf_samples_temporary"
            "$BCFTOOLS" query -l "$local_input" > "$bcf_samples_temporary"
            [[ -s "$bcf_samples_temporary" ]] || die \
                "bcftools found no samples in $local_input"
            mv -f -- "$bcf_samples_temporary" "$bcf_samples"
            mapfile -t indexed_contigs < <("$BCFTOOLS" index -s "$local_input")
            [[ "${#indexed_contigs[@]}" -eq 1 ]] || die \
                "expected one indexed contig in $local_input; found ${#indexed_contigs[@]}"
            IFS=$'\t ' read -r input_contig sequence_length n_records <<< \
                "${indexed_contigs[0]}"
            [[ "$sequence_length" =~ ^[1-9][0-9]*$ ]] || die \
                "BCF index did not report a positive contig length: ${indexed_contigs[0]}"
            normalized_contig="${input_contig#chr}"
            normalized_contig="${normalized_contig#CHR}"
            [[ "$normalized_contig" == "$chromosome" ]] || die \
                "BCF contig $input_contig does not match requested chromosome $chromosome"
            if [[ "$MASK_ENABLED" -eq 1 ]]; then
                "$AOU" workbench-mask \
                    --hardmask "$local_mask_source" \
                    --contig "$input_contig" \
                    --sequence-length "$sequence_length" \
                    --output "$local_mask" \
                    --audit-output "$mask_audit"
                if [[ "$UPLOAD" -eq 1 ]]; then
                    rsync_files "$OUTPUT_PREFIX/shared/masks" \
                        "$local_mask" "$mask_audit"
                fi
            fi
            input_staged=1
        fi

        safe_clear_run "$population_summary_dir" \
            "$summary" "$run_json" "$pairs_manifest" "$sample_list" \
            "$sample_audit" "$completion" "$decode_log" "$bitmatrix" \
            "${bitmatrix}.meta" "$candidate_regions" "$candidate_positions"
        safe_clear_candidate_directory "$candidate_directory"
        rm -f -- "$candidate_raw" "${candidate_raw}.meta"
        "$AOU" workbench-samples \
            --population "$population" \
            --ancestry "$local_ancestry" \
            --qc-exclusions "$local_qc_exclusions" \
            --relatedness-exclusions "$local_relatedness_exclusions" \
            --bcf-samples "$bcf_samples" \
            --output "$sample_list" \
            --audit-output "$sample_audit"

        n_population_samples="$(wc -l < "$sample_list")"
        n_population_haplotypes=$((2 * n_population_samples))
        available_pairs=$((n_population_haplotypes * (n_population_haplotypes - 1) / 2))
        if [[ "$EXCLUDE_WITHIN" -eq 1 ]]; then
            available_pairs=$((available_pairs - n_population_samples))
        fi
        (( N_RANDOM_PAIRS <= available_pairs )) || die \
            "$population has only $available_pairs eligible haplotype pairs; requested $N_RANDOM_PAIRS"

        decode_args=(
            decode
            --executable "$REPO/bin/gamma_smc"
            --input "$local_input"
            --input-format vcf
            --samples "$sample_list"
            --output "$summary"
            --bitmatrix "$bitmatrix"
            --theta "$THETA"
            --rho-over-theta "$RHO_OVER_THETA"
            --mutation-rate "$MUTATION_RATE"
            --generation-time "$GENERATION_TIME"
            --threshold-years "$THRESHOLD_YEARS"
            --recent-call "$RECENT_CALL"
            --no-output-at-hets
            --output-at-stride "$OUTPUT_STRIDE"
            --threads "$THREADS"
            --cache-size "$CACHE_SIZE"
            --pair-block "$PAIR_BLOCK"
            --n-random-pairs "$N_RANDOM_PAIRS"
            --pairs-seed "$PAIRS_SEED"
            --pairs-manifest "$pairs_manifest"
            --exp10 "$EXP10"
            --backward-alignment "$BACKWARD_ALIGNMENT"
        )
        if [[ "$MASK_ENABLED" -eq 1 ]]; then
            decode_args+=(--mask "$local_mask")
        fi
        if [[ "$EXCLUDE_WITHIN" -eq 1 ]]; then
            decode_args+=(--exclude-within)
        fi

        echo "Decoding $population chromosome $chromosome with $THREADS threads..."
        if [[ -x /usr/bin/time ]]; then
            /usr/bin/time -v "$AOU" "${decode_args[@]}" 2>&1 | tee "$decode_log"
        else
            "$AOU" "${decode_args[@]}" 2>&1 | tee "$decode_log"
        fi
        "$AOU" workbench-regions \
            --population "$population" \
            --chromosome "$chromosome" \
            --summary "$summary" \
            --sequence-length "$sequence_length" \
            --output "$candidate_regions" \
            --positions-output "$candidate_positions" \
            --threshold-years "$THRESHOLD_YEARS" \
            --minimum-fraction "$SIGNAL_FRACTION" \
            --merge-gap "$MERGE_GAP" \
            --output-at-stride "$OUTPUT_STRIDE" \
            --profile-half-width "$PROFILE_HALF_WIDTH"

        mkdir -p "$candidate_directory" "$candidate_raw_directory"
        if [[ -s "$candidate_positions" ]]; then
            n_candidate_positions="$(wc -l < "$candidate_positions")"
            check_candidate_capacity "$n_candidate_positions" "$N_RANDOM_PAIRS"
            candidate_decode_args=(
                decode
                --executable "$REPO/bin/gamma_smc"
                --input "$local_input"
                --input-format vcf
                --samples "$sample_list"
                --output "$candidate_summary"
                --raw-output "$candidate_raw"
                --theta "$THETA"
                --rho-over-theta "$RHO_OVER_THETA"
                --mutation-rate "$MUTATION_RATE"
                --generation-time "$GENERATION_TIME"
                --threshold-years "$THRESHOLD_YEARS"
                --recent-call "$RECENT_CALL"
                --no-output-at-hets
                --output-at-stride -1
                --output-positions-file "$candidate_positions"
                --threads "$THREADS"
                --cache-size "$CACHE_SIZE"
                --pair-block "$PAIR_BLOCK"
                --pairs-file "$pairs_manifest"
                --exp10 "$EXP10"
                --backward-alignment "$BACKWARD_ALIGNMENT"
            )
            if [[ "$MASK_ENABLED" -eq 1 ]]; then
                candidate_decode_args+=(--mask "$local_mask")
            fi
            echo "Replaying the exact $N_RANDOM_PAIRS pairs at $n_candidate_positions candidate positions..."
            "$AOU" "${candidate_decode_args[@]}"
            "$AOU" workbench-candidates \
                --population "$population" \
                --chromosome "$chromosome" \
                --summary "$summary" \
                --regions "$candidate_regions" \
                --raw-posteriors "$candidate_raw" \
                --pairs-manifest "$pairs_manifest" \
                --input "$local_input" \
                --sample-list "$sample_list" \
                --bcftools "$BCFTOOLS" \
                --contig "$input_contig" \
                --output-dir "$candidate_directory" \
                --mutation-rate "$MUTATION_RATE" \
                --threshold-years "$THRESHOLD_YEARS" \
                --profile-half-width "$PROFILE_HALF_WIDTH" \
                --variant-half-width "$VARIANT_HALF_WIDTH" \
                --minimum-genotype-pairs "$MIN_GENOTYPE_PAIRS" \
                --minimum-fraction "$SIGNAL_FRACTION" \
                --merge-gap "$MERGE_GAP"
            rm -f -- "$candidate_raw" "${candidate_raw}.meta"
        else
            "$AOU" workbench-candidates-empty \
                --population "$population" \
                --chromosome "$chromosome" \
                --regions "$candidate_regions" \
                --pairs-manifest "$pairs_manifest" \
                --output-dir "$candidate_directory"
        fi

        "$AOU" "${validation_args[@]}"

        if [[ "$UPLOAD" -eq 1 ]]; then
            upload_completed_chromosome "$remote_directory" "$summary" \
                "$run_json" "$pairs_manifest" "$sample_list" "$sample_audit" \
                "$decode_log" "$bitmatrix" "$candidate_regions" \
                "$candidate_positions" "$candidate_directory" \
                "$completion"
        fi
    done

    if [[ "$KEEP_INPUTS" -eq 0 ]]; then
        rm -f -- "$local_input" "$local_index" "$bcf_samples" \
            "${local_input}.cloud-source" "${local_index}.cloud-source"
    fi
done

if [[ "${CHR_SPEC,,}" == "all" ]]; then
    plot_scope="all"
else
    printf -v plot_scope 'chr%s_' "${CHROMOSOMES[@]}"
    plot_scope="${plot_scope%_}"
fi

for population in "${POPULATIONS[@]}"; do
    population_summary_dir="$LOCAL_ROOT/results/$population/chromosomes"
    population_plot_dir="$LOCAL_ROOT/results/$population/plots/$plot_scope"
    mkdir -p "$population_plot_dir"
    plot_args=(
        workbench-plot
        --population "$population"
        --summary-dir "$population_summary_dir"
        --chromosomes "${CHROMOSOMES[@]}"
        --output-dir "$population_plot_dir"
        --threshold-years "$THRESHOLD_YEARS"
        --top-n "$TOP_N"
        --signal-fraction "$SIGNAL_FRACTION"
    )
    if [[ "${CHR_SPEC,,}" == "all" ]]; then
        plot_args+=(--whole-genome)
    fi
    "$AOU" "${plot_args[@]}"

    if [[ "$UPLOAD" -eq 1 ]]; then
        remote_plot_dir="$OUTPUT_PREFIX/$population/plots/$plot_scope"
        rsync_directory "$population_plot_dir" "$remote_plot_dir" \
            --exclude='^plot_manifest\.json$'
        # The manifest commits the complete plot artifact set.
        rsync_single_file "$population_plot_dir/plot_manifest.json" \
            "$remote_plot_dir"
    fi
done

report_dir="$LOCAL_ROOT/results/summary/$plot_scope"
report_args=(
    workbench-report
    --results-root "$LOCAL_ROOT/results"
    --populations "${POPULATIONS[@]}"
    --chromosomes "${CHROMOSOMES[@]}"
    --output-dir "$report_dir"
    --threshold-years "$THRESHOLD_YEARS"
    --signal-fraction "$SIGNAL_FRACTION"
    --merge-gap "$MERGE_GAP"
    --top-n "$TOP_N"
    --gene-annotation "$local_gene_annotation"
    --hit-bin-size "$HIT_BIN_SIZE"
    --gene-context-flank "$GENE_CONTEXT_FLANK"
    --zoom-ymax "$ZOOM_YMAX"
    --hit-label-min-fraction "$HIT_LABEL_MIN_FRACTION"
)
if [[ "${CHR_SPEC,,}" == "all" ]]; then
    report_args+=(--whole-genome)
fi
"$AOU" "${report_args[@]}"

if [[ "$UPLOAD" -eq 1 ]]; then
    remote_report_dir="$OUTPUT_PREFIX/summary/$plot_scope"
    rsync_directory "$report_dir" "$remote_report_dir" \
        --exclude='^run_report_manifest\.json$'
    # The report manifest is synchronized only after all report artifacts.
    rsync_single_file "$report_dir/run_report_manifest.json" \
        "$remote_report_dir"
fi

echo
echo "All requested Gamma-SMC scans, population plots, and the combined report are complete."
