# chr1 Truvari 70% merge pilot

This workflow stages the full AoU chr1 panel, extracts structural alleles, runs Truvari 5.5.0 with **both sequence and size similarity at 0.7**, and compares the existing panel alleles with the new representatives. It writes data, summaries, logs, and figures only to the VM. There is no upload or Git push in the workflow.

The user confirmed this input on 2026-11-02, with outputs retained on the VM:

```text
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf
gs://rw-long-reads-transfer-2026-06-17/v9/lrWGS/panel/panel/panel_bubble_split_vcf/aou_lr_phase2_v1.chr1.bubble.split.bcf.csi
```

The source is recorded in this repository's `scripts/run_aou_workbench.sh`, `scripts/run_singleton_chr.py`, and `WORKBENCH_RUNNER.md`. The `hmmix_aou_postdecode` repository uses the same full-panel path in `src/aou_hmmix/postdecode.py` and its SV analysis code. Its README also mentions a different `panel_bubble_split_leaveout_vcf` callset; that is not this pilot's input. No HiPhase BCF, alignment files, reference FASTA, or ancestry masks are required for this pilot.

## Run in the authorized Workbench VM

Run this command cell there, not on the development machine. The workflow uses existing gcloud authentication and configuration without requiring `GOOGLE_PROJECT` or adding a billing-project flag. The default output directory is outside the Git checkout. The first invocation installs public tools and runs synthetic tests before staging controlled data.

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
  git -C "$REPO" switch AOU_run_opt
  git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
cd "$REPO"
export SV_MERGE_ROOT="$HOME/sv_merge_chr1_p70"
export SV_THREADS=4
export SV_MEMORY_GB=26
bash scripts/run_chr1_sv_merge_workbench.sh all
```

For separate phases, replace `all` with `plan`, `stage`, `prepare`, `merge`, or `compare`. `plan` does not access the bucket or read controlled data, though the shell launcher bootstraps public tooling. The Python `--phase plan` does neither. Run the phases in that order. `prepare`, `merge`, and `compare` use local files and do not contact GCS. Repeating a phase reuses checksum-verified artifacts when relevant inputs and parameters match. A new repository commit alone does not invalidate the cache.

`stage` uses `gcloud storage rsync SOURCE_DIRECTORY LOCAL_DIRECTORY`, with an anchored exclusion expression that selects only this BCF and its `.csi`. Neither rsync nor the metadata queries receive a billing-project flag by default. The Python entry point retains `--billing-project PROJECT` as an optional explicit override; it does not infer this override from `GOOGLE_PROJECT`. It never uses deletion of unmatched files. It verifies remote size, generation, CRC32C, and local checksums; it checks the source generation again after transfer. If an index is genuinely absent (404), it builds one locally; authorization failures are errors. Staging requires enough disk for a replacement copy plus a configurable 64 GiB scratch floor. This is a floor, not a reliable prediction of cohort-size-dependent scratch needs.

The standalone Python entry point accepts `--input-bcf /absolute/path/to/already-local.bcf` to use an indexed local BCF or VCF without any cloud access. It also accepts `--bcftools` and `--truvari` paths. `uv run --no-project --python 3.12 --with truvari==5.5.0 scripts/run_chr1_sv_merge.py --help` lists all options. Use a separate output root for a different threshold so the two results remain available.

## Allele definition and merge settings

The original panel contains small variants as well as SVs. The baseline here is **the panel's structural ALT alleles after multiallelic splitting**, not all panel rows and not a fresh 95% Truvari run. `prepare/summary.json` separately reports candidate source rows, source rows containing an SV, and SV allele records. `prepare/source_alleles.tsv.gz` maps unique pilot IDs back to source candidate row, original ID, position, and REF/ALT sequence hashes. Duplicate original IDs and same-position distinct alleles are retained.

`bcftools norm -m -any --multi-overlaps 0` decomposes multiallelic records. A target allele is coded as 1 and another ALT as 0; this is a per-allele non-target indicator, not proof of a true reference haplotype. For example, `1|2` becomes `1|0` and `0|1`. No reference-based left alignment or haplotype reconstruction is performed. Both catalogues use the same decomposition. Literal alleles qualify by net length change of at least 50 bp; equal-length substitutions and complex replacements with net change below 50 bp are outside this definition.

The initial eligibility range is **50–100,000 bp, inclusive**, explicitly extending Truvari's 50,000 bp maximum default. Larger variants remain unchanged and are counted. Symbolic, ambiguous-sequence, breakend, and complex replacement alleles also remain in the catalogues with explicit pass-through reasons. Symbolic variants of unknown size are retained and reported as unknown; known symbolic variants below 50 bp are excluded. Only pure sequence-resolved A/C/G/T INS and DEL alleles are merged. This prevents an unresolved `<INS>` from silently bypassing the sequence-similarity gate. A future symbolic-resolution analysis would require the appropriate reference and insertion sequences.

The two independent type-specific commands are equivalent to:

```bash
truvari collapse -i INS.vcf.gz -o INS.kept.vcf -c INS.removed.vcf \
  -p 0.7 -P 0.7 -r 500 -O 0 -s 50 -S 100000 \
  -k first --gt off --median-info
```

The DEL command has the same settings. Chaining and type mixing are off. The ordinary genotype consolidation behavior is enabled. No PASS-only filter is applied, so the original and new comparisons cover the same input FILTER categories. Unique source IDs and old Truvari cluster annotations are replaced for this pilot; source mapping is preserved. SVTYPE and SVLEN are made allele-specific. The defaults and full tool commands are recorded in local logs and manifests.

INS and DEL collapses can run concurrently because cross-type matching is disabled. Within a type, the entire chromosome is processed together, avoiding artificial genomic chunk boundaries. Truvari's collapse core is serial; four threads are used where supported for I/O, with up to two independent collapse processes. The default 26 GiB is a resource budget, not a hard OS memory limit; each concurrent `bcftools sort` uses a quarter of that budget. Request that memory from Slurm if using a scheduler. The existing seed 1729 is recorded; the analysis does not perform random sampling.

## Outputs and interpretation

All paths below are relative to `$SV_MERGE_ROOT`:

| Output | Contents |
| --- | --- |
| `prepare/original.vcf.gz` + `.tbi` | Original decomposed SV allele catalogue with unique pilot IDs and recomputed AC/AN/AF. |
| `prepare/source_alleles.tsv.gz`, `samples.txt`, `summary.json` | Mapping to original records, sample manifest, extraction/eligibility counts. |
| `merge/new.vcf.gz` + `.tbi` | Merged INS/DEL representatives plus unchanged pass-through SV alleles; coordinate sorted. |
| `merge/{INS,DEL}.removed.vcf.gz` + `.tbi` | All removed records and Truvari similarity annotations. |
| `compare/counts.tsv` | Original/new/removed counts, percent reduction, and eligible original count for ALL, INS, DEL, OTHER. |
| `compare/statistics.tsv` | Number of observations, quantiles, mean, zero counts, and distance counts at ≤50/500/1000 bp; size, AF, and missingness distributions. |
| `compare/size_counts.tsv` | Counts in 50–99, 100–499, 500–999, 1k–9.9k, 10k–99.9k, ≥100k, unknown bins. |
| `compare/density.tsv` | Occupied 1 Mb bins by callset/type; absent bins have zero records. Final bin endpoints are not clipped to contig length. |
| `compare/{original,new}.records.tsv.gz` | Per-allele positions, changed-reference intervals, lengths, carrier counts, AC/AN/AF, and missingness. |
| `compare/membership.tsv.gz` | Every original allele mapped exactly once to its representative; removed allele sequence/size similarities. |
| `compare/clusters.tsv.gz` | Cluster size, coordinate/size spans, same-sample co-carriage, carrier retention, representative AC and GT changes. |
| `compare/summary.json` | Key count/conservation and genotype audit results. |
| `compare/counts.{png,pdf}`, `spacing_{ALL,INS,DEL}.{png,pdf}` | Large-text letter-landscape figures; PDFs contain vector artwork and editable text. No inline display. |
| `logs/`, phase `done.json` files | Parameters, executed commands, tool output, checksums, sample-manifest identity, and elapsed times. |

**Inter-SV distance is defined explicitly.** Adjacent distance is `POS[i+1] - POS[i]` after coordinate sorting, including zero for same-position alleles. Nearest-neighbor POS distance uses the smaller of the left/right gaps, or the only available gap at a chromosome end. A group with fewer than two alleles has no distance observations. ALL mixes types; INS and DEL calculate within their respective type. These are descriptive coordinate spacings, not Truvari matching distances or estimates of unmerged duplicates.

The additional interval gap measures the distance from a changed-reference interval's start to the rightmost endpoint of the union of preceding intervals, clipped at zero. It handles nested deletions; insertions are points. Literal intervals remove shared REF/ALT prefix and suffix; symbolic spans use the available POS/END representation. The strict-overlap count uses start < preceding endpoint; adjacent touching intervals have zero gap but are not counted as overlapping. This interval metric is distinct from POS spacing.

The comparison verifies **original = retained + removed**, unique membership, same-type clusters, and agreement with Truvari's NumCollapsed. It also reports samples carrying multiple member alleles, any lost/gained carriers relative to the member union, and changes in each representative's GT and AC. AC/AN/AF are recalculated from final GT, including called alleles in partially missing GTs; missing-sample fraction counts any missing allele. The other inherited INFO/FORMAT annotations remain representative-specific and may no longer describe the whole cluster.

**This produces an exploratory merged callset, not validated association dosages.** Truvari consolidation does not in general add two member dosages or resolve different alleles on opposite haplotypes. Co-carriage counts highlight that issue but do not resolve it. Use the membership table for further haplotype-aware grouping or re-genotyping before association analysis. The 70% pilot and a reduction in counts alone do not establish improved biological accuracy. Merging an already merged panel cannot restore discarded original sequences.

## Validation and recovery

### Separate INS/DEL counts, raw lengths, and unique affected bases

When the four `merge/{INS,DEL}.{kept,removed}.vcf.gz` files and their indexes are complete, the following comparison can run even while the original workflow is still building the combined `merge/new.vcf.gz`. It writes into a separate `compare_ins_del/` directory. It does not stage inputs, rerun preparation or Truvari, read the combined VCF, or parse cohort genotypes. Do not rerun `all` to obtain these summaries.

```bash
%%bash
set -euo pipefail
REPO="$HOME/gamma_smc_ts"
if [[ ! -d "$REPO/.git" ]]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$REPO"
else
  git -C "$REPO" switch AOU_run_opt
  git -C "$REPO" pull --ff-only origin AOU_run_opt
fi
cd "$REPO"
"$REPO/.tools/uv-bin/uv" run --no-project --python 3.12 \
  scripts/compare_chr1_sv_alleles.py \
  --root "$HOME/sv_merge_chr1_p70" --threads 4
```

The existing run has already installed the private uv executable used by this cell. This standalone script needs Matplotlib, but not Truvari, bcftools, or gcloud. It verifies the prepared allele metadata and sample manifest against their saved checksums, reads each kept/removed VCF to EOF, checks REF/ALT hashes and raw lengths against the original manifest, and verifies disjoint retained/removed membership and cluster counts. The comparison fails if a source file changes while it is being read. This validates these separate outputs without requiring the overall merge completion marker; it does not certify the still-running combined-file phase. Existing cached preparation/merging artifacts remain unchanged.

Two scopes are reported separately, each split into DEL and INS:

- `merge_eligible`: alleles that entered Truvari, compared with their retained representatives.
- `all_ins_del`: all original INS/DEL alleles, with the same removed IDs excluded after merging. This includes unchanged pass-through calls reconstructed from `prepare/source_alleles.tsv.gz`, without rereading the large pass-through genotype VCF.

The source manifest's schema-1 literal `size` was computed directly from `abs(len(ALT)-len(REF))`; the script uses that for unchanged literal alleles and independently verifies it against actual REF/ALT strings for the merged subset. No INFO/SVLEN, allele frequency, genotype count, or median cluster size substitutes for raw length.

| Metric | Definition |
| --- | --- |
| `sv_records` | Number of biallelic SV records, including exact duplicate rows. |
| `unique_alleles` | Distinct `(CHROM, POS, REF, ALT)` representations, using saved sequence SHA-256 hashes. Same sequence at another position is a different allele. This is exact deduplication, not another similarity merge. |
| `raw_bp_sum_records` | Sum of `len(ALT)-len(REF)` for INS or `len(REF)-len(ALT)` for DEL across literal records. |
| `raw_bp_sum_unique_alleles` | The same sum, counting each distinct literal allele once. |
| `unique_inserted_bp` | INS raw-length sum once per distinct literal insertion allele. Different sequences at the same position each count. This is not a union of reference coordinates. |
| `unique_deleted_reference_bp` | Union length of simple literal deletion intervals, counting each deleted reference base once across overlaps, nesting, and duplicate records. Shared VCF anchors/padding are excluded. |

Symbolic calls remain in record/allele counts but are explicitly excluded from raw-base totals: the length of the text `<DEL>` or `<INS>` is not a biological allele length. Complex DEL replacements contribute their raw net lengths to sums, but are excluded from the deletion union because net length alone does not identify the deleted reference bases. Their exclusion counts are included in the table. Literal sequences containing N still have raw lengths; they contribute to deletion coverage only when their trimmed reference span equals their net deletion length.

Outputs in `compare_ins_del/`:

- `before_after.tsv`: before/after counts and all base totals for both scopes and types.
- `changes.tsv`: absolute and percentage reductions for each metric.
- `inter_sv_distance_summary.tsv`: adjacent POS distance summaries, separately by type, scope, and record/unique-allele basis. Zero distances are retained and no distance is calculated across contigs.
- `{scope}.{before,after}.deleted_union.bed`: the 0-based half-open deletion unions underlying the coverage totals.
- `{scope}.sv_records.{png,pdf}` and `{scope}.unique_affected_bp.{png,pdf}`: separate DEL/INS panels; deleted-reference coverage and inserted-allele bp have distinct axis labels.
- `{scope}.{DEL,INS}.inter_sv_distance.{png,pdf}`: before/after ECDFs of adjacent POS distances, using the record basis, a symlog distance axis, and zero-distance observations. Each figure is letter-landscape with large text and vector PDF artwork/editable text.
- `summary.json`, `comparison.log`, `done.json`: definitions, input paths/checksums, conservation results, runtime provenance and cache completion. Nothing is displayed inline or uploaded.

Six focused synthetic tests cover overlapping/nested deletions, exact duplicate records, distinct insertion sequences at one position, repeated sequences at different positions, raw lengths despite conflicting INFO/SVLEN, complex/symbolic exclusions, zero distances, corrupted gzip, invalid membership, and full comparison/cache reuse without a combined VCF. Run with `uv run --no-project --python 3.12 --with matplotlib==3.10.8 python -m unittest discover -s tests -p test_compare_chr1_sv_alleles.py -v`.

### Original full-workflow checks

Portable tests cover allele boundaries, symbolic exclusion from merging, partial and haploid GT counts, tied positions, nested deletions, namespaced cluster IDs, conservation failures, exact rsync selection, content corruption, comparison outputs, and cache reuse. A real-tool synthetic integration test covers multiallelic decomposition, eligible and pass-through calls, INS/DEL collapse, genotype consolidation, and idempotence. The shell launcher runs it before `merge` or `all`; a failure stops the run.

Development validation on 2026-11-02: ten portable tests passed, shell syntax and Python compilation passed, and synthetic PNG figures were visually inspected. The real-tool integration test could not run locally because both WSL distributions failed to start with socket error `0x80072747`; it is therefore an explicit preflight requirement in the VM launcher. No controlled input was downloaded or analyzed on the development machine.

To run just these tests in a configured Linux environment:

```bash
SV_PILOT_INTEGRATION=1 uv run --no-project --python 3.12 \
  --with truvari==5.5.0 --with matplotlib==3.10.8 --with google-crc32c==1.8.0 \
  python -m unittest discover -s tests -p test_chr1_sv_merge.py -v
```

Phase completion manifests are written only after successful work and checksum verification. If a phase fails, correct the logged issue and rerun that phase, followed by dependent phases. Existing partial outputs are not treated as completed. Full artifact SHA-256 checks detect changes even when file sizes and modification times are preserved. Cache keys use relevant inputs, sample manifests, tool versions and scientific parameters, not the repository commit. A per-root process lock prevents concurrent writers; after a forced termination, verify that the process has stopped before removing a stale `.running.lock`.

Implementation references: [Truvari 5.5.0 collapse](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py), [record comparison and sequence gates](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py), and the repository's [published parameter review](truvari_published_merging_parameters_20261102.md).
