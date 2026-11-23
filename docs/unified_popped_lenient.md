# Unified popped and lenient panel

Run `bash scripts/run_unified_sv_panel_workbench.sh` from a Linux Workbench checkout.
Default inputs are the completed chr1–21 popped merge and genotype repair under HOME
and `~/ashg_revised/v3_lenient_recollapse_same_start_90_70`. Default output is
`~/unified_popped_lenient_p70`. The launcher can take chromosome names for a trial.
Use a separate UNIFIED_ROOT for a subset so the aggregate tables always describe
the chromosomes requested in that invocation. JOBS defaults to four.

The repaired popped VCF supplies all small variants, and its retained SV
representatives through exactly 10,000 bp. For literal alleles the cutoff is
absolute net REF/ALT length. This includes long balanced replacements from
popped when their net change is small. An original member can cross the boundary
when its representative stays on the other side. Cluster membership, representative
alleles and phased GTs come directly from the saved run. Size-unknown popped SVs
are retained. Popped records above the cutoff have a null panel key in the map.
The independent ultralong lenient run supplies its p70 retained records including
symbolic events, using the audited source length (>10,000 bp). This is concatenation
of the two tiers in coordinate order without a new biological similarity match
between tiers. Two nearby records from different tiers can describe related events.

The builder checks input receipts, sample identity and any defined target contig
lengths. Long sample columns are reordered if necessary. The same reference build is
assumed when contig lengths are absent. SV genotypes, phase and missingness are
carried through. Incompatible non-reserved INFO/FORMAT declarations in long records
get a UL_ prefix. SVs carry INFO/UP_KEY for joining to the mapping. Small non-SV
rows are copied directly. Output and interim data stay on the VM. A fresh build
needs disk space for the full output VCFs. A completed chromosome is reusable via
its checksum manifest. The run log and record accounting show input and output counts.

Each chromosome directory contains an indexed VCF and `sv_mapping.jsonl.gz`.
Group mapping rows by non-null panel_key to recover the original members. Popped
simple clusters carry the original merge's exact old and new representation fields.
Other popped SVs map individually. For lenient records, old.audit_ID is the unique
RC identifier from the long run (the original BCF row ordinal); old.original_ID
can be nonunique. The long crosswalk starts at the supplied v3 callset. Original
discovery calls predating v3 are outside that crosswalk.

`sv_summary.tsv` reports each chromosome and ALL, by source and SV class.
SV statistics include literal alleles with at least 50 residual REF or ALT bases
after trimming identical flanks, and symbolic SVs. Simple INS/DEL have one empty
residual side. Complex replacements have bases on both sides. Ambiguous sequence
and symbolic types have separate rows. The INS/DEL literal_and_symbolic rows
pool unambiguous simple literals and the matching symbolic type. Use the combined/ALL
row for distinct SV positions across all classes. Subtotals overlap.

- positions counts distinct output CHROM/POS, allele_records counts VCF rows, and
  unique_allele_keys counts distinct trimmed breakpoint/REF/ALT identities. Symbolic
  identities use reported position, REF, ALT, END and size; identical symbolic
  representations do not establish identical biological sequence.
- unique_reference_span_bp is the reference interval union within the row: deleted
  bases for simple DEL, replaced reference span for complex alleles, or a verified
  END span for symbolic DEL/DUP/INV. Overlaps count once. A complex replacement can
  contain internal matches, so its span is not an edit distance.
- unique_alt_payload_bp sums residual ALT length once per distinct breakpoint and
  ALT sequence. For simple INS this is site-specific inserted sequence. Symbolic
  insertion payloads are unknown; symbolic_reported_length_sum is a separate
  sum of reported lengths, without sequence deduplication.
- allele_frequency_spectrum.tsv gives counts of SV records in called_AF bins and
  observed_AC bins. AF = observed ALT copies / called alleles, including called
  halves of partially missing GTs. AN=0 has its own bin. called_allele_fraction
  is summed AN divided by possible diploid autosomal copies. Discovery-only long
  calls can have low call rate, making their called AF unsuitable as a population
  prevalence estimate. Compare the popped and lenient rows separately.
