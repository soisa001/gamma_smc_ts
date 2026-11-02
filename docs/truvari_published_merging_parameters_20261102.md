# Truvari merging parameters for an association callset

Research date: 2026-11-02. Scope: primary human population studies, with association and QTL applications prioritized. Parameters below describe `collapse` unless explicitly labeled otherwise. This is a literature comparison, not a validated recipe for this project's callset. No controlled data were accessed.

Recommendation: assuming sequence-resolved long-read INS/DELs, test 90% sequence/size similarity at 500 bp as a conservative relaxation, and 80%/80% at 500 bp as the stronger association-callset candidate. Retain 95%/95% as the control. Use 70%/70% and 1000 bp as sensitivity settings rather than immediately adopting both relaxations. Choose the final setting from phenotype-blind validation; the current command, version, and call representation have not been supplied.

## Original Truvari paper and author guidance

English et al. (Genome Biology 2022) distinguish `bench` defaults (70% sequence/size) from `collapse` defaults (95% sequence/size); both use 500 bp reference distance and zero required reciprocal overlap. Their v3.1 intra-individual experiment compared strict merging with a loose `--hap --pctsim 0.70 --pctsize 0.70 --refdist 1000` merge. They selected strict merging for their per-person calls; intersample merging used defaults plus `--chain`. Thus the loose experiment is not their recommended population setting. They explicitly discuss association-power loss from undermerging and allele-frequency distortion from overmerging. [Original paper](https://link.springer.com/article/10.1186/s13059-022-02840-6).

The maintainer later addressed the specific hypothesis that similar SVs may have similar effects: this is plausible but unproven, and the default is intended mainly to resolve technical representations. He retained 95% as his preferred collapse default while acknowledging that future evidence could favor 90%. [Author discussion, July 2024](https://github.com/ACEnglish/truvari/discussions/217).

## Parameters in other primary studies

| Study and status | Version | Merging criteria actually reported | Context and limits |
| --- | --- | --- | --- |
| [Using the linear references from the pangenome to discover missing autism variants](https://www.nature.com/articles/s41467-026-68378-4), Nature Communications 2026 | 4.3.1 | Intersample: `--sizemin 0 --sizemax 1000000 --pctseq 0.90 --pctsize 0.90 --keep common --gt all`; distance omitted. | Merges autism families and control genomes after within-person caller validation. The separate within-person command uses `--intra --gt het -k maxqual --refdist 500`, also at 90% sequence/size. These two operations must not be conflated. |
| [Billingsley et al., diverse brain long reads](https://pmc.ncbi.nlm.nih.gov/articles/PMC11702628/), bioRxiv 2024 preprint | 4.2.1 in author-code installation comment | Initial assembly-call merging: `--pctseq 0.75 --pctsize 0.75 --passonly -f reference`; cross-caller merging adds `--keep maxqual --chain`; cross-cohort merging also sets `--sizemax 60000`. Distance omitted in these commands. | Direct eQTL/mQTL precedent in 351 brains. The final resource also underwent a much more permissive targeted multiallelic pass; it is not exclusively a 75% callset. [Author code, fixed revision](https://github.com/NIH-CARD/CARDlongread_NABEC_HBCC_QTL_manuscript/blob/81c02bfa25e472cd22fed8680bf0afad390abb00/variants/Sv_merging.py). |
| [Fan et al., Han-SV haplotype panel](https://doi.org/10.21203/rs.3.rs-6312842/v1), Research Square 2025 preprint | 4.1.0 | Sequence and size similarity both 0.5; other comparison thresholds default. The prose corresponds to `-p 0.5 -P 0.5`, not a fully printed command. | 943-person panel used for imputation and GWAS. Crucially, merged sites were re-genotyped with LRCaller before statistical phasing. This is an aggressive published precedent, not evidence that 50% is generally optimal. [Methods PDF, pp. 9–10](https://assets-eu.researchsquare.com/files/rs-6312842/v1/ccc4f5b8-f63d-4e0e-8e63-937c50bbcef7.pdf?c=1744263184). |
| [Wang et al., The 1000 Chinese Pangenome](https://pmc.ncbi.nlm.nih.gov/articles/PMC13233627/), Nature 2026 | 4.2.2 explicitly named for benchmarking elsewhere; collapse version not separately specified | `-k common -r 0 -p 0.7 -P 0.7` | Similar alleles within original graph-derived SV sites were collapsed before biallelic decomposition. The study includes eQTL analysis. The zero-distance, same-site design is a major restriction: this is not evidence for applying 70% across an unrestricted 500–1000 bp neighborhood. |
| [Lee et al., ADSP structural variation](https://link.springer.com/article/10.1186/s13059-026-04206-8), Genome Biology 2026 accepted manuscript | 3.1 | Same SV type; >95% sequence and size similarity; within 1000 bp; `--keep maxqual`. | Association study using BioGraph in 16,841 people. Replicate comparisons informed the stringent threshold: >97% of consistent replicate calls had at least 95% sequence/size similarity. Separate 70% replicate comparison criteria are not the intersample merging settings. [Accepted PDF, Methods 5.2 and 5.4](https://link.springer.com/content/pdf/10.1186/s13059-026-04206-8_reference.pdf). |
| [Liao et al., A draft human pangenome reference](https://www.nature.com/articles/s41586-023-05896-x), Nature 2023 | 3.1.0 | `-r 500 -p 0.95 -P 0.95 -s 50 -S 100000` | Collapsed highly similar alleles to compare graph-derived call counts with other callsets. This should not be described as an optimized association-analysis threshold. |
| [K-mer analysis of long-read alignment pileups for structural variant genotyping](https://www.nature.com/articles/s41467-025-58577-w), Nature Communications 2025 | 4.2.2 | `--sizemin 0 --sizemax 1000000 --keep common --gt all`; text specifies same type, 500 bp, and ≥95% size/sequence similarity. | Merged Sniffles discovery calls from 47 HPRC samples, followed by evaluation of genotyping. Supports the discovery-merge-then-genotype architecture, not relaxing thresholds to improve association power. |

## Important implementation detail in the brain QTL precedent

The author script has additional targeted passes after the 75% merge: first `--pctsize 0 --pctseq 0 --passonly --keep maxqual --chain --sizemax 60000`, then `--pctsize 0 --pctseq 0 --typeignore --keep maxqual --chain --sizemax 180000 --refdist 2000`. Those operate on extracted multiallelic loci. They group biologically distinct alleles much more broadly than ordinary redundancy removal. The script mixes notebook syntax and historical filenames, so it is evidence of the chosen settings, not a directly reusable production workflow. [Author code](https://github.com/NIH-CARD/CARDlongread_NABEC_HBCC_QTL_manuscript/blob/81c02bfa25e472cd22fed8680bf0afad390abb00/variants/Sv_merging.py).

The brain paper describes assigning reference genotypes only where assembly coverage supports doing so. The Han-SV paper explicitly re-genotypes the merged set. Both details matter: relaxing clusters alone does not establish accurate per-person allele dosage.

## Interpretation for a new association callset

The literature spans 50–95% similarity, but the studies differ in technology, representation, genotype handling, and scientific objective. There is no universal association-specific Truvari threshold established by these papers. A reasonable analysis design is to compare a modestly relaxed allele-definition candidate against the existing stringent resource, then use a separate locus/functional grouping layer when distinct alleles are hypothesized to share effects. Do not select merging thresholds by whichever produces the smallest association P values.

Before choosing a candidate, inspect pairs that failed the current merge and identify whether sequence similarity, size similarity, distance, type, size limits, or genotype restrictions caused the split. Evaluate representative loci, technical replicate or trio concordance where available, allele-count and genotype changes, and sensitivity across repeat/nonrepeat regions. These are recommendations inferred from the differences between published workflows, not claims that a particular parameter set has already been validated for this dataset.

## Source-code audit

Inspected release v5.5.0, peeled commit `fa97db45fd8a892953f9196408b3966427bf4b1e`, as well as historical v4.1.0. The installed version must be checked before adopting flags.

- `collapse` explicitly defaults to `-p 0.95 -P 0.95 -r 500 -O 0`, `--keep first`, `--gt off`, and size eligibility 50–50,000 bp. Filtered-size records pass through instead of being merged. Ordinary consolidation fills reference/missing genotypes from matching alternate records but does not generally add two heterozygous dosages. `--no-consolidate` disables this operation. The algorithm matches remaining records to a representative; `--chain` expands matching transitively. [Pinned collapse implementation](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py).
- `refdist` tests overlap against an expanded interval; it is not a guarantee that both breakpoint differences are below the supplied number. Sequence thresholds are checked only when both records are sequence-resolved. Therefore, symbolic `<INS>` records without insertion sequence do not become sequence-aware merely by setting `-p`. [Pinned record comparison](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py).
- Size similarity is `min(length1,length2)/max(length1,length2)`. A 100-bp versus 110-bp pair passes 90% but fails 95%; 100 versus 125 passes 80%. Truvari sequence similarity is an edit-distance score normalized by the sum of sequence lengths, with rolling/unrolling comparisons, not ordinary percent identity. [Pinned comparison functions](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/comparisons.py).
- `--gt all` blocks collapsing records carried by the same individual. That can protect distinct alleles but also block genuine duplicate calls; it is not a universal solution to genotype conflicts. `--hap` assumes one person's phased haplotypes, and `--intra` consolidates sample columns for one person. Neither is a general cohort-merging option. [Official collapse documentation](https://github.com/ACEnglish/truvari/wiki/collapse).
- Version changes matter: v5 changed symbolic-variant handling, and v5.5.0 fixed a `--no-consolidate` bug. Reproduce the old baseline before comparing an upgrade plus new thresholds. [Release notes](https://github.com/ACEnglish/truvari/wiki/Updates).

## Proposed pilot and acceptance criteria

These are proposed settings, not claims of an experimentally established optimum:

| Role | `--pctseq` | `--pctsize` | `--refdist` |
| --- | --- | --- | --- |
| Baseline | 0.95 | 0.95 | 500 |
| Conservative relaxation | 0.90 | 0.90 | 500 |
| Stronger candidate for confirmed undermerging | 0.80 | 0.80 | 500 |
| Lower-similarity sensitivity | 0.70 | 0.70 | 500 |
| Distance-only sensitivity | 0.95 | 0.95 | 1000 |
| Combined permissive sensitivity | 0.70 | 0.70 | 1000 |

Keep `--typeignore` and `--chain` off initially. Keep `-O 0` for the initial sequence-resolved INS/DEL comparison; a blanket overlap requirement can block shifted repeat representations. For large deletions, separately inspect interval/gene overlap and consider a DEL-specific overlap gate. Use identical size eligibility and representative selection across pilot arms; rerun the original settings exactly as a separate reproduction baseline. If testing 50–100,000-bp events, explicitly set `-s 50 -S 100000`, and account separately for events outside that range.

Change one threshold at a time for a small diagnostic subset when many candidate pairs still fail: lowering sequence similarity alone cannot rescue a pair that fails the size gate. Increasing distance alone cannot rescue either similarity failure. Revisit original per-sample calls or retained removed records if available: an already collapsed VCF may have lost representative sequences and genotype provenance that a second collapse cannot recover.

For each pilot, retain original IDs, representatives, cluster membership, similarity metrics, and all removed records. Measure within-cluster size and coordinate ranges, remaining near-duplicate pairs, same-person co-carriage, missingness, replicate/trio concordance, and genotype/allele-count changes. Stratify by INS/DEL, size, repeats, ancestry, and sequencing batch. Use ancestry-aware HWE diagnostics rather than pooled HWE as the sole criterion. Increased AF and fewer records alone are not evidence of improvement. A deterministic locus sample (seed 0 if sampling is needed) can support manual inspection without viewing controlled data outside its approved environment.

For association testing, distinguish two products: (1) an allele-resolved catalogue after technical deduplication and (2) explicit locus/functional groups whose members are hypothesized to have similar effects. Preserve membership so that separate-allele and grouped results can be compared. For phased groups, count affected haplotypes without counting two same-haplotype representations twice; for unphased overlapping alleles, do not simply sum and cap dosage. Re-genotype or resolve ambiguous groups. Recompute AC/AN/AF from the final genotypes and handle missingness explicitly. Repeat expansion/contraction lengths may warrant a length-dosage or multiallelic model rather than a single binary event.

No merging or association analysis was executed for this review. The code inspection is static; recommendations require validation against the actual callset.

## Provenance and verification

- Primary article Methods were checked for the brain, Han-SV, 1KCP, and ADSP studies. The ADSP accepted manuscript (publisher PDF) confirms the preprint's merging description.
- Brain author code was retrieved read-only from its public GitHub repository, latest file revision `81c02bfa25e472cd22fed8680bf0afad390abb00` at review time. No source-data locations in the script were accessed or executed.
- HPRC and kanpig parameter strings were confirmed in publisher-indexed primary-methods text; publisher page retrieval itself redirected to an inaccessible authentication endpoint during this session.
- Version omissions are explicit. Defaults are not silently substituted for unreported settings, and benchmarking thresholds are kept separate from population merging thresholds.
