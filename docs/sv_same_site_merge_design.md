# Same-site alternate-allele consolidation: feasibility and proposed design

Research date: 2026-11-06. This is a design note, not a production implementation or a claim that the controlled panel has been reprocessed.

## Two different operations

`bcftools norm -m +any` joins compatible biallelic records into multiallelic records. It does not cluster different ALT strings by a 0.7 or 0.9 similarity threshold. A site with 100 distinct ALTs may become one VCF record while still containing 100 alternate alleles. Count ALT alleles, as well as records, when evaluating consolidation. Bcftools' `merge` command instead primarily combines files with different sample sets. [Official bcftools documentation](https://samtools.github.io/bcftools/bcftools.html#norm)

Similarity consolidation at exactly the same `(CHROM, POS, REF)` is feasible and avoids cross-position haplotype reconstruction. POS alone is insufficient: different REF strings can describe different reference intervals. The proposed method below is custom allele grouping, not a native Truvari command.

## Proposed matching contract

1. Group literal alleles by exact CHROM/POS/REF, retaining the original allele index and provenance. Keep symbolic, spanning-deletion `*`, ambiguous-sequence and breakend alleles outside this sequence-comparison method unless their semantics are explicitly supported.
2. Compare ALT sequences directly for both net INS and net DEL. Truvari 5.5.0 instead chooses REF for a DEL's sequence comparison; a shared REF therefore fails to discriminate different replacement ALTs at that site. [Pinned Truvari source](https://github.com/ACEnglish/truvari/blob/v5.5.0/truvari/variant_record.py#L651-L684)
3. Define the metric explicitly. A candidate is normalized edit similarity `1 - edit_distance(ALT1, ALT2)/max(len(ALT1), len(ALT2))`. This is **not** the native Truvari score, which uses the sum of sequence lengths as denominator. Consequently 0.7/0.9 would be new settings, not directly equivalent to the prior Truvari runs. Use a documented convention for two empty comparison strings. [Native similarity implementation](https://github.com/ACEnglish/truvari/blob/v5.5.0/truvari/comparisons.py#L133-L148)
4. Keep net INS and DEL separate if retaining the previous type policy. A separate net-size similarity requirement can prevent long shared sequence from masking a large difference in net length. Do not prefilter at 50 bp; apply that requirement to final representatives as requested. Equal-length replacements need their own explicit policy.
5. Select representatives by decreasing called allele count with deterministic ties. Compare each candidate to the chosen representative; avoid connected-component chaining that can join mutually dissimilar endpoints. Representative-based matching does not guarantee every pair of members exceeds the threshold; complete-link clustering would be required for that stronger contract.

These are proposed analysis choices, not claims that sequence similarity establishes equivalent association effects.

## Shared padding

For a common coordinate frame, remove only a common prefix/suffix shared by REF and **every** ALT in the group when constructing comparison strings. Do not independently trim each REF/ALT pair and then compare their residual ALTs as if they covered the same interval. Different per-allele trimming boundaries describe different intervals and can erase the positional distinctions being tested.

Keep original padded representative strings in the emitted VCF, or normalize the whole output record consistently with valid VCF anchoring. Removing common flanks is optional for record representation but material to normalized sequence scores: identical padding increases similarity without changing the edit distance. Group-wide trimming will not remove every shared internal base, and one divergent ALT may prevent trimming a long flank; inspect score sensitivity rather than treating the metric as representation independent.

## Phased genotypes and carrier policy

Prefer original multiallelic records when available. VCF GT values index REF/ALT, and `|` plus PS encodes phase relationships. A many-to-one allele-index map can preserve each haplotype's position: if alleles 1 and 2 map to cluster 1, `2|0` becomes `1|0`. Preserve missing values, ploidy, separators and phase-set membership. [VCF specification, GT and PS](https://samtools.github.io/hts-specs/VCFv4.5.pdf)

Under unrestricted cluster recoding, `1|2` becomes `1|1`: both chromosomes carry the cluster, but their original sequence difference is lost. This is valid **cluster dosage**, not evidence of sequence homozygosity. It also changes the prior `--gt all` policy. To retain that policy, reject a merge whenever any sample carries both candidate allele groups; update the carrier union after every accepted member. Do not silently allow `1|2` to collapse.

The current bucket filename includes `bubble.split`; do not assume it still contains original multiallelic GTs. If only split records exist, reconstruct a per-haplotype assignment and validate phase sets, missingness and conflicting calls. Joining records is not proof that the original GT was recovered.

## Annotations and outputs

VCF Number=A, R and G fields depend on alternate alleles, all alleles, and genotypes respectively. After changing ALT indices, remap or recompute them with field-specific semantics, or remove them. Recompute AC/AN/AF from recoded GT; do not blindly copy PL/GL, AD, GQ or allele-specific annotations from the representative. [VCF specification, header Number definitions and genotype fields](https://samtools.github.io/hts-specs/VCFv4.5.pdf)

Retain an allele-to-cluster map and the original input, because approximation loses allele identity. Report original versus clustered ALT counts, positions, complex/simple categories, cluster sizes, carrier-conflict rejections, genotype missingness and final >=50 bp counts. Validate every haplotype's mapped index, unchanged phase sets and deterministic output. The same-site boundary removes distance and reciprocal overlap as tuning parameters, but it does not make the custom similarity and genotype contracts optional.
