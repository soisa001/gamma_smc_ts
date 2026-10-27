# Adding singleton SNVs to an immutable SV/SNV scaffold

Research date: 2026-10-27. Sources below are primary documentation, source code, or the methods paper. No controlled data were accessed. All scale estimates are planning scenarios, not measurements of the user's files.

**Final implementation decision:** the user selected direct physical phase-block lookup. The delivered workflow has **no SHAPEIT stage**. See [singleton_phasing.md](singleton_phasing.md) for the actual interface, tests, and operating limits. The SHAPEIT details below document the investigated alternative and why it does not retain target physical linkage; they are not claims that a SHAPEIT pipeline was delivered or benchmarked.

## Recommendation

Use direct physical phase transfer. Preserve original scaffold records as authoritative records and add new singleton records around them. Leave conflicting or unanchored physical blocks explicitly unresolved. Do not pass HiPhase genotypes through `phase_rare` expecting it to preserve their physical phase, or deliver SHAPEIT's whole output as a replacement for the scaffold.

## Verified tool constraints

The SHAPEIT documentation recommends `phase_rare` for more than 2,000 samples and describes a phased scaffold with smaller target intervals inside larger scaffold intervals. Its example uses 2.5 Mb targets with 0.5 Mb buffers on both sides. The fixed default seed is 15052011. Thus, 12,000 samples is within the intended large-cohort use case. [phase_rare documentation](https://odelaneau.github.io/shapeit/docs/documentation/phase_rare/)

Source inspection was verified against upstream commit `c34d4db3e99a2f7e23deb727671ae260901a5886` (current main at research time). These are implementation details of that commit, not a guarantee about every installed release:

- The target reader sends every heterozygous target genotype through `pushRareUnphasedHet`, without preserving its supplied phase. The reader explicitly rejects missing scaffold genotypes. Consequently, a physically phased HiPhase target is treated as unphased by this code path. [Target/scaffold reader](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/io/genotype_reader/genotype_reader_reading.cpp)
- The CLI has no `--use-PS` option. `--score-singletons` is described as experimental. The code warns that multithreading prevents reproduction by seed alone. For reproducible stochastic results, use one thread per SHAPEIT process and parallelize independent chunks with deterministic seeds and stable inputs. [CLI implementation](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/phaser/phaser_parameters.cpp)
- Samples must be identical and in identical order: the scanner compares the sample vectors directly. Variants shared between target and scaffold use the scaffold. Empty rare-target or empty scaffold intervals are errors. [Scanner implementation](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/io/genotype_reader/genotype_reader_scanning.cpp)
- The writer creates a new header and new records with GT, optional PP, and AC/AN. Original QUAL, FILTER, and arbitrary INFO/FORMAT annotations are not preserved by that writer. It also writes phased genotypes and computes AN from sample count. Therefore, copying the entire output would violate the user's preservation requirement. [Writer implementation](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/io/haplotype_writer.cpp)

HiPhase encodes physical phase using GT and PS. PS is the one-based position of the first variant in a block; PS values are shared across simultaneously phased small-variant and SV VCFs. Some heterozygotes are intentionally left unphased without PS. Consequently, a `0|1` orientation is meaningful within the sample's physical block, but cannot be assumed to equal the scaffold's chromosome-wide `0|1` orientation. [HiPhase user guide](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md)

SHAPEIT5 uses a distinct singleton model: it assigns the minor allele to the target haplotype with the shorter maximal shared segment. The UK Biobank paper reported singleton switch error rates of 35.1% and 36.6% in WGS validation, versus 50% chance expectation. These are not a promised accuracy for this cohort and should not be confused with accuracy for MAC>1 variants. That study's approximately GBP3,258 genome-wide estimate was for the full common-plus-rare phasing workflow in 147,754 WGS samples on 2022 UKB RAP pricing; it is not an appropriate direct quote for this singleton-only task or current GCP pricing. [SHAPEIT5 methods and results](https://www.nature.com/articles/s41588-023-01415-w)

## Implementation options evaluated (inferred from the above; SHAPEIT options were not adopted)

1. Validate the sample manifest, genome build/contigs, exact allele identities, sorting/indexes, genotype ploidy, and duplicate records. Determine singleton status from genotype counts in the final cohort rather than trusting stale INFO/AC. MAC is `min(AC, AN-AC)` for biallelic records; AC=1 alone does not cover reference-minor singletons. Missing genotypes must remain missing in the delivered records.
2. For each sample/chromosome/HiPhase-PS block, collect exact shared heterozygous SNV anchors. An anchor establishes whether the two physical haplotypes are in the same or reversed orientation as the scaffold. Exact shared SVs can also support alignment if representation and genotype agree, but approximate SV matching is unsuitable for a strict preservation workflow.
3. Transfer singleton phase only when supported anchors agree. Conflicting orientations imply a switch error, representation problem, or inconsistent scaffold block. Do not silently majority-vote those into apparent high-confidence phase. Unanchored and discordant cases should remain explicitly unresolved unless the user opts into statistical estimates.
4. A singleton tied to an independently phased HiPhase block requires a block orientation. Physical transfer must not infer that orientation from a different physical block or sample. If the final scaffold still has independent PS blocks, do not interpret all its haplotype indices as chromosome-wide phase; constrain transfer to a consistent scaffold block or require a global-phasing assertion.
5. Build a temporary SHAPEIT scaffold from complete, phased, biallelic common SNVs (e.g. MAF>=0.001) already in the authoritative callset. At 24,000 called haplotypes, MAF 0.001 corresponds to MAC>=24. Excluding SVs, multiallelic records, and missing data from this temporary statistical scaffold does not remove them from the final output.
6. Supply unresolved singleton genotypes with the exact same sample order. Do not rerun `phase_common` or ligation, which would risk changing established phase. Require a build-matched genetic map, useful anchor density, and disjoint target chunks with scaffold buffers. Skip chunks with no pending targets.
7. Copy only the inferred singleton carrier GT orientation (and separately named provenance/confidence if used) from statistical output into the original singleton record. Do not copy imputed missing genotypes or replace original annotations. Keep physical and statistical phase provenance distinguishable. Preserve original physical PS evidence separately or in a sidecar when changing output PS semantics.
8. Merge the new singleton records with the original scaffold. Verify that removing added singleton records reconstructs every original scaffold record, all sample fields, and all original metadata definitions. Compressed file bytes necessarily differ when adding records; the preservation invariant is original record content and ordering, not identical whole-file bytes.
9. Cache using manifest/input identity, relevant parameters, and verified output integrity, with atomic completion markers. Log versions, commands, seeds, input/output paths, counts, rejected records, physical anchor agreement, unresolved counts, wall time, and peak RSS. Avoid logging individual controlled genotypes into public logs.

For scope control, mixed ploidy/non-PAR X and multiallelic singleton records need explicit policies; an autosomal diploid biallelic SNV implementation should reject or separately report unsupported records instead of silently dropping them.

## Cost model for 12,000 samples

The reported 3 TB HiPhase input and 50 GB scaffold do not establish a runtime. Compressed representation, number of sites, format fields, joint versus per-sample layout, number of physical anchors, unresolved fraction, decompression throughput, and number of full scans matter. A joint cohort VCF can make Python genotype parsing the bottleneck; compressed disk size is not a memory estimate.

Verified US-dollar on-demand rates for Iowa (`us-central1`), retrieved 2026-10-27:

| Resource | Rate |
|---|---:|
| n2-highmem-4, 4 vCPU / 32 GiB | $0.262028/hour |
| n2-standard-32, 32 vCPU / 128 GiB | $1.553888/hour |

These are list VM prices before discounts, platform overhead, and taxes. A 32 GiB VM accommodates a 26 GB application budget plus overhead; it does not prove the pipeline will fit. [Google Compute Engine pricing](https://cloud.google.com/products/compute/pricing/general-purpose?hl=en)

Balanced persistent disk is $0.000136986/GiB-hour (approximately $0.10/GiB-month). Four TiB costs approximately $0.5611/hour or $13.47/day; eight TiB approximately $26.93/day. Provisioned capacity is charged even when unused, until the disk is deleted. Four TiB leaves limited expansion room after staging 3 TB; final retained FORMAT fields and intermediates could require substantially more. [Google disk pricing](https://cloud.google.com/compute/disks-image-pricing)

The following totals are arithmetic scenarios, **not runtime forecasts**. They assume one VM plus a 4 TiB balanced disk retained for the entire stated duration:

| Elapsed duration | 4 vCPU / 32 GiB + 4 TiB disk | 32 vCPU / 128 GiB + 4 TiB disk |
|---|---:|---:|
| 24 hours | $19.75 | $50.76 |
| 72 hours | $59.26 | $152.28 |
| 168 hours | $138.28 | $355.32 |

An 8 TiB disk adds approximately $13.47 per day to either column. The 32-core scenario is an optional scale-up, not the default allocation. Default to four workers and a 26 GB total application budget; locally never exceed the authorized 24 threads. VM memory includes operating-system overhead. Resource caps must apply across concurrent processes, not independently per worker.

Regional Standard Cloud Storage in us-central1 is approximately $0.02/GiB-month. Same-location transfer from a regional bucket to a VM is free of network transfer charges. A multi-region is not considered the same location as a region inside it; other North America transfers are $0.02/GiB, approximately $55.88 for 3 decimal TB. Storage operations, retrieval from non-Standard classes, and platform charges can add costs. Retaining a new 1,000 GiB output costs approximately $20/month. Existing source storage is a sunk ongoing cost unless duplicated. [Google Cloud Storage pricing](https://cloud.google.com/storage/pricing)

For 3 decimal TB, a single complete pass takes approximately 83.3 hours at 10 MB/s, 27.8 hours at 30 MB/s, or 8.3 hours at 100 MB/s. These are simple bytes/throughput calculations, not performance claims. Two or three parsing passes multiply the input-processing component. A sequential scan once or twice is preferable to repeatedly rescanning the same 3 TB for small intervals; indexed chromosome partitions or per-sample staging can improve locality. Avoid materializing full dense genotype matrices in Python.

The conditional 3-7 day scenarios cost roughly $60-$140 on a four-core VM with 4 TiB scratch, or $150-$355 on the 32-core scenario. This is arithmetic conditional on that runtime and disk allocation, not evidence that either configuration will finish within a week. Use a contingency for reruns and larger scratch; storage and I/O can dominate the transfer cost. The implemented workflow's measured synthetic throughput is not a compressed-file throughput estimate; a real pilot is required.

## How to replace the scenarios with a cohort estimate

Run a representative chromosome or several intervals at the full 12,000-sample count, entirely inside the approved controlled environment. Keep both ordinary and high-variant-density regions, and measure extraction/physical transfer and SHAPEIT fallback separately. A tiny synthetic correctness test cannot estimate production runtime or biological accuracy.

Record compressed bytes processed, site counts, singleton counts, unresolved fraction, shared-anchor counts, scan wall time, fallback wall time per chunk, peak RSS per worker, and output bytes. Then estimate:

- extraction hours = total relevant compressed bytes / measured aggregate parsing throughput;
- fallback core-hours = sum over chunk classes of measured core-hours per chunk times number of chunks (or per unresolved variant after checking fixed scan overhead);
- elapsed compute hours = core-hours / effective worker utilization, plus serial merges/validation/staging;
- peak memory = measured per-worker peak times concurrent workers, plus shared-process/system overhead;
- total cost = sum(VM hourly rate times VM hours) + disk GiB times $0.000136986 times retained hours + incremental bucket storage/operations/transfer/platform fees.

Keep a measured versus assumed column in any machine-generated estimate. For a production biological validation, hide phase at a subset of known physical singleton calls and evaluate recovery using held-out read or pedigree evidence, keeping leakage between anchor and test sets under control. Physical anchor agreement is a consistency check, not an independent truth estimate.

## Additional source audit: allele matching and sparse chunks

At the same pinned commit, the shared XCF reader explicitly sets HTSlib `COLLAPSE_NONE`; its `nextRecord` only admits biallelic records and requires scalar INFO/AC and INFO/AN in both VCF inputs. Recompute these tags from the selected sample genotypes in temporary inputs. The scanner reads the union of target and scaffold; a target containing only pending rare variants is valid and need not repeat scaffold SNPs. No duplicate-key or duplicate-position rejection is present in the inspected scanner. [Shared XCF reader](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/xcftools/common/src/utils/xcf.h)

HTSlib defines COLLAPSE_NONE as exact allele matching and documents reordering duplicate-position records while synchronizing files. Therefore two records at the same position with different ALT alleles are not automatically collapsed into one shared record. Exact key matching remains mandatory in our overlay. Reject identical duplicate keys in inputs; excluding all repeated-position SNP loci from the temporary statistical scaffold is an optional conservative restriction, not a SHAPEIT requirement. Never remove such loci from the authoritative final scaffold. [HTSlib synced reader contract](https://github.com/samtools/htslib/blob/1.20/htslib/synced_bcf_reader.h)

**Reference-minor singleton safeguard:** The inspected `phaseCoalescentViterbi` assigns literal ALT allele 1 to the shorter matching segment and does not consult `major_alleles`. Its transpose/merge path copies that result without an allele correction. Thus reference-minor singletons must be excluded from direct statistical fallback. Physical phase transfer remains valid for them. A temporary REF/ALT flip with genotype recoding and inverse mapping could work, but needs a dedicated allele-complement regression test; the original record must never be rewritten by that operation. The same function sets unscored singleton PP to exactly 0.5. Experimental singleton scoring uses `max(w0,w1)/(w0+w1)`; when both lengths are zero the expression has no explicit denominator guard. Do not present either score as calibrated confidence. [Singleton phasing implementation](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/containers/genotype_set/genotype_set_phasing.cpp)

Skip a chunk when there are no pending targets. Do not run SHAPEIT when its buffered interval contains zero eligible scaffold SNPs: the scanner errors. Merely having nonzero scaffold records is insufficient as a quality guarantee. PBWT evaluation requires MAC >= the configured threshold (default 2), and state-selection code assumes scaffold grouping and neighbor arrays exist. Preflight that there are multiple informative SNPs with distinct positions and nonzero genetic-map span; sparse or zero-span chunks should remain unresolved, be expanded/merged, or fail with an explicit diagnostic. The source does not specify an evidence-based minimum SNP count that guarantees acceptable singleton accuracy, so any stronger configurable threshold is our policy rather than a documented guarantee. [Conditioning initialization](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/containers/conditioning_set/conditioning_set_managment.cpp), [condition-state selection](https://github.com/odelaneau/shapeit/blob/c34d4db3e99a2f7e23deb727671ae260901a5886/phase_rare/src/containers/conditioning_set/conditioning_set_selection.cpp)
