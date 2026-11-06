# Truvari 5.5.0: complex literal alleles and a second merge

Reviewed 2026-11-06. No controlled data were accessed. Source is pinned to the
peeled v5.5.0 commit `fa97db45fd8a892953f9196408b3966427bf4b1e` (the annotated tag
object is `7d66247f7824ce7a1cfe6ec1db7584f38f0c7471`). Line anchors below refer to
the actual source files, not web extractors' reformatted line numbers.

## Literal replacements follow the ordinary INS/DEL path

For resolved, non-BND records without INFO/SVTYPE, `var_type()` labels a longer
ALT as INS, a shorter ALT as DEL, and equal-length multi-base alleles as UNK.
INFO/SVTYPE takes precedence. This does not establish whether the biological
event is a pure insertion/deletion or a replacement. INFO/SVLEN takes precedence
for size; otherwise unequal literal lengths use the absolute net length change.
An equal-length multi-base replacement uses REF length as size, not zero.
[Type and size source](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py#L835-L940).

`get_ref()` and `get_alt()` return the full stored literal alleles; they do not
trim shared prefix/suffix. DEL sequence comparison uses REF; other ordinary
types use ALT. At identical positions the selected sequences are compared
directly. With different boundaries, rolling/unrolling can increase similarity.
The direct score is `1 - editDistance / (len(sequence1) + len(sequence2))`, not
ordinary alignment percent identity.
[Allele access](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py#L396-L410),
[sequence selection](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py#L716-L752),
[score and size ratio](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/comparisons.py#L143-L179).

Consequences inferred from those operations:

- Two shorter ALTs at the same POS and REF have DEL sequence similarity 1,
  even if the ALTs encode different retained/replaced sequence. Same ALT length
  gives size similarity 1 as well. Raising both thresholds to 0.9 does not
  distinguish that case.
- For INS, ALT differences are examined, but shared retained sequence inside
  a long padded ALT contributes to the score. A high score can reflect shared
  context rather than similarity restricted to inserted sequence.
- As a constructed illustration, two 1,050 bp ALTs sharing 1,000 bases and
  differing at all 50 other positions have edit distance at most 50 and direct
  sequence score at least `1 - 50/2100 = 0.9762`. If both replace the same
  1,000 bp REF, their net sizes are both 50 and their size score is 1.
- Net length alone can understate a replacement's altered sequence. REF80/ALT70
  is net10 even if many bases change, whereas REF80/ALT20 is net60. Those two
  have size similarity 1/6 and fail 0.7 despite identical REF sequence scores.
- A same-POS restriction or positive reciprocal-overlap threshold does not
  distinguish different ALTs sharing the same POS and REF interval.

Truvari's `cpx_match()` terminology refers to matching a BND to a non-BND via
breakend decomposition. Two non-BND literal replacements go directly to
`var_match()`; enabling that decomposition is not an ALT-to-REF edit-script
analysis of each literal replacement.
[Dispatch](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py#L580-L597),
[BND decomposition](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/variant_record.py#L277-L370).

## What this repository sends to Truvari

`scripts/run_sv_phased_batch.py` currently streams chromosome extraction,
`bcftools norm -m -any --multi-overlaps 0`, then compression. It does not request
`norm -a` atomization or provide a reference with `norm -f` for left alignment.
Multiallelic splitting is not decomposition of every replacement into elementary
edits. The native merge uses `-s 0`, `-S -1`, `--keep common`, `--gt all`, `-r 500`,
`-O 0`, and no `--chain`; sequence and size thresholds are both 0.7 or both 0.9.
The final literal INS/DEL subset uses absolute net length at least 50 bp.

The comparison script's simple/complex classification is an audit category,
not a switch in Truvari. A deletion is simple for the coverage union only when
trimming shared prefix/suffix leaves no residual ALT. The earlier pure-indel
restriction therefore selected a different, much smaller record population.
The supplied chr8 aggregates indicate that 95.5% of original DEL records were
complex by this definition, and 98.9% of removed DEL/representative pairs shared
POS and REF. These observations explain why adding the complex population can
dramatically change counts, but do not by themselves establish the biological
correctness of individual merges.

## Merging an already merged callset

`collapse_chunk()` sorts candidate records, takes one representative, and
accepts matches to that representative. Without `--chain`, it does not expand
through matches to absorbed members. Members are not required to match one
another. The genotype mask is updated on each accepted member: `--gt all`
rejects a candidate if any sample carries ALT in both that candidate and the
current cluster. This is a sample-level compatibility condition, not a test of
equivalent molecular effects or phase-block agreement.
[Matching and clustering](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py#L122-L228),
[carrier mask](https://github.com/ACEnglish/truvari/blob/fa97db45fd8a892953f9196408b3966427bf4b1e/truvari/collapse.py#L63-L86).

A second pass over retained representatives cannot recover original alleles
removed in an earlier pass. It is generally not equivalent to applying the new
threshold once to all original alleles. The direction is not universally
overmerging or undermerging:

- **Remaining distinctions:** a member of an old cluster might match another
  cluster although the retained representatives do not. Representative-only
  matching loses that possible connection.
- **Broader effective grouping:** two representatives may match under relaxed
  thresholds while some hidden members of their original clusters would fail
  against the new representative. A match between representatives is not a
  guarantee about every original member.
- **Genotype dependence:** if prior aggregation correctly unions carriers,
  `--gt all` can prevent further merges that overlapping sample sets would
  otherwise allow. If prior genotype consolidation lost carriers, that protection
  may be weaker. The actual upstream provenance and GT behavior must be checked
  before attributing this to this panel.
- **No automatic non-idempotence claim:** a repeated run with identical
  parameters, ordering, alleles, and effective carrier constraints can be stable.
  Merely being a second pass is not evidence of a bug or additional collapse.

Keep the requested 0.7 and 0.9 runs independent from the same prepared panel,
not p90 applied to the p70 output. Compare simple and complex categories
separately, and distinguish same-POS/REF matches from shifted matches. These are
sensitivity analyses of a representative panel, not recovery of original
per-sample allele resolution or truth-set precision/recall benchmarks.

## Executed synthetic checks

`scripts/probe_truvari_complex_alleles.py` runs the installed, pinned PyPI
Truvari 5.5.0 through both its record API and 14 real `collapse` commands.
Each pair has disjoint sample carriers and phased GTs, with `--keep common
--gt all -r 500 -O 0 -s 0 -S -1`. No controlled sequences are used. Results
and installed-source hashes are in `truvari_complex_synthetic_results.json`.

| Constructed pair | Native sequence score | Size score | p70 / p90 retained records (input = 2) |
|---|---:|---:|---|
| Same REF80; different ALT20 replacement sequences | 1.0000 | 1.0000 | 1 / 1 |
| Same REF80; ALT20 versus ALT10 replacements | 1.0000 | 0.8571 | 1 / 2 |
| Simple INS: anchor plus 100 C versus anchor plus 100 G | 0.5050 | 1.0000 | 2 / 2 |
| Complex INS: REF = 1000 A; ALT = 999 A + 101 C or 101 G | 0.9541 | 1.0000 | 1 / 1 |
| Same complex INS after removing common 999-base prefix and shifting POS | 0.5000 | 1.0000 | 2 / 2 |
| Equal-length REF80/ALT80 replacements without SVTYPE | 0.5000 | 1.0000 | 2 / 2; native type UNK |

The padded/trimmed comparison preserves the represented haplotypes while
changing the native match outcome. This establishes representation sensitivity
in these synthetic cases, not the fraction of real panel merges that are wrong.
The complex DEL pair establishes that a score of 1 is not proof of identical
alternate sequence. Increasing to 0.9 therefore provides useful sensitivity
information but cannot by itself validate complex-allele equivalence.

Reproduce without cloud access:

```bash
uv run --no-project --python 3.12 scripts/probe_truvari_complex_alleles.py \
  --output docs/truvari_complex_synthetic_results.json
```
