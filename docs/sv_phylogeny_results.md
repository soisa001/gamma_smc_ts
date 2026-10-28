# Phylogenetic SV mutation-age benchmark with a chimp outgroup

All three proposed methods were run: IQ-TREE 2.4.0 followed by LSD2 dating,
BEAST 2.7.8, and Relate. Relate recovered a compatible SV-carrier clade at all
12 pilot loci. Direct local phylogenetic inference was less consistent,
particularly for recent mutations. None of these runs establishes a calibrated
mutation-age estimator. Chimp rooting solves an orientation problem; it does
not determine where along a long human ancestral branch the SV mutation arose.

The runnable implementation is `scripts/sv_phylogeny_benchmark.py`. Results,
editable vector PDFs, PNGs, and a compact synthetic input bundle are in
`docs/results/sv_phylogeny_20261028/`. All data in that directory are synthetic.

## Experiment

This is a small, balanced diagnostic panel from the earlier neutral
introgression simulations: four old imported mutations, four mutations arising
after introgression, and four ancestral-retention controls. The first four
eligible events in each history category were selected by replicate order from
the held-out regions, without filtering on inference performance. There are
eight independent source regions; the central imported-versus-new comparison
contains eight events in only four independent regions. The ancestral controls
are not all guaranteed to exhibit a specific ILS-discordant topology.

The original demography has a 2% pulse at 50 kya, a modern/archaic split at
500 kya, the EAS PHLASH population-size history, and heterogeneous recombination.
For each event, 32 modern haplotypes were selected using a fixed seed, with up
to 16 carriers and the remaining samples drawn from noncarriers. Sampling is
therefore enriched for carriers rather than random population sampling.

A contemporary chimp haplotype was added conditionally to each retained human
genealogy. The human roots were verified to be younger than a stipulated 6-Ma
species split. Ancestry was continued from that split with msprime, an ancestral
diploid Ne of 50,000, and recombination at 1e-8 per base per generation. Existing
human genealogies were preserved. A fresh, independent JC69 SNP overlay used
1.25e-8 mutations per base per generation and 25-year generations. This is a
simplified outgroup model, not a detailed reconstruction of human/chimp evolution.
There is no sampled archaic genome in the inference alignment.

The new SNP overlay and sample subset mean these results are **not a matched
head-to-head comparison with the earlier Gamma-SMC results**. Every method in
this new comparison receives the same human haplotypes and observed SNP states,
but uses a different amount of surrounding sequence.

The SV itself was withheld from sequence-tree inference, then mapped onto the
inferred trees as a single-origin, correctly polarized character. An exact
carrier clade is required for a one-branch placement. This specifically tests
the proposed infer-tree-then-place-mutation workflow; it does not test joint
inference with an SV monophyly constraint. The latter can recover information
that this diagnostic design deliberately leaves out of tree inference.

## Methods actually executed

**IQ-TREE + LSD2:** JC maximum likelihood, chimp outgroup, all 33 sequences
retained even when identical, 100 standard nonparametric bootstrap replicates,
then separate LSD2 dating with all tips contemporary and the mutation rate
fixed. Windows were 1, 10, and 25 kb and the true focal *human* genealogy
interval. The latter is an oracle boundary control unavailable in real data;
the chimp ancestral branches can still change within it. The 10-kb windows
contain 4–17 true human genealogies. Oracle intervals contain only 1–12 human
variable sites. The root was not calibrated to the simulated species split.

LSD2's `.timetree.nwk` file retains substitution units in this version, whereas
`.timetree.nex` contains dated edges and node-date annotations. The implementation
uses NEXUS and checks that root dates agree with root-to-tip distances. Tiny
rounding discrepancies in printed edges are handled with monotonic bottom-up
node heights. SNP invariant sites are included, and the focal SV base is masked.

**BEAST 2:** the oracle and 10-kb alignments, JC69, fixed strict clock, a constant
coalescent population size with a broad lognormal prior (log median 30,000,
log SD 2), and modern-human monophyly relative to chimp. This is a deliberately
simple, misspecified tree prior for an introgressed human/chimp alignment; it
does not model the species split or migration. IQ-TREE supplies an initial tree;
topology and node ages then vary in MCMC. Each analysis has two independent
seeds, initially two million steps per chain, sampling every 2,000 steps and
discarding 25% burn-in. One analysis (`r0023_1`, oracle window) was rerun with
four million steps per chain after failing the ESS threshold.

All 24 final analyses pass the implemented checks: minimum per-chain ESS above
200, split R-hat below 1.05 for posterior, likelihood, prior, population size,
root height, carrier MRCA, and its parent, plus between-chain agreement on
carrier-monophyly frequency within 0.05. The smallest ESS is 211.3 and largest
split R-hat is 1.0048. These checks do not prove complete topology mixing.
A placement is reported as usable only when carrier monophyly has posterior
probability at least 0.95 and convergence checks pass.

**Relate:** the entire 500-kb human haplotype segment, SNP ancestral states
polarized using the simulated chimp, mutation rate fixed, constant haploid
Ne=30,000, and a constant 1-cM/Mb map. Human sites with a third chimp allele are
excluded. Chimp is an external polarization source, not a haplotype inside
Relate's human coalescent model. The tree attached to the preceding SNP anchors
the focal coordinate; anchor distances are retained in the output. This leaves
breakpoint uncertainty between SNPs. Relate therefore borrows information over
a longer, recombining region than either conventional phylogenetic analysis.

WSL could not launch in this session. IQ-TREE and BEAST ran natively on Windows.
Relate was built from upstream commit
`b54ede259cbb0be095bc9c9a8bd18cdaf7e88b74`, with a Windows `getrusage` reporting
shim and a directory-cleanup change that retains scratch files when Windows
prevents deleting an open file. Inference algorithms were unchanged. The build
record and binary SHA256 are archived. This port was not cross-validated against
a Linux executable in this session. The Linux bootstrap installs the upstream
static 1.2.4 release, so exact cross-platform numerical identity is not claimed.

## Results

| Method / input | Available placements | True age within reported bounds | Median mutation-age fold error, available loci only |
| --- | ---: | ---: | ---: |
| IQ-TREE + LSD2, oracle human interval | 6/12 | 3/6 | 2.05 |
| IQ-TREE + LSD2, 1 kb | 5/12 | 4/5 | 1.54 |
| IQ-TREE + LSD2, 10 kb | 9/12 | 5/9 | 1.95 |
| IQ-TREE + LSD2, 25 kb | 6/12 | 3/6 | 2.05 |
| BEAST, oracle human interval | 6/12 | 4/6 | 1.67 |
| BEAST, 10 kb | 7/12 | 4/7 | 1.99 |
| Relate, 500 kb | 12/12 | 7/12 | 1.96 |
| True local tree | 12/12 | 12/12 | 1.73 |

The denominators and bounds have different meanings. IQ-TREE and Relate use a
compatible clade in a point tree and report its child/parent ages as conditional
branch bounds; these are **not confidence intervals**. Only six of the nine
10-kb IQ-TREE carrier clades have at least 95% bootstrap support. Relate branch
or topology posterior sampling was not run. BEAST mixes uniform mutation
placements over compatible posterior branches and reports conditional 90%
equal-tailed age intervals, subject to the criteria above. These intervals are
not empirically calibrated. Counts of bounds lying wholly on one side of
50 kya in the tables must not be read as equally confident statistical calls.

The point estimate for a fixed branch is its midpoint, equivalent to a uniform
placement prior on that branch. BEAST uses the median of the mixture over
posterior branches. Even the true-tree midpoint is not the true mutation age.
Conditional fold errors across different available subsets are not an unbiased
method ranking, and history-balanced ascertainment further affects coverage.

The recent-mutation group is particularly revealing:

| Method | Recent-mutation placements | Median fold error on those placements |
| --- | ---: | ---: |
| IQ-TREE + LSD2, 10 kb | 2/4 | 10.06 |
| BEAST, 10 kb | 1/4 | 17.01 |
| Relate | 4/4 | 6.03 |
| True tree, uniform branch midpoint | 4/4 | 3.71 |

All four true old-imported mutation branches span the 50-kya introgression
event in the 32-haplotype subset. Their timing remains ambiguous even with
perfect local trees. One of the four new-mutation true branches also spans the
pulse. No tested method should be interpreted as resolving that structural
ambiguity merely by assigning a midpoint or posterior median.

The practical result is that Relate is a promising existing tool for recovering
the relevant local clade in this setting. The evidence is insufficient to claim
that it accurately dates recent SV mutations or outperforms every calibrated
alternative. Longer windows do not consistently improve conventional tree
inference: they add SNPs while mixing genealogies. The next useful comparison
would condition tree inference on the SV's single-origin state pattern, test
population-history sensitivity, and evaluate more independent regions with
archaic SV observations where available.

## Artifacts and validation

- [Age comparison, vector PDF](results/sv_phylogeny_20261028/age_comparison.pdf)
- [All loci, including failures, vector PDF](results/sv_phylogeny_20261028/all_locus_outcomes.pdf)
- [Per-event estimates](results/sv_phylogeny_20261028/estimates.tsv)
- [Results by history](results/sv_phylogeny_20261028/by_history.tsv)
- [BEAST diagnostics](results/sv_phylogeny_20261028/beast_diagnostics.tsv)
- [Synthetic inference inputs](results/sv_phylogeny_20261028/benchmark_inputs.zip)
- [Commands, parameters, source and binary hashes](results/sv_phylogeny_20261028/provenance.json)

Seven targeted tests pass, covering mutation placement, incompatibility,
chimp-preserved human ancestry, corruption detection, time units, and convergence
diagnostics, and reuse of completed longer chains. All 12 data sets and 168 output/input receipts were audited. BEAST
tree-state indices match trace-state indices and all traces finish at the
requested chain lengths. PNGs were visually inspected. PDFs contain embedded
editable TrueType fonts and no raster image objects. Figures are saved to disk.

## Reproduce from the checked synthetic input bundle

Java 17+ and uv are required. The following Linux x86-64 workflow uses four
parallel jobs by default. BEAST is capped at a 5-GB Java heap per job; the default
four jobs fit within the requested 26-GB working budget. Local benchmarking used
up to 24 concurrent single-threaded inference jobs. Raw outputs are cached by
input hashes and parameters, not by the repository commit.

```bash
%%bash
set -euo pipefail
repo="$HOME/gamma_smc_ts"
if [ ! -d "$repo/.git" ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git "$repo"
fi
git -C "$repo" switch AOU_run_opt
git -C "$repo" pull --ff-only git@github.com:soisa001/gamma_smc_ts.git AOU_run_opt
cd "$repo"
uv run --no-project --python 3.12 python scripts/bootstrap_sv_phylogeny.py
out="$HOME/sim/sv_phylogeny_20261028"
for phase in restore iqtree beast relate; do
  bash scripts/run_sv_phylogeny.sh --out "$out" --phase "$phase"
done
bash scripts/run_sv_phylogeny.sh --out "$out" --phase beast \
  --event r0023_1 --window oracle --chain-length 4000000
bash scripts/run_sv_phylogeny.sh --out "$out" --phase analyze
bash scripts/run_sv_phylogeny.sh --out "$out" --phase audit
```

The input bundle eliminates the need to regenerate the earlier 160-region study.
`prepare` remains available to rebuild the panel from those original synthetic
tree sequences. The portable restore phase verifies archive and individual-file
hashes and refuses to overwrite different existing inputs. Repeating an inference
phase checks output hashes before reusing it. A completed longer BEAST run is
reused when its model, inputs, seed, binaries, and output hashes still match.

Method references: [IQ-TREE dating](https://www.iqtree.org/doc/Dating),
[LSD2](https://github.com/tothuhien/lsd2),
[BEAST 2](https://www.beast2.org/what-can-beast-do/),
[Relate paper](https://www.nature.com/articles/s41588-019-0484-x),
[Relate input and output formats](https://myersgroup.github.io/relate/input_data.html).
