# Mutation-context inference from assembled genomes: evidence and benchmark design

Research date: 2026-11-12. This note separates published empirical evidence from proposed simulation designs. Numerical pilot results belong in the experiment report, not in this literature note.

## Recommendation

Use **offspring base probabilities at supplied parental scaffold coordinates** as the primary task. Supply neutral calibration ancestor/descendant assemblies with known exposure, a new parental assembly, and query coordinates whose outcomes are held out. Ask for P(A), P(C), P(G), P(T) at each coordinate after the specified exposure. Aggregate these probabilities into expected mutation burden and directional counts for secondary scoring. Include targets with shifted sequence composition so that coarse rates cannot succeed simply by matching the calibration average. The prompt need not mention sequence context, k-mers, or any specific fitting method.

Simply asking for the mean mutation rate in the same fully observed ancestor/descendant pair is a weak benchmark: at low divergence, the number of differences divided by callable exposure already estimates the realized aggregate rate. Context becomes necessary when transporting a rate to a different composition, correcting unequal ascertainment, predicting local burden, or distinguishing mutation input from a downstream process. This is a mathematical property of the estimand, not a claim that any particular empirical context model must improve every task.

For offspring interpretation, specify a transmitted haplotype or a haploid mutation-accumulation lineage. Otherwise Mendelian inheritance at parental heterozygous sites, recombination, and phasing can dominate the requested base probabilities. Missing target outcomes must be held out from every training, alignment-derived feature, and hyperparameter-tuning stage. The known parental sequence remains available at those coordinates.

## Scoring the coordinate-probability task

For parental base a at site i, write `r_i = 1 - p_i(a)` and `q_i(b) = p_i(b)/r_i` for b different from a. In a one-generation rare-event model, `p_i(b) = mu_i(a>b)` and `p_i(a) = 1 - sum_b mu_i(a>b)`. With a longer supplied exposure, define the exact event model instead of silently treating an arbitrary branch probability as a one-generation rate.

Score full four-category log loss on held-out offspring bases. Its expectation under simulator truth is cross entropy; subtracting the oracle entropy gives KL divergence and permits stable evaluation even when one replicate contains almost no mutations. The binary mutation/no-mutation log score measures rate calibration, and the conditional alternate-base log score on mutation events measures spectrum discrimination. These are complementary: conditional scoring alone cannot test total mutation rate. The full log loss decomposes into the binary event loss plus the alternate-base loss for mutated observations.

Report scores in readable units, such as nats per million queried sites, alongside aggregate expected counts and observed/predicted ratios by prespecified probability bins. At human one-generation rates, almost every base is unchanged, so raw classification accuracy is essentially uninformative. A larger evaluation exposure or many independent descendants can provide enough events, but document that scaling. Do not select evaluation coordinates only because they mutated unless the benchmark explicitly conditions on mutation and provides the ascertainment rule. For a realistic single-generation release, exact expected scoring against hidden simulator probabilities is substantially less noisy than realized-event scoring.

A literature-derived 7-mer generator with 1-/3-/5-mer inference is a useful harder test: five bases may improve predictions without becoming the oracle. Include the source-derived 7-mer oracle only in evaluator diagnostics, never in solver-visible answers or features.

## Primary evidence

- **Chandra and Gao (2026)** is closely aligned with the motivating talk. It estimates CpG mutation rates using gnomAD polymorphisms, methylation, and a model correcting recurrent mutation. It studies CpG-centered 4- and 6-base windows, which are not identical to single-base-centered 3- and 5-mers. Both methylation and flanking sequence matter, including interactions. Expanded flanks improve fit, but a model separating the two flanking sides approximates the saturated model with far fewer parameters. Consequently, a sequence-only benchmark should hold the methylation distribution fixed or explicitly define its target as the rate averaged over that distribution. The paper was published June 1, 2026; the lab's publication list still points to its earlier preprint. [Published paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011957).
- **Carlson et al. (2018)** provides particularly direct support for five bases beyond three: every trinucleotide mutation subtype showed rate heterogeneity among its 16 constituent five-base contexts, with all reported tests having P < 10^-231. Their model uses approximately 36 million singleton variants from 3,560 genomes and validates predictions against de novo mutations. Their supplementary data contain context-specific relative rates. Very small significance values establish detectable heterogeneity in this large dataset; they do not quantify a guaranteed improvement in a smaller downstream task. [Paper](https://www.nature.com/articles/s41467-018-05936-5), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/30218074/).
- **Aggarwala and Voight (2016)** uses 1000 Genomes variation to show that seven-base context explains more than 81% of variation in substitution probabilities in its analysis, and applies the model to functional intolerance. The result concerns polymorphism/substitution probabilities, not an unconditional claim that sequence explains 81% of all germline mutation-rate variation. It supports considering larger context but does not establish five bases as biologically sufficient. [Paper record](https://pubmed.ncbi.nlm.nih.gov/26878723/).
- **Jónsson et al. (2017)** provides a human pedigree calibration of approximately 1.29e-8 mutations per base per generation, with parental-age and sex effects. A computationally accelerated simulation may use larger branch mutation probabilities, but should record the exposure explicitly and avoid calling those probabilities biological per-generation rates. [Paper](https://www.nature.com/articles/nature24018).
- **Glémin et al. (2015)** infers gBGC strength from derived allele-frequency spectra while accounting for mutation bias, polarization error, and spatial heterogeneity. Its design underscores why a single pairwise assembly comparison does not supply all information required for biological gBGC estimation. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4510005/).
- **Williams et al. (2015)** directly studies non-crossover gene conversion using pedigrees. It estimates both an affected-site rate (5.9e-6 per base per generation) and GC transmission probability (68%, with uncertainty). These are distinct quantities, and neither is automatically the same as an effective GC-favoring substitution coefficient. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4404656/).

## Candidate outputs and what each measures

Let c denote ancestral sequence context, a its central base, b a different base, E_c the number of callable base-generations in calibration, Y_cb the observed number of a-to-b events, and mu_cb their rate. For a target assembly, let N_c be its callable opportunity count and t its known exposure. W means A/T; S means G/C. Equations below define proposed benchmark metrics and are our derivations rather than published performance claims.

| Output | Definition | Why context can improve it | Priority |
| --- | --- | --- | --- |
| Target mutation burden | Lambda = t sum_c N_c sum_b mu_cb | Calibration and target can differ in motif frequencies despite similar GC content | Primary |
| Target directional burden | Lambda_WS and Lambda_SW from the corresponding terms above | Different mutation directions can respond differently to outer flanks | Primary |
| GC flux | (Lambda_WS - Lambda_SW)/(t sum_c N_c) | Predicts immediate expected change in GC fraction under the fitted mutation process | Primary |
| Standardized mutation spectrum | Apply each group's rates to the same supplied target sequence/opportunity distribution | Separates rate-process change from sequence-composition change | Primary |
| Relative rate multiplier alpha | Fit target counts with expectation alpha times its calibration-based expectation | Coarse models can misattribute a composition shift to altered rate | Primary |
| Effective GC-favoring log odds B | Excess log directional substitution odds after correcting neutral mutation opportunities | Coarse neutral expectations can mimic a GC-favoring process | Secondary, with explicit definition |
| Expected functional burden | Sum mutation rates over changes annotated as nonsense, missense, or a specified motif disruption | Codon and regulatory composition change the available mutation routes | Useful extension; requires annotation |
| Regional observed/expected depletion | Observed substitutions divided by the neutral predicted burden for each region | Prevents low-mutability sequence from automatically appearing constrained | Useful extension; not a unique selection coefficient |
| Long-run equilibrium GC fraction | Stationary composition under a fully specified sequence-evolution process | Motif frequencies change as the sequence evolves | Advanced extension |

For the effective bias task, a clear synthetic parameterization is

`lambda_WS = alpha * Lambda_WS * exp(B/2)`

`lambda_SW = alpha * Lambda_SW * exp(-B/2)`.

With nonzero directional counts, the direct point estimate is

`B_hat = log(Y_WS / Lambda_WS_hat) - log(Y_SW / Lambda_SW_hat)`.

Under this definition, B is a log ratio of directional multipliers. It is not a conversion-event probability. A physical gBGC model would separately specify demography, allele frequencies, recombination/conversion initiation, tract lengths, transmission bias, and their relationship to fixation. Zero counts require likelihood inference or a stated prior instead of an arbitrary undisclosed pseudocount.

Likewise, the two-state formula `GC* = u/(u+v)` for W-to-S rate u and S-to-W rate v describes a context-independent approximation. Under context dependence, the motif distribution evolves, so rates averaged over the starting genome do not automatically give the exact long-run equilibrium. Immediate GC flux is a cleaner first target.

## Essential identifiability assumptions

1. **Direction:** use a known ancestor, or enough outgroups and a documented polarization model. Arbitrarily labeling one extant assembly ancestral changes the interpretation of directional rates.
2. **Exposure:** mutations per generation require a generation or branch-exposure calibration. Without it, the identifiable object is rate times exposure.
3. **Mutation versus substitution:** neutral mutation-accumulation pairs can calibrate mutation input. Differences between natural species also reflect ancestry, fixation, selection, gene conversion, recurrent hits, and error.
4. **Callable opportunities:** every fitted rate needs the number of sites that could have undergone that mutation. Counting only mutated 5-mers estimates composition among events, not mutation probability.
5. **Ancestral context:** take flanks from the ancestor and require a callable full window. Descendant context can already contain additional mutations; near an indel, coordinate translation and alignment uncertainty matter.
6. **Shared process:** transport assumes the same mutation mechanism, or a specified multiplier, across calibration and target. Changing methylation or repair biology at the same time invalidates a pure composition-transport interpretation.

The first pilot can intentionally use error-free, collinear haploid assemblies with sparse substitutions and known ancestry/exposure. This isolates the statistical question. It should not be described as validation of a general de novo assembly variant caller. Fragmentation, repeats, rearrangements, indels, uncertain bases, assembly errors, and ancestral uncertainty are later stress tests. A single whole-genome global alignment is not a realistic universal solution for those cases.

## Experimental design that specifically tests five bases

Compare hierarchical 1-, 3-, and 5-base models using the same training variants, opportunities, callable mask, and exposure. With reverse complements collapsed, these have 6, 96, and 1,536 possible substitution channels. A saturated 5-mer model can lose to a 3-mer model in sparse data; shrinkage and held-out evaluation are therefore essential. Tune shrinkage on separate calibration validation data, not on hidden target truth.

Use a sequence-composition shift that leaves shorter contexts matched. One proposed generator comprises four independent stationary Markov chains, interleaved by position modulo four. Each chain has uniform marginal base composition and persistence `Pr(X_i = X_{i-4}) = q`. Consecutive four-base windows are uniform in expectation, while the outer bases of a centered five-base window match with probability q. Changing q between calibration and target therefore changes a fifth-base interaction without changing the expected 1-/3-/4-mer distribution. Finite samples are only approximately matched: report observed motif distributions and distances. Alternatively, explicitly construct or resample sequences with fixed shorter-word counts.

A planted five-base rate interaction makes the causal demonstration transparent, but it does not measure human effect size. Include all of these controls:

- A context-independent generator: complex models should have no systematic advantage.
- A 3-mer-only generator with the same composition shift: tests whether apparent 5-mer benefit is merely CpG modeling.
- A 5-mer generator with matched calibration/target composition: total-burden improvement can shrink even when local predictive likelihood improves.
- A 5-mer generator with shifted outer flanks: the intended transport test.
- Weak and strong outer-flank effects, plus multiple training sample sizes: expose the bias/variance tradeoff.
- A literature-derived relative-rate kernel: tests whether the conclusion survives removal of the planted outer-flank rule.

Use fixed paired seeds across models, report per-seed errors, and distinguish truth from realized counts. Recommended scores are relative error of expected burden, absolute error/RMSE of alpha and B, predictive Poisson or multinomial deviance on independent descendants, and interval coverage where intervals are implemented. For a zero GC-flux truth, absolute error is preferable to relative error. Do not score only in-sample fit or declare success from one favorable seed.

Keep per-site mutation probabilities low. If the simulator permits only one event per site using ancestral flanks, identify it as a one-step rare-event approximation. Longer evolutionary exposures require explicit recurrent events and evolving neighborhoods; simply dividing a high-divergence difference count by time is biased.

## Accessible empirical rate tables

The Carlson project's [official Mutation Rate Browser](https://mutation.sph.umich.edu/) links to the author's [Mr. Eel repository](https://github.com/carjed/mr-eel). Its [ERV_7bp_rates.txt](https://raw.githubusercontent.com/carjed/mr-eel/master/ERV_7bp_rates.txt) was fetched successfully in memory on the research date: 675,933 characters; 8,192 context rows plus a header. It has columns `Sequence`, `AT_CG`, `AT_GC`, `AT_TA`, `GC_AT`, `GC_CG`, `GC_TA`. A sequence is stored with its reverse complement, e.g. `AAAAAAA(TTTTTTT)`. These are relative singleton-based rates, not absolute per-generation rates. The author's [annotation script](https://github.com/carjed/mr-eel/blob/master/mr_eel.pl) confirms the following coding:

| Column | Mutations |
| --- | --- |
| AT_CG | A>C or T>G |
| AT_GC | A>G or T>C |
| AT_TA | A>T or T>A |
| GC_AT | G>A or C>T |
| GC_CG | G>C or C>G |
| GC_TA | G>T or C>A |

Nonapplicable central-base channels are zero. Expand both sequence orientations and preserve the matching destination allele when constructing a directed rate tensor. A practical second pilot can uniformly average the outermost pair of bases to produce a 5-mer kernel and normalize its mean under a stated reference composition. This is an **empirically shaped synthetic kernel**, not a reconstruction of the published human 5-mer rate estimate: original human opportunity-weighted marginalization is different. Preserve the original downloaded file, source URL, checksum, normalization, and averaging rule. An unmodified 7-mer kernel is also useful as a deliberate misspecification test for a fitted 5-mer model.

Chandra and Gao's [Supplementary Table S1](https://doi.org/10.1371/journal.pgen.1011957.s014) contains scaled rates and standard errors by CpG context, species, and methylation status. It is attractive for a more focused CpG experiment; its scaled rates need their own exposure calibration, and its CpG-centered even-length windows must not be silently relabeled as mutation-centered odd-length windows.

## Example task prompt without method hints

> Calibration data contain ancestral and descendant DNA assemblies from independent mutation-accumulation lineages, elapsed generations, and callable-region definitions. A separate file contains parental scaffolds and query coordinates for new lineages. For each query, predict the probability that the corresponding descendant base is A, C, G, or T after the specified number of generations. The parental sequence describes the transmitted haplotype, and only single-base substitutions are scored. Return one row per query in the supplied schema; each row must contain four nonnegative probabilities summing to one. Develop and validate your inference method using the calibration data. Target descendant bases are withheld for evaluation.

An aggregate-output alternative is:

> You are given ancestral and descendant assemblies from calibration mutation-accumulation lines, their elapsed generations, and callable-region definitions. You are also given ancestral assemblies for a new set of lines maintained under the same conditions. Estimate the expected total number of single-base substitutions, the expected A/T-to-G/C and G/C-to-A/T counts, and the expected change in GC fraction after the stated number of generations. Write the requested estimates and uncertainty measures to the supplied output schema. Develop and validate an inference method using the provided data; do not use hidden evaluation files.

For an alpha/B task, add target descendants and define those parameters mathematically in the prompt. Defining the requested estimand is necessary for a fair task and does not reveal the context-dependent solution. Keep evaluator truth, generator kernels, rate tables, analysis reports, suggestive filenames, and helper functions outside the solver-visible package. Neutral filenames alone do not prevent leakage if the solver can inspect the simulation source.
