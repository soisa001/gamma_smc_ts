# Chronological decision log

## 1. Initial inspection

Read only the six supplied files in this case directory. The task specifies error-free homologous haploid assemblies and single nucleotide changes over one transmission. Metadata provides no process parameters. Calibration FASTA files are about 4 MB each, the new parent about 1 MB. Queries start at position 4, suggesting that local sequence context may matter; this is a hypothesis, not an assumed generating rule.

Plan: quantify change counts, substitution spectra, sequence-context effects, position/contig heterogeneity, and neighboring-change dependence. Compare parsimonious models using held-out contigs. All random operations use seed 0; numerical thread limits are four. No external scientific references or other project directories will be inspected.

## 2. Initial empirical findings

Calibration has 80 contigs of 50,000 bases, with 20,138 differences among 4,000,000 bases (0.0050345). The new parent has 20 contigs of the same length. C followed by G and G preceded by C have elevated change rates near 0.011; corresponding non-CpG C/G rates are about 0.0032. A/T changes are near 0.0049. More distant single offsets show no obvious large effects. Neighboring changes do not show strong clustering at tested lags.

Some 5-mer contexts have more changes than the simple CpG model would suggest, so a CpG-only model is premature. Alternatives now being evaluated: base-only, reverse-complement pooled base and CpG models, and 3/5/7-mer models with a wide grid of shrinkage toward CpG probabilities. Validation is seeded five-fold cross-validation by entire contigs. The first exploration suggested candidate models using all calibration outcomes; cross-validation will measure predictive performance but is not an untouched confirmatory test of those initial observations.

## 3. A compact longer-context hypothesis

The 5-mer Pearson residual dispersion against the CpG model is 3.85 after reverse-complement pooling, compared with 1.17 for 3-mers. Inspection of the largest 5-mer residuals shows matching bases at offsets -2 and +2 (examples GCATG, CCAGC, TAACT, AAAAA). Direct counts show this equality increases A/T change probability from 0.00318 to 0.00996. In A-oriented coordinates, A-to-C and A-to-G increase, while A-to-T stays near 0.00075. Equality at offsets +/-1 or +/-3 has no analogous strong effect. C/G context rates do not show a comparable equality effect. I will compare a four-group motif model, a six-group model permitting equality effects for C/G, unpooled versions, and richer context models shrunk toward this motif model. I will also test a lower-parameter model sharing background transversion/transition rates across classes. These hypotheses were formed by inspecting the full supplied calibration, a limitation to report explicitly.

## 4. Validation and final model choice

The initial context-table grid favored reverse-complement pooled 5-mers with shrinkage 2,000 toward CpG (held-out log loss 0.03523016). A compact five-parameter motif model improves this to 0.03508965. A six-parameter model separating the two background transversion rates has the lowest tested loss (0.03508920), but its advantage is only 4.50e-7 per site, versus paired-contig standard error 5.73e-7. A 5-mer residual model with strong shrinkage is similarly close; richer and less regularized tables are worse. I select the five-parameter model as the parsimonious option within one paired standard error of the best. The complete four-group model does not significantly improve fit over it (likelihood-ratio p=0.131). The six-parameter contrast has nominal p=0.020 and is retained as uncertainty about background transversion equality, rather than asserted absent.

The first numerical optimization attempt terminated abnormally at very tight finite-difference tolerances. Replacing finite-difference gradients with exact analytic gradients resolved the issue; all final model fits require reported optimizer success. No failed fit is used. Final prediction generation uses explicit row-wise allele indexing and includes round-trip checks.

## 5. Final fitting and output verification

Fit the selected model on all calibration sites with a complete seven-base window (3,999,520 sites, 20,138 changes). The wider window is used consistently for all model comparisons although the final rule only needs two flanking positions. Generate all requested rows in original order. Estimate parameter uncertainty with 500 seeded whole-contig bootstrap replicates and four workers. Verify input hashes, query identities, nonnegative finite probabilities, unit row sums, and expected-change totals after serialization. Detailed results are recorded in summary.json and REPORT.md.

## 6. Bootstrap numerical refinement and context-shift check

Some bootstrap replicates still reached line-search precision limits despite analytic gradients. I added an exact-Hessian Newton fallback in scaled probability coordinates, requiring a maximum absolute likelihood gradient below 1e-10. The complete 500-replicate run then succeeded; no failed replicate was dropped. The primary final fit and all original cross-validation fits had already succeeded and are unchanged.

New-parent motif exposure is markedly different: 398,545/498,027 A/T sites have equal +/-2 flanks (80.025%), versus 500,124/1,998,519 in calibration (25.025%). This explains why the final expected total is 6,897.167986 rather than about 5,034 from simply extrapolating the overall calibration rate. I reweighted cross-validation to the new four-group mixture. The six-parameter model's tiny advantage remains smaller than its paired standard error; exact neighboring-base splits and the full four-group model are worse. This confirms the original compact-model choice under the observed mixture shift. These additional diagnostics inspect only the provided new parent and calibration outcomes.

All 999,880 prediction rows passed exact query-order/coordinate checks, finite/nonnegative checks, and row normalization (maximum serialization deviation 4.53e-13). Expected count from the serialized TSV agrees within 4.3e-7. Input hashes remain unchanged. Bootstrap SE of the expectation is 54.5311, percentile interval [6795.6964, 7004.1816]. The report explicitly retains uncertainty from residual five-base dispersion and adaptive model discovery. No plots were created.
