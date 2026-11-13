# Chronological decision log

## 1. Input inspection

Read only the six supplied inputs. The task is one-transmission haploid allele prediction with error-free homologous assemblies and circular scaffolds. Metadata supplies no numerical mutation model. Coordinates are one-based; all four probabilities must include the unchanged allele.

Initial hypotheses: a uniform source-base mutation matrix; dependence on local sequence context; longer sequence motifs or compositional domains; scaffold-specific heterogeneity; circular-boundary effects. I will estimate these from calibration data and compare on held-out scaffolds before selecting a model. New-parent outcomes are unavailable. All stochastic operations will use seed 0, four threads, and memory below 26 GB. Input files will not be edited.

## 2. Basic calibration audit

There are 80 calibration scaffold pairs of 50,000 bases (4,000,000 opportunities), 19,901 changes (0.00497525), and 20 new-parent scaffolds of the same length. All sequences contain only A/C/G/T and paired lengths match. Mutation counts and destination counts look approximately uniform at this coarse scale. No numerical generating parameters are inferred merely from their round appearance. Input hashes were recorded in input_audit.json.

Before detailed feature discovery, partition scaffolds with numpy default_rng(0).permutation(80): 40 discovery, 20 tuning/validation, 20 final evaluation. Detailed contextual, positional, and offset diagnostics will use only discovery. Whole-molecule circular offsets allow us to check for templating/copying signals that short-context tables could miss. Any resulting candidate must improve held-out predictive log loss.

## 3. Discovery findings and first validation

Discovery rate is 0.004951. Single-base and centered 3-mer mutation tests show no evidence of heterogeneity. Centered 7-mers have excess dispersion (chi-square 19,137.95 on 16,383 df; sparse cells make the nominal tail probability only diagnostic). Circular all-offset scans show no exceptional effect after accounting for tens of thousands of comparisons. Position and scaffold tests also do not indicate structure.

Smoothed context tables predict the separate validation scaffolds better than uniform: centered 7-mer gains 99.76 nats and centered 9-mer gains 104.28 nats in binary mutation log likelihood. This confirms a real, weakly localized contextual signal. Individual high-count contexts are sparse and some fail to replicate; no single top sequence will be treated as a certain rule. Next compare asymmetric contiguous windows and empirical-Bayes hotspot mixtures. A mixture can distinguish a small set of strongly mutable motifs from diffuse small differences more efficiently than a fixed pseudocount. These are candidate models, not assumptions about the hidden generator.

## 4. Context localization and destination structure

Scanning contiguous windows of 3–10 bases that contain the target localizes the transferable signal to the 8-mer at offsets -4 through +3. Its two-component rate model gains 333.67 validation nats, compared with roughly 100 nats for its best 7-mer subwindows and 106 for centered 9-mers. Increasing resolution beyond the 8-mer reduces sharing and predictive value. This supports an 8-mer lookup effect rather than broad GC or a simple single-position motif effect.

The eight-base context also predicts which alternative allele occurs. Among 1,689 discovery pairs of mutations sharing a context, 1,165 (68.98%) have the same destination, versus one third under uniform alternatives. A context-specific Dirichlet-smoothed destination model with total pseudocount 1 gains another 356.56 validation nats. A fitted beta prior for rate has shape parameters 0.9031 and 181.569; a flexible discrete rate prior gains 350.49 validation nats. I will compare a few distributional alternatives for destination probabilities, lock the specification using validation only, then evaluate once on the 20 untouched test scaffolds. The initial global count audit included all scaffolds; detailed model development has excluded the test partition.

## 5. Model specification locked before final test

Removing any of the eight positions reduces validation rate gains to 89–107 nats, confirming that all eight positions contribute to lookup identity. A smooth two-parameter logistic-normal rate prior gains 350.04 validation nats, essentially equal to the much more flexible discrete prior's 350.49. Choose the logistic-normal prior for parsimony and stability rather than its negligible validation difference. Its discovery-fit logit mean is -5.70128 and standard deviation 0.89794. The conditional alternative distribution is a symmetric Dirichlet with fitted total concentration 0.977995, gaining 356.72 validation nats; fixed-preferred-alternative and simplex-grid priors gain 353.27 and 355.55, respectively. Choose the simpler Dirichlet formulation.

Locked final model: circular parent 8-mer offsets [-4,+3]; exchangeable context-specific change rates under a fitted logistic-normal prior; exchangeable context-specific alternative-allele probabilities under a symmetric Dirichlet prior; posterior predictive allele probabilities. Fit hyperparameters and context counts on the combined 60 discovery/validation scaffolds, measure final performance once on the 20 test scaffolds, then refit this same specification to all calibration data for delivery. No choice will be revised based on test outcomes.

## 6. Final evaluation and delivery

The locked model improved final-test allele log likelihood by 889.58896 nats over the uniform baseline across 1,000,000 sites. A scaffold bootstrap (5,000 draws, seed 0) gave a gain interval of 792.936–991.619 nats. The test set had 4,972 changes; the selected model expected 4,973.801. Its mean allele negative log likelihood was 0.035903488 versus 0.036793077 for uniform, and mutation AUC was 0.584708. The model specification was retained without modification.

Refitting all 4,000,000 calibration sites gives logit-rate mu=-5.69395284739004, sigma=0.8946945036380947, and alternative Dirichlet concentration tau=1.048019815816202. The delivered million-row predictions.tsv expects 4,974.728737409099 changes. Every new-parent base is queried exactly once. A separate verification script checked 4,000 circular context encodings, including scaffold ends, and all million output rows. Maximum row-sum error is 5.50e-13; exact query order and all six original input hashes pass.

The fixed solve pipeline was rerun unchanged while recording final source hashes and verifying reproducibility; this recomputes the same held-out metrics without selecting or adapting a model. The exploratory two-point-mixture helper was guarded against extinct-component numerical divisions; that helper is not used by the final model. REPORT.md explains all hypotheses, selection, formulas, held-out results, reproducibility commands, and limits of conditional uncertainty estimates. No plot or sequence was displayed inline, and no repository or external outcome was accessed.
