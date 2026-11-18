# Cumulative posterior mass at intermediate selection strengths

Each of the 18 EAS treatments contains 100 independently retained replicates with a focal allele segregating in the present-day population. The simulations cross two origins (de novo and a 1% Neanderthal pulse at 50 kya), three selection windows (10 kya–present, 50 kya–present and 50–30 kya), and three additive selection coefficients (0.002, 0.003 and 0.004). Each replicate uses 100 sampled diploids and the same 1,000 haplotype pairs.

The completed run retained 900 de novo histories after 231,042 losses and 900 introgressed histories after 5,606 losses. No proposals fixed before or at the present, so these 18 treatments needed no tag fallback.

At the focal site (500 kb), the statistic is the mean pairwise posterior probability of a TMRCA below each cutoff from 1 to 50 kya in steps of 1 kya. The one-sided pointwise rank p-value is (1 + number of reference statistics at least as large as the query) / 101. A replicate is called once if any cutoff within the specified window, including its endpoints on the grid, passes the indicated strict threshold.

`rates.tsv` gives the TPR among selected replicates and the FPR among 100 held-out neutral replicates scored by the same rule. The reference and held-out neutral banks come from the preceding cumulative analysis. The de novo banks are matched to the allele's age at introduction. Both windows starting at 50 kya use the same de novo neutral banks.

`mass.npz` contains focal cumulative curves, source identities and present-day frequencies. `pointwise_p.npz` contains the replicate-by-cutoff p-values. `simulation_manifest.json` and `inputs.json` record the model, binary and input checksums, seeds, source paths and score checksums. The report phase of `scripts/run_intermediate_cumulative_grid.py` regenerates the tables and figures from these portable results.
