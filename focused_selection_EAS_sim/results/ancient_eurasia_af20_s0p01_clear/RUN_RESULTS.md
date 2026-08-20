# AncientEurasia 9K19: clear focal results

This additive bundle reused completed simulations and Gamma-SMC decodes; it did not simulate or decode any unit.

## Statistical contract

- Primary: selected overall focal P(TMRCA < x) versus the new single-copy, complete-simulated-Q5-present-Han-census-segregating neutral bank (100 units, s=0; unscaled Ne=6300, scaled census=1260 diploids/2520 genomes; no terminal AF target/band, selected-AF matching, pool-frequency gate, or panel-detection gate beyond the strict census 0<AF<1 condition). Because selected panels are exact-AF20 conditioned and this neutral bank is not, these are simulation-reference tails combining selection and endpoint-AF enrichment—not exchangeability-valid selection-only p-values.
- Selection-isolating sensitivity: selected alt/alt-minus-ref/ref versus the existing final-AF20-matched neutral contrast.
- Genotype-class limitation: each selected panel has only 3-6 alt/alt pairs (median 4.0), versus 63-66 ref/ref pairs (median 64.0). The alt/alt curves are therefore noisy descriptive summaries, not the primary test.
- All seven thresholds were prespecified. 50 kya is a pulse-aligned display cutoff; the primary simulation-reference minP retains whole seven-threshold vectors.

## Pulse-aligned display summary (50 kya)

| Source | Median selected score | Median neutral score | Median reference p | p<=0.05 fraction |
|---|---:|---:|---:|---:|
| Tree truth | 0.14 | 0.135 | 0.50495 | 0.000 |
| Gamma-SMC | 0.117684 | 0.130359 | 0.613861 | 0.000 |

## Corrected two-epoch comparison

The attached two-epoch neutral-null histogram is tree truth: its focal upper-tail result was 1/101. The corresponding Gamma-SMC focal decode of that linked low-AF selected base was 0.0305024 versus a neutral mean of 0.0234616, with 13/100 neutral exceedances and p=14/101=0.1386. It was the +/-500 kb positive-window density statistic—not the focal Gamma-SMC p-value—that reached 1/101.

This is a genealogy-linked diagnostic, not a fully identical mutation-layer comparison: the truth folder used mu=1.29e-9, while the corresponding official Gamma-SMC v0.2 container decode used mu=1.25e-8. It is not a native-optimized decode.

The two designs are also materially different: s=0.05 versus s=0.01, a 4.5-kya de novo allele versus 60-kya archaic standing variation, 2,000 within-diploid pairs versus 100, and 1-kb versus 10-kb Gamma output stride.

## Frequency-path limitation

No historical selected population-frequency trajectory was persisted. The timing figure is a schematic and the frequency panel shows observed uniform-pool/sample endpoints only; neither is a reconstructed trajectory.

## Output tables

`unit_scores.tsv` contains the complete seven-threshold truth/Gamma score grid. `selected_pointwise_pvalues.tsv` and `selected_minp.tsv` retain primary and sensitivity estimands separately. `input_cache_audit.tsv` lists every verified input receipt/artifact hash.
