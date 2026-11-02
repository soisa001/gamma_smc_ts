# Time-resolved selection: literature and design implications

Research date: 2026-11-02. Scope: primary-source context for an EAS all-pairs TMRCA-excess scan; no new simulation results.

## Akbari analogy

The likely reference is Akbari et al., *Ancient DNA reveals pervasive directional selection across West Eurasia*, Nature (2026), DOI [10.1038/s41586-026-10358-1](https://doi.org/10.1038/s41586-026-10358-1). Its 2024 preprint was titled *Pervasive findings of directional selection realize the promise of ancient DNA to elucidate human adaptation*.

The study tests allele-frequency trends against sample dates using a generalized linear mixed model that accounts for genetic relatedness/population structure. Extended Data Figure 9 estimates changing polygenic selection intensity in overlapping 2,000-year windows, advancing 100 years, requiring at least 2,000 people per window. This is a direct precedent for time-resolved windows. It uses dated ancient allele observations rather than pairwise coalescence times. [Primary paper, especially Extended Data Figure 9](https://doi.org/10.1038/s41586-026-10358-1).

## Closer genealogy-based comparison

Stern, Wilton and Nielsen's CLUES method models allele-frequency trajectories and selection through genealogies, integrating uncertainty in local trees and trajectories. It can test changes in selection strength and estimate sweep onset. This supplies a closer conceptual reference for the distinction between observed genealogical events and the latent evolutionary process producing them. [Stern et al., PLOS Genetics (2019)](https://doi.org/10.1371/journal.pgen.1008384).

## Implications for the proposed experiment

The following are design deductions, not claims that the cited papers validated this exact statistic.

- A moving upper cutoff gives the cumulative fraction F(t) = P(TMRCA <= t). A local peak in coalescence density requires interval mass F(b) - F(a), a smoothed derivative, or a kernel density. A cumulative excess can remain elevated long after its generating coalescence pulse.
- Keep the all-pairs denominator fixed across time; dividing interval mass by window width yields a density. A conditional coalescence hazard is a different estimand because it divides by the fraction of pairs not yet coalesced.
- Pairwise observations share haplotypes and tree branches. Calibrate the entire time scan with independent matched simulation replicates; treating pairs or adjacent time bins as independent would overstate information.
- Initially call the recovered object a **coalescence-excess interval**. Selection active from 30 to 20 kya need not produce excess limited to those same bounds. The relation depends on starting allele frequency, standing versus new variation, sweep strength, fixation/loss, recombination, demography, and genealogy reconstruction. Infer selection timing only after measuring the mapping in simulations.
- Compare true-tree and decoded-TMRCA performance on the same selected and neutral replicates. Fine scan steps increase sampling density but cannot create temporal resolution absent from the decoder.
- For a strict localization endpoint, score interval overlap with [20, 30] kya separately from genome/locus detection. Also report peak-time error and interval intersection-over-union; an arbitrarily broad interval otherwise earns an easy overlap success.
- Distinguish the probability of a completely disjoint interval call in selected replicates (temporal mislocalization) from the probability of any interval call in neutral replicates (neutral false-positive rate). Time-level false-positive duration and event-level false discoveries answer additional questions; F1 requires explicitly matched units for precision and recall.

These points motivate a pilot comparing cumulative cutoffs with interval-mass scans, followed by held-out calibration of full-curve thresholds and temporal accuracy before applying the procedure to observed data.
