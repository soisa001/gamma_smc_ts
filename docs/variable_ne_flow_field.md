# Gamma-SMC flow fields and a possible variable-Ne extension

Research and source audit, 9 September 2026. Production decoding remains unchanged: fixed Ne, mutation rate, and recombination rate. This note proposes a future decoder extension; it does not report an implemented variable-Ne decoder or an improvement in selection power.

**A variable Ne(t) flow field is feasible. The transition table can remain two-dimensional for a fixed demographic history, but a correct extension must also change the initial prior, forward/backward combination, posterior summaries, and clipping rules. Replacing only the prebuilt table is insufficient.** A particularly useful design is a two-parameter exponential tilt of the demographic coalescent prior. It preserves exact mutation updates and forward/backward combination within that family, and becomes the ordinary gamma family when Ne is constant.

## What neutral calibration already addresses

The user's point about simulation calibration is correct. If independent neutral regions come from the intended null simulation model and pass through the identical decoder and detector, empirical rank calibration can control the neutral-region call probability under that model even when the decoder is misspecified. A claim that constant Ne necessarily invalidates those calibrated calls would be wrong.

Calibration does not reconstruct information that the decoder discarded or mixed across locations. It sets an appropriate threshold for the resulting statistic. Selected simulations measure the power remaining at that threshold. Thus a demographic decoder would be an attempt to improve discrimination and age inference, not a prerequisite for the existing matched-simulation null. Rule selection, cutoffs, run searches, and any tuning must remain inside the training/calibration procedure. These are statistical implications of the experiment's design, not a claim made by the Gamma-SMC paper.

The original paper explicitly acknowledges constant population size in flow construction and describes demographic flow fields as a possible extension. Its demographic benchmark does not establish that a demographic extension improves this EAS introgression scan. [Schweiger and Durbin, 2023, main paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538485/)

## How the current flow field works

Let physical TMRCA in generations be g, choose a reference diploid population size N0, and write t = g/(2N0). The current parameters are theta = 4N0*mu and rho = 4N0*r. These are unit conversions; a reference N0 can still be used when the actual population size varies with age. [Current rate documentation](../README.md)

For each haplotype pair and genomic position, the decoder stores only two numbers. They describe a gamma approximation to the TMRCA distribution,

\[
f(t;\alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}t^{\alpha-1}e^{-\beta t},
\qquad E[t]=\alpha/\beta,\quad CV=1/\sqrt\alpha.
\]

Here beta is a **rate**, not a scale. The working coordinates are log10(mean) and log10(CV). A flow-field arrow says how those two coordinates change when moving one small recombination-distance step along the genome. The arrow is indexed by the current uncertain TMRCA distribution; it is not indexed by genomic coordinate or by historical date. [FlowField::at and coordinate conversion](../src/flow_field.h)

At each observation there are two conceptually distinct updates:

1. **Transition:** Allow recombination to change the local genealogy. Apply the SMC-prime transition kernel to the current TMRCA density, then approximate the result within the chosen two-parameter family.
2. **Emission:** Incorporate whether the two sequences match. Under the decoder's Poisson mutation approximation, a homozygous base multiplies the density by exp(-theta*t); a heterozygous base multiplies it by t*exp(-theta*t), up to a constant. The gamma update is therefore alpha += y and beta += theta, where y is zero or one. Missing bases omit the emission but still allow recombination.

The implementation precomputes compositions of transitions and homozygous/missing stretches, then interpolates cached maps to skip many bases at once. It also caches whether an informative site comes before or after a stretch, because forward and reverse operations have different ordering. This explains both the speed and why a changed demographic field must invalidate the derived cache. [FlowFieldCache::initialize_cache](../src/flow_field.h), [CachedPairwiseGammaSMC forward/backward passes](../src/gamma_smc.h)

### Constructing the arrows

The generator does **not** simply match the mean and variance of a recombined density. At each grid point, it computes the infinitesimal change A*f in the density and projects this onto the two tangent directions of the gamma family by weighted least squares. In ordinary alpha/beta coordinates, the tangent functions are

\[
\partial_\alpha f=f[-\psi(\alpha)+\log\beta+\log t],\qquad
\partial_\beta f=f[\alpha/\beta-t].
\]

For log10 coordinates, the columns receive the chain-rule multipliers alpha*ln(10) and beta*ln(10). The implementation evaluates the columns and target on a time grid, weights each row by sqrt(delta_t), and solves the two-column least-squares problem with Eigen's SVD. If the resulting velocities in log10(alpha), log10(beta) are u and v, the stored velocities are u-v for log10(mean), and -u/2 for log10(CV). [calculate_updated_gamma_parameters](../src/generate_canonical_flow_field.cpp)

This is a local approximation. Projection error, grid interpolation, repeated stepping, clipping, and restriction to a two-parameter density all remain even with a perfectly specified demographic history. A denser field only addresses numerical approximation of that same family.

The built-in grid is 51 mean values from 1e-5 to 100 in scaled time, and 50 CV values from 0.01 to 1. It is embedded in src/io.h. The README's reference to resources/default_flow_field.txt does not describe a file present in this checkout; the embedded table and the optional --flow_field loader are the operative sources. [load_default_flow_field and file loader](../src/io.h)

### Combining both directions

The reverse pass is a filter from the other end of the sequence, rather than an unnormalized backward likelihood. At a site, the combined posterior is proportional to

\[
p(t\mid \mathrm{all\ data})\ \propto\
\frac{q_F(t)q_R(t)}{\pi(t)},
\]

where the two messages cover nonoverlapping observations and pi is the marginal TMRCA prior. For constant Ne, pi(t)=exp(-t). Products and this prior division leave a gamma density:

\[
\alpha_{FR}=\alpha_F+\alpha_R-1,\qquad
\beta_{FR}=\beta_F+\beta_R-1.
\]

Both ends initialize at Gamma(1,1), implemented by zeroing log10(mean) and log10(CV). The minus-one arithmetic is explicit in mean_cv_to_alpha_beta_log10_vec_backward. This derivation depends on the **exponential prior**, so an alternative transition table cannot by itself make the whole decoder demographic. [Current implementation](../src/gamma_smc.h)

The original supplemental methods describe the projection and reverse-filter construction in sections C-D. They also make explicit the SMC-prime kernel underlying the table. [Original supplemental methods, author repository](https://api.repository.cam.ac.uk/server/api/core/bitstreams/48ddd6e8-c82a-4e78-a4a2-3f62bcbf28f1/content)

## What changes when population size varies

Keep the same fixed N0 for units and define

\[
\lambda(t)=\frac{N_0}{N_e(2N_0t)},\qquad
H(t)=\int_0^t\lambda(u)\,du,\qquad
\pi_N(t)=\lambda(t)e^{-H(t)}.
\]

Lambda is the pair coalescence hazard in scaled physical time. For positive piecewise-constant sizes with a finite positive ancient size, pi_N is normalized and the ancient tail is exponential. Lambda=1 recovers the existing Exp(1) prior. Variable population size changes both the marginal prior and the chances that a lineage pruned by recombination rejoins at different ages. A scalar Ne evaluated at the posterior mean cannot represent the relevant integrals over all earlier ages. [Variable population-size coalescent construction](https://pmc.ncbi.nlm.nih.gov/articles/PMC3697970/)

The relevant pairwise SMC-prime transition is already available analytically as integrals for arbitrary demographic histories; Carmi et al. give it in equation 57. The following is an implementation-oriented reduction of that kernel, derived here. [Carmi et al., 2014, section 5.2](https://arxiv.org/html/1403.1325#S5.SS2)

Define

\[
J(t)=\int_0^t e^{-2[H(t)-H(u)]}\,du.
\]

It can be computed stably from J'(t)=1-2*lambda(t)*J(t), J(0)=0. Within an epoch of constant lambda=k and width d,

\[
J(t+d)=J(t)e^{-2kd}+\frac{1-e^{-2kd}}{2k}.
\]

For a current TMRCA s, the density of a changed TMRCA t conditional on a recombination event is

\[
q_{\rm change}(t\mid s)=\frac{\lambda(t)}{s}J(\min(s,t))
e^{-[H(t)-H(s)]_+},\quad t\ne s.
\]

The remaining conditional mass is an invisible recombination that leaves TMRCA unchanged. For a total recombination probability rho*s per base, the infinitesimal off-diagonal kernel is s*q_change. To reproduce the existing local canonical generator, define instead

\[
K_N(s,t)=2\lambda(t)J(\min(s,t))e^{-[H(t)-H(s)]_+}.
\]

The factor two is intentional here: it matches the existing source, and is an audit issue discussed below. The total off-diagonal rate of this K_N is R_N(s)=s+J(s). With the standard rho*s convention, both K_N and R_N would be halved.

For any density f, its change under the source-compatible generator can be evaluated using only cumulative integrals:

\[
\begin{aligned}
(A_Nf)(t)={}&2\lambda(t)e^{-H(t)}\int_0^t e^{H(s)}J(s)f(s)\,ds\\
&+2\lambda(t)J(t)\int_t^\infty f(s)\,ds
-[t+J(t)]f(t).
\end{aligned}
\]

Thus table construction need not multiply a dense time-by-time transition matrix at every parameter point. After preparing H and J, a forward cumulative integral and a tail integral suffice on a quadrature grid. Production code should use scaled accumulators to avoid overflow in e^H, and split quadrature at demographic boundaries.

There are two useful analytic checks:

- **Constant-size reduction:** lambda=1 gives J(t)=(1-exp(-2t))/2. Substitution gives exactly the three terms in distribution_difference_pdf_gsl, including its diagonal coefficient exp(-2t)/2-0.5-t.
- **Detailed balance:** pi_N(s)*K_N(s,t)=pi_N(t)*K_N(t,s). Consequently A_N*pi_N=0, and the same transition kernel can be used in both genomic directions with the correct prior.

These checks identify the required transition law and stationary prior; they do not validate gamma projection, cached interpolation, or statistical power.

## Three implementation routes

| Route | What remains simple | Main difficulty | Assessment |
|---|---|---|---|
| Plain gamma in physical time | Current mutation update; two-parameter transition table | pi_N is generally not gamma; dividing two gamma messages by pi_N is not gamma; initialization and output projection need new approximations | Feasible prototype, but merely swapping the table is inconsistent |
| Gamma in cumulative-hazard time u=H(t) | Prior becomes Exp(1) | Mutations and recombination depend on physical branch time H^-1(u), so the original emission and transition formulas no longer apply | No shortcut; trades one nonconjugacy for another |
| A two-parameter tilt of pi_N in physical time | Exact prior, exact emission update, exact prior-corrected message combination within the family | New normalizers, moments, CDFs, projection and caches; still only two parameters per message | Preferred research prototype |

These assessments follow from the algebra above and below, rather than an existing implemented Gamma-SMC feature.

### Preferred prototype: preserve the demographic prior explicitly

Use the family

\[
q_N(t;a,b)=\frac{1}{Z_N(a,b)}\pi_N(t)t^a e^{-bt},\qquad
Z_N(a,b)=\int_0^\infty\pi_N(t)t^a e^{-bt}\,dt.
\]

This is an exponential family with two sufficient statistics, log(t) and t. It requires only two message parameters even when the base density contains many demographic epochs.

- The unobserved prior is exactly a=0, b=0.
- The Poisson emission update is exactly a += y, b += theta.
- The two filters combine as a_FR=a_F+a_R and b_FR=b_F+b_R, because dividing q_F*q_R by pi_N cancels one copy of the base density.
- With lambda=1, this is Gamma(alpha=a+1, beta=b+1), so the current algebra is recovered exactly.

Recombination does not preserve this family, just as it does not preserve the ordinary gamma family. It still requires projection. Its tangent directions are

\[
\partial_a q_N=q_N[\log t-E_q\log t],\qquad
\partial_b q_N=q_N[E_qt-t].
\]

These can be substituted into the same two-column least-squares construction with target A_N*q_N. An alternative is an information projection matching E[t] and E[log t]; that would be a deliberate method change requiring its own benchmark, not the original gamma least-squares rule.

For a piecewise-constant history, Z and truncated moments are sums of one-dimensional epoch integrals of powers times exponentials. They can be calculated offline and cached; CDFs at the six existing time cutoffs can also be cached. Some finite epochs can have b+lambda_k <= 0, which is not automatically invalid because their intervals are bounded. The infinite tail requires b+lambda_oldest > 0. With finite positive lambda at the present, a > -1 ensures integrability at zero. These constraints and stability near their boundaries must be explicit.

No online demography estimation is necessary. A fixed EAS history would produce a reusable demographic table. The number of epochs increases construction and validation cost; it need not add an epoch loop to every pair/site update once the maps are cached. Whether accuracy at practical grid sizes and speed remain competitive is untested.

There is a subtle limitation: closure of multiplication is algebraic, but two arbitrary individually integrable approximate messages can in principle produce a nonintegrable combined tilt. Admissible state ranges, the sign/size of b, and boundary behavior must therefore be checked during the prototype, rather than assuming every interpolated state is valid.

### Why cumulative-hazard rescaling is not enough

If u=H(t), then the prior on u is Exp(1). However, a match emits exp[-theta*H^-1(u)], and a mismatch emits H^-1(u)*exp[-theta*H^-1(u)] up to constants. These are not the gamma-conjugate factors exp(-theta*u) and u*exp(-theta*u), unless H is linear. Recombination pruning is uniform in physical branch time, not in u; its density gains a factor dt/du. Simply transforming the age axis while keeping the prebuilt field and mutation update would implement the wrong model. This is a direct change-of-variables calculation.

## Source-audit issues to resolve before a new demographic decoder

**Recombination factor of two.** With t=g/(2N0) and rho=4N0*r, total two-lineage branch length gives a per-base recombination probability approximately rho*t. The original supplement states this convention, but its section E replaces it by 2*rho*t during the flow derivation. The local generator follows that latter expression, and FlowFieldCache multiplies the resulting arrows directly by the supplied rho. This is a source/documentation inconsistency with a potential physical scaling consequence. It must be checked against a small independent SMC-prime reference before choosing the convention for a new decoder. It is not silently corrected here, and it has not been shown to explain the selected-versus-neutral performance gap. [Supplement sections B and E](https://api.repository.cam.ac.uk/server/api/core/bitstreams/48ddd6e8-c82a-4e78-a4a2-3f62bcbf28f1/content), [generator](../src/generate_canonical_flow_field.cpp), [cache construction](../src/flow_field.h)

The convention was checked beyond the generator formula:

- The supplement defines Ne as diploid size, scaled time as g/(2Ne), and rho as 4Ne*r on its printed page 22. Printed page 24 states rho*s; printed page 29 states it again and then introduces 2*rho*s in the next displayed transition expression. This is not a switch from haploid to diploid population-size units.
- In the current experiment, theta=0.00075 and rho/theta=0.8 reach the command-line parser unchanged, producing rho=0.0006. The same parser converts times using 2Ne=theta/(2mu)=30,000 generations. No compensating rate factor is present there. [Wrapper](../python/gamma_smc_aou/decoder.py), [rate and threshold construction](../src/gamma_smc.cpp)
- Five nontrivial embedded-grid entries were recomputed independently with SciPy using the source projection formula and its default 1,000-step quadrature settings. Embedded/source velocity ratios ranged from 0.9999991 to 1.0000181, consistent with the stored precision and numerical integration; they were approximately one, not one half. Thus the embedded table has not globally halved the generator. The probed zero-based (mean,CV) indices were (25,49), (30,42), (35,35), (40,49), and (40,42). This is a convention check, not full-grid accuracy validation. [Embedded table](../src/io.h)
- The cache constructor keeps rho unchanged; preprocess adds rho times each velocity once per base. Segment counts are physical base differences. Flattening copies the cached values, and the SIMD path takes ordinary bilinear weighted sums. The variable name `halfdots` refers to the two 128-bit halves of an AVX register and introduces no factor of one half. Likewise, the -0.5 factor converting log(alpha) velocity to log(CV) is the coordinate derivative, not a rate correction. [Cache/runtime code](../src/flow_field.h), [segment construction](../src/data_processor.h)

These checks support a factor-of-two discrepancy in the current source/table rate convention; they do not quantify its effect on the approximate final posterior or on calibrated scan performance. Resolving the convention in a separate controlled benchmark is preferable to changing rates in the ongoing mean/median/mass experiment.

**Entropy clipping is an inherited heuristic.** The existing cache clips gamma differential entropy at 1, the entropy of Exp(1), when enabled. This constant belongs to the current units and prior. Moreover, the general inequality that conditioning reduces entropy is an average over observations; it does not require every realized posterior to have lower entropy than its prior. A direct counterexample within the current emission model is one heterozygous observation with no recombination: Exp(1) updates to Gamma(2,1+theta), whose entropy is 1+EulerGamma-log(1+theta), approximately 1.576 at theta=0.00075. It is a valid exact posterior with entropy greater than 1. A new implementation should test positivity, normalization, projection residuals and boundary error directly; it should not carry over this pointwise bound as a mathematical necessity. [Existing entropy-clipping implementation](../src/flow_field.h)

The actual call order confirms that this heuristic can act on emissions, not just transition approximations. main constructs FlowFieldCache with entropy clipping explicitly true. In the forward cache, preprocess adds theta to beta, clips the homozygous candidate, then adds one to alpha and clips the heterozygous candidate before storing it. In the reverse cache, it constructs and clips the emission candidate before applying the cached stretch. The clipping routine preserves mean while lowering CV to its entropy boundary; the comment saying to increase CV is inconsistent with the operative table and with the supplement's description. The Gamma(2,1+theta) calculation is an exact-model counterexample, not a measured claim that the interpolated production state equals that gamma. Its intermediate mean is about two, within the active entropy-clipping branch. No clipping rates have been measured on the selected simulations, and no resulting power loss is asserted. [Forward and reverse cache assembly](../src/flow_field.h), [enabled constructor argument](../src/gamma_smc.cpp)

**The CV ceiling is not demographic.** A gamma CV limit of 1 is a chosen family/grid restriction. A valid variable-size prior can have CV above 1, even before seeing data. The illustrative numerical history below has prior CV=1.13694. Keeping the existing mean/CV bounds would exclude that exact prior. New coordinates and grids require explicit investigation.

**Caching and units need stronger metadata.** A demographic field must record the history, interpolation convention, units, N0, generator normalization, prior/family definition, projection method, quadrature bounds, grid limits, tolerances, and hashes. Derived skip caches also depend on theta, rho, clipping and implementation parameters. Output posterior metadata must distinguish ordinary gamma parameters from tilt parameters; existing gamma readers cannot interpret them unchanged. Cache reuse should compare these scientific inputs and relevant outputs, not only a git commit.

**A scalar Ne(t) does not represent the entire introgression model.** The EAS simulation contains lineage population states and an archaic pulse. A panmictic time-dependent coalescence hazard can represent a size history, or match a marginal TMRCA distribution, without reproducing the two-locus transition dependence of a structured population. Additional lineage-state information would be required for a faithful structured-coalescent decoder. The proposed extension is therefore a demographic approximation, while neutral/selected simulation calibration remains useful.

## Numerical feasibility checks completed

A deterministic, small quadrature calculation is saved as [check_variable_ne_flow.py](../scripts/check_variable_ne_flow.py), with machine-readable [results](results/variable_ne_flow_checks.json). It does not read simulation data, change decoding, build the production table, or use random numbers. It checks one constant history and an illustrative four-epoch history with starts [0,0.05,0.15,0.4] and hazards [0.25,4,0.5,1] in scaled time. The latter is a mathematical example, not an estimated EAS history.

The following checks passed with absolute errors below 2e-15, against a required tolerance of 1e-8:

- Prior normalization; J from the epoch recurrence versus direct integration.
- Integrated off-diagonal rates versus s+J(s).
- Detailed balance and stationarity of pi_N under A_N.
- Constant-history kernel and generator versus the expression in the existing C++ source.
- Constant-history tilted density versus scipy's ordinary gamma density.
- Variable-history emission and prior-corrected filter combination versus separately normalized products.

These checks verify the formulas at the evaluated points. They do not establish numerical convergence of a full flow grid, correctness of a decoder, or improved power.

## A bounded implementation and validation plan

1. **Reference model:** Implement a small, accurate discretized SMC-prime decoder for constant and piecewise-constant size histories, with physical-rate units explicit. Resolve the factor-of-two convention using the transition kernel and no-emission segment behavior.
2. **Demographic state family:** Implement pi_N, the tilted-family normalizer, mean, CDF and quantiles, and test analytic constant-size reduction, emission/merge closure, and tail integrability. Preserve the existing definitions of frac_recent_T as fractions of pairs passing mean, median or a posterior-mass criterion.
3. **Offline construction:** Generate a modest two-dimensional demographic table, including posterior states near recent and ancient epochs. Compare projected transitions with direct numerical transitions, refine the time grid and parameter grid, and report projection error separately from interpolation error.
4. **Forward/backward and caches:** Compare uncached steps, cached steps, and reverse-filter combination against the reference. Include zero observations, missing intervals, a single mismatch, long matches, observations around every demographic boundary, and exact output-coordinate alignment.
5. **Simulation comparison:** Hold the simulated genealogies, sequences, haplotypes, pair manifest, mu, r and age cutoffs fixed. Compare constant and demographic decoding on paired inputs. Evaluate pair-age errors and posterior calibration, then held-out selected-region versus neutral-region call fractions with identical complete region searches. Any detector choice must be trained separately from final evaluation.

This separates a potentially useful decoder improvement from changes to the biological experiment. A variable-Ne implementation should only become the production decoder after both the reference checks and the held-out regional detection comparison support it.

To rerun the standalone feasibility check in the repository's pinned uv environment:

```bash
%%bash
if [ ! -d gamma_smc_ts/.git ]; then
  git clone --branch AOU_run_opt git@github.com:soisa001/gamma_smc_ts.git
fi
cd gamma_smc_ts
git pull --ff-only
uv run python scripts/check_variable_ne_flow.py \
  --output docs/results/variable_ne_flow_checks.json
```

Source audit reference: repository commit 916a829a513bee6c09e20c1cf649d7d9cb61c503. Primary article and supplemental methods were read from the author-hosted Cambridge copies when publisher/PMC full-page access failed. The original paper's general motivation, the supplement's projection recipe, and Carmi's demographic kernel are distinguished above from the new algebra, proposed family, numerical checks, and implementation recommendations.
