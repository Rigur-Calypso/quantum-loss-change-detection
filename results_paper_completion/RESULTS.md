# Paper-specific implementation additions

Source: [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).
The previous robustness extensions and original source arrays remain unchanged.

- BPSK Gaussian-mixture densities, stable exact LLRs, sampling, KL and uniform-input
  mutual information are implemented. Unknown signs are marginalized, not supplied
  by an oracle. Physical raw homodyne units are used consistently before/after loss.
- Figure 3 now includes 16 new BPSK runs on the authors' eight energy points.
  Source scalar thresholds are fixed and newly checked using two independent
  250,000-cycle ensembles per law. No threshold is fitted to
  latency. All non-Gaussian detection checks pass: **True**.
  Mean absolute relative difference from the eight published BPSK delay ratios:
  **6.40%**. Author errors and thresholds are unknown, so disagreement
  remains descriptive; this does not establish exact published-latency reproduction.
- Figure 4's fixed-squeezer analytical family is generated for n=2 through 256.
  The chosen seed range is 0–20 dB; these are formula-derived curves, not digitized data.
- Eq.22 one-photon homodyne densities, Eq.23 CDFs and exact mixture sampling are
  implemented. The sampler mixes vacuum and a signed chi-square(3) radius; it does
  not need numerical inverse-CDF interpolation. Independent tests check its CDF.
- The DV benchmark KL is **0.206340 nats**, matching the paper's rounded
  0.2063 value. Its squeezing crossover is **1.812 dB**
  at eta1=0.9/eta2=0.85. The independently calibrated DV mean conditional delay is
  **54.640 ± 0.129 pulses**
  (±1 Monte Carlo SE), with source-style change at 1000 and horizon 5000.
- Figure 11(a) and Figure 12 model families are evaluated. Figure 11(a) uses the
  stated rounded tap=0.944; the separate benchmark uses exactly 0.9/0.85. The weak
  entangled example uses n=8; Kennedy residuals use PNR Eq.7 and an explicit residual
  amplitude convention. These choices are saved rather than silently equated.
- Appendix A Eqs.10–11 phase-mixture homodyne and PNR Kennedy receivers are
  implemented, including phase quadrature and thermal occupation in Eq.11. The
  saved LO sweep uses chosen N=100, Na=1, tap=0.9, fractional amplitude error=0.01,
  and phase SDs 0,0.01,0.03. It is not an exact Figure 5 parameter reproduction.
  Zero-noise limits and quadrature refinement are tested independently.
- Figure 11(b) now has a model-comparison plot using the earlier unmodulated
  simulations, the new BPSK simulations and independently calibrated DV result.
  The author DV reference is the plotting script's signed mean; its two negative
  latencies are retained. It is explicitly distinguished from our conditional
  first-alarm delay. Exact numerical agreement is not claimed.
- Appendix C's Hadamard-network covariance is propagated independently and
  agrees with the analytic entangled covariance to machine precision. Appendix
  B's energy derivative is evaluated by differentiating Eq.3 directly, and its
  high-energy squeezing/classical crossover is computed without dropping the
  small variance-entropy term. The approximate Appendix B.3 threshold is retained
  separately. Appendix F's shared-path multiple-threshold first-alarm estimator
  is implemented and tested, including same-pulse overshoots and censoring.
- Bernoulli photon counting and Poisson likelihoods support the discrete receiver
  models. Ideal Kennedy nulling has infinite KL and infinite no-change ARL; its
  first post-change click still has a geometric waiting time, not zero expected
  latency at finite pulse energy.

The weak-probe and receiver-noise curves are formula/model evaluations with
explicit conventions and chosen grids. Their qualitative ranking is not imposed
in code; it can depend on loss range, residual-error convention and PNR versus
on/off observations. They are not presented as exact digitizations of Figures
5 or 12. The open coding/communications tradeoff problems are not solved here.

Every changed simulation stops at its first strict Page crossing; early alarms
and censoring are explicit. Independent moderate-ARL direct controls check the
generic non-Gaussian renewal engine. Successful iid weights are bounded by exp(-h),
with overshoot retained. No iid engine is applied to entangled correlated increments.

See the root PAPER_IMPLEMENTATION.md for complete figure-by-figure coverage and
remaining reproduction work. Full optical hardware and capacity-achieving codes
are schematics/discussion in the source, not claimed software implementations.
