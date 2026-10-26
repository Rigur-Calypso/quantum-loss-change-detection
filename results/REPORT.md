# Combined results and implementation report

This project implements the paper's principal computational models, validates
their numerical behavior, and investigates robustness beyond the plotted
benchmarks. It does not claim full published-curve agreement, new quantum
hardware, a new quantum advantage theorem or a literature-established novelty.

## Reproduced or independently checked

| Item | Evidence | Limit |
| --- | --- | --- |
| Figure 2 | Coherent/squeezed and n=2–256 entangled relative-entropy curves; Gaussian KL and dense covariance crosschecks | Formula reproduction; original numerical figure array unavailable |
| Detection model | Page reflection, strict first crossing, conditional likelihoods, mixed-loss blocks and one-based timing | Known model assumptions; no claim of a measured optical implementation |
| ARL calibration | 61/61 original author-grid settings internally confirmed near two million pulses | Approximate Monte Carlo evidence, not exact author thresholds |
| Eq.1 trend | 14 scalar settings, target ARL 100–100 million; high-ARL normalized slopes within 0.014% of one | Finite cold-start delays have offsets; not a worst-case theorem |
| Single-photon entropy | 0.206340 nats vs quoted 0.2063; exact sampler and independent integration | Numerical model validation |
| Figure 12 | Seven sparse comparisons at 1 dB within 2% | Digitization uncertainty and sparse coverage; no original full arrays |

## Additional work relative to the paper's reported benchmarks

These are project extensions and implementation/validation contributions.
They should not be described as novel methods without a broader literature review.

1. **Independent rare-event ARL validation.** Regenerative cycles and a
   history-dependent mixture proposal make extremely rare false alarms
   computationally accessible. The full mixture density corrects each stopped
   path. Direct moderate-ARL checks, independent holdouts, effective sample size,
   weight concentration and repeatability diagnostics are retained. The adaptive
   follow-ups confirm both previously unresolved settings, supporting 25/25 on
   that exploratory grid. Derivation and evidence.
2. **Sensitivity to the initial channel.** With a detector assuming η₁=0.9,
   changing actual transmission from 0.899 to 0.901 changes classical ARL from
   about 1.31 million to 3.01 million, against a two-million target. Initial-channel
   calibration matters even over this narrow interval. Scalar calculations use
   independent grid refinement and actual-null path weights; correlated blocks
   have a separate validation. Scalar study,
   correlated study.
3. **Conservative thresholds under channel uncertainty.** Frozen thresholds
   pass independent validation on 25 sampled detector/channel settings.
   Adjusted approximate lower ARLs exceed two million, at delay costs of about
   3.61–4.39 pulses for the classical detector and 1.12–2.33 for squeezed/entangled
   examples. These are sampled-grid Monte Carlo results, not a guarantee over
   a continuous interval. Full protocol.
4. **Explicit ideal dual-rail channel and timing.** Arbitrary qubit states,
   Kraus loss maps, heralded survival, and changes between the two time bins are
   simulated. Uniform loss preserves a surviving qubit; unequal rail loss can
   distort it. Two independent holdouts support the two-million-physical-bin ARL
   target. This develops the Appendix E.2 logical model; it does not simulate
   cavity hardware or prove QKD security. Evidence.
5. **Statistical/provenance checks and teaching material.** Saved first alarms,
   explicit censoring, independent seeds, checksum inventories, portable evidence
   and an executable lesson make the results inspectable. The DV re-analysis
   distinguishes signed author latencies from a conditional first-alarm mean.
   DV uncertainty study.

## Differences that are not reproduced results

| Comparison | Our observed difference | Interpretation |
| --- | --- | --- |
| Figure 3 | 7.49% unmodulated / 6.40% BPSK mean absolute relative delay-ratio difference | Analytical KL coordinates agree; original thresholds and raw latency records unavailable |
| Figure 7 | Mean delays about 13.6–34.4% higher by family | Internal ARL calibration passes; published latency agreement incomplete |
| Figure 4 | Fixed-seed equations imply n=1 ≥ n=2 ≥ n=4, unlike the published inset | Dense Gaussian calculations and the A/n+B identity agree; original plotting calculation needed |
| Figure 5 | PNR Kennedy KL about 11.9–12.6 times four digitized values | Verified legend parameters, independent quadrature and samples support our calculation; cause unresolved |
| Figure 11(b) | Our conditional DV delay 54.640 vs source signed mean 48.034 | Source includes two negative latencies; nonnegative re-summary is 50.181, with event-interpretation assumptions |

The DV conditional comparison leaves a 4.459-pulse gap; its empirical bootstrap
interval is approximately [2.665, 6.206] under independence and missing-event
assumptions. This does not identify a source error. A roughly ten-photon Kennedy
counterfactual and a fitted-loss diagnostic are explicitly exploratory; neither
changes the paper's stated model or counts as a successful reproduction.

Hardware schematics, an optical experiment, capacity-achieving communication
codes, cavity dynamics and the paper's open optimal communication/detection
tradeoff are outside this computational release. The remaining numerical
questions are also documented in specific clarification questions.

## Figure-by-figure implementation

| Paper item | Implemented computation | Evidence and remaining limits |
| --- | --- | --- |
| Figure 1: transmitter families | Coherent, independent squeezing and correlated entangled probe statistics | `models.py`, `entangled.py`; optical-network covariance independently checked. The drawing represents hardware, not a software benchmark. |
| Figure 2: per-pulse relative entropy | All coherent, squeezed and n=2,4,…,256 entangled curves | Analytical figure and data. Formula reproduction with independent Gaussian/covariance checks; not pixel digitization. |
| Figure 3: detection improvement, with/without BPSK | Both unmodulated and BPSK simulations on all eight saved author grid points, with ARL near two million | Author comparison, BPSK addition. KL coordinates agree to machine precision. Mean absolute relative delay-ratio differences: 7.49% unmodulated, 6.40% BPSK. Exact latency agreement remains unresolved. |
| Figure 4: fixed seed squeezing | Formula-derived n=1,2,…,256 family, with Na=sinh²(s)/n, over the verified 0–18 dB axis | Equation check. Independent dense Gaussian KL and A/n+B identity agree. At fixed seed, equations give n=1 ≥ n=2 ≥ n=4; the published inset places n=2 above n=1. This discrepancy is unresolved. Earlier n=2-only ranking excluded n=1 and does not establish published curve agreement. |
| Figure 5: imperfect local oscillators | Eqs.7,10,11: phase-mixed homodyne and PNR Kennedy laws, plus on/off coarse-graining | Verified-setting comparison. PDF confirms N=100, Na=1, tap=0.9, epsilon=0.01, phase variance=0.01 (SD=0.1), loss 0–6 dB. PNR is explicitly specified by Appendix A Eq.10. Calculated Kennedy entropy is approximately 11.9–12.6 times the plotted values at four digitized points. Binary coarse-graining does not close the gap. |
| Figure 6: receiver schematic | Statistical models for homodyne, PNR and on/off receivers | `receivers.py`, `paper_models.py`. No optical/electronic hardware implementation claimed. |
| Figure 7: correlated detection latency | Conditional-LLR Page detector for all five detector types and nine saved author energies; exact mixed-loss blocks and pulse-level first alarms | Full comparison. All 45 ARL checks support the internal model, but mean delays differ from author data by approximately 13.6–34.4% across families. Detector phase, conventions and threshold diagnostics remain preserved. |
| Figures 8–10: construction and entanglement network | Displaced squeezed states, BPSK mixtures, normalized Hadamard network and resulting covariance | Covariance propagation agrees with the analytic model within 2.8×10⁻¹⁷. A capacity-achieving coded transceiver or physical Green Machine is not supplied by this simulation. |
| Figure 11(a): discrete-variable augmentation | Eq.22 density, Eq.23 CDF, exact vacuum/Fock-1 mixture sampler, KL integration, CV/DV comparison | Plot. Uses stated rounded tap=0.944; the separate benchmark uses exactly eta1=0.9, eta2=0.85. Computed DV KL 0.206340 matches quoted 0.2063. |
| Figure 11(b): CV/DV latency | Independently calibrated DV detector, existing CV simulations and new BPSK simulations; all retain first-alarm/censoring information | [Comparison](figure11.png). Our DV conditional delay is 54.640±0.129 pulses (one SE); the source plotting script uses a signed mean containing two negative latencies. Those statistics are labeled separately. Published latency matching is incomplete. |
| Figure 12: weak one-photon probes | Coherent homodyne, squeezed vacuum, entangled vacuum, Fock homodyne, Bernoulli photon counting, residual-nulling Kennedy models | [Verified-setting curves](figure12.png). PDF verifies N=0, Na=1, n=8, tap=0.9, residual photons 0.0001/0.00001 and loss axis 0–5 dB. Seven approximate raster readings at 1 dB agree within 2%; this is a sparse pixel check, not full author-data agreement. Zero-loss photon-counting KL is infinite and recorded separately. |


## Algorithms and appendices

- **Page detection and change-point estimate:** strict one-based first alarm,
  reflected evidence, and the main-text argmin-prefix estimate of the change
  point. Ties use the earliest prefix. No estimate is fabricated without an alarm.
- **BPSK communication/detection:** unknown signs are marginalized in the actual
  Gaussian-mixture likelihood. No transmitted sign or unknown true channel is
  provided to the detector. Uniform-input BPSK mutual information is calculated
  in bits/pulse. The illustrative capacity sweep includes low-energy cases beyond
  the paper's bright-pulse regime; it is not a capacity-achieving coding simulation.
- **Appendix A:** phase Gaussian truncated to [-pi,pi], normalized quadrature,
  amplitude/phase-sensitive PNR nulling, homodyne squeezed/antisqueezed quadratures
  and optional thermal occupation. The strong-LO homodyne model has no fractional
  amplitude-error term. Absolute residual amplitude (Eq.7) and fractional LO
  amplitude error (Eq.10) are kept distinct.
- **Appendix B:** energy derivatives computed directly from Eq.3; exact positive
  high-energy squeezing/classical crossover 930.513 photons, versus the Appendix
  B.3 large-energy approximation 900 photons for N=100, eta1=0.9. Conditional
  block likelihoods follow Eqs.16–18; Page reflection never erases block history.
- **Appendix C:** normalized Hadamard mixing of one squeezed input plus vacuum
  inputs reproduces the output covariance; uniform loss is then applied. This is
  a statistical network calculation, not a hardware pulse-control simulator.
- **Appendix D:** exact non-Gaussian photon quadrature sampler, density, CDF and
  likelihoods. The entropy crossover with squeezing is r=0.208618, Na=0.044157,
  or 1.812 dB, close to the text's rounded 1.83 dB.
- **Appendix E:** Bernoulli survival/heralding and Poisson counting likelihoods,
  including boundary/infinite-KL cases. An ideal null has no false alarms, but
  still waits for a post-change photon click at finite energy. For a dual-rail
  qubit, one heralding observation covers two time bins; pulse-time conversion
  must account for that when using the Bernoulli API.
  The new explicit dual-rail model propagates
  arbitrary 2x2 qubit states into vacuum/early/late output states using three
  trace-preserving Kraus operators. Equal rail loss preserves the heralded
  state. The detector uses only total-survival heralds, with readout after the
  two-bin pair; a mixed-loss boundary pair has its actual state-dependent
  probability. Two independent ARL holdouts, 240,000 changed trajectories and
  1,500 state samples are retained. Ideal receiver phase correction is assumed;
  cavity dynamics, QKD security and optimal code tradeoffs are not simulated.
- **Appendix F:** shared, unreset null trajectories for several increasing alarm
  thresholds, including multiple crossings on the same overshooting pulse. The
  saved low-SNR BPSK demonstration agrees with independent renewal estimates at
  all three thresholds. Those shared-path ARL estimates are correlated. The exact
  Fock sampler replaces numerical inverse-CDF sampling while preserving Eq.23.


## Selected numerical comparisons

### Relative entropy, Figure 2

![Relative entropy, Figure 2](figure2.png)

### Delay improvement, Figure 3

![Delay improvement, Figure 3](figure3.png)

### Detection delay, Figure 7

![Detection delay, Figure 7](figure7.png)

### Discrete-variable comparison, Figure 11

![Discrete-variable comparison, Figure 11](figure11.png)

### Weak-probe comparison, Figure 12

![Weak-probe comparison, Figure 12](figure12.png)

### Independent discrepancy checks

![Independent discrepancy checks](figure5.png)

### Finite-ARL scaling

![Finite-ARL scaling](delay_scaling.png)

### Conservative thresholds and delay costs

![Conservative thresholds and delay costs](robustness.png)

### Ideal dual-rail channel

![Ideal dual-rail channel](dual_rail.png)

## Evidence and interpretation

All current and historical study records are packaged in `data/metadata.tar.gz`;
raw ensembles are the two GitHub release downloads. Run the README's evidence
command to materialize records under `data/studies`. Their original scientific
arrays and tables are preserved. Old failed runs and model diagnostics retain
their labels; they are not additional successful paper reproductions.

The main paper check and four extension auditors recalculate saved evidence.
Historical exact-regeneration checks are recorded separately. Delay means
E[T−ν | T≥ν]; early alarms and censoring are separate. Approximate Monte Carlo
intervals are not certified bounds or Lorden worst-case delay results.

The remaining author clarifications concern Figure 4's fixed-seed calculation,
Figure 5's probe-energy/PNR plotting calculation, and Figures 3/7/11 thresholds,
false-alarm calibration and first-alarm event conventions. The manuscript and
original figure archives are linked in the records, not redistributed here.

See [validation summary](validation.json), [numeric summary](summary.csv), and
[the repository README](../README.md) for commands and attribution.

## Regenerate computational experiments

Use fresh output folders. Regenerating the full ensembles can take much longer
than checking the saved evidence. Settings and seed records are retained in the
archived studies. The author dataset is parsed as data, never executed as MATLAB.

```sh
python -m lossqcd.author_comparison --output results/author_new
python -m lossqcd.paper_completion --source results/author_new --output results/paper_new
python -m lossqcd.audit_paper_completion results/paper_new --regenerate --regenerate-calibration
python -m lossqcd.verified_figures --output results/figures_new
python -m lossqcd.equation_validation --source results/figures_new --output results/equations_new
python -m lossqcd.audit_equation_validation results/equations_new --regenerate
python -m lossqcd.delay_scaling --output results/scaling_new
python -m lossqcd.audit_delay_scaling results/scaling_new --regenerate
python -m lossqcd.dual_rail_experiment --output results/dual_rail_new
python -m lossqcd.audit_dual_rail results/dual_rail_new --regenerate
```

The independent ancillary analyses are `python -m lossqcd.dv_uncertainty` and
`python -m lossqcd.source_audit`. The latter requires the original arXiv source
archives from the URLs/checksums in the restored source-audit records; those
manuscript archives are not redistributed. Hardware operation, cavity dynamics,
QKD security and the optimal communication/detection tradeoff remain outside
this computational implementation.
