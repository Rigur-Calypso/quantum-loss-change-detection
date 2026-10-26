# Findings and limits

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
   that exploratory grid. [Derivation and evidence](results_adaptive_validation/RESULTS.md).
2. **Sensitivity to the initial channel.** With a detector assuming η₁=0.9,
   changing actual transmission from 0.899 to 0.901 changes classical ARL from
   about 1.31 million to 3.01 million, against a two-million target. Initial-channel
   calibration matters even over this narrow interval. Scalar calculations use
   independent grid refinement and actual-null path weights; correlated blocks
   have a separate validation. [Scalar study](results_null_channel_refined/RESULTS.md),
   [correlated study](results_entangled_null_diagnostics/RESULTS.md).
3. **Conservative thresholds under channel uncertainty.** Frozen thresholds
   pass independent validation on 25 sampled detector/channel settings.
   Adjusted approximate lower ARLs exceed two million, at delay costs of about
   3.61–4.39 pulses for the classical detector and 1.12–2.33 for squeezed/entangled
   examples. These are sampled-grid Monte Carlo results, not a guarantee over
   a continuous interval. [Full protocol](results_robust_thresholds/RESULTS.md).
4. **Explicit ideal dual-rail channel and timing.** Arbitrary qubit states,
   Kraus loss maps, heralded survival, and changes between the two time bins are
   simulated. Uniform loss preserves a surviving qubit; unequal rail loss can
   distort it. Two independent holdouts support the two-million-physical-bin ARL
   target. This develops the Appendix E.2 logical model; it does not simulate
   cavity hardware or prove QKD security. [Evidence](results_dual_rail/RESULTS.md).
5. **Statistical/provenance checks and teaching material.** Saved first alarms,
   explicit censoring, independent seeds, checksum inventories, portable evidence
   and an executable lesson make the results inspectable. The DV re-analysis
   distinguishes signed author latencies from a conditional first-alarm mean.
   [DV uncertainty study](results_dv_uncertainty_final/RESULTS.md).

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
questions are also documented in [specific clarification questions](source_audit/CLARIFICATION.md).
