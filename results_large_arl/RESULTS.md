# Large-ARL numerical calibration and detection experiment

Source: [arXiv:2503.12276v2, Figure 3 and Appendix F](https://arxiv.org/html/2503.12276v2#S0.F3).

Target ARL=2e+06 pulses, N=100, eta1=0.9, eta2=0.85. Na grid=[0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0]. Equal total transmitted energy N+Na within each pair. Change at pulse 1000, horizon 5000; change at pulse 1 is an additional initialization check.

## What was verified

- All 24 comparisons with the prior moderate-ARL Monte Carlo validation passed the declared criterion: relative deviation <= max(2%, 3 Monte Carlo SE / numerical ARL). These were previously saved independent no-change ensembles at target ARLs 500, 1000 and 5000, not new two-million-pulse ensembles.
- High-ARL thresholds were calibrated on 1024 positive cells, then evaluated at the SAME thresholds on half and twice as many cells. Every final base-to-fine change and fine-grid deviation from target was <= 0.1%. Grid refinement is an empirical numerical check, not a rigorous error bound.
- Raw changed-run stopping times/event flags were saved. Censored or unavailable unrestricted delay estimates cause the experiment to stop, rather than silently discarding runs.

| Na | Fine-grid ARL ordinary / squeezed | Ordinary / squeezed delay (nu=1000) | Delay ratio (approx. 95% interval) | KL ratio | Early alarms ordinary / squeezed |
|---|---:|---:|---:|---:|---:|
| 0 | 2000648 / 2000648 | 76.80 / 76.78 | 1.000 [0.995, 1.005] | 1.000 | 0.030% / 0.030% |
| 0.02 | 2000648 / 2000524 | 76.35 / 60.83 | 1.255 [1.249, 1.262] | 1.284 | 0.094% / 0.056% |
| 0.05 | 2000648 / 2000465 | 76.65 / 53.39 | 1.436 [1.428, 1.443] | 1.476 | 0.054% / 0.044% |
| 0.1 | 2000648 / 2000408 | 76.90 / 46.39 | 1.658 [1.649, 1.666] | 1.716 | 0.056% / 0.048% |
| 0.25 | 2000647 / 2000319 | 76.39 / 36.10 | 2.116 [2.105, 2.127] | 2.256 | 0.050% / 0.046% |
| 0.5 | 2000646 / 2000251 | 76.29 / 28.02 | 2.723 [2.709, 2.737] | 2.935 | 0.060% / 0.042% |
| 1 | 2000643 / 2000190 | 76.04 / 21.12 | 3.601 [3.582, 3.619] | 3.933 | 0.058% / 0.034% |

Largest base-to-fine ARL change: 0.03241%.

## Method and interpretation

For Page statistic w, let g(w) be expected remaining pulses to the first strict threshold crossing and F the pre-change LLR CDF. First-step conditioning gives g(w)=1+F(-w)g(0)+integral_(0,h] g(y)dF(y-w). We retain the point mass at zero and integrate exact Gaussian/quadratic-Gaussian LLR probabilities across each positive cell, approximating g by its midpoint value. A linear solve gives g(0). Threshold search sets this numerical estimate to the desired ARL. The CDF handles unequal variances, so squeezed LLRs retain their quadratic term.

This differs from Appendix F’s shared-trajectory Monte Carlo calibration. It is a sampled, Figure-3-style numerical reproduction with an explicitly different ARL method. The threshold results have been checked against moderate-ARL simulations and mesh refinement; the actual rare-event ARLs have NOT been independently validated by no-change Monte Carlo at two million pulses.

Reported detection delays are E[T-nu | T>=nu], not worst-case delays. False alarms before the change are reported separately; no restart is performed. nu=1 consumes T pulses but has paper-convention delay T-1. Zero squeezing uses a common threshold for identical models and independent observation streams. Uncertainty intervals are pointwise approximate Monte Carlo intervals conditional on the selected thresholds and exclude discretization/model uncertainty. No multiple-comparison correction is applied.

Finite-threshold delay gains need not equal KL ratios. A closer approach to the diagonal at large ARL supports the paper’s asymptotic relationship but does not prove the limit or establish identical numerical results to its published scatter plot. The experiment does not include entangled probes.

## Evidence and next work

report.json contains settings, versions, seeds, crosschecks, threshold-search histories and all grid diagnostics. Per-condition report.json and raw_runs.npz hold the corresponding delay ensembles. summary.csv includes both change times and uncertainties. The original moderate-ARL evidence remains unchanged.

Next: an independent rare-event ARL check (for example a regenerative estimator) and the correlated n-mode entanglement model. All numerical ARL claims should retain their discretization qualification until independently verified.
