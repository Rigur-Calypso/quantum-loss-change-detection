# Correlated conditional detector at large false-alarm ARL

Source model and detector: [paper Appendix B.4](https://arxiv.org/html/2503.12276v2#A2.SS4).
The regenerative extension and its validation here are our implementation work; they are not a method reported in the optical paper. General background: [Calvin, Glynn & Nakayama (2001)](https://web.stanford.edu/~glynn/papers/2001/CalvinGNakayama01.pdf).

N=100, Na=0.1, eta1=0.9, eta2=0.85. Target ARL=2,000,000 pulses. All checks passed: **True**.

| Detector | Threshold | ARL ± SE (pulses) | Matched | Repeat z | Moderate check z |
|---|---:|---:|---|---:|---:|
| Coherent | 11.942005 | 1,996,125 ± 3,451 | True | 0.13 | 1.51 |
| Squeezed n=1 | 12.253225 | 1,995,490 ± 2,965 | True | 0.89 | 1.95 |
| Entangled n=2 | 12.338530 | 1,999,570 ± 2,267 | True | 1.11 | 0.50 |
| Entangled n=4 | 12.411838 | 1,994,058 ± 3,768 | True | 0.99 | 1.05 |
| Entangled n=8 | 12.441420 | 1,993,296 ± 9,047 | True | 1.42 | 0.19 |

## Detection delays

Each run stops at its first alarm. Delay is E[T−nu | T>=nu]; early alarms are reported separately and never counted as detections. Pulse nu is the first changed pulse; block phase is (nu−1) mod n + 1. No changed-run delay in this report is censored.

| Detector | Change pulse | Phase | Delay ± SE | Early alarms | Coherent/detector delay |
|---|---:|---:|---:|---:|---:|
| Coherent | 1 | 1 | 79.721 ± 0.139 | 0.000% | 1.000 |
| Coherent | 1005 | 1 | 76.724 ± 0.140 | 0.052% | 1.000 |
| Squeezed n=1 | 1 | 1 | 48.242 ± 0.085 | 0.000% | 1.653 |
| Squeezed n=1 | 1005 | 1 | 46.548 ± 0.087 | 0.038% | 1.648 |
| Entangled n=2 | 1 | 1 | 40.158 ± 0.071 | 0.000% | 1.985 |
| Entangled n=2 | 1005 | 1 | 38.804 ± 0.072 | 0.026% | 1.977 |
| Entangled n=4 | 1 | 1 | 32.127 ± 0.057 | 0.000% | 2.481 |
| Entangled n=4 | 1005 | 1 | 31.050 ± 0.057 | 0.050% | 2.471 |
| Entangled n=8 | 1 | 1 | 25.401 ± 0.044 | 0.000% | 3.139 |
| Entangled n=8 | 1005 | 5 | 25.870 ± 0.045 | 0.056% | 2.966 |

## Why the rare-event calculation is different

A cycle starts with zero Page evidence and an empty block prefix. It ends either at a strict alarm on any pulse, including inside a block, or at a block end with Page evidence exactly zero. Only that block-end zero renews the full state. Reflections inside a block retain the measured prefix. Blocks are independent under the homogeneous no-change law, so renewal gives ARL = mean pre-change cycle length / pre-change cycle alarm probability.

The numerator uses pre-change cycles. An independent post-change proposal estimates the denominator with ordinary weights exp(−sum of ALL raw conditional LLRs), multiplied by the cycle-alarm indicator. The sum includes increments before Page reflections and stops at the actual alarm pulse. The weight is not exp(−Page), not a complete-block likelihood at a partial-block alarm, and not normalized by the sum of weights. The ratio estimator is consistent, not exactly unbiased for finite ensembles.

Post-proposal weights can be heavy-tailed; no iid exp(−h) bound is assumed. We report ESS, largest weight fraction, relative SE and two independent repetitions. Screening requires ESS>=1000, largest contribution<=1%, and relative SE<=2%. Repeats must agree within 4 combined SE; the new method must also agree within 4 combined SE with an independent direct simulation at the earlier approximately 1000-pulse setting. These diagnostics do not prove finite variance or rule out unobserved tail contributions.

Thresholds use 500,000 cycles per ensemble per search evaluation (three times as many for n>=8), with four fixed log-ARL corrections. They freeze before two independent holdout ensembles, each using 1,000,000 cycles per law (three times as many for n>=8). Validation never retunes thresholds. Target and inter-detector approximate 95% intervals must lie fully inside ±5%. All intervals are pointwise Monte Carlo approximations conditional on thresholds, not rigorous bounds or model uncertainty.

Raw renewal-cycle lengths, raw sums, final Page scores and event/completion flags are saved separately for both validation repetitions and both laws. Changed-run stopping times and event flags, calibration histories, random seeds, model settings, versions, CSV and plots are saved alongside. Output folders are never overwritten.

This is a single Na=0.1 slice with the five detector types in Figure 7, not the full figure or a pixel-level reproduction. The physical model remains ideal homodyne with known losses and no excess receiver noise; the conditional detector is not claimed optimal. One-based indexing explicitly fixes the phase; pulse 1005 is phase 1 for n=2 and n=4 and phase 5 for n=8.

Next: sweep Na and change phase using these matched-ARL checks. Larger blocks (such as n=16) need improved proposal design or substantially more samples: exploratory post-proposal pilots had unstable weights, so they are excluded from this confirmed comparison.
