# Parameter sweep: independent squeezing versus coherent light

Grid: Na=[0.0, 0.01, 0.1, 1.0]; target ARL=[500.0, 1000.0, 5000.0]. Fixed N=100, eta1=0.9, eta2=0.85. Same total energy N+Na within each comparison.

**12/12 settings passed the prespecified approximate ARL equivalence criterion:** the entire pointwise 95% ARL-ratio interval lies within ±5%, and both validated ARL estimates lie within that tolerance of target. Unmatched settings, if any, remain visible and are exploratory only.

| Na | Target ARL | Validated ordinary / squeezed ARL | Delay ordinary / squeezed (nu=100) | Delay gain ratio | Approx. 95% ratio interval | Reduction | Matched |
|---|---:|---:|---:|---:|---:|---:|---|
| 0 | 500 | 499.1 / 498.3 | 19.99 / 19.96 | 1.001 | [0.989, 1.014] | 0.1% | True |
| 0.01 | 500 | 496.8 / 499.8 | 20.06 / 17.41 | 1.152 | [1.138, 1.166] | 13.2% | True |
| 0.1 | 500 | 501.7 / 505.2 | 19.90 / 13.32 | 1.495 | [1.476, 1.513] | 33.1% | True |
| 1 | 500 | 495.5 / 498.9 | 19.75 / 6.64 | 2.976 | [2.939, 3.014] | 66.4% | True |
| 0 | 1000 | 998.3 / 1003.1 | 24.31 / 24.36 | 0.998 | [0.987, 1.009] | -0.2% | True |
| 0.01 | 1000 | 1000.1 / 1010.7 | 24.33 / 21.06 | 1.155 | [1.143, 1.168] | 13.4% | True |
| 0.1 | 1000 | 1002.3 / 1003.0 | 24.16 / 15.87 | 1.523 | [1.506, 1.539] | 34.3% | True |
| 1 | 1000 | 1002.5 / 998.6 | 24.23 / 7.82 | 3.099 | [3.065, 3.133] | 67.7% | True |
| 0 | 5000 | 5014.8 / 4984.0 | 35.16 / 34.76 | 1.011 | [1.002, 1.021] | 1.1% | True |
| 0.01 | 5000 | 4980.9 / 4987.0 | 34.94 / 30.09 | 1.161 | [1.150, 1.172] | 13.9% | True |
| 0.1 | 5000 | 4981.1 / 4949.7 | 34.92 / 22.10 | 1.580 | [1.565, 1.595] | 36.7% | True |
| 1 | 5000 | 5058.4 / 4964.4 | 34.73 / 10.61 | 3.275 | [3.244, 3.306] | 69.5% | True |

## Control and interpretation

At Na=0, both physical models and thresholds are identical. Only the observations are drawn independently. The delay ratio should be near 1; the 95% intervals need not all contain 1 by chance. This control checks that the pipeline does not manufacture a squeezing advantage.

Na=0.01, 0.1 and 1 test progressively stronger squeezing. Each condition includes changes at pulse 1 and pulse 100. The CSV and JSON report conditional delays, survivor counts and early-alarm probabilities for both. nu=1 has no pre-change alarms; nu=100 conditions on no earlier alarm. This is not a worst-case delay calculation.

The plot uses equally spaced categories for the tested Na values. The dashed KL ratio is an asymptotic reference, not an exact finite-ARL delay prediction. A growing delay gain as ARL increases is suggestive of that limit, not a proof of convergence. Exact equality of true ARLs is not established by finite Monte Carlo samples.

Calibration uses 30,000 trials/evaluation, separate validation 50,000 per method/setting, and 30,000 delay trials per method/change point. The broader sweep uses a 5% equivalence band by default; the earlier high-precision single-setting comparison retained its separate 2% band.

All intervals are pointwise approximate Monte Carlo intervals; there is no multiple-comparison correction and no simultaneous confidence claim for the whole grid. Threshold calibration uncertainty is not included in the fixed-threshold delay intervals. Every condition was prespecified; no failing validation was discarded or used to retune its thresholds.

## Saved evidence and scope

Each condition folder contains report.json with model parameters, stream seeds and calibration history, plus raw_runs.npz containing all validation/delay stopping times and event flags. The aggregate report.json, summary.csv and figures use those records. A fresh output folder is required. No relevant calibration/validation/delay ensemble was censored; the runner refuses unrestricted comparisons otherwise.

This extends the moderate-ARL numerical validation. It does not reproduce Figure 3 at ARL=2,000,000, implement entangled blocks, or validate a real optical experiment. Next: efficient large-ARL calibration and the correlated entangled model.
