# Coherent versus independently squeezed detection

Same total energy per pulse; separately calibrated thresholds; independent ARL validation and delay ensembles.

Parameters: N=100, Na=0.1, eta1=0.9, eta2=0.85; seed 20261022. Target ARL 1000 pulses.

| Method | Threshold | Validated ARL ± MC SE |
|---|---:|---:|
| Ordinary light | 4.377441 | 998.87 ± 3.12 |
| Independent squeezing | 4.680176 | 993.88 ± 3.10 |

Approximate ARL ratio 95% interval: [0.9963, 1.0137]. Entire interval inside ±2% equivalence band: **True**. Each validated ARL point estimate is also required within that tolerance of target. Validation was not reused to retune thresholds.

| Change at pulse | Ordinary delay ± SE | Squeezed delay ± SE | Delay ratio (approx. 95% interval) | Delay reduction | Early alarms: ordinary / squeezed |
|---|---:|---:|---:|---:|---:|
| 1 | 27.113 ± 0.072 | 17.449 ± 0.046 | 1.554 [1.542, 1.565] | 35.6% | 0.00% / 0.00% |
| 100 | 24.331 ± 0.075 | 15.797 ± 0.048 | 1.540 [1.527, 1.553] | 35.1% | 8.03% / 8.43% |
| 1000 | 24.307 ± 0.118 | 15.790 ± 0.076 | 1.539 [1.519, 1.560] | 35.0% | 63.19% / 62.97% |

## Interpretation and limits

Delay means E[T-nu | T>=nu], with zero delay for an alarm on the first changed pulse. Early alarms are excluded from this conditional mean and reported separately; this is not a worst-case delay. nu=1 is a fresh detector with no pre-change false alarms. Longer pre-change periods can produce substantial early-alarm rates even at the same ARL.

There were no censored runs in calibration, validation or the delay ensembles (the script refuses unrestricted comparisons if any are censored). Calibration uses 100,000 trials/evaluation; validation 100,000 per method; delays 50,000 per method/change point. Log-likelihoods include unequal variances for squeezed measurements.

The analytical KL ratio is 1.716219. It is a large-ARL reference, not an exact finite-threshold delay prediction. Uncertainty intervals describe Monte Carlo sampling at the fixed chosen thresholds; they do not include a hypothetical recalibration's variability or optical modeling errors. Exact equality of true ARLs is not established by finite simulations.

This is a moderate-ARL, single-Na computational comparison. It does not reproduce Figure 3 at ARL=2,000,000, implement entangled probes, or constitute a hardware experiment. Next: validate a range of Na and ARL settings, then implement efficient large-ARL calibration.

## Reproducibility

The JSON report records settings, stream seeds, calibration histories, distribution parameters, package versions and all summaries. Raw per-run stopping times and event flags are in raw_runs.npz. Streams 0/1 calibrate coherent/squeezed, 2/3 validate, and 4+2j/5+2j simulate the jth change point. All ensembles are independent. Each calibration stream spawns independent streams for its threshold evaluations.
