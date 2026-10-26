# Independent rare-event ARL validation

The existing numerical-calibration thresholds were frozen. This experiment estimates their no-change ARLs using an independent regenerative importance-sampling algorithm, with no grid solver and no retuning.

**Moderate ARL checks:** 24/24 passed against historical direct Monte Carlo. **Large ARL checks:** 14/14 passed against the prior fine-grid numerical estimates.

Agreement means absolute relative difference <= max(2%, 3 combined Monte Carlo SE / reference). The numerical reference has no assigned Monte Carlo SE; its discretization diagnostics are in the preceding experiment. All settings and thresholds were retained, including any failures.

| Na | Method | Fixed threshold | Independent ARL ± SE | Approx. 95% interval | Numerical reference | Agrees |
|---|---|---:|---:|---:|---:|---|
| 0 | coherent | 11.941600 | 1995557 ± 6909 | [1982062, 2009144] | 2000648 | True |
| 0 | squeezed | 11.941600 | 2003149 ± 6928 | [1989616, 2016774] | 2000648 | True |
| 0.02 | coherent | 11.941738 | 2006019 ± 6959 | [1992426, 2019705] | 2000648 | True |
| 0.02 | squeezed | 12.096169 | 1997706 ± 6438 | [1985128, 2010364] | 2000524 | True |
| 0.05 | coherent | 11.941944 | 2008205 ± 6979 | [1994571, 2021931] | 2000648 | True |
| 0.05 | squeezed | 12.176754 | 1996155 ± 6161 | [1984115, 2008268] | 2000465 | True |
| 0.1 | coherent | 11.942289 | 2004738 ± 6930 | [1991201, 2018367] | 2000648 | True |
| 0.1 | squeezed | 12.258676 | 2002168 ± 5933 | [1990574, 2013830] | 2000408 | True |
| 0.25 | coherent | 11.943321 | 1993559 ± 6869 | [1980141, 2007068] | 2000647 | True |
| 0.25 | squeezed | 12.392862 | 2008524 ± 5550 | [1997676, 2019431] | 2000319 | True |
| 0.5 | coherent | 11.945038 | 1994822 ± 6905 | [1981335, 2008401] | 2000646 | True |
| 0.5 | squeezed | 12.502344 | 2006890 ± 5150 | [1996823, 2017009] | 2000251 | True |
| 1 | coherent | 11.948455 | 1998854 ± 6896 | [1985384, 2012415] | 2000643 | True |
| 1 | squeezed | 12.598543 | 1993131 ± 4749 | [1983845, 2002460] | 2000190 | True |

## False-alarm matching within each pair

The entire pointwise approximate 95% ordinary/squeezed ARL-ratio interval must fit inside a prespecified ±3% band, and each independently estimated ARL must lie within 3% of target. This is approximate equivalence at finite simulation precision, not proof of exact equality.

| Na | Independent ARL ratio | Approx. 95% interval | Equivalence passed |
|---|---:|---:|---|
| 0 | 0.99621 | [0.98665, 1.00576] | True |
| 0.02 | 1.00416 | [0.99484, 1.01348] | True |
| 0.05 | 1.00604 | [0.99687, 1.01520] | True |
| 0.1 | 1.00128 | [0.99235, 1.01022] | True |
| 0.25 | 0.99255 | [0.98396, 1.00114] | True |
| 0.5 | 0.99399 | [0.98559, 1.00238] | True |
| 1 | 1.00287 | [0.99463, 1.01111] | True |

## Exact identities and estimator

Let ell=log(p2(X)/p1(X)). Starting at zero, run an unreflected sum S until it returns to/below zero or strictly crosses h. This is one cycle, with length tau and upper-crossing event A. A failed cycle resets Page CUSUM to zero, giving m=E_pre[tau]+P_pre(A-complement)*m; hence ARL=m=E_pre[tau]/P_pre(A).

Estimate E_pre[tau] with ordinary pre-change cycles. Independently sample cycles under the post-change law. For a stopped path the original/proposal density ratio is product(p1/p2)=exp(-S_tau), so P_pre(A)=E_post[exp(-S_tau) 1_A]. Every failed proposal cycle contributes ZERO, and the estimate divides by ALL proposal cycles. It is ordinary importance sampling, not self-normalized weighting. The terminal overshoot S_tau-h is retained. Successful weights are bounded by exp(-h). The two sample means are unbiased individually; their ratio is consistent but not exactly unbiased at finite sample size.

This is a specialization of measure-specific regenerative importance sampling; general background: [Calvin, Glynn & Nakayama (2001), section 3.2](https://web.stanford.edu/~glynn/papers/2001/CalvinGNakayama01.pdf). The source discusses estimating a regenerative numerator under the original measure and the rare-event denominator under an independent changed measure.

## Evidence, uncertainty and scope

Per model/setting: 100,000 pre and post cycles for each moderate-ARL check, and 500,000 each for the large-ARL checks. Seeds are independently spawned for every setting and measure. No cycles were unfinished; the estimator refuses to discard unfinished cycles or count them as failed alarms.

Per-check JSON files contain streams, cycle moments, proposal-success counts, effective weight sample size, actual simulated pulse counts and uncertainty. Compressed raw length/terminal-sum/upper-event/completion arrays are saved for both measures, so every estimate can be recomputed. Source report hashes identify the exact earlier thresholds. Code uses no new dependencies.

Uncertainties use an independent-means delta-method SE and a positive log-scale approximate 95% interval. They do not include physical-model uncertainty or provide simultaneous coverage across all settings. Agreement checks are diagnostics, not exact hypothesis-free guarantees. This is weighted Monte Carlo validation; it does not mean that full no-change sequences were simulated until millions of pulses per alarm.

The preceding detection-delay estimates are unchanged. The new check addresses the previously unverified rare-event false-alarm calibration under the same ideal Gaussian models. Entangled probes remain unimplemented; the next model extension is correlated n-mode blocks.
