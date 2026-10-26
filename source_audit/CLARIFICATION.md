# Questions needed to finish numerical reproduction

Paper: [Quantum-enhanced quickest change detection of transmission loss,
arXiv:2503.12276v2](https://arxiv.org/abs/2503.12276v2),
[Physical Review Letters 135, 210801](https://doi.org/10.1103/czg5-3y3d).

The Python implementation follows the stated observation laws, strict Page
first-alarm rule and conditional Gaussian block likelihoods. Independent dense
covariance calculations, adaptive phase integration, Monte Carlo experiments
and rare-event ARL estimates support its internal calculations. The following
questions concern exact agreement with the published numerical figures, rather
than the qualitative advantage of squeezing or entanglement.

## 1. Figure 4: fixed-seed conversion and covariance

At fixed seed squeezing s, the stated energy convention gives
Na=sinh²(s)/n. The collective variance under channel j is
v_j=[1−eta_j+eta_j exp(−2s)]/4, independent of n. Eqs.4–6 then reduce the
per-pulse Gaussian KL to B+A/n, with

```
q = v2/v1
A = (q − 1 − log(q))/2 >= 0
B = N*(sqrt(eta2) − sqrt(eta1))²/(2*v1)
```

For N=100, eta1=0.9, eta2=0.85 and seed squeezing 18 dB, this gives
1.28675 nats/pulse at n=1 and 1.26862 at n=2. The published inset instead
places n=2 above n=1. Its PNG is identical in both arXiv versions.

**Requested clarification:** the original Figure 4 calculation, including
conversion from seed dB to per-pulse energy, displacement energy, covariance
matrix and division by block length. Is a different energy constraint intended?

## 2. Figure 5: Kennedy energy and noise convention

Under the legend's N=100, Na=1, tap=0.9, epsilon=0.01 and sigma²=0.01,
use sigma=0.1 and energy-matched alpha²=101. At zero initial loss,
Eq.10 gives KL 0.0589971 nats/pulse; the published black marker is approximately
0.00486 (manual pixel uncertainty approximately ±5.6%). The independent
adaptive integration and phase/Poisson simulation agree with the former value.

The same model with 10 photons yields approximately 0.004828 at that point.
A one-point energy fit of 10.057 photons also comes within approximate pixel
uncertainty at three held-out losses. This is only a non-unique diagnostic;
it does not identify what was actually used to create the figure.

**Requested clarification:** the numerical alpha/energy and epsilon/sigma
values passed to the original calculation; whether the phase was averaged
before the KL calculation as in Eq.10; count cutoff/normalization; and the
original plotting array or script. Please also confirm any intended distinction
between absolute residual amplitude in Eq.7 and fractional LO error in Eq.10.

## 3. Figures 3/7/11: original stopping-time experiment

The released Zenodo tables omit the simulation thresholds and raw first-alarm
arrays for Figures 3/7. The supplied plotting script reads previously saved
results. Our comparisons calibrate thresholds to no-change ARL near two
million and count the first strict crossing only. Detection delay is
E[T−nu | T>=nu]; early first alarms and censored paths are separate.

For example, at Figure 7 Na=0.1, our coherent mean delay is 76.434±0.182
pulses (one SE) versus the released mean 64.624. This computation assumes
eta1=0.9, eta2=0.85, N=100 and change pulse 1005; the Figure 7 data file does
not separately save eta1/eta2. The independent scalar first-step equations
support our mean. Correlated cases have their own conditional-LLR calculations
and validation rather than inheriting scalar thresholds.

The supplemental single-photon latency array contains two negative values.
Its plotting script uses their signed mean, while our conditional first-alarm
delay is 54.640±0.129 pulses. Those quantities are not silently equated.

**Requested clarification:** original per-scheme thresholds and ARL lookup
tables; true and assumed channel parameters; whether an early alarm ends a
run or is reset; treatment of censored paths; the exact change-pulse convention
and block alignment; trial counts; and raw first-alarm/latency arrays with
uncertainty estimates. Does the shared-threshold ARL routine assign the same
pulse to every threshold crossed by one overshooting observation?

## Local evidence available for review

- [Full coverage and remaining limits](../PAPER_IMPLEMENTATION.md)
- [Author-grid comparison](../results_author_comparison/RESULTS.md)
- [PDF-verified parameters and pixel comparisons](../results_verified_figures_v2/RESULTS.md)
- [Independent integration and dense simulation](../results_equation_validation/RESULTS.md)
- [Public-source availability audit](REPORT.md)

This note has been prepared locally for review. No message or file has been
sent to the authors, and no original threshold has been inferred by fitting
published latency means.
