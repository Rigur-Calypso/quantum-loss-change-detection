# Independent checks of the Figure 4/5 discrepancies

The stated receiver equations survive independent integration and simulation.
A lower-energy Kennedy model is consistent with the few published pixels
checked, but this does not establish the authors' actual plotting parameters.

Source: [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2),
Eqs.4–6 and Appendix A Eq.10. The source settings/pixel readings are the
unchanged [PDF-verified study](../results_verified_figures_v2/RESULTS.md).
Its report checksum is saved in this study's `report.json` and checked again
after calculation. No original study outputs are replaced.

![Independent calculation and energy comparison](independent_checks.png)

## Kennedy equation: three independent calculations

For initial loss 0, 2, 4 and 6 dB we compare:

1. The existing production phase-quadrature implementation of Eq.10.
2. Direct adaptive integration of the complex-amplitude Poisson law for each
   photon count 0–128, followed by the marginal relative-entropy sum.
3. Direct random Gaussian phase draws, complex-amplitude nulling, Poisson
   counts and average marginal log-likelihood ratio: 500,000 draws per setting.

All use 101 probe photons, eta2=0.9 eta1, fractional amplitude error 0.01 and
phase standard deviation 0.1. Phase integration omits at most the |Z|>12
Gaussian mass (<4×10⁻³³); physical truncation is [-pi,pi].

The two integration methods differ by at most **2.27×10⁻¹⁴ relative**.
All four random estimates lie within **1.94 standard errors** of the integrated
value. All count-mean checks agree with independently derived mixed-Poisson
moments. The largest integrated probability above count 128 is **1.52×10⁻²⁶**;
extending the KL sum from 64 to 128 changes it by at most **1.09×10⁻¹³ nats**.
This is strong numerical evidence for the stated calculation, not a claim of a
formal interval bound on every integration error.

## A concrete energy diagnostic, with held-out points

The published Kennedy curve can be approximated by the same model using
**10 photons**. Its relative differences from the four pixel readings are
−0.61%, +1.39%, +4.26% and +2.17%. This counterfactual energy is kept separate
from the paper's stated bright pulse and from all earlier simulations.

As a separate inverse-model diagnostic, fitting energy to **only the 0 dB
point** gives **10.0570 photons**. With this energy then fixed:

| Initial loss | Role | Difference from published pixel reading |
| --- | --- | --- |
| 0 dB | Fitted | Approximately zero by construction |
| 2 dB | Held out | +2.00% |
| 4 dB | Held out | +4.87% |
| 6 dB | Held out | +2.76% |

All three held-out differences are within the prior study's approximate ±5.6%
pixel uncertainty. That uncertainty is a manual digitization estimate, not a
statistical confidence interval. This is a **non-unique diagnostic explanation**:
it does not prove a particular source-code error, identify the author energy,
justify changing N=100 to N=10, or resolve the published-setting reproduction.
No detection latency or ARL threshold is fitted here.

## Fixed-seed entanglement: independent dense simulations

Six additional cases use seed squeezing 6/18 dB and n=1/2/8. They build the full
covariance directly from Eqs.5–6, sample via dense Cholesky factors, and calculate
Gaussian likelihoods using dense solves. Each has 500,000 fresh draws, with
independent seeds. Every KL sample average is within **1.14 standard errors**
of the equation-based prediction. The calculation supports the previously
derived ordering n=1 ≥ n=2 ≥ n=8 at fixed seed; it does not reproduce the
published inset's n=2 > n=1 ordering. No covariance modification is introduced
to force agreement.

## Files, validation and reproduction

There are **five million independent draws** across ten cases. Raw count draws,
integrated probability tables and dense Gaussian per-pulse likelihood scores
are stored in `independent_samples.npz`; settings, seeds and standard errors
are in `report.json`. `energy_diagnostic.csv` records the stated setting, fixed
10-photon counterfactual and one-point fit separately. Samples concern relative
entropy, not new CUSUM stopping-time or false-alarm estimates.

The full suite passes **373 tests**. The new saved evidence passes **114 audit
checks**, including exact regeneration of all ten random ensembles and all
eight adaptive probability tables. Earlier source results are unchanged.

```sh
python -m lossqcd.equation_validation --output results_equation_validation_repeat
python -m lossqcd.audit_equation_validation results_equation_validation_repeat --regenerate
python -m pytest -q
```

The output directory must be new. Reproduction uses the existing environment
and cached source report; it needs no network access. The raw draw seed is
2026102401, expanded into ten independent child streams.

## Specific remaining source questions

- For Figure 5, which numerical probe amplitude/energy was passed to the Eq.10
  evaluation? Can its original plotting values and photon-count cutoff be
  supplied? The supplied legend settings give 0.0589971 nats at zero initial
  loss, whereas the black plotted marker is approximately 0.00486.
- For Figure 4, what transformation was used from fixed seed s to per-pulse
  squeezing energy and covariance? With Na=sinh²(s)/n, Eqs.4–6 give A/n+B,
  A≥0. At 18 dB they give n=1: 1.28675 and n=2: 1.26862 nats per pulse.
- For the remaining Figure 3/7/11 latency gaps, can the original thresholds,
  false-alarm calibration and first-alarm arrays be supplied? The existing
  released data do not identify those choices.

These are clarification questions for a future source comparison. No message
has been sent to the authors.
