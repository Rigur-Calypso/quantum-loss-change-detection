# Entangled initial-channel robustness

True initial transmission is 89.9%, 90% or 90.1%. The detector still assumes
90% before the change and 85% after it. True post-change transmission stays
85%; N=100, Na=0.1. Every threshold is frozen at the source calibration.
No author latency is fitted. Fresh estimates are required because the no-change
law changes. ARLs below are in pulses, ±1 Monte Carlo SE:

| Entangled block n | True eta1=0.899 | True eta1=0.900 | True eta1=0.901 |
| --- | --- | --- | --- |
| 2 | 1,253,419 ± 6,223 | 1,984,484 ± 10,098 | 3,150,930 ± 16,525 |
| 4 | 1,240,181 ± 4,489 | 1,977,829 ± 7,308 | 3,202,433 ± 12,227 |
| 8 | 1,218,423 ± 3,321 | 1,995,763 ± 5,554 | 3,289,063 ± 9,443 |

Relative to each fresh matched-channel estimate, the lower-transmission cases
reduce ARL by 36.8–38.9%, while higher-transmission
cases increase it by 58.8–64.8%.
The largest observed delay span across the three initial channels for any block
size is 0.200 pulse. Those are Monte Carlo point estimates, with SEs
saved separately; they are not deterministic bounds on true delay changes.

Nine settings each have two independent 250,000-cycle
ensembles per sampling law and a 30,000-trial changed
ensemble. Numerator and denominator moments are pooled before forming their
ratio; ARL ratios are never averaged. Six h=3 controls compare fresh renewal
estimates against independent 10,000-trial direct no-change simulations.
All diagnostics, repetitions, controls and censoring checks pass:
**True**. Approximate intervals and diagnostic thresholds
are not finite-variance proofs or certified coverage guarantees.

Each renewal begins with an empty block and W=0, ending at a first alarm or a
block boundary with W=0. Reflection within a block retains its observation
history. The adaptive proposal mixes the actual null conditional law with the
assumed post conditional law. The exact path weight uses that whole mixture,
not its sampled component. Page evidence still uses the assumed pre/post
score. The weight and score therefore use different likelihood ratios under
mismatch and are saved separately. No clipping or self-normalization is used.

Changed observations retain the exact physical covariance for mixed-loss blocks.
The change is at pulse 1005; delay is E[T−1005 | T>=1005]. Early first alarms
are counted separately and never restarted; censored survivors make unrestricted
mean delay unavailable. Dense-covariance oracle tests check the physical sampler
and conditional likelihood weights; scalar reductions match numerical equations.
No iid solver is applied to correlated increments.

These are robustness results, not recovered author settings or a reconciliation
of the published latency gap. Original source evidence is checksum-protected.
Full cycle arrays, path weights, changed/control runs, seeds, CSV and plots are saved.

Next: use these measured sensitivities to design a detector that tolerates initial
channel uncertainty, with fresh false-alarm calibration and held-out validation.
