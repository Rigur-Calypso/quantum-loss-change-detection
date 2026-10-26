# Channel-parameter diagnostic

A single true post-change transmissivity was fitted to the published coherent
mean at Na=0.1: **eta2=0.84519511**, instead of the
detector's assumed 0.85. Pre-change eta1 stays 0.9, N stays 100. Every detector
keeps its original likelihoods and ARL-calibrated threshold. This changes the
physical post-change observations, including squeezed-state variance.

The fitted point is training data; **8/17 held-out scalar means**
are predicted within ±5%. This descriptive tolerance is not a statistical test:
author uncertainty estimates are unavailable. Correlated entangled blocks are
outside this experiment's scope. The fit does not establish the authors' loss.
Figure 7 source data do not separately save eta; Figure 3 explicitly specifies
0.9/0.85 and is not refitted here.

All 18 predictions use 512/1024-cell deterministic calculations and fresh
30,000-trial simulations, with independent seeds per setting.
All numerical/MC and mesh checks pass: **True**.
Early first alarms are recorded separately; delay is E[T−1005 | T>=1005].
Censored delays are never reported as unrestricted means.

The no-change law is mathematically identical to the source study because
only the true post-change law changes. Existing ARL evidence (near two million
pulses) is reused unchanged; no new ARL estimate or calibration is claimed.
Original report, raw scalar runs and author files are checksum-protected.
See predictions.csv, report.json, raw runs and channel_predictions.png.

The common change reduces the gap but does not bring all held-out scalar means
within the chosen tolerance. Next: extend the physical-mismatch diagnostic to
correlated blocks, preserving their exact mixed-loss covariance and prefix
likelihoods. Requesting the authors' original
thresholds and simulator would resolve assumptions more directly; no message
has been sent.
