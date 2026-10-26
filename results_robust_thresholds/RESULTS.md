# Conservative thresholds for uncertain initial transmission

This extends the paper implementation with a robustness experiment. Score laws
remain at eta1=0.9/eta2=0.85, N=100, Na=0.1; true post-change eta2 stays 0.85.
Five initial transmissions, 0.899 through 0.901 in steps of 0.0005, are sampled.
Original results remain unchanged. Only the Page alarm threshold is raised.

Four fixed calibration rounds use independent 100,000-cycle ensembles per law
at every grid point, targeting a minimum sampled ARL of 2.15 million. These
calibration streams are separate from all holdout and changed-run streams.
Thresholds are then frozen. No holdout results or published delays are fitted.

Independent validation passes at **25/25 sampled settings**;
overall status: **True**. Each setting has two independent
250,000-cycle ensembles per sampling law. Raw
cycle moments are pooled before forming the ARL ratio. Each adjusted approximate
lower ARL must exceed two million, with ESS, largest-weight, relative-SE and
repeatability checks also passing.

| Detector | Original h | Conservative h | Smallest adjusted lower ARL | Observed delay cost (pulses) |
| --- | --- | --- | --- | --- |
| Coherent | 11.93340 | 12.44615 | 2,073,242 | 3.61–4.39 |
| Squeezed n=1 | 12.24644 | 12.79556 | 2,111,868 | 1.63–2.33 |
| Entangled n=2 | 12.32943 | 12.88009 | 2,117,885 | 1.46–1.94 |
| Entangled n=4 | 12.40393 | 12.97099 | 2,099,537 | 1.31–1.61 |
| Entangled n=8 | 12.45083 | 13.04262 | 2,129,459 | 1.12–1.31 |

The lower bounds use one-sided Bonferroni allocation of nominal 5% error over
25 points and asymptotic log-scale delta SEs. This is approximate sampled-grid
evidence, not certified coverage or a guarantee for every point in the continuous
interval. No monotonicity in true transmission is assumed for correlated scores.
Renewal boundaries and exact actual-null conditional mixture corrections are
retained; no iid solver is substituted for correlated data.

Delay costs use independent fresh 20,000-trial
ensembles at the original and conservative thresholds, at each true initial
channel. The change is at pulse 1005; delay is E[T−1005 | T>=1005]. SEs for the
differences combine the two independent errors. Early first alarms stay separate,
and censored survivors make unrestricted means unavailable. The plotted delay
error bars are pointwise ±1.96 Monte Carlo SE, not simultaneous bounds.

The source evidence is checksum-protected. Calibration histories and seeds,
all holdout cycle arrays and path weights, baseline/conservative first-alarm runs,
CSV and plots are saved. Calibration arrays can be regenerated from recorded
seeds. This does not reconcile the publication's latency gap or identify author
parameters. The conservative thresholds have deliberately unequal ARLs away
from the least favorable sampled channel, so this is not exact ARL matching
across detectors at every channel.

Next: test channels between the calibration points with fresh samples, then
study estimating the initial channel from a separate training period.
