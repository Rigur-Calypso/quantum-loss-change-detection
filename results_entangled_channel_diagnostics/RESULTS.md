# Entangled channel-loss diagnostic

We reused **eta2=0.84519511**, previously fitted using
only the coherent Na=0.1 mean. No additional author means or thresholds were
fitted. Detector likelihoods still assume eta1=0.9 and eta2=0.85; N=100.
The true pre-change eta1 remains 0.9. The exact mixed-loss covariance is retained
when the physical change at pulse 1005 falls inside a block. Conditional history
resets at block boundaries; Page evidence persists, including through reflections.

**0/27 entangled means** fall within a descriptive ±5%
of published values. Including the earlier scalar diagnostic gives
**8/44 held-out means** within ±5%; the one
fitted classical point is excluded from this total.

- n=2: 0/9 within ±5%.
- n=4: 0/9 within ±5%.
- n=8: 0/9 within ±5%.

Mean absolute relative differences across each nine-point energy grid:

| Block size | Original matched simulation | Common-loss diagnostic |
| --- | --- | --- |
| 2 | 34.39% | 13.26% |
| 4 | 30.60% | 10.34% |
| 8 | 28.63% | 11.04% |

Each correlated setting has two independent 30,000-trial
ensembles, pooled at the raw stopping-time level. All repeatability/censoring
checks pass: **True**. Reported uncertainty is our Monte Carlo
SE only; author uncertainties are unavailable. ±5% is not a significance test.
Delay means E[T−1005 | T>=1005]; early first alarms are counted separately and
never restarted. A censored survivor makes unrestricted mean delay unavailable.

The entire no-change stopping law is unchanged, so original near-two-million
pulse ARL evidence applies unchanged; it is not re-estimated here. Unit tests
compare the streaming implementation with an independent dense-covariance
conditional-Gaussian calculation, including within-block changes and first alarms.
No iid scalar delay solver is applied to correlated blocks.

The fitted loss improves agreement but is a diagnostic assumption, not an
estimate of the authors' actual channel. Figure 7 does not save eta separately.
Figure 3's specified channel is not refitted. Original source files and scalar
results are protected by checksums; new raw ensembles, seeds and summaries are saved.
See combined_predictions.csv and entangled_channel_predictions.png.

Next: consolidate the remaining discrepancies and test sensitivity to the true
pre-change channel, which can also change false-alarm performance. Such a study
needs new ARL checks; the current unchanged-null argument would no longer apply.
