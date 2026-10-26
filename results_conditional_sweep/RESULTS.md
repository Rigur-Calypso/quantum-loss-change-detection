# Energy and block-phase sweep

Model and pulse-by-pulse detector: [paper Figure 7 and Appendix B.4](https://arxiv.org/html/2503.12276v2#A2.SS4).
This is an independent, sampled Figure 7-style study with explicit ARL matching and one-based indexing. It is not a digitized or supplemental-data error comparison.

**23/25 settings confirmed.** Zero-energy delay-equivalence control: **True**.
N=100, eta1=0.9, eta2=0.85; target ARL=2,000,000 pulses. Na values: [0.0, 0.01, 0.1, 1.0, 5.0]. Each detector uses offsets [1, 2, 3, 4, 5, 6, 7, 8] after 1000 pre-change pulses.

## Change on pulse 1005

| Na | Detector | ARL ± SE | Delay ± SE | Coherent / detector delay | Confirmed |
|---:|---|---:|---:|---:|---|
| 0 | Coherent | 1,998,737 ± 3,460 | 76.299 ± 0.222 | 1.000 | True |
| 0 | Squeezed n=1 | 1,998,737 ± 3,460 | 77.031 ± 0.224 | 0.991 | True |
| 0 | Entangled n=2 | 1,998,737 ± 3,460 | 76.864 ± 0.224 | 0.993 | True |
| 0 | Entangled n=4 | 1,998,737 ± 3,460 | 76.907 ± 0.225 | 0.992 | True |
| 0 | Entangled n=8 | 1,998,737 ± 3,460 | 76.933 ± 0.225 | 0.992 | True |
| 0.01 | Coherent | 1,999,633 ± 3,466 | 76.508 ± 0.222 | 1.000 | True |
| 0.01 | Squeezed n=1 | 2,002,246 ± 3,296 | 64.956 ± 0.188 | 1.178 | True |
| 0.01 | Entangled n=2 | 2,000,113 ± 2,512 | 61.077 ± 0.178 | 1.253 | True |
| 0.01 | Entangled n=4 | 1,999,120 ± 2,117 | 55.727 ± 0.162 | 1.373 | True |
| 0.01 | Entangled n=8 | 2,068,263 ± 12,942 | 49.998 ± 0.142 | 1.530 | True |
| 0.1 | Coherent | 2,003,744 ± 3,463 | 76.940 ± 0.223 | 1.000 | True |
| 0.1 | Squeezed n=1 | 2,006,419 ± 2,981 | 46.515 ± 0.134 | 1.654 | True |
| 0.1 | Entangled n=2 | 1,990,401 ± 2,256 | 38.683 ± 0.114 | 1.989 | True |
| 0.1 | Entangled n=4 | 2,004,672 ± 1,944 | 31.188 ± 0.091 | 2.467 | True |
| 0.1 | Entangled n=8 | 2,056,089 ± 20,014 | 25.944 ± 0.072 | 2.966 | False |
| 1 | Coherent | 1,994,537 ± 3,439 | 75.594 ± 0.220 | 1.000 | True |
| 1 | Squeezed n=1 | 2,003,116 ± 2,390 | 21.027 ± 0.064 | 3.595 | True |
| 1 | Entangled n=2 | 2,004,274 ± 2,016 | 16.449 ± 0.051 | 4.596 | True |
| 1 | Entangled n=4 | 1,995,582 ± 2,443 | 13.694 ± 0.041 | 5.520 | True |
| 1 | Entangled n=8 | 1,979,027 ± 11,829 | 14.672 ± 0.037 | 5.152 | True |
| 5 | Coherent | 2,004,018 ± 3,413 | 73.194 ± 0.213 | 1.000 | True |
| 5 | Squeezed n=1 | 2,000,692 ± 2,149 | 11.877 ± 0.039 | 6.162 | True |
| 5 | Entangled n=2 | 1,998,783 ± 2,054 | 10.683 ± 0.035 | 6.851 | True |
| 5 | Entangled n=4 | 2,002,355 ± 2,327 | 10.524 ± 0.033 | 6.955 | True |
| 5 | Entangled n=8 | 2,003,434 ± 18,511 | 13.241 ± 0.032 | 5.528 | False |

At change pulse 1005, the largest confirmed speedup on this sampled grid is 6.955 times (Entangled n=4, Na=5; delay 10.524 pulses). This is a grid result, not a global optimum.

Unconfirmed settings and their validation screens:
- Entangled n=8, Na=0.1: ARL match=True; repeat/weight screen=False; moderate direct crosscheck=True.
- Entangled n=8, Na=5: ARL match=True; repeat/weight screen=False; moderate direct crosscheck=True.

## Phase sensitivity for n=8

Ranges below are descriptive minima/maxima over the eight sampled change times, not worst-case Lorden delay or confidence intervals for extrema. Phase estimates share a burn-in and are correlated.

| Na | Minimum delay | Maximum delay | Range (pulses) | Confirmed |
|---:|---:|---:|---:|---|
| 0 | 76.217 | 76.986 | 0.769 | True |
| 0.01 | 49.693 | 50.076 | 0.383 | True |
| 0.1 | 24.943 | 25.944 | 1.001 | False |
| 1 | 12.894 | 14.672 | 1.778 | True |
| 5 | 11.234 | 13.241 | 2.007 | False |

## Calibration, validation and control rules

Each positive-energy setting gets a fresh calibration, two independent large-ARL holdouts, and a fresh direct no-change check at h=4.5. Four log-ARL calibration corrections are followed by one independent refinement using twice the per-evaluation sample budget. All thresholds freeze before holdouts; no holdout can retune them. The direct check must agree with a separate importance-sampling estimate within four combined SE, and that weighted estimate must pass the diagnostic screen.

The default budgets are 500,000 cycles per law per calibration evaluation, 1,000,000 per law per holdout repetition, and 500,000 per law for the moderate check. All these budgets are multiplied by three for n>=4. Direct moderate checks use 10,000 first-alarm runs with horizon 100,000; changed-run phase branches use 20,000 trials with horizon 10,000. All settings and seeds are in report.json.

Confirmation requires full pointwise approximate 95% ARL intervals within ±5% of target and ARL-ratio intervals within ±5% of the same-energy coherent baseline. Two holdouts must agree within four combined SE. Weighted ESS>=1000, largest weight fraction<=1%, and relative SE<=2% must hold in both repetitions, pooled holdout and moderate estimate. A detector also requires its coherent comparator to pass repeat and moderate checks. Failing settings remain in tables, with crosses in the energy plot and asterisks in the phase heatmap. Their delay gains are exploratory, not confirmed fair comparisons.

At Na=0 the covariance is I/4 and every model has the identical coherent pulse law. One scalar calibration/ARL validation is therefore shared explicitly across all five labels. No independent-ratio SE is applied to this shared ARL. Each label still gets an independent changed-run experiment. Every zero-control pointwise approximate 95% delay-ratio interval must fit inside ±3% of one.

The batched estimator pools the numerator and denominator first/second moments, then forms their ratio; it never averages batch ARL ratios or normalizes likelihood weights. Each batch retains its seed, full moment/weight diagnostics and counts. One 100,000-cycle audit batch per holdout retains raw lengths, raw LLR sums, Page scores, success/completion flags for both sampling laws. The remaining cycle arrays can be regenerated exactly from recorded seeds, sizes and pinned dependencies. Changed stopping times and event flags are all retained. This reduces storage without discarding failed or extreme weighted batches.

Phase branches share a common pre-change prefix ending at a block boundary, retaining Page evidence and first-alarm flags. Each continuation has an independent random stream and the exact mixed-loss within-block covariance. No early alarm is restarted. Delay is E[T−nu | T>=nu]; early alarms are reported separately. Relevant censoring refuses unrestricted estimates. Pulse 1005 is phase 1 for n=2,4 and phase 5 for n=8 under the explicit one-based convention.

Intervals describe Monte Carlo sampling uncertainty conditional on fixed thresholds, not simultaneous coverage across all settings or optical-model uncertainty. Importance-sampling diagnostics do not prove finite variance or exclude unseen tail contributions. Zero-control labels share ARL evidence, and phase branches share burn-in; they are not independent ARL confirmations or independent phase measurements.

The supplemental Zenodo dataset linked by the paper could not be retrieved in this session. No numerical agreement with its raw Figure 7 data is claimed. The selected Na grid spans [0,5] with five samples; it does not resolve the entire continuous curve. The model remains ideal homodyne with known losses. Next: resolve any unconfirmed settings and compare directly with the authors’ supplemental data when accessible.
