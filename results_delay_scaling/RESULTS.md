# Finite-false-alarm check of the paper's detection-delay scaling

The new study supports the logarithmic delay growth predicted by
[the paper's Eq.1](https://arxiv.org/html/2503.12276v2#S0.E1), while showing
that its first-order expression is not an exact finite-target delay formula.
This is an additional equation check, not a reproduction of a new published
figure or a solution of an open research question.

![Finite-ARL delay scaling](delay_scaling.png)

## Setup

Use coherent and independently squeezed probes, N=100, Na=0.1, eta1=0.9 and
eta2=0.85. Seven no-change ARL targets are 100, 1,000, 10,000, 100,000,
2,000,000, 10,000,000 and 100,000,000 pulses. Thresholds are calibrated to
ARL using numerical first-step equations, never fitted to delay or author data.
All logarithms are natural and entropy is in nats per pulse.

For each of the 14 settings:

- Solve ARL at 512 cells, checking the same threshold at 256/1024 cells.
- Solve post-change stopping time at 512/1024 cells.
- Independently check ARL with 150,000 pre-law and 150,000 post-proposal
  renewal cycles, retaining ordinary likelihood weights and overshoots.
- Directly simulate 40,000 changed trajectories starting on pulse 1.

At target ARL two million, each family also has 40,000 independent runs with
change on pulse 1,000. Their early first alarms are counted separately; delay
is conditional on surviving to the change. No unrestricted mean is reported
from a censored ensemble. All changed ensembles here finish without censoring.

## Findings

| Target no-change ARL | Coherent cold-start delay | Eq.1 coherent prediction | Squeezed cold-start delay | Eq.1 squeezed prediction |
| --- | --- | --- | --- | --- |
| 100 | 12.711 | 32.197 | 8.740 | 18.761 |
| 2 million | 79.892 | 101.447 | 48.259 | 59.109 |
| 100 million | 107.260 | 128.811 | 64.202 | 75.050 |

Cold-start delay means **E[T−1]**, so an alarm on the first changed pulse has
zero delay. Expected post-change samples are E[T], exactly one larger; both
are saved. The first-order comparison is log(numerical ARL)/D(post||pre).
Using the independently checked numerical ARL rather than its nominal target
keeps finite-grid calibration error out of the horizontal coordinate.

Between ARL two million and 100 million, normalize the observed increase in
delay by the predicted increase `log(ARL_high/ARL_low)/D`. The resulting slopes
are **1.000138 coherent** and **1.000084 squeezed**, deviations of about
0.014% and 0.0084% from one. At the same time, the absolute delay divided by
the first-order expression at 100 million is only **0.8327** and **0.8554**.
These results support the logarithmic growth rate while exposing a finite
offset that still matters even at very large ARLs.

The scalar cold-start experiment is deliberately distinguished from the earlier
change-at-1,000 experiments. Fresh latter means at target ARL two million are
**76.473 coherent** and **46.530 squeezed** pulses. Their ensembles contain
24 and 22 early first alarms respectively, out of 40,000 runs; those alarms
are excluded from conditional delay and explicitly saved. First-step solutions
and the independent simulations agree within 1.01 standard errors in these
two comparisons.

This does not explain the remaining author-data latency gaps: the source's
simulation thresholds and event conventions are still unavailable. Nor does
this finite experiment establish a Lorden worst-case theorem, extend that
theorem to correlated pulse-level increments, or make an empirical small-grid
refinement criterion into a certified numerical error bound.

## Validation and retained evidence

All **14/14 settings pass**. The largest ARL discrepancy from independent
renewal sampling is **1.62 SE**, and the largest cold-start delay discrepancy
from direct simulation is **1.67 SE**. Maximum relative 512→1024-cell changes
are **0.311% for ARL** and **0.000646% for post-change stopping time**.
The predeclared ARL mesh tolerance is 0.5%; the post-time criterion is 0.1%.
Importance-weight effective sample sizes exceed 1,000 and all renewal cycles
complete. Uncertainty intervals are approximate Monte Carlo diagnostics.

The study saves 4.2 million renewal cycles, 560,000 cold-start runs and 80,000
change-at-1,000 runs, with independent seeds. `report.json` contains all model
settings, calibration history, grids, seeds, summary statistics and criteria.
Each case directory retains raw renewal arrays and first-alarm arrays.
`scaling.csv` provides the main numerical table.

The full project suite passes **379 tests**. The new study passes **267 saved-data
checks**, including exact regeneration of every cycle and changed ensemble.

```sh
python -m lossqcd.delay_scaling --output results_delay_scaling_repeat
python -m lossqcd.audit_delay_scaling results_delay_scaling_repeat --regenerate
```

Use a fresh output directory. These commands use the existing numerical
environment and require no network access or additional packages. Earlier
study outputs and first-alarm conventions are unchanged.
