# Controlled diagnostics of the author-latency gap

Reference: [authors’ supplemental data](https://doi.org/10.5281/zenodo.17460928) and [paper Appendix F](https://arxiv.org/html/2503.12276v2#A6). The source comparison report and raw runs are preserved; this folder contains separate diagnostic experiments.

## Findings

- Across all 45 Figure 7 settings, including signed early alarms changes our saved conditional mean by at most **0.355 pulses**. A one-pulse index shift or other listed convention brings any point within ±5% of its author mean: **False**. These checks concern our observed first alarms; the authors’ event processing remains unknown.
- An independent deterministic calculation agrees with the scalar coherent/squeezed Monte Carlo delays at both meshes for all 18 settings: **True**. The method does not use importance sampling or Monte Carlo observations. It shares the Gaussian LLR CDF and physical assumptions with the implementation, so it is not an independent check of every optical assumption.
- For six representative scalar cases, fitting a LOWER diagnostic threshold to the published latency mean yields independently estimated ARLs of **187,147–734,575 pulses**, rather than two million. All diagnostic ARL intervals lie below 95% of the original target: **True**. All fresh diagnostic validation checks pass: **True**.
- The n=8 probe varies the first changed pulse across all eight block positions at one unchanged threshold. Its phase ranges and the author reference are shown separately below. Phase estimates share a burn-in; their range is descriptive, not a worst-case Lorden bound.

These results rule out the tested convention changes as a sufficient explanation within our fixed model. A lower threshold can imitate a published mean, but sacrifices the validated false-alarm interval. This does not establish that the authors used a lower threshold or that their implementation is incorrect. Their thresholds, event arrays, calibration tables and uncertainty estimates are not supplied.

## Deterministic scalar method

Finite-volume transition matrices include the exact Gaussian LLR CDF cell probabilities, the Page atom at zero and first-alarm absorption. Starting at W=0, propagate the pre-change distribution for nu−1 pulses, renormalizing survivors and accumulating their log probability. Solve (I−Q_post)g=1 for the expected number of post-change samples to the first alarm. The conditional delay is survivor_distribution @ g − 1, giving zero delay for detection on the first changed pulse. No restart after an early alarm is allowed. This method applies to scalar iid coherent/squeezed laws; it is not applied to correlated n>1 increments.

The saved source thresholds are checked with 512 and 1024 cells by default. Agreement requires |MC−fine| <= 4 MC SE + |fine−coarse| and mesh change below 0.05 pulse. Grid differences are empirical refinement checks, not rigorous discretization bounds. Source pre-change false-alarm fractions are also compared against deterministic survival probabilities in deterministic_checks.csv.

## Counterfactual thresholds: diagnostics only

Each inverse calculation fixes N=100, eta1=0.9, eta2=0.85, the probe model and change pulse 1005. A root finder chooses h to match the author point mean on a 256-cell delay grid; 512- and 1024-cell evaluations then inspect that SAME h. The author mean has unknown uncertainty, so this inverse fit has no claimed confidence interval for h. These fitted thresholds are deliberately not calibrated to two-million-pulse ARL and are never written into production/source results.

Each fitted threshold is tested on 50,000 fresh changed runs and two independent adaptive-importance-sampling ensembles of 500,000 cycles per law. Raw audit cycles, ordinary importance-weight moments, raw changed stopping times and seeds are saved. Agreement of numerical and MC delay is checked using four MC SE plus observed mesh differences; ARL is checked against a scalar numerical first-step solution and between independent repetitions. Weight diagnostics retain ESS>=1000, largest contribution<=1%, relative SE<=2%. All intervals are pointwise approximate Monte Carlo intervals conditional on fixed fitted h, not simultaneous bounds or uncertainty in author means.

| Probe | Na | Original h | Diagnostic h | Author mean | New MC delay ± SE | Diagnostic ARL ± SE |
|---|---:|---:|---:|---:|---:|---:|
| Coherent | 0.1 | 11.933398 | 10.225939 | 64.6240 | 64.6985 ± 0.1280 | 357,298 ± 1,960 |
| Squeezed | 0.1 | 12.246435 | 9.888712 | 36.8382 | 36.8875 ± 0.0762 | 187,147 ± 943 |
| Coherent | 0.3 | 11.954340 | 10.449967 | 66.0598 | 65.8667 ± 0.1306 | 447,190 ± 2,451 |
| Squeezed | 0.3 | 12.426099 | 10.200655 | 27.3914 | 27.3747 ± 0.0567 | 216,563 ± 1,042 |
| Coherent | 0.5 | 11.943558 | 10.940937 | 69.3478 | 69.3616 ± 0.1331 | 734,575 ± 4,041 |
| Squeezed | 0.5 | 12.491721 | 10.258629 | 22.6968 | 22.8218 ± 0.0477 | 212,185 ± 999 |

## Phase probe

Na=0.1, n=8, frozen h=12.450826934. Across change pulses 1001..1008, mean conditional delay ranges from 25.004 to 25.917 pulses. The author mean at pulse 1005 is 19.9756. This probe changes the physical change time and mixed block accordingly, while retaining the same ARL-calibrated threshold. It does not assume the paper used the same one-based indexing convention.

The convention table compares different estimands: conditional delay, inclusive (+1), shifted (−1), signed all-first-alarm times, and clipping early alarms to zero. Signed/clipped means are diagnostics and are not replacements for conditional delay. Event censoring and no-survivor cases are refused.

Next: obtain the original simulator or per-setting threshold/ARL tables, or test explicitly documented alternative optical/receiver assumptions. No messages to authors are sent by this workflow.
