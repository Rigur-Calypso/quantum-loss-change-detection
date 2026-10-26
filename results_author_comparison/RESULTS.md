# Direct comparison with the authors’ supplemental dataset

Source: Guha, Saikat and Basu, Prithwish, [Zenodo record 17460928](https://doi.org/10.5281/zenodo.17460928), published 2025-10-28, version v1, CC BY 4.0. The four downloaded ASCII/Octave plotting files match the published MD5 checksums; SHA-256 hashes and original metadata are preserved in data/zenodo_17460928. The MATLAB/Octave script is read as source material, never executed.

## What matches and what does not

All eight published Figure 3 analytical KL ratios match our Gaussian formulas to maximum absolute error **3.95e-14**. This independently verifies the optical/statistical formulas at the published coordinates. It is not a validation of detection delays.

The published latency means and our independently calibrated simulations differ. Values below are descriptive numerical differences, not combined statistical tests: author per-point standard errors and raw Figure 3/7 stopping times are not supplied. We do not tune thresholds, optical parameters or time offsets to match the reported means.

| Figure 7 detector | Mean relative difference (ours/author−1) | Maximum absolute difference | Points within ±5% | Our confirmed ARL settings |
|---|---:|---:|---:|---:|
| Coherent | 13.65% | 18.28% | 0/9 | 9/9 |
| Squeezed n=1 | 25.89% | 32.51% | 0/9 | 9/9 |
| Entangled n=2 | 34.39% | 41.64% | 0/9 | 9/9 |
| Entangled n=4 | 30.60% | 40.27% | 0/9 | 9/9 |
| Entangled n=8 | 28.63% | 34.34% | 0/9 | 9/9 |

Figure 3 unmodulated delay ratios: 3/8 within ±5% numerical difference; mean absolute relative difference 7.49%; 8/8 of our pairs pass ARL confirmation.

## Figure 7 comparison at Na=0.1

| Detector | Authors’ delay | Our conditional delay ± SE | Relative difference |
|---|---:|---:|---:|
| Coherent | 64.6240 | 76.4343 ± 0.1817 | 18.28% |
| Squeezed n=1 | 36.8382 | 46.2874 ± 0.1101 | 25.65% |
| Entangled n=2 | 30.7980 | 38.7986 ± 0.0930 | 25.98% |
| Entangled n=4 | 24.8668 | 31.1004 ± 0.0732 | 25.07% |
| Entangled n=8 | 19.9756 | 25.9746 ± 0.0588 | 30.03% |

## Matched inputs and explicit assumptions

- Figure 7 Na=0.1:0.05:0.5 is read directly from the Octave double_range; the column-to-block mapping comes from blocksizes=[1,2,4,8] and the supplied plot script. N=100, change pulse 1005 and ARL0=2,000,000 are in the data file.
- Figure 3 contains eight KL/latency ratios, N=100, r_min=0.1, N_a=1, cp=1000 and ARL0=2,000,000. Its r coordinates are inferred from the shared supplemental r_range (0.1:0.1:asinh(sqrt(1)), ending at 0.8). This mapping is checked by reproducing every supplied KL coordinate to machine precision; Figure 3 does not explicitly save its energy array.
- eta1=0.9, eta2=0.85 are explicitly stated for Figure 3 in the paper and used throughout our comparison. They are not saved in the Figure 7 file. Figure 7 uses these same paper-context values as an explicit assumption, not independently verified per-file metadata.
- The saved n=5000 is treated as an observation horizon, consistent with the Figure 3 paper caption and plot titles; it is not assumed to be the number of Monte Carlo trials. Our new simulations use horizon 5000 and report their separate trial count.
- Our detector makes strict one-based Page crossings, retains prefix history through reflections, and generates the exact mixed-loss covariance for a straddling block. Delay is E[T−nu | T>=nu]; early first alarms are excluded and counted separately. All unrestricted estimates refuse relevant censoring.
- BPSK data are retained and shown for context, but this implementation compares unmodulated homodyne models. No BPSK reproduction is claimed. The Supplemental Figure 4b experiment is archived but is not simulated in this milestone.

## Independent simulation protocol

Each of the 45 Figure 7 settings and 16 Figure 3 scalar settings receives five fixed log-ARL calibration updates using 200,000 cycles per law. Two fresh holdout repetitions each use 500,000 cycles per law, batched at 100,000. Adaptive proposal offset b=2 is fixed. Each changed-run ensemble has 30,000 trials, a separate stream, and horizon 5000. Thresholds are frozen before holdouts and are never fitted to author latency values.

Confirmation requires both holdouts and their pool to pass ESS>=1000, largest weighted contribution<=1%, relative SE<=2%, repeat agreement within four combined SE, plus approximate target and same-energy coherent ARL-ratio 95% intervals fully inside ±5%. Any unconfirmed cases remain visible in CSV and JSON. Intervals are pointwise Monte Carlo uncertainty conditional on thresholds; they do not include model uncertainty or uncertainty in published author means. All raw changed-run arrays, audit cycle arrays, calibration histories, weighted batch moments, seeds and versions are retained.

## Unresolved differences and information needed

The dataset supplies plotting data and a plotting script, not the simulator or ARL lookup tables. It does not provide the detector thresholds, per-point ARL validation, raw first-alarm/event flags, seeds, author per-point SEs, or explicit Figure 7 transmissivities. These missing details prevent attribution of the latency gap to one cause. Agreement of KL formulas and our ARL diagnostics does not prove that either implementation is wrong.
A one-sample indexing convention changes a delay by one pulse; the observed coherent discrepancy is larger, and cannot be explained solely by block phase because that detector has block size one. No index correction is applied to reduce the gap.
The separate Supplemental Figure 4b raw latency vector contains 2 negative entries among 500 values, and its plotting script applies mean(lat). That is a different experiment; it does not establish how Figure 3 or 7 early alarms were processed. We record it only as a reason to request explicit first-alarm handling rather than infer it from summary means.
Next research step: obtain or reconstruct the original per-setting thresholds and simulator conventions, then test candidate explanations with independent ARL checks. Contacting the authors is not performed by this script.
