# Entangled probe model and conditional CUSUM implementation

Source: [paper Eqs. (4)–(6), Appendix C.1–C.2 and Appendix B.4](https://arxiv.org/html/2503.12276v2#A2.SS4).

## Reproduced and implemented

- Full Figure 2 analytical curve family: ordinary coherent light, independent squeezing, and entangled blocks n=2,4,8,...,256, for N=100, eta1=0.9, eta2=0.85 and Na in [0,5]. The seed has n*Na photons. This is a formula-based reproduction, not a pixel/digitized-error comparison.
- Correlated full-block Gaussian sampling, covariance, joint log density, log-likelihood ratio and per-pulse KL. Dense multivariate calculations check the reduction to the collective statistic Y=sum(X)/sqrt(n). n=1 reduces to independent squeezing; zero squeezing gives the ordinary coherent information per pulse.
- The paper's Appendix B.4 pulse-by-pulse conditional-LLR detector. Conditioning history resets at block boundaries; Page evidence persists across blocks and reflects at zero after each pulse. A reflection does not erase the measured prefix. This variant is not claimed to satisfy iid-CUSUM optimality.
- Within-block physical changes use the exact mixed-loss covariance, not independent pre/post samples within an entangled block. Full-block decisions are also provided as a separate reference implementation; the benchmark below uses pulse-by-pulse decisions.

## Initial conditional-detector benchmark

Na=0.1; target ARL=1000 PULSES. Calibration uses fresh direct no-change simulations for each detector, followed by independent holdout ARL validation and independent changed-run ensembles. **6/6 detector settings passed** the predeclared 5% approximate ARL equivalence criterion against coherent light.

The entire pointwise approximate 95% ARL-ratio interval must fit within the equivalence band, and both point estimates must be close to target. Thresholds are never retuned on validation. Any unmatched detector remains visible and its comparison is exploratory only. No calibrated iid scalar two-million-pulse threshold is reused here: conditional increments are correlated, so those ARLs do not transfer.

| Detector | Threshold | Validated ARL ± SE (pulses) | Match passed |
|---|---:|---:|---|
| Coherent | 4.377441 | 987.58 ± 4.37 | True |
| Squeezed n=1 | 4.689941 | 1000.15 ± 4.41 | True |
| Entangled n=2 | 4.763184 | 999.91 ± 4.44 | True |
| Entangled n=4 | 4.831543 | 1005.90 ± 4.48 | True |
| Entangled n=8 | 4.870605 | 1005.55 ± 4.45 | True |
| Entangled n=16 | 4.860840 | 992.99 ± 4.39 | True |

### Change on pulse 1

| Detector | Phase in block | Delay ± SE (pulses) | Early alarms | Surviving runs | Coherent / detector delay |
|---|---:|---:|---:|---:|---:|
| Coherent | 1 | 27.116 ± 0.092 | 0.00% | 30000 | 1.000 |
| Squeezed n=1 | 1 | 17.462 ± 0.060 | 0.00% | 30000 | 1.553 |
| Entangled n=2 | 1 | 14.866 ± 0.051 | 0.00% | 30000 | 1.824 |
| Entangled n=4 | 1 | 12.425 ± 0.042 | 0.00% | 30000 | 2.182 |
| Entangled n=8 | 1 | 10.631 ± 0.032 | 0.00% | 30000 | 2.551 |
| Entangled n=16 | 1 | 11.696 ± 0.022 | 0.00% | 30000 | 2.318 |

### Change on pulse 101

| Detector | Phase in block | Delay ± SE (pulses) | Early alarms | Surviving runs | Coherent / detector delay |
|---|---:|---:|---:|---:|---:|
| Coherent | 1 | 24.229 ± 0.096 | 8.20% | 27541 | 1.000 |
| Squeezed n=1 | 1 | 15.816 ± 0.063 | 8.72% | 27383 | 1.532 |
| Entangled n=2 | 1 | 13.767 ± 0.054 | 8.91% | 27327 | 1.760 |
| Entangled n=4 | 1 | 11.534 ± 0.044 | 8.73% | 27381 | 2.101 |
| Entangled n=8 | 5 | 10.975 ± 0.037 | 9.01% | 27297 | 2.208 |
| Entangled n=16 | 5 | 10.793 ± 0.034 | 8.83% | 27350 | 2.245 |

### Change on pulse 1005

| Detector | Phase in block | Delay ± SE (pulses) | Early alarms | Surviving runs | Coherent / detector delay |
|---|---:|---:|---:|---:|---:|
| Coherent | 1 | 24.305 ± 0.158 | 63.78% | 10865 | 1.000 |
| Squeezed n=1 | 1 | 15.954 ± 0.101 | 63.38% | 10986 | 1.523 |
| Entangled n=2 | 1 | 13.695 ± 0.085 | 63.57% | 10929 | 1.775 |
| Entangled n=4 | 1 | 11.543 ± 0.068 | 63.28% | 11017 | 2.106 |
| Entangled n=8 | 5 | 11.020 ± 0.058 | 62.81% | 11157 | 2.206 |
| Entangled n=16 | 13 | 12.614 ± 0.054 | 63.69% | 10892 | 1.927 |

## Derivation and time conventions

K_eta=I/4+(v_eta-1/4)uu^T, where u=ones(n)/sqrt(n) and v_eta=[1-eta+eta*exp(-2s)]/4. Orthogonal modes have the same variance 1/4 under both hypotheses, so they cancel from the homogeneous joint LLR. The per-pulse CRE is [(v2/v1)-1-log(v2/v1)]/(2n)+N*(sqrt(eta2)-sqrt(eta1))^2/(2v1).

For k previously measured pulses, c=(v_eta-1/4)/n gives conditional mean mu+c/(1/4+k*c)*(sum(previous X)-k*mu) and conditional variance (1/4)*(1/4+(k+1)c)/(1/4+k*c). These match the general Gaussian conditioning equations. Unreflected conditional LLRs sum to the full block LLR; the per-pulse Page reflections make the conditional detector distinct from a complete-block detector.

If transmissivities differ within a block, the covariance is I/4+q*w*w^T, where q=[exp(-2s)-1]/(4n) and w_i=sqrt(eta_i). Sequential generation uses Gaussian conditioning for that actual covariance. The detector score compares homogeneous pre/post conditional hypotheses; it is not a maximization over all within-block change hypotheses.

Pulses are one based: block 1 is pulses 1..n, change nu makes pulse nu the first changed pulse, and delay is T-nu. The table states the phase explicitly. Indexing and ARL settings are therefore explicit rather than inferred from the caption of the published Figure 7. Delay conditions on T>=nu; early alarms are reported, not counted as successful detections. nu=1 checks initialization without early alarms. A partially observed final block is allowed for the conditional detector.

## Evidence and limitations

Calibration uses 20,000 runs/evaluation, holdout validation 50,000, and each delay ensemble 30,000. The horizon is 100,000 pulses. No calibration/validation/surviving changed runs were censored; unavailable unrestricted estimates cause refusal rather than silent omission.

The covariance_check.npz contains 100,000 independently sampled n=8 blocks, expected and empirical covariances. report.json records moment errors and dense/conditional likelihood checks. Separate model folders contain calibration histories, seed streams and raw validation/delay stopping times with event flags. CSV files retain plot coordinates and summaries. All uncertainty is pointwise approximate Monte Carlo uncertainty, conditional on chosen thresholds, not simultaneous coverage or optical-model uncertainty.

This completes the entangled Gaussian model, full Figure 2 analytical family, and an initial Appendix B.4 implementation. It is not a reproduction of Figure 7 at all of its settings, nor an independently calibrated two-million-pulse conditional-detector study. The earlier iid two-million-pulse results remain valid for their original detectors. Next: calibrate and independently validate the correlated conditional detector at the larger ARL, then sweep change phase and Na.
