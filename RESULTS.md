# Delay versus false-alarm scaling

**Dual-rail logical-channel implementation:** [new Appendix E.2 results](results_dual_rail/RESULTS.md)
now include explicit state propagation, ideal heralding and two-bin timing.
Uniform loss preserves surviving-state coherence; a mid-pair change is treated
as differential rail loss. The study passes 92 checks with exact regeneration,
and the full suite passes 406 tests. The [current report](paper_reproduction_v5/REPORT.md)
includes this seventh core study. These logical-channel checks do not resolve
the earlier Figure 4/5 or published-latency discrepancies.

**DV uncertainty follow-up:** [new results](results_dv_uncertainty_final/RESULTS.md)
show a remaining 4.459-pulse conditional gap after separately accounting for
two negative source entries, under the complete-first-alarm interpretation.
Its bootstrap percentile interval is [2.665, 6.206]. This does not identify a
model or source-code error. The full suite passes 382 tests. The current
[consolidated report](paper_reproduction_v4/REPORT.md) links this ancillary
analysis and the public-source availability audit while preserving older reports.

**Guided implementation lesson:** the [executed notebook](walkthrough_run_v2/01_paper_walkthrough.ipynb)
now demonstrates measurements, Figure 2's formulas, fair thresholds, first-alarm
paths, independent repeated trials and censoring. It uses ARL 1,000 as a teaching
example, separately from the paper studies. The [current consolidated report](paper_reproduction_v3/REPORT.md)
links the lesson while retaining the previous scientific findings and gaps.

The [new Eq.1 study](results_delay_scaling/RESULTS.md) validates 14 scalar
settings at ARL targets from 100 to 100 million. High-ARL delay growth agrees
with the 1/KL slope within 0.014%, but finite absolute delays remain below the
first-order prediction. All new ensembles pass independent checks and 267
audit checks with exact regeneration; the project suite remains 379 tests.
The [updated consolidated report](paper_reproduction_v2/REPORT.md) includes this
sixth study. Earlier studies and the first report remain preserved.

# Independent equation validation

**Consolidated review:** [visual report](paper_reproduction/report.html),
[complete comparison](paper_reproduction/REPORT.md) and
[single-command verification guide](REPRODUCING.md) now collect the five core
paper studies. The checker preserves historical audit files and hashes the
recorded code/data/report evidence. Complete published numerical agreement
remains unresolved; the report keeps those gaps explicit.

Five million fresh draws and an independent adaptive integrator support the
existing Figure 4/5 equations. A roughly ten-photon Kennedy counterfactual
approximates the published pixels, including three points held out from a
one-point energy fit. This is a diagnostic clue; the stated-energy discrepancy
remains unresolved. The full suite passes 373 tests, and the new study passes
114 checks with exact sample regeneration. See [findings and limits](results_equation_validation/RESULTS.md).

# PDF-verified figure comparison

The subsequent [verified-setting study](results_verified_figures_v2/RESULTS.md)
reruns Figures 5/12 from their original PDF legends. Figure 12 has sparse pixel
agreement within 2% across seven curves at 1 dB. Figure 4's fixed-seed ordering
and an approximately twelve-fold Figure 5 Kennedy discrepancy remain unresolved.
The current suite has 368 passing tests, and the new saved tables pass 958 checks.
Historical milestone results below remain intact.

# Initial reproduction results

**Paper-specific completion milestone:** BPSK mixture likelihoods, mutual information, displaced single-photon homodyne detection, Bernoulli/Poisson receivers, phase-noisy homodyne/PNR Kennedy models, fixed-seed entanglement, independent Hadamard covariance propagation, ML change-point estimates and the Appendix F shared-path ARL estimator are implemented. All **17/17 new large-ARL cases** (16 BPSK plus DV) pass fresh independent confirmations; generic non-Gaussian and shared-path controls also pass. **359 tests and 382 data checks pass**, including regeneration of all new calibration/validation evidence. BPSK delay ratios differ from published values by **6.40% mean absolute relative error**; their agreement with our unmodulated ratios supports the high-SNR modulation comparison. The DV KL is **0.206340 nats** and its theoretical squeezing crossover **1.812 dB**. Our independently calibrated DV conditional delay is **54.640±0.129 pulses**, so absolute published-latency agreement is still incomplete. Figures 4/11a/12 model families and an explicit Eq.10–11 LO sweep are saved with CSV/PNG/SVG outputs. Exact Figure 5/12 plotting conventions remain to be verified, and published timing discrepancies stay visible. See [the complete coverage map](PAPER_IMPLEMENTATION.md), [new study results](results_paper_completion/RESULTS.md) and [Figure 11(b) model comparison](results_paper_completion/figure11b_comparison.png). Earlier robustness studies are separate extensions and do not substitute for reproduction.

**Conservative threshold milestone:** all five detector types have thresholds calibrated on five sampled true initial transmissions from 0.899 to 0.901 at Na=0.1, while assumed score laws remain 0.9/0.85 and true post-change transmission stays 0.85. Four fixed selection rounds target a minimum sampled ARL of 2.15 million. Independent frozen-threshold validation passes at **25/25 settings**: all multiplicity-adjusted approximate lower ARLs exceed two million, with the smallest **2,073,242 pulses**. Fresh delay comparisons show added mean delays of **3.61–4.39 pulses (coherent)**, **1.63–2.33 (independent squeezed)**, **1.46–1.94 (n=2)**, **1.31–1.61 (n=4)** and **1.12–1.31 (n=8)**. Thresholds, calibration histories and all independent validation/changed-run seeds are saved. **324 tests and 740 checks pass**, including exact regeneration of calibration summaries, weighted cycle arrays and both sets of first-alarm ensembles. See [the results and statistical limits](results_robust_thresholds/RESULTS.md), [protection/cost plot](results_robust_thresholds/protection_and_cost.png) and [validation CSV](results_robust_thresholds/validation.csv). Approximate Bonferroni-adjusted log-scale bounds give sampled-grid evidence; they do not prove continuum robustness or exact simultaneous coverage. This extension preserves the earlier reproduction and discrepancy studies.

**Entangled initial-channel robustness milestone:** at Na=0.1 and fixed source thresholds, fresh ARLs for true initial transmission 89.9%, 90%, 90.1% are **1.253, 1.984, 3.151 million pulses (n=2)**, **1.240, 1.978, 3.202 million (n=4)** and **1.218, 1.996, 3.289 million (n=8)**. Likelihoods still assume 90%/85%, and true post-change transmission stays 85%. Endpoint ARL changes relative to each fresh matched-null estimate are approximately **−37% to −39%** and **+59% to +65%**. The largest observed conditional mean-delay span across this grid is **0.200 pulse**; sampling errors are retained, so this is not a deterministic bound. Nine settings each have two independent 250,000-cycle ensembles per law and 30,000 changed runs. Six direct-simulation controls support the actual-null adaptive weighting. The correct path correction uses the actual null conditional density, while Page evidence still uses the assumed score; block histories persist through reflection and cycles renew only at block boundaries. **312 tests and 230 saved-data checks pass**, including exact regeneration. See [results and uncertainty](results_entangled_null_diagnostics/RESULTS.md), [plot](results_entangled_null_diagnostics/entangled_null_sensitivity.png) and [CSV](results_entangled_null_diagnostics/sensitivity.csv). This is a new robustness result, not a fitted reconciliation of published latencies.

**Initial-channel robustness milestone:** with original scalar likelihoods and thresholds fixed at assumed eta1=0.9/eta2=0.85, true eta1 values 0.899–0.901 produce fresh ARLs of **1,309,760–3,010,713 pulses for coherent light** and **1,268,083–3,089,618 for independent squeezed light**, at Na=0.1. True post-change eta2 stays 0.85. Relative to each matched-null numerical ARL, this is approximately −34%/+52% and −36%/+56%, respectively. Numerical conditional delay spans only 0.412 and 0.241 pulse across the grid. Thus this small initial-channel perturbation changes false-alarm performance strongly while barely changing latency. No author mean is fitted. Fresh actual-null renewal estimation uses the full path density ratio independently of the detector score; two repetitions per setting agree with a refined numerical calculation. Four moderate-threshold controls also agree with direct no-change simulation. **291 tests and 264 data checks pass**, including exact regeneration. The initial 512/1024-cell ARL mesh failed its fixed 0.1% tolerance; the preserved follow-up changes only the ARL grid to 1024/2048, retaining identical raw samples and achieving maximum mesh change 0.0337%. See [results and limits](results_null_channel_refined/RESULTS.md), [plot](results_null_channel_refined/null_channel_sensitivity.png) and [CSV](results_null_channel_refined/sensitivity.csv). Correlated blocks require a separate extension.

**Entangled channel-sensitivity milestone:** the scalar-fitted physical eta2=0.84519511 is frozen across all 27 correlated Figure 7 settings, without fitting any additional mean or changing any likelihood or threshold. Mean absolute relative differences fall from **34.39% to 13.26% (n=2)**, **30.60% to 10.34% (n=4)** and **28.63% to 11.04% (n=8)**. Nevertheless **0/27** entangled means falls within a descriptive ±5%; combined with the scalar experiment, **8/44 held-out means** do. This common loss change improves agreement without reconciling the curves. Each setting has two independent 30,000-trial ensembles; all repeats agree within four combined standard errors (maximum observed z=1.59), and no post-change survivors are censored. Exact within-block physical covariance is retained. Original ARL evidence is reused because the entire no-change stopping law stays fixed; no new ARL validation is claimed. **283 tests and 511 audit checks pass**, including dense-covariance oracle tests and exact regeneration of all 54 ensembles. See [the report](results_entangled_channel_diagnostics/RESULTS.md), [plot](results_entangled_channel_diagnostics/entangled_channel_predictions.png) and [combined comparison CSV](results_entangled_channel_diagnostics/combined_predictions.csv). Author errors are unavailable; these descriptive comparisons do not prove an author error or identify their actual channel.

**Channel-sensitivity milestone:** a one-point diagnostic fit changes only the true post-change transmissivity from 0.85 to **0.84519511**; detector likelihoods, thresholds and true pre-change laws remain fixed. The fitted coherent Na=0.1 mean is training data. **8/17 held-out scalar means** fall within a descriptive ±5% tolerance, so the common channel change reduces the gap without fully explaining it. All 18 numerical predictions agree with fresh 30,000-trial ensembles; no post-change survivors were censored. Maximum 512/1024-cell discrepancy is below 0.001 pulse. Original no-change ARL evidence is reused because its entire stopping law is unchanged, rather than claiming a new estimate. **256 tests and 246 audit checks pass**, including exact seeded regeneration. See [results and caveats](results_channel_diagnostics/RESULTS.md), [plot](results_channel_diagnostics/channel_predictions.png) and [prediction table](results_channel_diagnostics/predictions.csv). This scalar diagnostic does not establish the authors' parameters or cover correlated blocks.

**Gap-diagnostic milestone:** Timing-convention changes do not close the observed gap on any of the 45 Figure 7 points within ±5%; including signed early alarms changes our means by at most 0.355 pulse. A separate deterministic scalar delay calculation agrees with simulation for all 18 checked settings. In six scalar diagnostic fits, matching published means within our model lowers the independently validated ARL to 187,147–734,575 pulses instead of two million. A same-threshold block-phase probe at Na=0.1 also leaves a gap. All 246 tests and 90 diagnostic checks pass. See [the controlled diagnostics and limits](results_gap_diagnostics/RESULTS.md). This narrows tested explanations; it does not establish the original authors' thresholds or prove a source error.

**Direct author-data comparison (2026-10-22):** The Zenodo v1 dataset is retrieved and all four published source checksums verified. All eight published analytical KL ratios match to 3.95×10⁻¹⁴. Fresh simulations on the exact author grids pass ARL confirmation in 61/61 settings, but numerical latency agreement remains unresolved: Figure 7 mean relative differences range from 13.6% to 34.4% across detector families; Figure 3 unmodulated latency ratios have 7.49% mean absolute relative difference (3/8 within ±5%). All 229 tests and 185 saved-data checks pass. See [the direct comparison, metadata gaps and assumptions](results_author_comparison/RESULTS.md). Historical internal-validation milestones below remain valid for their own samples; they are not claims of exact author-latency reproduction.

**Adaptive follow-up milestone:** Both previously flagged n=8 settings pass independent validation at their unchanged sweep thresholds. New ARLs are 2,040,296 ± 2,841 pulses at Na=0.1 and 1,985,577 ± 3,399 at Na=5 (±1 Monte Carlo SE). Relative SE falls to 0.139% and 0.171%, respectively. With the original 23 confirmations, 25/25 sampled grid settings now have supporting evidence. All 209 tests and 31 follow-up data checks pass; original failed runs remain unchanged. See [the derivation, independent checks and limits](results_adaptive_validation/RESULTS.md) and [the combined energy plot](results_adaptive_validation/combined_energy_sweep.png).

**Energy/phase milestone:** The subsequent Figure 7-style sweep covers five Na values from 0 to 5, five detector types and eight change positions. 23/25 settings pass all confirmation checks; the zero-energy equivalence control passes. At Na=5 and change pulse 1005, confirmed n=4 delay is 10.524 pulses versus 73.194 for coherent light (6.955 times faster). Two n=8 settings remain flagged for unstable weight contributions in one holdout repetition each. All 185 tests pass. See [the sweep findings and limitations](results_conditional_sweep/RESULTS.md). Earlier experiments below remain historical evidence for their own thresholds and sample streams.

**Latest milestone (2026-10-21):** The correlated Appendix B.4 detector now has independently validated ARL calibration at two million pulses, for the five Figure 7 detector types at Na=0.1. All matching and repeat checks passed. At change pulse 1005, delays are 76.724 (coherent), 46.548 (independent squeezing), 38.804 (n=2), 31.050 (n=4), and 25.870 (n=8) pulses. All 153 tests pass. See [full results and limitations](results_conditional_large_arl/RESULTS.md). The older reports below retain their original scope; the full Na and block-phase sweep is next.

**Update:** The subsequent equal-ARL coherent/squeezed comparison is complete. See [the new comparison results](results_equal_arl/RESULTS.md). The numbers below document the original baseline experiment.

Source: [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2), Eqs. (2)–(3), Figure 2, and Appendix F.

## Analytical Figure 2 curves

- Classical CRE: 0.1428863093 nats/pulse at Na=0, rising to 0.1500306247 at Na=5.
- Independently squeezed CRE at Na=5: 1.0368931674 nats/pulse.
- At Na=0.1, squeezed/classical CRE ratio: 1.716219, consistent with the paper’s approximately 1.72.
- The analytical Gaussian KL was checked independently by integrating the log density ratio under the post-change density, including unequal variances.

## Initial classical Monte Carlo study

N=100, Na=0, eta1=0.9, eta2=0.85, seed=20261020. Each threshold uses 4,000 no-change runs and a separate 4,000 changed runs. Change at pulse 100; horizon 100,000. Uncertainties below are ±1 Monte Carlo standard error.

| h | No-change ARL (pulses) | Conditional delay (pulses) | Pre-change false alarms | Surviving changed runs |
|---|---:|---:|---:|---:|
| 2 | 71.32 ± 1.05 | 9.19 ± 0.24 | 75.33% (3013/4000) | 987 |
| 3 | 232.38 ± 3.50 | 15.64 ± 0.22 | 33.60% (1344/4000) | 2656 |
| 4 | 663.50 ± 10.03 | 21.71 ± 0.25 | 11.60% (464/4000) | 3536 |

No no-change or surviving changed runs were censored in these ensembles. Delays exclude pre-change alarms and condition on survival to pulse 100. At h=2, most trials already alarm before the change; its small conditional delay must therefore be read with the 75.33% early false-alarm rate. Increasing h trades slower detection for longer time to false alarm.

## Validation and limits

- 35 automated tests passed. Coverage includes independent Gaussian density/integration checks, Eq. (2)/(3) and energy allocation, exact CUSUM versus exhaustive suffix search, strict equality at the threshold, one-based time indices, replay of seeded observations, and false-alarm/censoring edge cases.
- Saved PNG figures were visually inspected. CSV coordinates, raw run/event arrays, and JSON metadata accompany them.
- Reproduced: the classical and independently squeezed Figure 2 analytical curves. This is not a full reproduction of Figure 2’s entanglement family.
- The Monte Carlo plots are an initial baseline validation, not Figure 3. The later equal-ARL comparison at target ARL=1000 is documented separately; calibration to ARL=2,000,000 remains planned.
- The moderate-ARL calibration and squeezed comparison are now complete for Na=0.1. Next: sweep parameters, extend the large-ARL solver and implement the entangled block model.


**Further update:** The 12-setting parameter sweep is complete. See [the sweep results](results_sweep/RESULTS.md). All settings passed the sweep’s 5% approximate ARL equivalence criterion, and the zero-squeezing control stayed close to unit gain.


**Large-ARL update:** A fast numerical solver and sampled Figure-3-style experiment at target ARL=2,000,000 are complete. See [large-ARL results and validation limits](results_large_arl/RESULTS.md). The rare-event ARLs are numerically calibrated with mesh checks and moderate-ARL Monte Carlo crosschecks; they have not been independently checked by rare-event Monte Carlo.


**Independent validation update:** The rare-event ARL check is complete: 24/24 moderate comparisons and 14/14 large-ARL comparisons passed, with all 7 large-ARL pairs meeting the 3% approximate equivalence criterion. See [regenerative importance-sampling evidence](results_rare_validation/RESULTS.md). This provides grid-independent weighted Monte Carlo evidence for the prior calibration; original numerical-only reports remain unchanged.


**Entangled-model update:** Full Figure 2 analytical curves and the Appendix B.4 conditional detector are implemented. See [entangled results](results_entangled/RESULTS.md). Six initial settings passed new direct-Monte-Carlo ARL matching at 1000 pulses. This does not transfer the earlier two-million-pulse iid calibration to the correlated conditional detector.
