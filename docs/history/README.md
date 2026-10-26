# Archived local project history

These milestones describe earlier study snapshots. Use the repository README for current scope and status.

# Transmission-loss quickest change detection

**Portable core copy:** [bundle instructions](../../REPRODUCING.md) and
[its checked report](../../REPRODUCING.md) provide
the seven core studies with relative working paths. The export preserves all
315 scientific array/table files and retains original metadata plus a migration
record. Validation blocks reads from the original project, except the installed
Python environment. This improves reproducibility; it adds no new scientific
evidence and leaves the published-figure differences unresolved.

**Appendix E.2 now includes explicit qubit states:** [dual-rail model and results](../../results_dual_rail/RESULTS.md) add arbitrary qubit density matrices, pure-loss Kraus maps, heralded erasures and CUSUM readout once per two physical time bins. Uniform loss preserves the surviving state; a change between rails is modeled separately. Two independent ARL holdouts support the two-million-bin target, and 92 checks regenerate all new calibration, validation and changed samples. The full suite now passes 406 tests. This is an ideal logical-channel simulation, with no hardware or QKD-security claim.

**Latest single-photon comparison:** the [uncertainty re-analysis](../../results_dv_uncertainty_final/RESULTS.md) keeps the author's signed mean (48.034 pulses) separate from a nonnegative re-summary (50.181, conditional interpretation assumed). Compared with our conditional mean 54.640, a 4.459-pulse gap remains; its empirical-bootstrap interval is about 2.665–6.206 pulses under the stated assumptions. This narrows the descriptive discrepancy but does not identify its cause. The full suite now passes 382 tests. The [source audit](../../source_audit/REPORT.md) documents the missing original simulator/thresholds and [clarification questions](../../source_audit/CLARIFICATION.md); no author has been contacted.

**Learn the implementation step by step:** open the [guided notebook](../../notebooks/01_paper_walkthrough.ipynb) or its [executed copy with plots and results](../../walkthrough_run_v2/01_paper_walkthrough.ipynb). It generates measurements, reproduces Figure 2's formulas, calibrates fair thresholds, plots first-alarm evidence, compares independent trials and demonstrates censoring. To run it without installing a notebook server, use `python -m lossqcd.walkthrough --output walkthrough_repeat` in the existing environment. This is a teaching experiment at ARL 1,000, with early false alarms explicitly counted; it is separate from the paper's two-million-pulse studies.

**New Eq.1 check:** fourteen coherent/squeezed settings now test how delay grows as false alarms become rarer. Independent ARL and delay simulations support the logarithmic growth rate; at high ARLs its slope differs from the predicted value by less than 0.014%, while the absolute first-order formula still has a material finite-target offset. All 14 settings and 267 exact-regeneration audit checks pass; the full suite remains 379 tests. See [scaling findings](../../results_delay_scaling/RESULTS.md), [plot](../../results_delay_scaling/delay_scaling.png) and [updated consolidated report](../../paper_reproduction_v2/REPORT.md).

**Start with the consolidated report:** [visual report](../../paper_reproduction_v5/report.html), [plain report](../../paper_reproduction_v5/REPORT.md) and [reproduction guide](../../REPRODUCING.md). One command now checks the seven core paper studies and creates a fresh report: `python -m lossqcd.reproduction_report --output paper_reproduction_repeat`. It preserves historical evidence and records source/output checksums. Equation agreement, sparse pixel agreement and unresolved published-data differences remain explicitly distinguished.

**Independent discrepancy checks:** five million fresh samples and a separate adaptive integrator support the implemented Figure 4/5 equations. A counterfactual Kennedy model with about 10 photons matches the four checked published pixels; a one-point energy fit predicts the other three within approximate digitization uncertainty. This is a diagnostic clue, not a change to the stated 100-photon model or proof of the authors' parameters. **373 tests and 114 new evidence checks pass**, including exact sample regeneration. See [findings and source questions](../../results_equation_validation/RESULTS.md) and [comparison plot](../../results_equation_validation/independent_checks.png).

**Latest verification:** Figures 5 and 12 have now been rerun using settings read directly from the original PDF. Figure 12's seven curves agree with approximate published readings at 1 dB within 2%. Two gaps remain explicit: Figure 4's fixed-seed ordering and Figure 5's Kennedy curve (about twelve times below our equation-based calculation). The full suite has **368 passing tests**; the new tables pass **958 audit checks**. See [verified results and limits](../../results_verified_figures_v2/RESULTS.md) and [plots](../../results_verified_figures_v2/verified_receiver_curves.png). This milestone verifies settings and exposes discrepancies; it does not claim complete numerical reproduction.

**Paper implementation status:** the principal paper models and algorithms now include coherent/squeezed/entangled homodyne detection, BPSK mixtures and mutual information, displaced single-photon probes, photon counting, imperfect homodyne/Kennedy receivers, fixed-squeezer comparisons, the optical Hadamard covariance, change-point estimation and Appendix F shared-path ARL computation. **359 tests and 382 new-study data checks pass.** All 16 new BPSK settings plus the independently calibrated DV setting pass fresh ARL validation. DV entropy is **0.206340 nats**, agreeing with the quoted 0.2063; its entropy crossover is **1.812 dB**, close to the paper's rounded 1.83 dB. Exact reproduction of all published numerical latencies remains incomplete. Start with [the full figure-by-figure coverage and remaining work](../../PAPER_IMPLEMENTATION.md), [new results](../../results_paper_completion/RESULTS.md), [BPSK/fixed-seed plots](../../results_paper_completion/bpsk_and_fixed_seed.png), [single-photon comparisons](../../results_paper_completion/photon_probe_comparisons.png) and [LO/communications models](../../results_paper_completion/receiver_and_communications.png). Robust-threshold and uncertain-channel studies below are separate research extensions.

**Research extension — conservative threshold calibration:** all five detector types now have separately calibrated thresholds for the sampled initial-channel uncertainty grid 0.899, 0.8995, 0.9, 0.9005 and 0.901 at Na=0.1. **25/25 independent validation settings** pass diagnostic and repetition checks, with multiplicity-adjusted approximate lower ARL estimates above two million pulses; the smallest is **2,073,242 pulses**. Fresh original/conservative delay ensembles show a cost of **3.61–4.39 pulses for coherent light**, **1.63–2.33 for independent squeezing**, and **1.12–1.94 for entangled detectors**, across this grid. **324 tests and 740 audit checks pass**, including exact regeneration of all calibration, holdout and delay samples. See [calibration results and limits](../../results_robust_thresholds/RESULTS.md), [protection/cost plot](../../results_robust_thresholds/protection_and_cost.png) and [validation table](../../results_robust_thresholds/validation.csv). These are empirical sampled-grid robustness results, not a guarantee for every channel between grid points or a reproduction of a new paper figure.

**Earlier entangled initial-channel robustness milestone:** the same 89.9%–90.1% initial-transmission check now covers n=2, 4 and 8 at their original thresholds. Fresh no-change ARLs span **1.22–3.29 million pulses**. Relative to each fresh matched-channel estimate, the endpoint shifts are approximately **−37% to −39%** and **+59% to +65%**. The largest observed mean-delay span is only **0.200 pulse**, with Monte Carlo uncertainty recorded separately. Nine settings each use two independent 250,000-cycle ensembles per sampling law and 30,000 changed runs. All six direct no-change controls agree with the weighted method. **312 tests and 230 data-audit checks pass**, including exact regeneration of all saved ensembles. See [correlated robustness results](../../results_entangled_null_diagnostics/RESULTS.md), [sensitivity plot](../../results_entangled_null_diagnostics/entangled_null_sensitivity.png) and [CSV](../../results_entangled_null_diagnostics/sensitivity.csv). No author mean or detector threshold is fitted.

**Earlier scalar initial-channel robustness milestone:** at Na=0.1, changing true pre-change transmission from 89.9% to 90.1% while the scalar detectors assume 90% shifts fresh no-change ARLs to approximately **1.27–3.09 million pulses**. Relative to each matched-channel calculation, the extreme changes are roughly **−36% to +56%**. Numerical conditional delays move by less than half a pulse across this range, so this initial-channel perturbation has a much stronger effect on false alarms than on the latency gap. True post-change transmission stays 85%; no published mean is fitted. Fresh renewal estimates, refined numerical solutions and four direct-simulation controls support the findings. **291 tests and 264 data checks pass**. The initial coarse-grid failure is preserved; refinement keeps every random sample and threshold identical. See [robustness findings](../../results_null_channel_refined/RESULTS.md), [sensitivity plot](../../results_null_channel_refined/null_channel_sensitivity.png) and [numerical table](../../results_null_channel_refined/sensitivity.csv). This milestone covers coherent and independent squeezed light only.

**Earlier entangled channel-sensitivity milestone:** the same loss value fitted from one classical point has now been tested at all 27 correlated Figure 7 settings, with two independent 30,000-trial ensembles each. Average absolute relative differences shrink from **28.6–34.4% to 10.3–13.3%** across the three block sizes, but **0/27** entangled means falls within ±5% of the published means. Including the previous scalar checks gives **8/44 held-out means within ±5%**, excluding the one fitted point. No additional mean is fitted. The source thresholds, likelihoods and true pre-change laws remain fixed, so original ARL evidence applies unchanged. **283 tests and 511 audit checks pass**, including exact regeneration of all 54 raw ensembles. See [findings and caveats](../../results_entangled_channel_diagnostics/RESULTS.md), [entangled overlay](../../results_entangled_channel_diagnostics/entangled_channel_predictions.png) and [all 45 comparison points](../../results_entangled_channel_diagnostics/combined_predictions.csv).

**Earlier scalar channel-sensitivity milestone:** fitting only the coherent Na=0.1 published mean gives a diagnostic true post-change transmissivity of **0.845195**, while the detector still assumes 0.85. This single parameter predicts **8/17 held-out scalar means within ±5%**, improving agreement without fully explaining the gap. All thresholds, assumed likelihoods and true pre-change laws stay fixed, so the original near-two-million-pulse ARL evidence applies unchanged; ARL is not re-estimated here. All 18 deterministic predictions agree with fresh 30,000-trial simulations. **256 tests and 246 data-audit checks pass**, including exact regeneration of every new raw ensemble. See [channel findings and limits](../../results_channel_diagnostics/RESULTS.md) and [prediction plot](../../results_channel_diagnostics/channel_predictions.png). This is a diagnostic fit, not an inference of the authors' actual channel. Correlated blocks remain outside this milestone.

**Earlier gap-diagnostic milestone:** controlled gap diagnostics are complete. Changing early-alarm handling moves our means by at most **0.355 pulse** across 45 settings; none of the tested conventions closes the numerical gap within ±5%. A deterministic first-step calculation corroborates the scalar Monte Carlo delays in **18/18 settings**. Six diagnostic threshold fits reproduce published mean latencies within our fixed model, but their independently checked false-alarm intervals fall to **187,147–734,575 pulses**, below the required two million. The eight-position block-phase probe also leaves a gap at Na=0.1. Original thresholds and source runs remain unchanged. All **246 tests** and **90 diagnostic data checks** pass. See [diagnostic findings and caveats](../../results_gap_diagnostics/RESULTS.md) and [threshold/false-alarm plot](../../results_gap_diagnostics/counterfactual_thresholds.png). These counterfactual fits do not establish the authors' thresholds or identify an error in their implementation.

**Direct author-data comparison:** all four published file checksums pass, and all eight analytical KL ratios match to 3.95×10⁻¹⁴. New simulations on the authors' exact grids pass independent ARL confirmation for 61/61 settings, while Figure 7 delays average 13.6%–34.4% higher in our simulation. Figure 3 unmodulated delay ratios have 7.49% mean absolute relative difference. The source dataset lacks detector thresholds, raw Figure 3/7 first-alarm arrays and author uncertainty estimates. See [the original comparison](../../results_author_comparison/RESULTS.md), [Figure 7 overlay](../../results_author_comparison/figure7_overlay.png) and [Figure 3 overlay](../../results_author_comparison/figure3_overlay.png). That experiment and its original 229-test/185-check milestone remain preserved.

**Earlier adaptive-validation milestone:** both flagged n=8 settings pass fresh independent validation with an evidence-dependent mixture sampler at their unchanged original thresholds. Together with the earlier 23 confirmations, 25/25 sampled grid settings have supporting validation. The original failed repetitions remain preserved. Estimated relative Monte Carlo errors for the two follow-ups are 0.139% and 0.171%, versus about 1% before, using fewer cycles. All 209 tests and 31 follow-up data checks passed at that milestone. See [follow-up results and derivation](../../results_adaptive_validation/RESULTS.md) and [combined energy plot](../../results_adaptive_validation/combined_energy_sweep.png).

The sampled Figure 7-style study covers **Na=0, 0.01, 0.1, 1 and 5**, five detector types and all eight change positions after 1000 pre-change pulses. Its original [energy/phase report](../../results_conditional_sweep/RESULTS.md) retains its 23/25 outcome; the two new confirmations are separate fixed-threshold follow-ups. At Na=5 and change pulse 1005, n=4 takes **10.52 pulses versus 73.19** for coherent light (6.95 times faster), while n=8 takes 13.24 pulses (5.53 times faster). The zero-energy equivalence control passes. These are five sampled Na values and approximate ARL-matching checks, not a direct comparison with the authors' supplemental dataset.

**Earlier single-energy milestone:** all five detector types passed separate large-ARL checks at Na=0.1, with n=8 about 2.97 times faster. Those thresholds, samples and results remain preserved in [the earlier experiment](../../results_conditional_large_arl/RESULTS.md). The new sweep is a fresh experiment and retains its failed diagnostic repetitions rather than substituting the earlier favorable run.

This version-2 folder continues the Desktop SOP implementation. The original Desktop files are preserved; see [project origin](../../PROJECT_ORIGIN.md). A local `.venv` with the pinned dependencies is installed here.


**Previous milestone:** independent rare-event simulation now checks the two-million-pulse false-alarm calibration. All 24 moderate-ARL checks and 14 large-ARL checks passed; all 7 ordinary/squeezed pairs passed the 3% ARL equivalence diagnostic. All **86 tests pass**. See [the latest validation results](../../results_rare_validation/RESULTS.md). The entangled-pulse model is the next extension.


A computational SOP project under Prof. Tiju Cherian John, implementing the first reproducible steps of [Guha, John, Gong & Basu, *Quantum-enhanced quickest change detection of transmission loss*, arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).

## Start here

To use the already-installed environment in this version-2 folder:

```sh
source .venv/bin/activate
python -m pytest -q
```

To reproduce the new paper-specific additions and audit all saved simulations:

```sh
python -m lossqcd.paper_completion --output results_paper_completion_repeat
python -m lossqcd.audit_paper_completion results_paper_completion_repeat --regenerate --regenerate-calibration
```

This reads the preserved author-comparison input and checksummed author data. It runs all eight BPSK Figure 3 pairs, independently calibrates and validates the non-Gaussian single-photon detector, checks moderate-scale controls, saves fixed-seed/weak-probe/LO/capacity curves, and demonstrates the Appendix F shared-path estimator. It preserves unknown-sign marginalization, exact non-Gaussian likelihoods, first-alarm handling, censoring and ordinary likelihood weighting. It never fits a threshold to a published latency. Models/plots use explicit recorded parameters; the LO sweep and weak-probe comparisons are not claimed as exact digitizations of Figures 5/12. Existing latency differences remain visible. The output folder must be new. See PAPER_IMPLEMENTATION.md for the distinction between implemented models and verified numerical reproduction.

The generic iid detector supports BPSK, Fock homodyne, Bernoulli and Poisson observations. The phase-mixture receiver models also provide observation likelihoods/sampling, with quadrature convergence checks. Photon-number-resolving and binary on/off Kennedy outputs remain distinct. An ideal null can have infinite relative entropy while requiring a nonzero post-change click waiting time. The separate correlated detector remains responsible for entangled block increments; the iid engine is never substituted for it.

To reproduce conservative threshold calibration and its independent validation (research extension):

```sh
python -m lossqcd.robust_thresholds --output results_robust_thresholds_repeat
python -m lossqcd.audit_robust_thresholds results_robust_thresholds_repeat --regenerate --regenerate-calibration
```

Four fixed calibration rounds evaluate every grid point and aim for a minimum sampled ARL of 2.15 million pulses, providing a buffer above the two-million requirement. Calibration uses 100,000 cycles per sampling law per evaluation. After thresholds are frozen, independent holdouts use two 250,000-cycle ensembles per law at each of 25 model/channel settings. All calibration, holdout, original-delay and conservative-delay seed slots are separate; the audit checks that separation. No validation result or author delay is used to fit a threshold. The original score laws, source thresholds and historical data remain preserved.

Validation requires weight diagnostics, repeated estimates within four combined standard errors, no censored changed-run survivors, and an approximate lower ARL above two million. The lower bounds allocate nominal 5% error across the 25 points using a one-sided Bonferroni adjustment with asymptotic log-scale SEs. This is not certified coverage or a continuum guarantee. Delay costs use independent 20,000-trial ensembles at original and conservative thresholds; difference SEs account for both ensembles. Threshold histories, seeds, raw holdout cycles, exact conditional path corrections, baseline/conservative first alarms, CSV and plots are saved. Calibration arrays are reproducible from recorded seeds. The output folder must be new. Away from the least favorable sampled channel, the conservative detectors intentionally have longer and unequal ARLs, so the curves should not be described as exact equal-ARL comparisons at every channel.

To repeat the earlier entangled initial-channel robustness experiment:

```sh
python -m lossqcd.entangled_null_diagnostics --output results_entangled_null_repeat
python -m lossqcd.audit_entangled_null_diagnostics results_entangled_null_repeat --regenerate
```

The detector's assumed channels remain eta1=0.9/eta2=0.85, true post-change eta2 stays 0.85, and only the true pre-change eta1 changes. Its conditional Gaussian law includes the correct physical covariance. Cycles renew only at a block boundary with Page evidence zero, or end at the first pulse alarm. The adaptive proposal mixes the **actual** null conditional law with the assumed post law; its complete mixture likelihood corrects the path independently of the assumed Page score. Both score and weight are saved. Raw numerator/denominator moments are pooled before their ratio is formed. Six h=3 controls use independent 10,000-trial direct simulations. Outputs retain all raw cycles, changed and control runs, seeds, weights, diagnostic repetitions, CSV and plots; the output folder must be new. No iid numerical solver is applied to correlated data. Approximate error bars are Monte Carlo diagnostics, not certified bounds.

`simulate_conditional(..., sampling_eta1=...)` now supports true initial-channel mismatch, including mixed physical blocks after a change. `estimate_conditional_arl(..., sampling_eta1=...)` supports the corresponding actual-null renewal law. Default matched-model arrays remain exactly reproducible. Tests compare the sampler with dense covariance calculations, mixture corrections with full conditional densities, scalar reductions with numerical first-step equations, and correlated estimates with direct simulation.

To repeat the earlier initial-channel scalar robustness experiment on the finer grid:

```sh
python -m lossqcd.null_channel_diagnostics --output results_null_channel_repeat
python -m lossqcd.audit_null_channel_diagnostics results_null_channel_repeat --regenerate
```

It fixes original source thresholds and assumed channel laws at eta1=0.9/eta2=0.85 while varying the **true initial** channel over 0.899, 0.8995, 0.9, 0.9005 and 0.901. The true post-change channel stays 0.85; this is separate from the earlier post-change fit. Each of the ten scalar settings uses 1024/2048-cell ARL solutions, 512/1024-cell delay solutions, two independent 200,000-cycle ensembles per sampling law, and 30,000 changed runs. Four h=3 controls also use direct no-change simulation. Outputs include all raw cycles, path likelihood weights, first-alarm arrays, seeds, CSV, plots and audit results. Output folders must be new.

The general scalar renewal estimator stores the assumed detector score separately from the correct path density ratio between the actual null law and the proposal. It does not reuse the matched-model exp(−score) formula under mismatch, clip weights, discard unfinished cycles, or substitute an old ARL estimate. Diagnostics and intervals are approximate. Correlated-block initial-channel uncertainty is not yet implemented.

The saved `results_null_channel_diagnostics` run retains a failed 0.1% ARL grid-refinement check (observed changes 0.124–0.135% for coherent cases). `results_null_channel_refined` repeats the same seeds, samples, thresholds and physical parameters with a finer ARL grid, whose maximum refinement change is 0.0337%. Its audit verifies every raw array is unchanged. To reproduce that explicit provenance trail into a new folder:

```sh
python -m lossqcd.null_channel_diagnostics --output results_null_channel_refined_repeat --refine-from results_null_channel_diagnostics
python -m lossqcd.audit_null_channel_diagnostics results_null_channel_refined_repeat --regenerate
```

To reproduce the earlier entangled extension using the preserved scalar fit:

```sh
python -m lossqcd.entangled_channel_diagnostics --output results_entangled_channel_diagnostics_repeat
python -m lossqcd.audit_entangled_channel_diagnostics results_entangled_channel_diagnostics_repeat --regenerate
```

The true post-change transmissivity stays at the previous diagnostic fit, 0.84519511. All entangled points are held out. The simulator keeps the original assumed eta1=0.9/eta2=0.85 conditional likelihoods, and generates the exact physical covariance across a block that straddles the change. Dense-matrix oracle tests independently check the streaming conditional generation and first-alarm timing. Two independent ensembles are pooled using raw stopping times; early first alarms and censoring stay explicit. The study records source checksums, seeds, raw arrays, repeatability diagnostics, CSV, plots and notes. The output folder must be new. Original source and scalar artifacts remain unchanged. The source dataset lacks author uncertainty estimates, so these numerical gaps do not establish a source error or the authors' actual channel parameters.

To repeat the earlier one-point channel diagnostic and regenerate/audit its raw simulations:

```sh
python -m lossqcd.channel_diagnostics --output results_channel_diagnostics_repeat
python -m lossqcd.audit_channel_diagnostics results_channel_diagnostics_repeat --regenerate
```

This reads the preserved author-comparison source, fits only one coherent training mean, then freezes that physical post-change loss for 17 held-out scalar settings. It changes both the squeezed-state mean and variance consistently. Detector likelihoods remain at eta1=0.9, eta2=0.85, with unchanged source thresholds. The true pre-change channel also stays at eta1=0.9. The numerical solver and direct simulator now accept separate physical sampling laws, with defaults preserving the original matched-model behavior. Outputs include CSV/JSON, independent seeded first-alarm arrays, two numerical meshes, a plot and checksum checks on original evidence. The output folder must be new. Author errors are unknown, so ±5% is a descriptive comparison, not a significance test. Figure 3's explicitly specified channel is not refitted.

To reproduce the controlled diagnostic experiments and audit their saved data:

```sh
python -m lossqcd.gap_diagnostics --output results_gap_diagnostics_repeat
python -m lossqcd.audit_gap_diagnostics results_gap_diagnostics_repeat
```

The study reads the preserved `results_author_comparison` report and raw runs. It compares timing conventions across 45 settings, solves a scalar first-step delay equation on two meshes for 18 settings, and performs six clearly labeled inverse-threshold diagnostics at Na=0.1, 0.3 and 0.5. Those inverse fits intentionally match latency instead of ARL; they are not production thresholds and do not infer the authors' actual settings. Fresh simulation and independent ARL estimation quantify their false-alarm cost. A separate n=8 experiment tests every block position at a frozen source threshold. The audit checks source hashes, raw statistics, numerical results and diagnostic likelihood weights. Output folders must be new. See [information needed to reconcile the source study](../../results_gap_diagnostics/REPRODUCTION_NOTES.md).

To repeat the direct author-data comparison and audit its outputs:

```sh
python -m lossqcd.author_comparison --output results_author_comparison_repeat
python -m lossqcd.audit_author_comparison results_author_comparison_repeat
```

The included `data/zenodo_17460928` folder preserves four small original text files, the source metadata and checksum manifest under CC BY 4.0; see its attribution README. No network access or MATLAB installation is needed for the comparison. The Octave ASCII parser reads numeric values and nested cells without executing the supplied plotting script. Figure 7 uses the nine author energy values from 0.1 to 0.5, the five detector types and change pulse 1005. Figure 3 uses eight scalar coherent/squeezed pairs at change pulse 1000; its energy mapping is inferred from the shared supplemental squeezing range and verified against every original KL coordinate. BPSK reference ratios are retained for context but are not simulated.

The full run uses five fixed ARL calibration updates, two independent holdouts and a separate 30,000-trial changed ensemble per setting, with observation horizon 5000. The saved author `n=5000` is treated as this horizon, consistent with the paper caption, not as a replicate count. The output contains author/our numerical comparison CSVs, overlays, calibration histories, weighted batch moments and seeds, raw audit cycles and all new changed-run arrays. Output directories must be new. Neither thresholds nor optical parameters are fitted to author latency means. The paper's Figure 3 eta values 0.9 and 0.85 are used; Figure 7 files do not independently record eta, so this is an explicit paper-context assumption there. Author per-point errors and first-alarm handling are unknown. Next: investigate the original thresholds and simulator conventions behind the latency gap.

To reproduce the new fixed-threshold follow-up and verify its saved data:

```sh
python -m lossqcd.adaptive_validation --output results_adaptive_validation_repeat
python -m lossqcd.audit_adaptive_validation results_adaptive_validation_repeat --regenerate
```

It reads the original sweep's thresholds and failed diagnostic records, fixes the proposal offset at 2, and runs two new independent holdouts plus fresh direct/weighted checks at h=2 and h=4.5. Only the ARL sampling method changes. Existing detector thresholds, physical models and changed-run stopping times remain fixed. The new method keeps a separate path likelihood correction for the adaptive mixture, including reflections and partial blocks; ordinary exp(−raw LLR sum) is reserved for the original post-only proposal. The raw audit files retain these two distinct quantities. No weights are clipped or self-normalized.

The output includes a source-hashed copy of each original delay ensemble, every weighted batch's moments and seed, raw audit batches, direct-check runs, a comparison plot, and a combined CSV whose `validation_source` column distinguishes the original sweep from the adaptive follow-up. The audit verifies source hashes, unchanged thresholds, raw likelihood summaries, direct/changed-run summaries, pooled moments and the combined CSV; `--regenerate` also recreates audit arrays exactly from recorded seeds. Output folders must be new. Importance-sampling diagnostics and pointwise intervals remain approximate evidence, not rigorous finite-variance or simultaneous-coverage guarantees.

To repeat the latest energy and phase sweep into a new folder, then audit its saved data:

```sh
python -m lossqcd.conditional_sweep --output results_conditional_sweep_repeat
python -m lossqcd.audit_conditional_sweep results_conditional_sweep_repeat
```

This generates 25 detector/energy settings and 200 changed-run summaries. Zero-energy labels explicitly share one scalar ARL calibration; all positive-energy settings use separate calibration and two independent holdouts. The saved sweep uses 148 million validation cycles across the two sampling laws. It retains every batch's moments, diagnostics and seeds, raw audit batches, and all changed-run stopping times. Outputs include plots, CSV, JSON and a report that keeps unconfirmed settings visible. The output folder must not already exist. Remaining raw renewal cycles can be regenerated from their recorded batch seeds and sizes.

Implementation: `phase_simulation.py` shares a pre-change prefix across phase branches while preserving first alarms and Page evidence, and simulates the exact mixed-loss block after each change. `batched_rare_event.py` pools raw first/second moments before forming the ARL ratio; tests compare this with concatenating complete raw arrays. Shared-prefix phases are correlated, so the reported phase ranges are descriptive, not independent-means significance tests. The later adaptive follow-up addresses the two originally unconfirmed n=8 settings. Next: compare with the supplemental Figure 7 data when accessible and test uncertain channel parameters.

To reproduce the latest large-ARL conditional experiment into a **new** folder:

```sh
python -m lossqcd.conditional_large_arl --output results_conditional_large_arl_repeat
```

It reads the preserved moderate-ARL direct simulation in `results_entangled`, runs independent rare-event calibration and validation, and saves new plots, tables, raw cycle/run data and seeds. It refuses to overwrite an existing output directory. The default comparison uses Na=0.1, n=1,2,4,8, two million-pulse target ARL and changes on pulses 1 and 1005. It samples tens of millions of short cycles, rather than simulating every two-million-pulse waiting time directly.

For a fresh installation instead:

From this folder in a terminal:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m lossqcd.experiments --output results
```

Python 3.10+ is supported by the source and dependency ranges; the saved experiment was tested on Python 3.14 on Apple Silicon. For the exact tested dependency versions, install `requirements-lock.txt` first (these pins may require a newer Python), then `python -m pip install --no-build-isolation --no-deps -e .`. No quantum hardware or quantum SDK is needed. The console command `lossqcd-experiments` runs the same experiment. The default experiment takes seconds to tens of seconds after dependencies and fonts are ready.

## What is reproduced

- General Gaussian log-likelihood ratio `log p2(x)/p1(x)`, including the log-normalization and quadratic terms for unequal variances.
- Gaussian `D(P2 || P1)` with independent numerical integration checks in the tests.
- Classical coherent homodyne model: `mean_j = sqrt(eta_j)*sqrt(N+Na)`, `variance_j = 1/4`.
- Independent squeezed homodyne model: `mean_j = sqrt(eta_j)*sqrt(N)`, `variance_j = [eta_j*exp(-2r)+(1-eta_j)]/4`, `r=asinh(sqrt(Na))`.
- Figure 2 classical dashed curve (Eq. 2) and independent squeezed solid curve (Eq. 3), with `N=100`, `eta1=0.9`, `eta2=0.85`, `Na` from 0 to 5. This is an analytical reproduction of two curves, not a digitized comparison or a reproduction of the entangled curves.
- Initial classical Monte Carlo ARL and conditional detection-delay study with independent random streams and explicit false-alarm/censoring accounting.

The equal total energy comparison is essential: the classical probe puts all `N+Na` photons into displacement, while the squeezed probe uses `N` for displacement and `Na` for squeezing. Eq. (3) uses the PRE-change variance in the mean-shift denominator. All relative entropies and likelihood ratios use natural logarithms (nats).

## Detector and statistical conventions

`W_0=0`, `W_t=max(0,W_(t-1)+log[p2(X_t)/p1(X_t)])`. Stop at the first **strict** crossing `T=inf{t>=1: W_t>h}`. Thresholds are nonnegative. The reflection at zero allows an empty suffix; it has the same positive threshold crossings as the paper's maximum over nonempty suffixes. A cumulative sum without reflection is not this detector.

Time indices are one based: the change at `nu` makes samples 1 through `nu-1` pre-change and sample `nu` the first post-change sample. A detection on `nu` has delay zero, matching the paper's `T-nu` convention. Add one to obtain the alternate convention counting post-change samples consumed.

Each simulation stops at its **first** alarm and never restarts it. ARL uses an independent no-change ensemble and estimates `E_infinity[T]`. Changed runs estimate `E[T-nu | T>=nu]`, reporting `P(T<nu)` and the number of surviving runs separately. This conditional mean at a fixed change point is not a worst-case Lorden delay, nor an unconditional performance guarantee. Low thresholds can give a small conditional delay while generating many false alarms.

At the finite horizon `H`, stopping times and event flags are stored separately. A genuine alarm at `H` counts as an event; an unobserved alarm time is censored. If any relevant runs are censored, the unrestricted ARL/delay is reported as null, never as the mean of the observed detections. Restricted means `mean(min(T,H))` and `mean(min(T,H)-nu | T>=nu)` are still reported as lower-bound estimates. If all changed runs alarm early, conditional delay is unavailable. Standard errors are Monte Carlo sampling errors, not rigorous confidence bounds or physical modeling uncertainty.

The exact physical assumptions are ideal real-quadrature homodyne detection, independent pulses, known pre/post transmissivities, aligned squeezing, and no excess receiver noise. The simulator samples observations from the Gaussian model and computes their likelihood ratios; it does not replace observations with their mean likelihood increment.

That independence statement applies to the original coherent/independently squeezed simulators. In the entangled extension, blocks are independent and pulses within each block are correlated. The latest ARL estimator renews only at a block end with Page evidence zero. It weights successful post-proposal paths by the entire **raw conditional** likelihood sum, including increments before reflections and any partial final block. See `src/lossqcd/conditional_rare_event.py`; the experiment and its validation are in `src/lossqcd/conditional_large_arl.py`. Heavy importance weights limit larger blocks: n=16 pilots failed the diagnostic screen, and are recorded separately rather than included in the confirmed comparison. Approximate Monte Carlo intervals and diagnostic agreement do not constitute rigorous error bounds.

## Files and outputs

- `src/lossqcd/models.py`: Gaussian distributions and paper Eqs. (2)–(3).
- `src/lossqcd/cusum.py`: Page statistic and one-based stopping time.
- `src/lossqcd/simulation.py`: seeded ensembles and censor-aware summaries.
- `src/lossqcd/experiments.py`: runnable research script.
- `tests/`: density/integration checks, paper parameter checks, exhaustive suffix comparison, strict crossing, deterministic simulation replay, and edge cases for false alarms/censoring.
- `results/figure2_classical_squeezed.png` and `.svg`: reproduced curves.
- `results/figure2_data.csv`: all analytical plot coordinates.
- `results/classical_monte_carlo.png`: ARL, conditional delay and false-alarm fraction, with ±1 Monte Carlo SE bars.
- `results/classical_trace.png`: one seeded CUSUM trajectory. The statistic continues after crossing solely to illustrate its behavior; the actual detector stops.
- `results/monte_carlo_summary.csv`, `report.json`: estimates, parameters, seed, platform and dependency versions.
- `results/monte_carlo_runs.npz`: per-run stopping times and event flags, indexed by thresholds `[2,3,4]` as `arl_0_stop`, `arl_0_event`, etc.
- `RESULTS.md`: numerical findings and the scope of this initial reproduction.

Defaults: seed `20261020`, 4,000 runs per threshold per ensemble, `nu=100`, horizon 100,000, classical simulation `Na=0`. Streams come from `SeedSequence(seed).spawn(7)`: streams 0/1 serve h=2 ARL/delay, 2/3 h=3, 4/5 h=4, and 6 the trace. Reproducibility requires the same settings and compatible NumPy version. Changing trial count or early stopping changes random-number consumption. Running again into the same output folder intentionally replaces generated results; choose `--output results_new` to preserve a run.

For a larger study:

```sh
python -m lossqcd.experiments --trials 10000 --horizon 200000 --change-at 100 --seed 123 --output results_new
```

## Implementation milestones

1. **Completed for Na=0.1, target ARL=1000:** separate Monte Carlo threshold calibration and independent validation, with a 2% ARL equivalence criterion. See `results_equal_arl/RESULTS.md`.
2. **Completed at that setting:** coherent/squeezed conditional-delay comparisons for changes at pulses 1, 100 and 1000. The subsequent 12-setting sweep over Na and target ARL is complete; see `results_sweep/RESULTS.md`.
3. Implement Appendix F shared-trajectory threshold crossing or a Fredholm ARL solver for rare alarms. Validate at moderate ARLs before the paper's two-million-pulse setting.
4. Implement correlated n-mode entangled probes and their per-pulse CRE (the remaining Figure 2 curves), respecting block length and change placement within a block.
5. Reproduce Figure 3 only after threshold calibration and large-ARL validation; the current study does not claim this reproduction. Excess receiver noise and unknown change magnitudes are later extensions.

## A useful first session

Run the tests, open the two main plots, then read the paper's “Quickest change detection,” “Baseline problem setup,” and “Squeezing-augmented transmitter” sections alongside `models.py`. Change the seed and compare shifts to the reported Monte Carlo standard errors. Then raise the threshold and observe both the ARL increase and the decrease in early false alarms.


## Equal-ARL comparison: completed

The new experiment compares ordinary and independently squeezed light at the same total energy, N+Na=100.1 photons/pulse. It calibrates separate thresholds to target ARL=1000, validates them on independent no-change runs, then measures detection delays on fresh changed runs. The whole approximate 95% ARL-ratio interval must lie inside a predeclared ±2% equivalence band; matching point estimates alone is insufficient. This establishes approximate equivalence, not exact equality of the underlying ARLs.

**Result:** validated ARLs 998.87±3.12 and 993.88±3.10 pulses (±1 MC SE). With a change at pulse 100, conditional delay falls from 24.33 to 15.80 pulses, about 35% shorter with squeezing. Early-alarm rates are 8.03% and 8.43%, separately reported. The comparable gains at change points 1 and 1000 are in `results_equal_arl/RESULTS.md`. No calibration, validation or changed runs were censored.

All **69 tests** pass, including added calibration and squeezed-simulation checks. Existing baseline results remain in `results/`; new artifacts are in `results_equal_arl/`:

- `equal_arl_comparison.png` and `.svg`: ARL, conditional delays and gain plot.
- `RESULTS.md`: estimates, uncertainty and limitations.
- `report.json`: settings, streams, model parameters, calibration history and summaries.
- `raw_runs.npz`: per-run validation/delay stopping times and event flags.
- `delay_summary.csv` and separate calibration JSON files.

To reproduce with the installed environment, use a fresh output directory:

```sh
source .venv/bin/activate
python -m pytest -q
python -m lossqcd.comparison --output results_equal_arl_new
```

Defaults: seed 20261022, 100,000 calibration runs per threshold evaluation, 100,000 independent ARL validation runs per method, and 50,000 delay runs per method/change point. No new dependencies are required. This experiment refuses an existing output folder, preserving previous studies. For a quick exploratory run, reduce `--calibration-trials`, `--validation-trials` and `--delay-trials`; lower precision may fail the equivalence check and is reported as such.

The observed delay ratio is about 1.54, while the KL ratio is about 1.72. The KL ratio describes a large-ARL asymptotic relationship; it is not an exact delay prediction at ARL=1000. This completes a first moderate-ARL comparison, not the paper's Figure 3 at ARL=2,000,000 or the entangled-probe implementation.


## Parameter sweep: completed

The new `lossqcd.sweep` experiment tests Na in [0, 0.01, 0.1, 1] and target ARL in [500, 1000, 5000], with N=100, eta1=0.9, eta2=0.85. All 12 settings passed the prespecified approximate ARL matching check. The sweep uses a 5% equivalence band with independent validation; the earlier single-setting result retains its tighter 2% check. The same total energy N+Na is used within each ordinary/squeezed comparison.

Zero squeezing is a control: both models use exactly the same threshold, but independently generated observations. Its observed delay gain stays near 1. The nonzero settings show progressively larger gains: for changes at pulse 100, delays are about 13–14% shorter at Na=0.01, 33–37% shorter at Na=0.1, and 66–69% shorter at Na=1 across the three target ARLs. Changes at pulse 1 are also included. Delay means condition on no earlier alarm; false-alarm fractions are recorded separately.

All **69 tests pass**. New evidence is in `results_sweep/`: `RESULTS.md`, `summary.csv`, `report.json`, `parameter_sweep.png` and `.svg`, plus per-condition reports and raw validation/delay stopping times and event flags. Summaries were independently checked against the saved arrays. No relevant runs were censored. Original baseline and single-setting experiment outputs are preserved.

To reproduce into a fresh output directory:

```sh
source .venv/bin/activate
python -m pytest -q
python -m lossqcd.sweep --output results_sweep_new
```

Defaults: seed 20261023, 30,000 calibration trials/evaluation, 50,000 validation runs per method/setting, and 30,000 delay runs per method/change point. Calibration and validation streams are independent; validation is never used for retuning. Each condition has a reproducible seed derived from the master seed. No additional dependencies are required. This broader run takes longer than the initial experiments.

These are pointwise approximate Monte Carlo intervals, not simultaneous confidence intervals across all settings. The large-ARL trend is suggestive, not a proof of asymptotic convergence. The paper's ARL=2,000,000 experiment and entangled probes remain planned; the next implementation work is efficient large-ARL calibration and the correlated entangled model.


## Large-ARL numerical experiment: completed with validation limits

A new first-step integral-equation solver (`arl_solver.py`) avoids repeated long no-change simulations when calibrating rare alarms. It handles both equal-variance and quadratic unequal-variance Gaussian likelihood increments and includes the Page statistic's point mass at zero. Its numerical estimates passed all 24 checks against the saved moderate-ARL Monte Carlo validation ensembles.

We then ran a **sampled Figure-3-style experiment** at target ARL=2,000,000, with N=100, eta1=0.9, eta2=0.85, Na in [0, 0.02, 0.05, 0.1, 0.25, 0.5, 1], change at pulse 1000 and horizon 5000. Each method/setting uses 50,000 changed runs. An additional change-at-pulse-1 ensemble checks initialization. Same total energy N+Na is used within each comparison.

At Na=0.1, conditional delay falls from about 76.90 to 46.39 pulses, a gain of 1.658 (about 39.7% shorter delay), closer to the asymptotic KL ratio of 1.716. Zero squeezing gives unit gain within sampling error. Thresholds calibrated on 1024 cells were checked on 512 and 2048 cells at the SAME thresholds; the largest base-to-fine ARL change was 0.03241%. Grid refinement is an empirical check, not a rigorous error bound.

**The two-million-pulse ARLs are numerically calibrated, not independently verified by rare-event no-change Monte Carlo.** Appendix F uses a different Monte Carlo calibration method. We have not established identical numerical results to the full published scatter plot. Finite-threshold delay ratios remain below the asymptotic KL ratios, and their uncertainty bars exclude discretization uncertainty. This distinction is preserved in all result notes.

All **69 tests pass**. New outputs are in `results_large_arl/`: `RESULTS.md`, `report.json`, `summary.csv`, `large_arl_comparison.png`/`.svg`, `solver_validation.png`, historical crosschecks and per-condition raw delay runs. All 28 delay summaries were checked against their raw arrays, and no changed-run ensemble was censored. Prior results are preserved.

SciPy is needed for this solver. It is already installed in this project's environment via the development dependencies; minimal new installations can use `python -m pip install -e '.[numerics]'`.

To reproduce into a fresh folder (the default benchmark file is the saved parameter sweep):

```sh
source .venv/bin/activate
python -m pytest -q
python -m lossqcd.large_arl --output results_large_arl_new
```

The numerical-method derivation is in the new results notes and source docstrings. Threshold search histories, all mesh checks, software versions, seeds and benchmark-file hash are recorded. The RNG seed is 20261024. Next: independent rare-event ARL validation, followed by the entangled multi-pulse model. The older milestone list above records the original sequence; this section records the current large-ARL progress and its limits.


## Independent rare-event validation: completed

The new `rare_event.py` module estimates ARL without a numerical grid. It splits Page CUSUM into independent cycles ending at a return to zero or a strict alarm-threshold crossing. Mean cycle length is estimated under the pre-change law. A separate post-change proposal ensemble estimates the rare probability of an alarm in one cycle, weighted by exp(-terminal log-likelihood sum). The terminal overshoot is retained, failed proposal cycles contribute zero, and weights are averaged over ALL cycles. The exact renewal identity is ARL=mean cycle length / cycle alarm probability. The ratio estimator is consistent, not exactly unbiased at finite sample size.

`rare_validation.py` freezes all previous thresholds and uses new independent streams. All **24** checks against historical direct Monte Carlo and **14** checks against the large-ARL numerical references passed the declared tolerance max(2%, 3 combined Monte Carlo SE / reference). All **7 pairs** passed a separate 3% equivalence diagnostic based on their independent ARL estimates and pointwise approximate 95% ratio intervals. No failed checks were dropped, and thresholds were not retuned.

For Na=0.1, independent ARL estimates are 2,004,738±6,930 pulses for coherent light and 2,002,168±5,933 for independently squeezed light (±1 Monte Carlo SE). This strengthens the earlier delay comparison of 76.90 versus 46.39 pulses. It does not assert exact equality of true ARLs or validate physical-model assumptions.

The earlier large-ARL results describe their original numerical-only evidence. The new validation supersedes that earlier limitation through weighted regenerative Monte Carlo, not through brute-force full no-change sequences. Original reports, flags, plots and raw experiment outputs remain preserved for provenance.

New artifacts: `results_rare_validation/RESULTS.md`, `report.json`, `summary.csv`, `rare_event_validation.png`/`.svg`, and per-check JSON plus compressed raw cycle arrays. All 38 estimates were independently recomputed from those raw arrays. No cycles were unfinished; the estimator refuses to discard unfinished cycles or silently count them as failed alarms. Raw evidence occupies about 152 MB.

Defaults: seed 20261025; 100,000 pre and post cycles per moderate-ARL check; 500,000 each per large-ARL check. Every setting uses independent pre/post proposal streams. Cycle counts are fixed in advance. The estimate records actual pulse counts, rare probability, successful proposal count, effective weight sample size and uncertainty. This module adds no dependencies.

To reproduce into a fresh output folder:

```sh
source .venv/bin/activate
python -m pytest -q
python -m lossqcd.rare_validation --output results_rare_validation_new
```

All **86 tests pass**, including exact cycle boundaries, overshoot weights, zero-squeezing-related baseline laws, geometric ARL at zero threshold, independent direct-simulation checks, unfinished-cycle rejection, underflow-resistant uncertainty and reproducibility. The full derivation and a primary method reference are in the new results notes. Confidence intervals are pointwise asymptotic approximations; they exclude model uncertainty and do not establish simultaneous coverage across all checks.

Next: implement the paper's correlated n-mode entangled probe model, starting with its covariance, multivariate likelihood and Figure 2 per-pulse relative entropy before a sequential block detector.


## Entangled model and conditional detector: implemented

`entangled.py` implements the paper's correlated homodyne measurements. For n pulses, the squeezed seed has n*Na photons, so s=asinh(sqrt(n*Na)). The covariance has one collective eigenmode with variance [1-eta+eta*exp(-2s)]/4 and n-1 orthogonal modes with variance 1/4. The mean of each pulse is sqrt(eta*N). This provides full Gaussian sampling, density, joint likelihood and KL/n, together with the exact collective-statistic reduction. Default pair-model energy Na=0 matches the existing pair APIs; experiments specify Na=0.1 explicitly.

All Figure 2 analytical curves are reproduced in `results_entangled/figure2_entangled.png` and `.svg`, with data in `figure2_all_curves.csv`: coherent, independently squeezed, and n=2,4,...,256 entangled probes. This is a formula-based reproduction rather than a digitized comparison to published pixels.

`conditional_detection.py` implements Appendix B.4: condition each pulse's likelihood on the measured prefix of its block. Reset that conditioning history at block boundaries while retaining the Page statistic across blocks. Reflect the Page statistic at zero after each pulse without erasing the measured prefix. The true generator uses the exact mixed-loss covariance if a physical transmission change occurs inside a block. This variant is not claimed to have the iid Page detector's optimality guarantee.

The initial benchmark uses Na=0.1 and **new direct-Monte-Carlo calibration** to target ARL=1000 pulses, with independent holdout validation. Six detector settings passed the predeclared 5% ARL equivalence check: coherent, independent squeezing, and entangled n=2,4,8,16. At a change on pulse 101, conditional delays are 24.23 pulses (coherent), 15.82 (independent squeezing), 13.77 (n=2), 11.53 (n=4), 10.98 (n=8), and 10.79 (n=16). Early-alarm rates are about 8–9% and are reported separately. Changes on pulses 1 and 1005 show initialization and block-phase effects; delay gain need not increase monotonically with n.

The complete-block detector is a separate reference in `block_detection.py`, with decisions only at block ends and explicit pulse-time conversion. Its minimum-ARL feasibility check uses stable Gaussian tails. It is not the detector used in the main benchmark and is not a substitute for the paper's per-pulse conditional version.

**The old scalar iid ARL calibration and rare-event checks cannot be transferred to conditional correlated increments.** The new conditional experiment is a moderate-ARL implementation check, not a reproduction of all Figure 7 settings or a two-million-pulse conditional-detector study. These distinctions are recorded in the results and metadata; all prior results remain intact.

All **129 tests pass**, including dense Gaussian covariance/likelihood/KL checks through n=256, conditional chain-rule identities, seeded covariance checks, exact scalar and zero-squeezing reductions, phase-aware sampling, full-block time conventions, calibration-stream independence and censor rejection. All 6 ARL and 18 delay summaries were recomputed from raw run/event arrays. No relevant runs were censored. Dense/reduced and conditional/joint likelihood calculations agree to about 1e-13 in the saved model check.

Saved evidence in `results_entangled/`: `RESULTS.md`, `report.json`, full Figure 2 CSV/PNG/SVG, `conditional_summary.csv`, conditional-detector plots, covariance plots and a full 100,000-block sample archive, plus per-detector calibration histories, seeds and raw stopping times/event flags. No new dependencies are needed for the main model and experiment. The optional complete-block feasibility helper uses the already-installed SciPy numerics extra.

To reproduce into a fresh folder:

```sh
source .venv/bin/activate
python -m pytest -q
python -m lossqcd.entangled_experiments --output results_entangled_new
```

Defaults: seed 20261021, Na=0.1, target ARL=1000 pulses, 20,000 calibration runs/evaluation, 50,000 holdout runs/detector, and 30,000 changed runs per detector/change point. Horizon 100,000 pulses; change points 1,101,1005; benchmark entangled sizes 2,4,8,16. The analytic Figure 2 family separately fixes N=100, eta1=0.9, eta2=0.85 and covers Na in [0,5]. All output paths must be fresh to preserve previous studies.

Next: develop and independently validate large-ARL calibration for the conditional, correlated detector, then sweep augmentation energy and the exact phase of a within-block change.
