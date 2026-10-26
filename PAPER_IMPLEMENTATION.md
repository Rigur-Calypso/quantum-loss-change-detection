# Paper implementation coverage

Source: [Guha, John, Gong and Basu, arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).

The [consolidated report](paper_reproduction_v5/REPORT.md) and
[visual version](paper_reproduction_v5/report.html) now collect the main-paper
comparisons, validation outcomes and remaining source questions. The
[reproduction guide](REPRODUCING.md) documents a single command for checking
the seven core studies while preserving their original audit records.

For learning the code, the [guided notebook](notebooks/01_paper_walkthrough.ipynb)
and [executed copy](walkthrough_run_v2/01_paper_walkthrough.ipynb) now demonstrate
the whole scalar workflow at ARL 1,000. It uses the existing models and detector,
retains first-alarm arrays and explicitly illustrates early alarms and censoring.
This teaching example is separate from the seven core paper studies.

The implementation now covers the principal observation models and detection
algorithms in the main text and Appendices A–F. **Exact numerical reproduction
of the entire paper is not complete.** Formula checks, independently validated
simulations, comparisons with author data and exploratory extensions are distinct
below. Hardware schematics and open research questions are not presented as
implemented optical hardware or solved coding problems.

The latest [independent discrepancy study](results_equation_validation/RESULTS.md)
adds adaptive count-by-count integration and five million fresh simulation
draws. These support the existing Figure 4/5 equations. A roughly 10-photon
Kennedy counterfactual agrees with sparse published pixel readings; that is an
energy-convention clue, not a verified author parameter or a replacement for
the stated setting. The unresolved reproduction gaps below therefore remain.

## Main figures

**Main-text Eq.1:** the [finite-ARL scaling study](results_delay_scaling/RESULTS.md)
adds 14 coherent/squeezed settings with independent ARL and cold-start delay
checks, spanning targets 100 to 100 million pulses. High-ARL slopes agree with
1/KL within 0.014%, while finite absolute delays remain below log(ARL)/KL.
This supports the asymptotic scaling without claiming a worst-case theorem
or resolving the published absolute-latency gaps. Two fresh nu=1000 ensembles
keep early first alarms separate. Its 267 audit checks include exact sample
regeneration. The earlier five-study consolidated report remains preserved;
the updated six-study report has a new output directory.

| Paper item | Implemented computation | Evidence and remaining limits |
| --- | --- | --- |
| Figure 1: transmitter families | Coherent, independent squeezing and correlated entangled probe statistics | `models.py`, `entangled.py`; optical-network covariance independently checked. The drawing represents hardware, not a software benchmark. |
| Figure 2: per-pulse relative entropy | All coherent, squeezed and n=2,4,…,256 entangled curves | [Analytical figure and data](results_entangled/RESULTS.md). Formula reproduction with independent Gaussian/covariance checks; not pixel digitization. |
| Figure 3: detection improvement, with/without BPSK | Both unmodulated and BPSK simulations on all eight saved author grid points, with ARL near two million | [Author comparison](results_author_comparison/RESULTS.md), [BPSK addition](results_paper_completion/RESULTS.md). KL coordinates agree to machine precision. Mean absolute relative delay-ratio differences: 7.49% unmodulated, 6.40% BPSK. Exact latency agreement remains unresolved. |
| Figure 4: fixed seed squeezing | Formula-derived n=1,2,…,256 family, with Na=sinh²(s)/n, over the verified 0–18 dB axis | [Equation check](results_verified_figures_v2/fixed_seed_equation_check.png). Independent dense Gaussian KL and A/n+B identity agree. At fixed seed, equations give n=1 ≥ n=2 ≥ n=4; the published inset places n=2 above n=1. This discrepancy is unresolved. Earlier n=2-only ranking excluded n=1 and does not establish published curve agreement. |
| Figure 5: imperfect local oscillators | Eqs.7,10,11: phase-mixed homodyne and PNR Kennedy laws, plus on/off coarse-graining | [Verified-setting comparison](results_verified_figures_v2/RESULTS.md). PDF confirms N=100, Na=1, tap=0.9, epsilon=0.01, phase variance=0.01 (SD=0.1), loss 0–6 dB. PNR is explicitly specified by Appendix A Eq.10. Calculated Kennedy entropy is approximately 11.9–12.6 times the plotted values at four digitized points. Binary coarse-graining does not close the gap. |
| Figure 6: receiver schematic | Statistical models for homodyne, PNR and on/off receivers | `receivers.py`, `paper_models.py`. No optical/electronic hardware implementation claimed. |
| Figure 7: correlated detection latency | Conditional-LLR Page detector for all five detector types and nine saved author energies; exact mixed-loss blocks and pulse-level first alarms | [Full comparison](results_author_comparison/RESULTS.md). All 45 ARL checks support the internal model, but mean delays differ from author data by approximately 13.6–34.4% across families. Detector phase, conventions and threshold diagnostics remain preserved. |
| Figures 8–10: construction and entanglement network | Displaced squeezed states, BPSK mixtures, normalized Hadamard network and resulting covariance | Covariance propagation agrees with the analytic model within 2.8×10⁻¹⁷. A capacity-achieving coded transceiver or physical Green Machine is not supplied by this simulation. |
| Figure 11(a): discrete-variable augmentation | Eq.22 density, Eq.23 CDF, exact vacuum/Fock-1 mixture sampler, KL integration, CV/DV comparison | [Plot](results_paper_completion/photon_probe_comparisons.png). Uses stated rounded tap=0.944; the separate benchmark uses exactly eta1=0.9, eta2=0.85. Computed DV KL 0.206340 matches quoted 0.2063. |
| Figure 11(b): CV/DV latency | Independently calibrated DV detector, existing CV simulations and new BPSK simulations; all retain first-alarm/censoring information | [Comparison](results_paper_completion/figure11b_comparison.png). Our DV conditional delay is 54.640±0.129 pulses (one SE); the source plotting script uses a signed mean containing two negative latencies. Those statistics are labeled separately. Published latency matching is incomplete. |
| Figure 12: weak one-photon probes | Coherent homodyne, squeezed vacuum, entangled vacuum, Fock homodyne, Bernoulli photon counting, residual-nulling Kennedy models | [Verified-setting curves](results_verified_figures_v2/verified_receiver_curves.png). PDF verifies N=0, Na=1, n=8, tap=0.9, residual photons 0.0001/0.00001 and loss axis 0–5 dB. Seven approximate raster readings at 1 dB agree within 2%; this is a sparse pixel check, not full author-data agreement. Zero-loss photon-counting KL is infinite and recorded separately. |

## Algorithms and appendices

- **Page detection and change-point estimate:** strict one-based first alarm,
  reflected evidence, and the main-text argmin-prefix estimate of the change
  point. Ties use the earliest prefix. No estimate is fabricated without an alarm.
- **BPSK communication/detection:** unknown signs are marginalized in the actual
  Gaussian-mixture likelihood. No transmitted sign or unknown true channel is
  provided to the detector. Uniform-input BPSK mutual information is calculated
  in bits/pulse. The illustrative capacity sweep includes low-energy cases beyond
  the paper's bright-pulse regime; it is not a capacity-achieving coding simulation.
- **Appendix A:** phase Gaussian truncated to [-pi,pi], normalized quadrature,
  amplitude/phase-sensitive PNR nulling, homodyne squeezed/antisqueezed quadratures
  and optional thermal occupation. The strong-LO homodyne model has no fractional
  amplitude-error term. Absolute residual amplitude (Eq.7) and fractional LO
  amplitude error (Eq.10) are kept distinct.
- **Appendix B:** energy derivatives computed directly from Eq.3; exact positive
  high-energy squeezing/classical crossover 930.513 photons, versus the Appendix
  B.3 large-energy approximation 900 photons for N=100, eta1=0.9. Conditional
  block likelihoods follow Eqs.16–18; Page reflection never erases block history.
- **Appendix C:** normalized Hadamard mixing of one squeezed input plus vacuum
  inputs reproduces the output covariance; uniform loss is then applied. This is
  a statistical network calculation, not a hardware pulse-control simulator.
- **Appendix D:** exact non-Gaussian photon quadrature sampler, density, CDF and
  likelihoods. The entropy crossover with squeezing is r=0.208618, Na=0.044157,
  or 1.812 dB, close to the text's rounded 1.83 dB.
- **Appendix E:** Bernoulli survival/heralding and Poisson counting likelihoods,
  including boundary/infinite-KL cases. An ideal null has no false alarms, but
  still waits for a post-change photon click at finite energy. For a dual-rail
  qubit, one heralding observation covers two time bins; pulse-time conversion
  must account for that when using the Bernoulli API.
  The [new explicit dual-rail model](results_dual_rail/RESULTS.md) propagates
  arbitrary 2x2 qubit states into vacuum/early/late output states using three
  trace-preserving Kraus operators. Equal rail loss preserves the heralded
  state. The detector uses only total-survival heralds, with readout after the
  two-bin pair; a mixed-loss boundary pair has its actual state-dependent
  probability. Two independent ARL holdouts, 240,000 changed trajectories and
  1,500 state samples are retained. Ideal receiver phase correction is assumed;
  cavity dynamics, QKD security and optimal code tradeoffs are not simulated.
- **Appendix F:** shared, unreset null trajectories for several increasing alarm
  thresholds, including multiple crossings on the same overshooting pulse. The
  saved low-SNR BPSK demonstration agrees with independent renewal estimates at
  all three thresholds. Those shared-path ARL estimates are correlated. The exact
  Fock sampler replaces numerical inverse-CDF sampling while preserving Eq.23.

## Validation and reproduction

The current full suite has **421 passing tests**, including six portable-export
checks. The paper-focused study
passes **382 saved-data checks**, including regeneration of all changed runs,
weighted renewal cycles, DV calibration summaries and shared-threshold paths.
All 16 BPSK settings and the DV setting have two independent ARL repetitions
whose approximate 95% intervals lie within ±5% of two million, passing weight
and repeatability diagnostics. These are approximate Monte Carlo checks, not
certified finite-sample bounds. Reported delay is conditional E[T−nu | T>=nu],
not Lorden's worst-case delay. Early first alarms are separate; censored means
are never silently reported as unrestricted means.

The current workspace already includes all required source results and the
checksummed CC BY 4.0 author dataset. With the installed `.venv` active:

```sh
python -m pytest -q
python -m lossqcd.paper_completion --output results_paper_completion_repeat
python -m lossqcd.audit_paper_completion results_paper_completion_repeat --regenerate --regenerate-calibration
```

The output directory must be new. No network or MATLAB execution is needed.
To regenerate the original author-comparison input first, run
`python -m lossqcd.author_comparison --output results_author_comparison_repeat`,
then pass `--source results_author_comparison_repeat` to the new study. Old
unmodulated/entangled analytic figure commands remain documented in the README.

The new files are `paper_models.py`, `receivers.py`, `iid_detection.py`,
`paper_completion.py` and `audit_paper_completion.py`. Existing coherent,
squeezed and correlated-detector results remain intact. Numerical tables,
PNG/SVG plots, seeds, full validation cycles and first-alarm arrays are saved in
`results_paper_completion`. Original source checksums are audited before/after.

The subsequent PDF-verified study is in `results_verified_figures_v2`; its
958 checks recalculate saved tables, verify quadrature refinement and physical
inequalities, and verify the original PDF hash when that file is supplied.
See its [results and reproduction instructions](results_verified_figures_v2/RESULTS.md).
The earlier receiver sweep remains an explicitly chosen evaluation and is not
silently relabeled as the Figure 5 reproduction.

The independent follow-up in `results_equation_validation` passes 114 saved-data
checks, including exact regeneration of ten 500,000-draw ensembles and adaptive
probability tables. Its one-point energy fit reserves three points as held out.
All earlier thresholds and detection experiments remain unchanged.

The dual-rail addition in `results_dual_rail` passes 92 audit checks with exact
regeneration, including six calibration rounds, two independent holdouts, a
direct no-change control, eight encoding/change-phase ensembles and state
propagation. Its ARL unit is explicitly pairs or physical time bins, never
silently interchanged.

## What remains before claiming complete numerical reproduction

The [ancillary DV uncertainty analysis](results_dv_uncertainty_final/RESULTS.md)
now distinguishes the author signed mean from a nonnegative re-summary. Under
the complete-first-alarm interpretation, the remaining conditional gap is
4.459 pulses, bootstrap percentile interval [2.665, 6.206]. This is conditional
on independence and missing-event assumptions, not a diagnosis of source error.
The [source audit](source_audit/REPORT.md) confirms no simulator in the inspected
arXiv archives or released Zenodo files, with explicit search/access limits.

1. Resolve the Figure 4 fixed-seed ordering discrepancy and the approximately
   twelve-fold Figure 5 Kennedy gap. Legend settings, axis ranges and Appendix A's
   PNR convention have now been verified, and Figure 12 has supporting sparse
   pixel agreement. Original Figure 4/5/12 numerical grids and plotting arrays
   remain unavailable; the new smooth grids are chosen here.
2. Reconcile absolute/relative latency differences in Figures 3,7,11 using the
   original simulator thresholds and event-handling conventions. These are not
   saved in the supplied Figure 3/7 data. Parameter-fitting diagnostics are not
   substitutes for identifying the authors' actual model and implementation.
3. Update the consolidated figure-by-figure report if original source arrays
   or thresholds become available. The current report and combined saved-data
   checker are complete; full numerical agreement is not. Author error bars
   are not inferred, negative source latencies are retained, and detector
   thresholds are never fitted to published delay values.

The conservative thresholds, uncertain-channel studies and fitted-loss probes
are **separate research extensions**, retained in their own result directories.
They do not establish reproduction of a published figure. The general optimal
communications/detection tradeoff and alternative encodings are explicitly open
research directions in the paper, not missing solved algorithms to implement.
