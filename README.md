# Quantum-enhanced transmission-loss change detection

An independent Python implementation and partial numerical reproduction of
[Guha, John, Gong and Basu, *Quantum-enhanced quickest change detection of transmission loss*](https://arxiv.org/html/2503.12276v2).

The problem is to notice a sudden increase in optical-channel loss as quickly
as possible. For example, transmission might fall from 90% to 85%, while noisy
measurements make that small change difficult to identify. The paper adds
squeezed light or entanglement to ordinary laser pulses so that a receiver can
collect evidence of the change faster.

This project simulates those measurement statistics and applies **Page CUSUM**:
a running evidence score that triggers an alarm after crossing a threshold.
Detection speed is compared at a matched false-alarm requirement. Everything
runs on a classical computer, including the quantum-state measurement models.

**Status: partial numerical reproduction.** Several analytical benchmarks agree;
some published detection delays and receiver curves remain unresolved. There
is no optical hardware experiment or claim of complete agreement with the paper.

[Compiled report](results/report.pdf) · [Detailed results](results/report.md) ·
[Source notebook](notebooks/lesson.ipynb) · [Executed notebook](notebooks/example.ipynb)

## What is implemented

| Component | Implementation |
| --- | --- |
| Coherent and squeezed light | Gaussian homodyne models, unequal-variance likelihood ratios and analytical KL divergence. |
| Entangled pulses | Correlated Gaussian blocks, per-pulse conditional likelihoods and changes occurring inside a block. |
| Sequential detection | Page CUSUM, first-alarm times, change-point estimates and explicit handling of early alarms and censoring. |
| Other probes and receivers | BPSK modulation, single-photon probes, photon counting and imperfect Kennedy receivers. |
| Performance estimation | Monte Carlo delay estimates, numerical ARL calibration, regenerative rare-event estimation and independent holdouts. |

Here **ARL** means the average number of pulses before a false alarm when no
change occurs. **KL divergence** measures how distinguishable the before- and
after-change measurement distributions are. Higher KL generally supports
faster detection in the asymptotic regime; it is not an exact finite-delay formula.

## Main results

| Paper comparison | Outcome |
| --- | --- |
| Figure 2 | Recalculated coherent, squeezed and entangled KL curves for N=100, transmission 0.9 → 0.85, added energy 0–5 photons and entangled groups of 2–256 pulses. This reproduces the formulas; no original numeric curve array was available in the inspected dataset. |
| Figure 3 | Analytical KL-ratio coordinates match the author data. Mean absolute relative differences in simulated delay ratios are 7.49% without modulation and 6.40% with BPSK. |
| Figure 7 | All 45 internal ARL settings pass validation. Our family-average detection delays are 13.6–34.4% higher than the authors' values. |
| Figure 11(a) | Single-photon KL is **0.206340 nats**, matching the paper's rounded **0.2063**. |
| Delay scaling | At high ARL, normalized delay slopes agree with inverse-KL scaling within **0.014%** across the checked coherent/squeezed cases. |

The Figure 4 ordering and Figure 5 Kennedy-receiver curve also remain unresolved.
For Figure 11(b), our conditional first-alarm delay and the authors' signed
latency average use different event summaries; the numerical gap is documented
rather than treated as a reproduced result. The [results report](results/report.md)
explains the assumptions, uncertainty and remaining questions for each comparison.

![Figure 2: relative entropy for coherent, squeezed and entangled probes](results/figure2.png)

## Additional studies

- **Rare false alarms:** independent validation of large ARLs using regenerative
  cycles and importance sampling, with weight and effective-sample-size checks.
- **Channel uncertainty:** measured how small errors in the assumed initial
  transmission affect false alarms, then evaluated more conservative thresholds
  and their detection-delay cost on a sampled grid.
- **Ideal dual-rail qubits:** explicit loss maps, heralded photon survival and
  two-bin timing, developing the paper's Appendix E.2 model.

These are additional computational and validation studies. Their research
novelty has not been established through a broader literature review.

## Installation and first run

Python **3.10+**. From the repository root, on macOS or Linux:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m lossqcd.lesson --output results/lesson_new
```

This creates an isolated environment, installs the package and development
dependencies, runs the tests, and executes the seven-cell introductory lesson.
The lesson saves an executed notebook, plots and simulation arrays in the
chosen output folder. **Tests and the lesson require no data download.**
Use a fresh output directory for each run.

NumPy and Matplotlib provide the core numerical and plotting dependencies.
The development installation also includes SciPy, pytest and Ruff.
`requirements.txt` records the original Python 3.14.2 dependency snapshot;
GitHub checks run on Python 3.11 and 3.14.

## Check the saved simulations

```sh
python -m lossqcd.evidence --scope core
python -m lossqcd.report --output results/review_new
python -m lossqcd.report --check results/review_new
```

The first command downloads and restores checksum-verified evidence. The next
two recalculate the core saved-data checks and verify the generated report's
file hashes. They **do not regenerate every random trajectory**.

The [v0.1.0 release](https://github.com/Rigur-Calypso/quantum-loss-change-detection/releases/tag/v0.1.0)
contains two raw simulation archives: core evidence is about **570 MiB**;
`python -m lossqcd.evidence --scope all` restores both archives, about **1.8 GiB**
in download size. Detailed metadata are bundled in `data/metadata.tar.gz`.
Restored studies live under `data/studies/` and are ignored by Git.
The [evidence manifest](data/evidence.json) lists file hashes and archive URLs.

For commands that generate new ensembles, see [Reproduce](results/report.md#reproduce).
The full evidence checks are separate from the unit tests.

**Recorded validation:** 428 tests, 1,869 core checks and 984 additional study
checks passed; all 739 raw simulation array files were checksum-verified.
Passing these checks establishes internal consistency, not complete agreement
with the published figures. See the [validation record](results/validation.json).

## Repository layout

| Folder | Contents |
| --- | --- |
| `src/lossqcd/` | Physical models, detectors, calibration, experiments and evidence checks. |
| `tests/` | Tests of formulas, stopping rules, estimators and saved-data handling. |
| `notebooks/` | Introductory lesson and its executed example. |
| `data/` | Original author data, attribution and evidence manifests/metadata. |
| `results/` | Compiled report, comparison plots and numerical summaries. |

Start with [models.py](src/lossqcd/models.py) and [cusum.py](src/lossqcd/cusum.py)
for the baseline, then [entangled.py](src/lossqcd/entangled.py) and
[corr_cusum.py](src/lossqcd/corr_cusum.py) for correlated pulses.
The [detailed report](results/report.md#follow-the-implementation) maps the
remaining implementation steps to their source files.

## References and attribution

- **Paper:** S. Guha, T. C. John, Z. Gong and P. Basu,
  [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).
- **Author data:** S. Guha and P. Basu,
  [Zenodo dataset, version v1](https://doi.org/10.5281/zenodo.17460928).
  The original files are preserved with [checksums and attribution](data/authors/README.md).

By **Bhagyesh Kamal**. This repository is independent of the paper's authors.
Code: **MIT**. Data and plots: **CC BY 4.0**, with original attribution retained.
See [LICENSE](LICENSE) for details and [CITATION.cff](CITATION.cff) to cite this project.
