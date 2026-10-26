# Quantum-enhanced transmission-loss change detection

Independent Python implementation and computational study of **Guha, John,
Gong and Basu, “Quantum-enhanced quickest change detection of transmission
loss”**, [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).

The principal statistical models and detection algorithms are implemented.
**This is a partial numerical reproduction, with documented disagreements;
it is not a claim to reproduce every published curve or build optical hardware.**
The project was developed with AI coding assistance. Its numerical claims are
supported by the saved tests, independent checks and explicitly stated limits.
It is independent of the paper's authors and is not their official code.

![Implemented Figure 2 relative-entropy curves](results_entangled/figure2_entangled.png)

## Start here

- [Findings: reproduced results, additional studies and unresolved gaps](FINDINGS.md)
- [Current consolidated numerical report](paper_reproduction_public/REPORT.md)
- [Full paper coverage](PAPER_IMPLEMENTATION.md)
- [Guided notebook](notebooks/01_paper_walkthrough.ipynb) and [executed lesson](walkthrough_run_v2/01_paper_walkthrough.ipynb)
- [Reproduction instructions and evidence downloads](REPRODUCING.md)
- [Archived study history](docs/history/README.md)

## Implemented

| Area | Computation |
| --- | --- |
| Classical baseline | Coherent homodyne Gaussian measurements, mean √η α and variance 1/4 |
| Quantum augmentation | Independent squeezing and correlated Gaussian entangled blocks |
| Sequential detection | Gaussian unequal-variance log likelihoods, strict Page CUSUM, first-alarm and change-point handling |
| Performance | ARL calibration, conditional delay, separate early alarms/censoring, independent holdouts |
| Receiver models | Phase-noisy homodyne, photon-resolving/on-off Kennedy, photon counting |
| Communication/probes | Marginalized BPSK mixtures, mutual information, exact single-photon quadrature sampler |
| Logical qubits | Ideal dual-rail pure-loss channel, heralded erasures and two-time-bin detector timing |

The entangled states and receivers are simulated mathematically on a classical
computer. A quantum computer or optical laboratory is not needed to run this code.

## Install and try

Python 3.10+ is required. The retained numerical environment uses Python 3.14.2;
its exact versions are recorded in `requirements-lock.txt`.

```sh
git clone https://github.com/Rigur-Calypso/quantum-loss-change-detection.git
cd quantum-loss-change-detection
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m lossqcd.walkthrough --output walkthrough_repeat
```

The tests and teaching example do not require the large release downloads.
Every experiment output should use a fresh folder. Run commands from the
repository root. On Windows, activate with `.venv\Scripts\activate` instead.
Exact historical sample regeneration requires the recorded dependency versions;
tests on a compatible newer environment are a different validation.

## Review the saved evidence

Plots, tables, summaries, source and tests are tracked in Git. The **739 large
raw NumPy ensembles are retained locally and are not publicly uploaded yet**.
Their planned two-archive inventory, exact sizes and SHA-256 checksums are in
[release/evidence.json](release/evidence.json). Public saved-data summaries and
plots can be reviewed immediately; rerunning the saved-array audits requires
those raw ensembles or newly regenerated experiments.

If you have received the matching archives from the maintainer, restore them
without network access:

```sh
python -m lossqcd.evidence --scope core --archive-dir /path/to/archives
python -m lossqcd.reproduction_report --output review_repeat
python -m lossqcd.reproduction_report --check review_repeat
```

Use `--scope all` for the additional studies. The loader verifies checksums and
refuses conflicting evidence. The downloader's default GitHub asset URLs are
reserved for a future approved upload and are not currently available.
See [regeneration instructions](REPRODUCING.md) for running new simulations.

## Validation snapshot

The publication copy passes **421 unit tests**, **1,869 core saved-evidence
checks**, and **984 additional extension checks**. Validation blocks reads of
the original local project except its installed dependency environment.
The public report rechecks saved arrays; it does not regenerate every historical
trajectory. See [validation records](release/VALIDATION.md).

## Results in brief

- Figure 2's coherent, squeezed and n=2–256 entangled formulas are reproduced.
- All 61 original-grid detector settings have supporting internal ARL checks.
- The single-photon KL is **0.206340 nats**, agreeing with the quoted **0.2063**.
- Figure 12 has seven sparse pixel comparisons within **2%** at 1 dB; this is
  not full curve reproduction.
- Figure 3 delay-ratio differences remain **7.49%** unmodulated and **6.40%** BPSK
  on average in absolute relative terms. Figure 7 means are approximately
  **13.6–34.4% higher** by detector family.
- Figure 4 ordering and the roughly **12-fold Figure 5 Kennedy gap** remain
  unresolved despite independent equation/integration checks.

See [FINDINGS.md](FINDINGS.md) for our additional experiments and interpretation.
Hardware absence does not explain away numerical disagreements. Missing original
thresholds, event conventions and plotting calculations are separate limitations.

## Attribution and use

Code: [MIT](LICENSE). Our generated numerical data/plots: CC BY 4.0; see
[license boundaries and attribution](THIRD_PARTY_NOTICES.md). The authors'
four imported Zenodo data/script files remain unchanged and retain their
[original CC BY 4.0 attribution](data/zenodo_17460928/README.md).
The paper PDF, TeX archives and original paper figure images are not redistributed.
The source-availability audit records their checksums and source URLs.
Use [CITATION.cff](CITATION.cff) to cite this implementation; cite the paper and
the [author dataset](https://doi.org/10.5281/zenodo.17460928) for their results.
