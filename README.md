# Quantum-enhanced transmission-loss change detection

Independent Python implementation and computational study of
[Guha, John, Gong and Basu, arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).
The principal statistical models and detection algorithms are implemented;
**complete numerical reproduction remains unresolved**. This is a classical
simulation, with no optical hardware or new quantum theorem claimed.

**[Combined results and implementation report](results/REPORT.md)** ·
[Source notebook](notebooks/01_paper_walkthrough.ipynb) ·
[Executed notebook](notebooks/02_executed_walkthrough.ipynb) ·
[Raw evidence downloads](https://github.com/Rigur-Calypso/quantum-loss-change-detection/releases/tag/v0.1.0)

![Figure 2 relative-entropy curves](results/figure2.png)

## Repository layout

```text
src/        Models, detectors, experiments and evidence checks
tests/     Unit and numerical-validation tests
data/      Author dataset, compressed study records and evidence manifest
notebooks/ Guided lesson and executed example
results/   One combined report, selected plots and numeric/validation summaries
```

The original proliferation of result folders has been consolidated. Detailed
historical records are compressed in `data/metadata.tar.gz`; raw ensembles live
in two release downloads. Restoring evidence creates only `data/studies/`.

## Install and run

Python 3.10+ is supported; the recorded numerical environment uses Python 3.14.2
and `requirements-lock.txt`. Run commands from the repository root.

```sh
git clone --depth 1 https://github.com/Rigur-Calypso/quantum-loss-change-detection.git
cd quantum-loss-change-detection
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m lossqcd.walkthrough --output results/walkthrough_new
```

On Windows, activate using `.venv\Scripts\activate`. Tests and the lesson need
no large downloads. Use a fresh output folder for each experiment. The pinned
lock file supports exact regeneration on compatible Python/platform versions;
compatible newer dependencies may change exact floating-point/random results.

## Recheck the saved numerical evidence

```sh
python -m lossqcd.evidence --scope core
python -m lossqcd.reproduction_report --output results/review_new
python -m lossqcd.reproduction_report --check results/review_new
```

The core archive is approximately 570 MiB. `--scope all` restores both archives
(about 1.8 GiB total), including the additional studies. Downloads and individual
files are SHA-256 verified; conflicting local evidence is never overwritten.
For offline use, add `--archive-dir /path/to/archives`. Exact sizes, hashes and
archive members are in [data/evidence.json](data/evidence.json).

The report recalculates saved formulas/summaries, not every historical random
trajectory. [The combined report](results/REPORT.md) explains exact-regeneration
commands, first-alarm/censoring conventions and remaining source questions.
To check the four additional studies after restoring all evidence:

```sh
python -m lossqcd.audit_adaptive_validation data/studies/results_adaptive_validation
python -m lossqcd.audit_null_channel_diagnostics data/studies/results_null_channel_refined
python -m lossqcd.audit_entangled_null_diagnostics data/studies/results_entangled_null_diagnostics
python -m lossqcd.audit_robust_thresholds data/studies/results_robust_thresholds
```

## Validation

The cleaned implementation passes **419 tests**, **1,869 core evidence checks**
and **984 additional-study checks**. All scientific tests are retained.

## Findings

Figure 2's analytical curves and the quoted single-photon entropy are checked.
Figures 3/7 latency comparisons and Figures 4/5 still have documented differences.
Seven sparse Figure 12 comparisons agree within 2%; this is not full-curve
reproduction. Missing author thresholds/plotting details are distinct from
hardware limitations.

Additional work covers independent rare-event ARL validation, channel-mismatch
sensitivity, conservative thresholds and an explicit ideal dual-rail channel.
These are implementation/validation extensions, not established novelty claims.
See [all findings](results/REPORT.md), [numeric summary](results/summary.csv) and
[validation record](results/validation.json).

## Attribution

Developed by Bhagyesh Kamal with AI coding assistance; independent of the paper's
authors. Code is MIT; generated data/plots and the unchanged imported Zenodo files
follow CC BY 4.0 attribution, as detailed in [LICENSE](LICENSE) and
[the author-data notice](data/zenodo_17460928/README.md).
Use [CITATION.cff](CITATION.cff) for this implementation and also cite the
[paper](https://arxiv.org/abs/2503.12276) and
[author dataset](https://doi.org/10.5281/zenodo.17460928).
Manuscript source/PDF and original paper figure files are not redistributed.
