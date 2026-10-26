# v0.1.0 — Computational implementation and evidence

Independent implementation of the principal statistical models and sequential
detectors in arXiv:2503.12276v2, including coherent/squeezed/entangled homodyne,
BPSK mixtures, photon/receiver models and an ideal dual-rail channel.

This is a partial numerical reproduction. Figure 2 formulas and the quoted
single-photon entropy are checked; complete published-figure matching remains
unresolved. See FINDINGS.md for measured differences and hardware/scope limits.

Additional studies cover rare-event ARL validation, initial-channel sensitivity
and conservative sampled-grid thresholds. These are project extensions, not
claims of a new quantum theorem or established research novelty.

The proposed local assets contain unchanged raw NumPy arrays, split into core and extensions. They are not attached to a public release yet.
The repository tracks plots, numerical tables and summaries. Archive/member
hashes and exact sizes are in release/evidence.json. Restore using
`python -m lossqcd.evidence --scope core --archive-dir /path/to/archives` or `--scope all` with the matching local archives.

The public report recalculates saved evidence; it does not regenerate every
historical random trajectory. Code is MIT; generated data/plots and imported
author data follow the CC BY 4.0 notices. Original manuscript sources are linked,
not redistributed. This release makes no claim of author endorsement.
