# Original-source availability audit

Checked 25 October 2026. This audit narrows the remaining reproduction questions; it does not claim that an unavailable original implementation does not exist.

## Public source files inspected

- [arXiv v1 source archive](https://arxiv.org/src/2503.12276v1): 11 ordinary files.
- [arXiv v2 source archive](https://arxiv.org/src/2503.12276v2): 13 ordinary files.
- [Zenodo record](https://zenodo.org/records/17460928) and [version list](https://zenodo.org/api/records/17460928/versions).
- [Published article landing page](https://journals.aps.org/prl/abstract/10.1103/czg5-3y3d).

Both archives contain manuscript/figure material, with no files classified as a simulator program. The TeX was read as text, never compiled or executed. The archives were inspected without extracting arbitrary paths.

The dataset version list has one public record. Its seven files are three numeric data tables, one plotting script and three generated PNG images. The four originally imported data/script MD5 values still match the released metadata. The plotting script loads saved arrays; it does not supply the detector simulations or thresholds.

## What changed between paper versions

| Figure asset | Evidence | Implication |
| --- | --- | --- |
| fixed_squeezing.png (Figure 4) | Identical SHA-256 in v1/v2 | The plotted ordering discrepancy is not resolved by selecting the other arXiv version. |
| Homodyne_Kennedy_comparison.png (v1) | Present only in v1 | It plots absolute residual-photon Kennedy cases, rather than the later phase/amplitude-noise comparison. |
| Imperfection_comparison.png (v2 Figure 5) | Present only in v2 | Its stated phase/amplitude model must be compared with v2 Eq.10; the old figure is not a like-for-like validation. |
| cusum_ARL.png | Present only in v2 | The original conditional detector code and thresholds are still absent from the inspected sources. |

The Figure 2, Figure 3, Figure 4, weak-probe comparison and CV/DV augmentation image files are unchanged byte-for-byte where the same filenames occur. The full per-file comparison and hashes are saved alongside this report. The v1/v2 manuscript texts differ, so unchanged image hashes alone do not establish equality of every surrounding mathematical statement.

## Bounded repository search

Queries included the exact paper title with code/repository terms, the arXiv identifier with GitHub, and GitHub-restricted transmission-loss/CUSUM searches. No original repository was identified in the returned results. One arXiv monthly-list excerpt mentioned GitHub for a different paper; that was not treated as evidence of code for this work. This is not an exhaustive search of every author account or private repository.

The APS public page currently marks supplemental material as subscription-required. That material was not retrieved, and no subscription or login was used. Its contents therefore remain an additional source uncertainty.

## Reproduction consequence

Keep the implemented physical models and independently validated calculations unchanged. The roughly ten-photon Kennedy model remains a non-unique diagnostic hypothesis; the source audit does not verify that energy. Exact matching still needs the Figure 4/5 calculation and the original Figure 3/7/11 threshold/event conventions.

See [the concrete clarification questions](CLARIFICATION.md). No author has been contacted.

## Recreate this availability audit

```sh
python scripts/audit_arxiv_sources.py --output source_audit_repeat
```

The cached archives and metadata are under `data/arxiv_source_audit`. The command uses only the Python standard library, requires an empty output folder, reads source material without executing it, and saves an inventory plus comparison table. The dated web-search/access observations are recorded research notes, not refreshed by this offline command. Earlier scientific arrays, audit records and code are unchanged.
