# Authors' supplemental data

Guha, Saikat and Basu, Prithwish (2025), *Dataset accompanying PRL publication on quantum-enhanced quickest change detection of transmission loss*, version v1, Zenodo. DOI: [10.5281/zenodo.17460928](https://doi.org/10.5281/zenodo.17460928). Retrieved 2026-10-22. License: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/), recorded as `cc-by-4.0` in the original API metadata.

The following original files are preserved unchanged:

- `prl_fig3_data`: eight analytical KL ratios and two latency-ratio series, plus settings.
- `prl_fig7_data`: nine energy values and five mean-latency series, plus settings.
- `prl_supplem_fig4b_data`: supplemental series, squeezing range and one raw latency vector.
- `plot_PRL_figs.m`: author-supplied plotting script. It is read to establish series mappings and is never executed by our Python implementation.

`manifest.json` records source URLs, file sizes, the four published MD5 checksums and locally computed SHA-256 hashes. All four MD5 checks passed. `record_metadata.json` preserves the original Zenodo API response, including author attribution and license. Python parses only a bounded subset of the Octave ASCII data format; it does not execute these files as code.

Our normalized coordinates, simulations, CSV comparisons and derived plots are separate in `results_author_comparison`. They are independent computational work and are identified separately from the authors' means. The Figure 3 energy mapping uses the supplemental squeezing range and is checked against all eight original analytical KL ratios. The saved `n=5000` is interpreted as an observation horizon in accordance with the paper's Figure 3 caption; it is not assumed to be a Monte Carlo replicate count. The original files do not contain Figure 3/7 detector thresholds, raw stopping-time/event arrays or per-point standard errors.
