# Information needed to reconcile the published latency means

These are local research notes. No outreach has been sent.

Reference: [Zenodo 17460928](https://doi.org/10.5281/zenodo.17460928), Figures 3 and 7, with the paper's Gaussian homodyne model and Appendix B.4 detector. Our source-file checksums and analytical KL ratios match. The published simulator and threshold tables are not included.

The most useful missing information is:

1. The actual threshold for each detector and Na in Figure 7, plus its numerical or measured no-change ARL. In particular, thresholds for the coherent and n=1 squeezed cases at Na=0.1 would provide a simple scalar comparison without block-phase ambiguity.
2. ARL calibration details: pre-change simulation horizon, number of independent runs, threshold grid, and handling of runs that never alarm before the horizon.
3. Delay aggregation details: first alarm only or detector restarts; exclusion, clipping, or inclusion of alarms before the change; definitions of the delay and change index; number of Monte Carlo trials and uncertainty estimates. The saved n=5000 is an observation horizon in the Figure 3 caption, not an independently specified trial count.
4. Explicit Figure 7 transmissivities and coherent displacement/quantum energy conventions. Our comparison assumes eta1=0.9, eta2=0.85 and the paper's equal-total-energy coherent comparator; the Figure 7 data file does not store these eta values.
5. For entangled blocks, whether conditional likelihood history persists after Page reflections, how the actual covariance changes inside a straddling block, and the indexing convention for pulse 1005 relative to block boundaries.

For a concrete comparison at Na=0.1 and change pulse 1005:

- The published coherent mean is 64.624 pulses. Our source run gives 76.434 pulses at h=11.933398 and independently estimated ARL about 1.981 million. A separate diagnostic fit in our fixed model yields h=10.225939 and ARL about 357,298 pulses.
- The published n=1 squeezed mean is 36.8382 pulses. Our source run gives 46.287 pulses at h=12.246435 and ARL about 1.985 million. A diagnostic fit yields h=9.888712 and ARL about 187,147 pulses.

The diagnostic fitted thresholds are not estimates of the authors' true thresholds and are not substituted into the validated implementation. Unknown author uncertainty, receiver conventions or physical parameters could alter the comparison. Matching the published point means under one fixed model does not identify the source of a discrepancy.
