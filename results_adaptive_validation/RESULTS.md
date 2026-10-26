# Adaptive-mixture validation at the original thresholds

**2/2 fixed-threshold follow-ups confirmed.** Combined with the original 23 confirmations, 25/25 grid settings now have supporting evidence. The original 23/25 report and its failed repetitions remain unchanged.

## Method and derivation

At each pulse let p1(x|prefix) and p2(x|prefix) be the two conditional Gaussian densities, W the Page evidence BEFORE the update, and ell=log(p2/p1). With a fixed b=2, sample from q=(p1+exp(W−b)*p2)/(1+exp(W−b)). Equivalently choose the post component with probability sigmoid(W−b), then sample that conditional Gaussian. The component choice is auxiliary: the likelihood correction always uses the ENTIRE mixture density.

The exact one-step correction is log(p1/q)=log(1+exp(W−b))−log(1+exp(W+ell−b)). Stable logaddexp computes it. Accumulate these corrections across every pulse of the cycle, including Page reflections and partial final blocks. Using the selected component density, exp(−sum ell), or the reflected W at the end would produce incorrect weights for this proposal. The raw sum is retained for auditing but is not used as the adaptive path weight.

The mixture is predictable from observed history and has full support. Multiplying these conditional density ratios gives the stopped path density ratio. Cycles still renew ONLY when W=0 at a block boundary, or end at the first strict pulse alarm W>h. Thus ARL remains E_pre[cycle length] / E_q[exp(path log weight)*cycle-alarm indicator]. The numerator and denominator use independent ensembles; the denominator averages over all cycles, including failures with zero success weight. Neither self-normalization nor clipping is used.

The general idea of retaining the target density in a proposal mixture is discussed in [Hesterberg (1995)](https://doi.org/10.1080/00401706.1995.10484303). This state-dependent sequential proposal and its Page-based schedule are our implementation, not a method claimed from the optical paper or a direct application of a static mixture variance theorem. There is no claimed global bound on stopped-path weights and no proof of finite variance.

## Frozen protocol and independent evidence

The exploratory pilots are in results_adaptive_pilots/pilots.json. Offset b=2 was fixed before these follow-up holdouts to reduce tail concentration while keeping successes common. No follow-up threshold is tuned: both are copied exactly from the earlier sweep. Each setting uses two fresh repetitions of 1,000,000 cycles per sampling law in batches of 100,000. Two fresh moderate checks use h=2 and h=4.5, each with 50,000 direct runs and 500,000 cycles per law.

Confirmation keeps the earlier criteria: each target and coherent ARL-ratio approximate 95% interval within ±5%; both repetitions and their pool pass ESS>=1000, largest contribution<=1%, relative SE<=2%; repeats agree within four combined SE; direct moderate comparisons agree within four combined SE and pass weighted diagnostics. Each first audit batch must also satisfy total path-weight mass one and weighted proposal cycle length equal to the independent pre-change cycle mean within four SE. These are pointwise approximate Monte Carlo checks, not rigorous bounds or simultaneous coverage.

| Na | Original threshold | New ARL ± SE | Repeat z | Matched to coherent | Confirmed |
|---:|---:|---:|---:|---|---|
| 0.1 | 12.471578048 | 2,040,296 ± 2,841 | 1.25 | True | True |
| 5 | 12.239478974 | 1,985,577 ± 3,399 | 0.22 | True | True |

## Sampling diagnostics

| Na | Method | Cycles per law (pooled) | Relative SE | Success-weight ESS | Largest weight contribution |
|---:|---|---:|---:|---:|---:|
| 0.1 | Original post proposal | 6,000,000 | 0.973% | 10,542 | 0.87023% |
| 0.1 | New adaptive mixture | 2,000,000 | 0.139% | 443,393 | 0.00139% |
| 5 | Original post proposal | 6,000,000 | 0.924% | 11,693 | 0.78924% |
| 5 | New adaptive mixture | 2,000,000 | 0.171% | 295,887 | 0.00135% |

Original failures occurred in an individual repetition, not necessarily in the pooled value in this table. The old and new sample sizes differ; this is a descriptive diagnostic comparison, not a controlled runtime benchmark.

## Direct moderate-ARL and likelihood-mass checks

| Na | h | Direct ARL ± SE | Adaptive ARL ± SE | Difference / combined SE | Passed |
|---:|---:|---:|---:|---:|---|
| 0.1 | 2 | 48.156 ± 0.203 | 47.704 ± 0.092 | 2.03 | True |
| 0.1 | 4.5 | 694.308 ± 3.077 | 688.059 ± 1.846 | 1.74 | True |
| 5 | 2 | 52.920 ± 0.226 | 52.962 ± 0.112 | 0.16 | True |
| 5 | 4.5 | 860.590 ± 3.842 | 858.902 ± 2.828 | 0.35 | True |

For each holdout audit batch, mass and weighted-cycle-length z scores are in report.json. Raw audit arrays retain the adaptive path log weight in addition to the raw conditional LLR sum, Page score, cycle length and success/completion flags. Batch summaries and random seeds retain every weighted batch, including extreme observations.

## Detection delays carried forward

Only the ARL estimator changed. The detector, its frozen thresholds and its Gaussian observation model did not change, so all eight phase-delay ensembles are carried forward from the sweep and copied into this result folder with source hashes. No new detection-speed improvement is inferred from changing the sampler.

| Na | Change pulse | Original n=8 delay ± SE | Coherent / n=8 delay | New ARL confirmation |
|---:|---:|---:|---:|---|
| 0.1 | 1005 | 25.944 ± 0.072 | 2.966 | True |
| 5 | 1005 | 13.241 ± 0.032 | 5.528 | True |

[Combined energy plot](combined_energy_sweep.png) and [combined CSV](combined_summary.csv) overlay the new ARL evidence on copies of the original changed-run summaries. The CSV validation_source column distinguishes original_sweep from adaptive_followup. Squares in the plot mark the two follow-up points. Original results are not edited.

Delay is E[T−nu | T>=nu]; first alarms before the change remain excluded and reported separately. The phase branches share burn-in, so their ranges are descriptive. This follow-up supports the two previously flagged comparisons; it does not rewrite the historical sweep as having passed its original protocol. The combined grid is still five Na samples, not a direct reproduction against the authors’ supplemental data. Next: compare the numerical curves with that supplemental dataset, then test sensitivity to imperfect or uncertain channel parameters.
