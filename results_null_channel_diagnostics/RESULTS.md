# Initial-channel sensitivity

Detector assumptions stay eta1=0.9, eta2=0.85; N=100, Na=0.1. True eta1 varies
from 0.899 to 0.901 (89.9% to 90.1% transmission). The true post-change eta2
stays 0.85. Original source thresholds stay fixed. No author mean is fitted,
and the earlier post-change diagnostic fit is not used in this separate study.

Fresh numerical no-change ARLs, in pulses:

| Detector | True eta1=0.899 | True eta1=0.900 | True eta1=0.901 |
| --- | --- | --- | --- |
| Coherent | 1,309,355 | 1,982,297 | 3,009,699 |
| Squeezed n=1 | 1,267,837 | 1,975,669 | 3,088,963 |

Each of 10 scalar settings uses 512/1024-cell first-step solutions, two independent
200,000-cycle ensembles per sampling law,
and a fresh 30,000-trial changed ensemble. Four moderate-threshold
controls also compare direct no-change simulation against renewal estimates
and the numerical solution. All checks pass: **False**.

The renewal numerator samples the actual null law. The independent denominator
uses the detector's assumed post law as proposal, with full path log weight
sum(log p_actual_null − log p_proposal). The assumed detector score is stored
separately. Substituting exp(−sum(assumed score)) would be wrong under mismatch.
No likelihood weights are clipped or self-normalized; unfinished cycles are
rejected. ESS, largest-weight fraction and relative-SE diagnostics are recorded.
Intervals and mesh refinements are approximate diagnostics, not rigorous bounds.

Changed delay is E[T−1005 | T>=1005], with early first alarms reported separately.
Censoring is explicit. Both initial-state survival and post-change timing are
checked numerically. The earlier matched-null ARL evidence is not substituted
for fresh estimates at the mismatched initial channel.

This scalar robustness study does not determine the authors' actual channel or
resolve the publication's latency discrepancy. Correlated entangled blocks are
outside its scope; they need distinct block-boundary renewals and corrected
conditional path weights under their actual null law. Original report and raw
source runs are checksum-protected. Seeds, full renewal arrays, changed runs,
control runs, CSV and plots are saved here.

Next: extend actual-null weighting to the correlated block detector and validate
it independently before drawing conclusions about entangled robustness.
