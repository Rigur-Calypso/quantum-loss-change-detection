"""Grid-independent regenerative importance-sampling ARL estimation.

A cycle starts at zero and ends at tau=inf{n>=1: S_n<=0 or S_n>h}, where
S is the UNREFLECTED sum of log(p2/p1) within that cycle. Let A={S_tau>h}.
First-step renewal gives E_pre[T] = E_pre[tau] / P_pre(A).

The numerator is sampled under p1. An INDEPENDENT ensemble sampled under p2
estimates the denominator with ordinary (not self-normalized) weights:
 P_pre(A) = E_post[ exp(-S_tau) * 1_A ].
Both Gaussian laws have full support; the path-density ratio is exp(-S_tau).
Stopping at a fixed boundary makes post-change successes common, while each
successful weight is bounded by exp(-h). Overshoot must be retained. The ratio
of sample means is consistent, not exactly unbiased at finite sample size.

General reference: Calvin, Glynn & Nakayama (2001), section 3.2:
https://web.stanford.edu/~glynn/papers/2001/CalvinGNakayama01.pdf
This module uses no numerical ARL solver or state-space grid.
"""

from dataclasses import dataclass
import numpy as np
from .models import Gaussian


@dataclass(frozen=True)
class Cycles:
    """Stopped iid likelihood paths under the pre or post sampling law.

    terminal_sum is unreflected evidence, not the Page statistic. Incomplete
    cycles remain explicit and cannot contribute an unrestricted ARL estimate.
    """

    length: np.ndarray
    terminal_sum: np.ndarray
    upper: np.ndarray
    complete: np.ndarray
    pre: Gaussian
    post: Gaussian
    threshold: float
    proposal: str
    max_steps: int

    def __post_init__(self):
        n = len(self.length)
        if (
            n < 2
            or self.length.ndim != 1
            or self.terminal_sum.shape != (n,)
            or self.upper.shape != (n,)
            or self.complete.shape != (n,)
            or self.upper.dtype != bool
            or self.complete.dtype != bool
            or not np.issubdtype(self.length.dtype, np.integer)
            or np.any(self.length < 1)
            or np.any(self.length > self.max_steps)
            or not np.all(np.isfinite(self.terminal_sum))
            or self.proposal not in ("pre", "post")
        ):
            raise ValueError("invalid cycle arrays or proposal")
        if np.any(self.upper & ~self.complete):
            raise ValueError("an upper crossing must be a completed cycle")
        if np.any(self.upper & (self.terminal_sum <= self.threshold)):
            raise ValueError("upper crossings must strictly exceed the threshold")
        if np.any(self.complete & ~self.upper & (self.terminal_sum > 0)):
            raise ValueError("completed failures must return to or below zero")


def sample_cycles(
    pre, post, threshold, cycles=500000, proposal="pre", seed=20261025, max_steps=100000
):
    """Sample unreflected iid paths until return to zero or strict upper crossing.

    The proposal changes sampling only; every path is scored by the same
    post/pre likelihood. Preserve overshoots for the stopped-path weight.
    """
    if not np.isfinite(threshold) or not 0 <= threshold <= 500:
        raise ValueError("threshold must be finite and in [0,500] for stable weighting")
    if post.kl(pre) <= 0:
        raise ValueError("distinct pre/post distributions required")
    if proposal not in ("pre", "post"):
        raise ValueError("proposal must be pre or post")
    for name, value in [("cycles", cycles), ("max_steps", max_steps)]:
        if not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError(f"{name} must be an integer >=2")
    rng = np.random.default_rng(seed)
    distribution = pre if proposal == "pre" else post
    sums = np.zeros(cycles)
    length = np.full(cycles, max_steps, dtype=np.int64)
    upper = np.zeros(cycles, dtype=bool)
    complete = np.zeros(cycles, dtype=bool)
    active = np.arange(cycles)
    for t in range(1, max_steps + 1):
        x = rng.normal(distribution.mean, np.sqrt(distribution.variance), active.size)
        sums[active] += post.llr(x, pre)
        hit_upper = sums[active] > threshold
        finished = hit_upper | (sums[active] <= 0)
        upper[active[hit_upper]] = True
        complete[active[finished]] = True
        length[active[finished]] = t
        active = active[~finished]
        if active.size == 0:
            break
    return Cycles(length, sums, upper, complete, pre, post, float(threshold), proposal, max_steps)


def success_weights(post_cycles):
    """Return ordinary post-to-pre weights, with zero for completed failed cycles.

    Reject unfinished cycles rather than silently treating them as failures.
    """
    if post_cycles.proposal != "post":
        raise ValueError("importance weights require post-proposal cycles")
    if not np.all(post_cycles.complete):
        raise RuntimeError("unfinished cycles cannot be discarded or assigned zero weight")
    weights = np.zeros(len(post_cycles.length))
    weights[post_cycles.upper] = np.exp(-post_cycles.terminal_sum[post_cycles.upper])
    return weights


def _mean_se(x):
    scale = float(np.max(np.abs(x)))
    if scale == 0:
        return 0.0, 0.0
    scaled = x / scale
    return float(np.mean(scaled) * scale), float(np.std(scaled, ddof=1) / np.sqrt(len(x)) * scale)


def summarize_cycles(pre_cycles, post_cycles):
    """Estimate ARL as mean pre-cycle length divided by weighted alarm probability.

    Numerator and denominator ensembles must be independent and complete.
    Return delta-method SEs and weight diagnostics, not certified bounds.
    """
    if pre_cycles.proposal != "pre" or post_cycles.proposal != "post":
        raise ValueError("need independent pre and post ensembles, in that order")
    if (
        pre_cycles.pre != post_cycles.pre
        or pre_cycles.post != post_cycles.post
        or pre_cycles.threshold != post_cycles.threshold
    ):
        raise ValueError("ensembles must use the same models and threshold")
    if not np.all(pre_cycles.complete) or not np.all(post_cycles.complete):
        raise RuntimeError("unfinished cycles: increase max_steps; no unrestricted ARL available")
    numerator, numerator_se = _mean_se(pre_cycles.length)
    weights = success_weights(post_cycles)
    probability, probability_se = _mean_se(weights)
    if not probability > 0:
        raise RuntimeError("no positive weighted successes; increase cycles or change proposal")
    scaled_weights = weights / weights.max()
    arl = numerator / probability
    se = arl * np.sqrt((numerator_se / numerator) ** 2 + (probability_se / probability) ** 2)
    # Log-scale delta interval stays positive; SE itself is on the original scale.
    interval = (arl * np.exp(-1.96 * se / arl), arl * np.exp(1.96 * se / arl))
    return {
        "arl": float(arl),
        "arl_se": float(se),
        "relative_se": float(se / arl),
        "approx_95_interval": [float(x) for x in interval],
        "mean_pre_cycle_length": numerator,
        "mean_pre_cycle_length_se": numerator_se,
        "cycle_alarm_probability": probability,
        "cycle_alarm_probability_se": probability_se,
        "pre_cycles": len(pre_cycles.length),
        "post_cycles": len(post_cycles.length),
        "pre_observed_alarm_cycles": int(pre_cycles.upper.sum()),
        "post_proposal_successes": int(post_cycles.upper.sum()),
        "success_weight_ess": float(scaled_weights.sum() ** 2 / np.sum(scaled_weights**2)),
        "largest_success_weight": float(weights.max()),
        "weight_upper_bound": float(np.exp(-pre_cycles.threshold)),
        "pre_pulses_simulated": int(pre_cycles.length.sum()),
        "post_pulses_simulated": int(post_cycles.length.sum()),
        "longest_pre_cycle": int(pre_cycles.length.max()),
        "longest_post_cycle": int(post_cycles.length.max()),
        "unfinished_pre_cycles": int((~pre_cycles.complete).sum()),
        "unfinished_post_cycles": int((~post_cycles.complete).sum()),
        "uncertainty_method": "independent-means delta SE; log-scale approximate 95% interval",
    }


def estimate_regenerative_arl(pre, post, threshold, cycles=500000, seed=20261025, max_steps=100000):
    """Use separate seed streams for renewal length and alarm-probability samples."""
    seeds = [
        int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(seed).spawn(2)
    ]
    numerator = sample_cycles(pre, post, threshold, cycles, "pre", seeds[0], max_steps)
    denominator = sample_cycles(pre, post, threshold, cycles, "post", seeds[1], max_steps)
    stats = summarize_cycles(numerator, denominator)
    stats.update(
        {
            "seed": seed,
            "stream_seeds": {"pre": seeds[0], "post": seeds[1]},
            "threshold": float(threshold),
            "max_steps": max_steps,
        }
    )
    return stats, numerator, denominator


def cycle_arrays(pre_cycles, post_cycles):
    """Export both laws with stable field names for exact saved-array audits."""
    return {
        f"{name}_{field}": getattr(batch, field)
        for name, batch in [("pre", pre_cycles), ("post", post_cycles)]
        for field in ["length", "terminal_sum", "upper", "complete"]
    }
