"""Regenerative ARL estimation for the Appendix B.4 conditional detector.

A renewal occurs ONLY at a block boundary with Page statistic zero. A zero
inside a block retains prefix information and is not a renewal. Cycles end at
the first pulse alarm W>h or such a renewal, and lengths are measured in pulses.
ARL = E_pre[length] / P_pre[alarm in cycle]. Independent pre/post ensembles
estimate the numerator and denominator. Under the post proposal, successful
weights are exp(-sum(raw conditional LLR)), including ALL reflections and the
partial final block. Substituting exp(-W) is wrong for correlated cycles.

These ordinary importance weights are not bounded by exp(-h) for n>1.
Diagnostics and independent checks are necessary; confidence intervals are
asymptotic Monte Carlo intervals, not certified error bounds.

An optional state-dependent proposal uses q=(p1+exp(W-b)*p2)/(1+exp(W-b)).
It mixes the two CONDITIONAL densities using the current Page score W and
fixed offset b. Its exact one-step log(p1/q) is
 logaddexp(0,W-b) - logaddexp(0,W+ell-b).
This uses the pre-update W, raw ell, and the whole mixture density, not the
selected component density. Its stopped path weight must be accumulated
separately from the raw conditional LLR sum. It is not globally bounded by
the standard iid success-weight bound or the static defensive-mixture bound.

With sampling_eta1, the actual null block may differ from the assumed pre
score model. The adaptive proposal mixes the ACTUAL null conditional law
and assumed post law, but Page still uses the assumed pre/post score. Its
weight uses ell_target=log(p_assumed_post/p_actual_null), not that Page score.
A fixed post proposal likewise stores its actual-null path weight explicitly.
The optional fifth parameters entry records actual eta1 for model matching.
"""

from dataclasses import dataclass
import numpy as np
from .conditional_detection import detector_models
from .entangled import conditional_parameters
from .cusum import validate_threshold
from .rare_event import _mean_se


@dataclass(frozen=True)
class ConditionalCycles:
    """Pulse-level renewal paths with separate raw score, Page state and path weight.

    A return to Page zero renews only at a block boundary. log_weight is
    required for adaptive proposals or mismatched actual-null sampling.
    """

    length: np.ndarray
    raw_sum: np.ndarray
    page: np.ndarray
    upper: np.ndarray
    complete: np.ndarray
    kind: str
    n: int
    threshold: float
    proposal: str
    parameters: tuple
    log_weight: np.ndarray | None = None
    proposal_offset: float | None = None


def sample_conditional_cycles(
    kind="entangled",
    n=2,
    threshold=10.0,
    cycles=100000,
    proposal="pre",
    seed=20261021,
    max_steps=100000,
    N=100.0,
    Na=0.1,
    eta1=0.9,
    eta2=0.85,
    proposal_offset=2.0,
    sampling_eta1=None,
):
    """Sample correlated renewal paths while keeping score and proposal laws distinct.

    Conditional prefixes survive Page reflections within a block. Each cycle
    ends at its first alarm or a zero-Page block boundary; incomplete cycles
    are retained for rejection by the summary routine.
    """
    validate_threshold(threshold)
    (pre, post), _, _ = detector_models(kind, n, N, Na, eta1, eta2)
    actual_eta1 = eta1 if sampling_eta1 is None else sampling_eta1
    (actual_pre, _), _, _ = detector_models(kind, n, N, Na, actual_eta1, eta2)
    if post.kl(pre) <= 0:
        raise ValueError("distinct pre/post models required")
    if proposal not in ("pre", "post", "adaptive"):
        raise ValueError("proposal must be pre, post or adaptive")
    if not np.isfinite(proposal_offset) or not 0 <= proposal_offset <= 10:
        raise ValueError("proposal_offset must be finite and in [0,10]")
    for name, value in [("cycles", cycles), ("max_steps", max_steps)]:
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or value < 2
        ):
            raise ValueError(f"{name} must be an integer >=2")
    rng = np.random.default_rng(seed)
    page = np.zeros(cycles)
    raw = np.zeros(cycles)
    past = np.zeros(cycles)
    log_weight = (
        np.zeros(cycles)
        if proposal == "adaptive" or (proposal == "post" and actual_eta1 != eta1)
        else None
    )
    length = np.full(cycles, max_steps, dtype=np.int64)
    upper = np.zeros(cycles, dtype=bool)
    complete = np.zeros(cycles, dtype=bool)
    active = np.arange(cycles)
    for t in range(1, max_steps + 1):
        k = (t - 1) % n
        if k == 0:
            past[active] = 0.0
        m1, v1 = conditional_parameters(pre, past[active], k)
        m2, v2 = conditional_parameters(post, past[active], k)
        m0, v0 = conditional_parameters(actual_pre, past[active], k)
        if proposal == "adaptive":
            beta = 1 / (1 + np.exp(proposal_offset - page[active]))
            choose_post = rng.random(len(active)) < beta
            mean = np.where(choose_post, m2, m0)
            var = np.where(choose_post, v2, v0)
        else:
            mean, var = (m0, v0) if proposal == "pre" else (m2, v2)
        x = rng.normal(mean, np.sqrt(var), len(active))
        ell = 0.5 * (np.log(v1 / v2) + (x - m1) ** 2 / v1 - (x - m2) ** 2 / v2)
        # Proposal correction follows the actual null, independently of the Page score.
        if log_weight is not None:
            target_ell = 0.5 * (np.log(v0 / v2) + (x - m0) ** 2 / v0 - (x - m2) ** 2 / v2)
            log_weight[active] += (
                adaptive_log_weight(page[active], target_ell, proposal_offset)
                if proposal == "adaptive"
                else -target_ell
            )
        raw[active] += ell
        page[active] = np.maximum(0.0, page[active] + ell)
        past[active] += x
        alarm = page[active] > threshold
        # A zero inside the block still carries a correlated prefix.
        renewal = (k == n - 1) & (page[active] == 0.0)
        finished = alarm | renewal
        upper[active[alarm]] = True
        complete[active[finished]] = True
        length[active[finished]] = t
        active = active[~finished]
        if not len(active):
            break
    parameters = (float(N), float(Na), float(eta1), float(eta2))
    if actual_eta1 != eta1:
        parameters = parameters + (float(actual_eta1),)
    return ConditionalCycles(
        length,
        raw,
        page,
        upper,
        complete,
        kind,
        n,
        float(threshold),
        proposal,
        parameters,
        log_weight,
        float(proposal_offset) if proposal == "adaptive" else None,
    )


def adaptive_log_weight(page, increment, offset=2.0):
    """Exact predictable-mixture log(p1/q) at a single pulse."""
    return np.logaddexp(0.0, np.asarray(page) - offset) - np.logaddexp(
        0.0, np.asarray(page) + increment - offset
    )


def summarize_conditional_cycles(pre, post):
    """Combine independent cycle means using ordinary stopped-path weights.

    Scale log weights before exponentiation to avoid overflow. Report ESS,
    weight concentration and SE checks alongside the ARL estimate.
    """
    if (
        pre.proposal != "pre"
        or post.proposal not in ("post", "adaptive")
        or (pre.kind, pre.n, pre.threshold, pre.parameters)
        != (post.kind, post.n, post.threshold, post.parameters)
    ):
        raise ValueError("need independent pre/post ensembles for the same detector")
    if not np.all(pre.complete) or not np.all(post.complete):
        raise RuntimeError("unfinished cycles: increase max_steps; no unrestricted ARL available")
    if not np.any(post.upper):
        raise RuntimeError("no weighted successes")
    numerator, numerator_se = _mean_se(pre.length)
    if post.proposal == "adaptive" or post.log_weight is not None:
        if (
            post.log_weight is None
            or post.log_weight.shape != post.raw_sum.shape
            or not np.all(np.isfinite(post.log_weight))
        ):
            raise ValueError("adaptive cycles require finite path log weights")
        log_weights = post.log_weight[post.upper]
    else:
        log_weights = -post.raw_sum[post.upper]
    offset = float(log_weights.max())
    scaled = np.zeros(len(post.length))
    scaled[post.upper] = np.exp(log_weights - offset)
    mean, se = _mean_se(scaled)
    log_probability = float(np.log(mean) + offset)
    log_arl = np.log(numerator) - log_probability
    if not -700 < log_arl < 700:
        raise RuntimeError("ARL outside stable floating point range")
    arl = float(np.exp(log_arl))
    relative_se = float(np.hypot(numerator_se / numerator, se / mean))
    ess = float(scaled.sum() ** 2 / np.sum(scaled**2))
    fraction = float(scaled.max() / scaled.sum())
    return {
        "arl": arl,
        "arl_se": arl * relative_se,
        "relative_se": relative_se,
        "approx_95_interval": [
            float(arl * np.exp(-1.96 * relative_se)),
            float(arl * np.exp(1.96 * relative_se)),
        ],
        "mean_pre_cycle_length": numerator,
        "mean_pre_cycle_length_se": numerator_se,
        "cycle_alarm_probability": float(np.exp(log_probability)),
        "cycle_alarm_probability_se": float(np.exp(log_probability) * se / mean),
        "success_weight_ess": ess,
        "largest_weight_fraction": fraction,
        "max_log_success_weight": offset,
        "diagnostics_passed": bool(ess >= 1000 and fraction <= 0.01 and relative_se <= 0.02),
        "pre_cycles": len(pre.length),
        "post_cycles": len(post.length),
        "pre_observed_alarm_cycles": int(pre.upper.sum()),
        "post_proposal_successes": int(post.upper.sum()),
        "pre_pulses_simulated": int(pre.length.sum()),
        "post_pulses_simulated": int(post.length.sum()),
        "longest_pre_cycle": int(pre.length.max()),
        "longest_post_cycle": int(post.length.max()),
        "unfinished_pre_cycles": 0,
        "unfinished_post_cycles": 0,
        "proposal": post.proposal,
        "proposal_offset": post.proposal_offset,
        "uncertainty_method": "independent-means delta SE; log-scale approximate 95% interval",
        "diagnostic_rule": "ESS >=1000, largest weight fraction <=1%, relative SE <=2%; necessary screening only",
    }


def estimate_conditional_arl(
    kind="entangled",
    n=2,
    threshold=10.0,
    cycles=100000,
    seed=20261021,
    max_steps=100000,
    proposal="post",
    proposal_offset=2.0,
    **parameters,
):
    """Estimate pulse ARL with independent null and denominator-proposal streams."""
    if proposal not in ("post", "adaptive"):
        raise ValueError("denominator proposal must be post or adaptive")
    seeds = [
        int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(seed).spawn(2)
    ]
    pre = sample_conditional_cycles(
        kind, n, threshold, cycles, "pre", seeds[0], max_steps, **parameters
    )
    post = sample_conditional_cycles(
        kind,
        n,
        threshold,
        cycles,
        proposal,
        seeds[1],
        max_steps,
        proposal_offset=proposal_offset,
        **parameters,
    )
    stats = summarize_conditional_cycles(pre, post)
    stats.update(
        {
            "seed": seed,
            "stream_seeds": seeds,
            "threshold": float(threshold),
            "max_steps": max_steps,
            "kind": kind,
            "n": n,
            "parameters": dict(zip(["N", "Na", "eta1", "eta2"], pre.parameters)),
        }
    )
    if len(pre.parameters) == 5:
        stats["parameters"]["sampling_eta1"] = pre.parameters[4]
    return stats, pre, post


def conditional_cycle_arrays(pre, post):
    """Export renewal arrays and any explicit proposal path weights without recomputing them."""
    result = {
        f"{name}_{field}": getattr(batch, field)
        for name, batch in [("pre", pre), ("post", post)]
        for field in ["length", "raw_sum", "page", "upper", "complete"]
    }
    if post.log_weight is not None:
        result["post_log_weight"] = post.log_weight
    return result
