"""Block-boundary CUSUM decisions with time measured in PULSES."""

import numpy as np
from scipy.special import log_ndtr
from .arl_solver import increment_coefficients
from .cusum import validate_threshold
from .entangled import entangled_pair, mixed_collective, validate_block_size
from .simulation import Runs


def minimum_arl_feasibility(pre, post, n, target_pulses):
    """At h=0, Page alarms on the first positive LLR: geometric in blocks.

    Positive thresholds only increase ARL. Large entangled blocks can already
    exceed the desired ARL at h=0, making exact matching unattainable within
    this nonnegative-threshold, full-block detector family. Log tails avoid
    subtracting an extremely small probability from one.
    """
    validate_block_size(n)
    if not np.isfinite(target_pulses) or target_pulses <= n:
        raise ValueError("target ARL must be finite and exceed one block")
    a, b, c = increment_coefficients(pre, post)
    if a == 0:
        logp = float(log_ndtr(c / abs(b))) if b != 0 else (-np.inf if c <= 0 else 0.0)
    else:
        vertex = c - b * b / (4 * a)
        if a > 0 and vertex >= 0:
            logp = 0.0
        elif a < 0 and vertex <= 0:
            logp = -np.inf
        else:
            radius = np.sqrt(-vertex / a)
            low, high = -b / (2 * a) - radius, -b / (2 * a) + radius
            if a > 0:
                logp = float(np.logaddexp(log_ndtr(low), log_ndtr(-high)))
            else:
                left, right = (
                    (log_ndtr(-low), log_ndtr(-high))
                    if low >= 0
                    else (log_ndtr(high), log_ndtr(low))
                )
                logp = float(left + np.log(-np.expm1(right - left)))
    log_min = float(np.log(n) - logp)
    min_value = float(np.exp(log_min)) if log_min < np.log(np.finfo(float).max) else None
    return {
        "minimum_arl_pulses": min_value,
        "log_minimum_arl_pulses": log_min,
        "log_first_positive_probability": logp,
        "target_arl_pulses": target_pulses,
        "attainable": bool(log_min < np.log(target_pulses)),
        "reason": "nonnegative strict-threshold full-block Page CUSUM",
    }


def simulate_entangled_blocks(
    N=100.0,
    Na=0.1,
    eta1=0.9,
    eta2=0.85,
    n=2,
    threshold=5.0,
    trials=30000,
    horizon_pulses=8192,
    change_at=None,
    seed=20261021,
):
    """Sample exact collective statistics of independent entangled blocks.

    Blocks cover [1,n], [n+1,2n], ... . Decisions only happen at their ends.
    The raw correlation model is implemented and tested separately; using its
    sufficient statistic avoids generating uninformative orthogonal modes.
    A straddling block is sampled with its mixed-transmission covariance.
    """
    validate_threshold(threshold)
    pre_block, post_block = entangled_pair(N, Na, eta1, eta2, n)
    pre, post = pre_block.collective, post_block.collective
    if not isinstance(trials, (int, np.integer)) or trials < 1:
        raise ValueError("trials must be a positive integer")
    if (
        not isinstance(horizon_pulses, (int, np.integer))
        or horizon_pulses < n
        or horizon_pulses % n
    ):
        raise ValueError("horizon must be a positive multiple of block size")
    if change_at is not None and (
        not isinstance(change_at, (int, np.integer)) or not 1 <= change_at <= horizon_pulses
    ):
        raise ValueError("change_at must be a pulse in the observation window")
    rng = np.random.default_rng(seed)
    stop = np.full(trials, horizon_pulses, dtype=np.int64)
    event = np.zeros(trials, dtype=bool)
    active = np.arange(trials)
    w = np.zeros(trials)
    for end in range(n, horizon_pulses + 1, n):
        start = end - n + 1
        if change_at is None or end < change_at:
            dist = pre
        elif start >= change_at:
            dist = post
        else:
            dist = mixed_collective(N, Na, n, end - change_at + 1, eta1, eta2)
        y = rng.normal(dist.mean, np.sqrt(dist.variance), active.size)
        w[active] = np.maximum(0.0, w[active] + post.llr(y, pre))
        crossed = w[active] > threshold
        stopped = active[crossed]
        stop[stopped] = end
        event[stopped] = True
        active = active[~crossed]
        if active.size == 0:
            break
    return Runs(stop, event, horizon_pulses, change_at)


def simulate_collective_blocks(pre, post, n, threshold, trials, horizon_pulses, change_at, seed):
    """Homogeneous independent blocks, requiring a block-aligned change.

    Used for ordinary-light/single-squeezed controls (n=1 in the experiment).
    """
    from .simulation import simulate_stopping_times

    validate_block_size(n)
    if horizon_pulses % n or (change_at is not None and (change_at - 1) % n):
        raise ValueError("horizon and change must align with block boundaries")
    block_change = None if change_at is None else (change_at - 1) // n + 1
    runs = simulate_stopping_times(
        pre, post, threshold, trials, horizon_pulses // n, block_change, seed
    )
    return Runs(runs.stop * n, runs.event, horizon_pulses, change_at)
