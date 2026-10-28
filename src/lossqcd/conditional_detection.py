"""Paper Appendix B.4 conditional-LLR CUSUM, with one decision per pulse.

Conditioning history resets at block boundaries. The Page statistic persists
across blocks. History is NOT reset when the Page statistic reflects to zero.
The actual signal covariance is mixed when a physical change straddles a block;
the score still compares the two homogeneous pre/post conditional hypotheses.
This variant is not covered by the iid-CUSUM optimality theorem.
"""

import numpy as np
from .models import coherent_pair
from .entangled import (
    entangled_pair,
    conditional_parameters,
    mixed_next_parameters,
    validate_block_size,
)
from .cusum import validate_threshold
from .simulation import Runs, summarize_arl


def detector_models(kind, n, N=100.0, Na=0.1, eta1=0.9, eta2=0.85):
    """Return assumed block laws, displacement amplitude and rank-one covariance strength.

    Coherent comparison uses n=1 and total energy N+Na. Entangled n=1 is
    independent squeezing, whose displacement energy remains N.
    """
    validate_block_size(n)
    if kind == "coherent":
        if n != 1:
            raise ValueError("ordinary-light comparator uses n=1")
        pre, post = coherent_pair(N, Na, eta1, eta2)
        from .entangled import SymmetricGaussianBlock

        blocks = (
            SymmetricGaussianBlock(1, pre.mean, 0.25),
            SymmetricGaussianBlock(1, post.mean, 0.25),
        )
        return blocks, float(np.sqrt(N + Na)), 0.0
    if kind != "entangled":
        raise ValueError(
            "kind must be coherent or entangled; n=1 entangled is independent squeezing"
        )
    blocks = entangled_pair(N, Na, eta1, eta2, n)
    q = float(np.expm1(-2 * np.arcsinh(np.sqrt(n * Na))) / (4 * n))
    return blocks, float(np.sqrt(N)), q


def simulate_conditional(
    kind="entangled",
    n=2,
    N=100.0,
    Na=0.1,
    eta1=0.9,
    eta2=0.85,
    threshold=5.0,
    trials=30000,
    horizon=100000,
    change_at=None,
    seed=20261021,
    sampling_eta2=None,
    sampling_eta1=None,
):
    """First-alarm simulation with optional true pre/post transmissivities.

    sampling_eta2 changes only the physical post-change observations. The
    likelihood models remain eta1/eta2. sampling_eta1 changes the true initial
    channel, including its block covariance, without changing those score laws.
    Mixed-block generation retains the exact physical covariance and history.
    """
    validate_threshold(threshold)
    (pre, post), alpha, q = detector_models(kind, n, N, Na, eta1, eta2)
    actual_eta2 = eta2 if sampling_eta2 is None else sampling_eta2
    actual_eta1 = eta1 if sampling_eta1 is None else sampling_eta1
    if not np.isfinite(actual_eta2) or not 0 <= actual_eta2 <= 1:
        raise ValueError("sampling_eta2 must be finite and in [0,1]")
    if not np.isfinite(actual_eta1) or not 0 <= actual_eta1 <= 1:
        raise ValueError("sampling_eta1 must be finite and in [0,1]")
    if not isinstance(trials, (int, np.integer)) or trials < 1:
        raise ValueError("trials must be a positive integer")
    if not isinstance(horizon, (int, np.integer)) or horizon < 1:
        raise ValueError("horizon must be a positive integer")
    if change_at is not None and (
        not isinstance(change_at, (int, np.integer)) or not 1 <= change_at <= horizon
    ):
        raise ValueError("change_at must lie in the observation window")
    if n == 1:
        from .simulation import simulate_stopping_times

        (actual_pre, actual_post), _, _ = detector_models(kind, n, N, Na, actual_eta1, actual_eta2)
        return simulate_stopping_times(
            pre.collective,
            post.collective,
            threshold,
            trials,
            horizon,
            change_at,
            seed,
            sampling_post=actual_post.collective,
            sampling_pre=actual_pre.collective,
        )
    rng = np.random.default_rng(seed)
    stop = np.full(trials, horizon, dtype=np.int64)
    event = np.zeros(trials, dtype=bool)
    active = np.arange(trials)
    w = np.zeros(trials)
    past = np.zeros(trials)
    weighted = np.zeros(trials)
    eta_sum = 0.0
    physical_history = change_at is not None or actual_eta1 != eta1
    for t in range(1, horizon + 1):
        k = (t - 1) % n
        # Block history resets here; Page evidence continues across blocks.
        if k == 0:
            past.fill(0.0)
            weighted.fill(0.0)
            eta_sum = 0.0
        eta = actual_eta2 if change_at is not None and t >= change_at else actual_eta1
        m1, v1 = conditional_parameters(pre, past[active], k)
        m2, v2 = conditional_parameters(post, past[active], k)
        if not physical_history:
            true_mean, true_var = m1, v1
        else:
            true_mean, true_var = mixed_next_parameters(alpha, q, eta, eta_sum, weighted[active])
        x = rng.normal(true_mean, np.sqrt(true_var), active.size)
        increment = 0.5 * (np.log(v1 / v2) + (x - m1) ** 2 / v1 - (x - m2) ** 2 / v2)
        w[active] = np.maximum(0.0, w[active] + increment)
        past[active] += x
        if physical_history:
            weighted[active] += np.sqrt(eta) * (x - alpha * np.sqrt(eta))
            eta_sum += eta
        hit = w[active] > threshold
        stopped = active[hit]
        stop[stopped] = t
        event[stopped] = True
        active = active[~hit]
        if not active.size:
            break
    return Runs(stop, event, horizon, change_at)


def calibrate_conditional(
    kind,
    n,
    target_arl=1000.0,
    trials=20000,
    horizon=100000,
    seed=20261021,
    iterations=10,
    N=100.0,
    Na=0.1,
    eta1=0.9,
    eta2=0.85,
):
    """Calibrate in PULSES by direct no-change simulation; fresh holdout required.

    Scalar iid grid/renewal calibration does not apply to this conditional,
    correlated process. This routine never substitutes per-block CRE for ARL.
    """
    if not np.isfinite(target_arl) or target_arl <= 1:
        raise ValueError("target ARL must be finite and >1")
    if not isinstance(iterations, (int, np.integer)) or not 2 <= iterations <= 30:
        raise ValueError("iterations must be an integer in [2,30]")
    streams = np.random.SeedSequence(seed).spawn(iterations + 24)
    history = []

    def evaluate(h):
        seed_i = int(streams[len(history)].generate_state(1, dtype=np.uint64)[0])
        runs = simulate_conditional(
            kind, n, N, Na, eta1, eta2, threshold=h, trials=trials, horizon=horizon, seed=seed_i
        )
        stats = summarize_arl(runs)
        if stats["arl"] is None:
            raise RuntimeError("censored conditional calibration; increase horizon")
        history.append({"threshold": h, "seed": seed_i, **stats})
        return stats["arl"]

    low, high = 0.0, 5.0
    if evaluate(low) >= target_arl:
        raise ValueError("target below nonnegative-threshold range")
    for _ in range(20):
        if evaluate(high) > target_arl:
            break
        high += 1.0
    else:
        raise RuntimeError("could not bracket conditional ARL")
    for _ in range(iterations):
        mid = (low + high) / 2
        if evaluate(mid) < target_arl:
            low = mid
        else:
            high = mid
    return {
        "threshold": (low + high) / 2,
        "target_arl_pulses": target_arl,
        "kind": kind,
        "n": n,
        "history": history,
        "seed": seed,
        "calibration_trials": trials,
        "parameters": {"N": N, "Na": Na, "eta1": eta1, "eta2": eta2},
        "method": "direct conditional-LLR Monte Carlo, noisy bisection; no retuning on validation",
    }
