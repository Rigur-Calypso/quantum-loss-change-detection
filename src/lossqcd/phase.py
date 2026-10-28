"""Branch a common no-change prefix into several changed-run ensembles.

Burn-in ends at a block boundary. Each branch retains the Page statistic,
first-alarm flags and stopping times, but starts a fresh block prefix. Branches
have independent continuation streams. They share the burn-in, so phase
estimates are correlated; their SEs must not be combined as independent means.
Within each branch the trials are independent. Early alarms are never restarted.
"""

import numpy as np
from .corr_cusum import detector_models
from .entangled import conditional_parameters, mixed_next_parameters
from .cusum import validate_threshold
from .simulation import Runs


def simulate_phases(
    kind="entangled",
    n=8,
    Na=0.1,
    threshold=12.0,
    trials=20000,
    burnin=1000,
    offsets=tuple(range(1, 9)),
    horizon=10000,
    seed=20261022,
):
    """Branch one retained burn-in into independently seeded change-phase continuations.

    Each trial keeps its original first alarm. Branches share the burn-in,
    so estimates across phases are correlated even though continuations differ.
    """
    validate_threshold(threshold)
    (pre, post), alpha, q = detector_models(kind, n, Na=Na)
    for name, value in [("trials", trials), ("burnin", burnin), ("horizon", horizon)]:
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or value < 1
        ):
            raise ValueError(f"{name} must be a positive integer")
    if burnin % n:
        raise ValueError("burnin must end at a block boundary")
    offsets = tuple(offsets)
    if (
        not offsets
        or len(set(offsets)) != len(offsets)
        or any(
            isinstance(o, (bool, np.bool_))
            or not isinstance(o, (int, np.integer))
            or not 1 <= o <= horizon - burnin
            for o in offsets
        )
    ):
        raise ValueError("unique positive offsets must fit inside horizon")
    seeds = [
        int(s.generate_state(1, dtype=np.uint64)[0])
        for s in np.random.SeedSequence(seed).spawn(1 + len(offsets))
    ]
    rng = np.random.default_rng(seeds[0])
    w = np.zeros(trials)
    past = np.zeros(trials)
    stop = np.full(trials, horizon, dtype=np.int64)
    event = np.zeros(trials, dtype=bool)
    active = np.arange(trials)
    for t in range(1, burnin + 1):
        k = (t - 1) % n
        if k == 0:
            past[active] = 0.0
        m1, v1 = conditional_parameters(pre, past[active], k)
        m2, v2 = conditional_parameters(post, past[active], k)
        x = rng.normal(m1, np.sqrt(v1), len(active))
        ell = 0.5 * (np.log(v1 / v2) + (x - m1) ** 2 / v1 - (x - m2) ** 2 / v2)
        w[active] = np.maximum(0.0, w[active] + ell)
        past[active] += x
        hit = w[active] > threshold
        stop[active[hit]] = t
        event[active[hit]] = True
        active = active[~hit]
        if not len(active):
            break
    runs = {}
    for i, offset in enumerate(offsets):
        nu = burnin + offset
        rng = np.random.default_rng(seeds[i + 1])
        bw = w.copy()
        bs = stop.copy()
        be = event.copy()
        active = np.flatnonzero(~be)
        past = np.zeros(trials)
        weighted = np.zeros(trials)
        eta_sum = 0.0
        for t in range(burnin + 1, horizon + 1):
            if not len(active):
                break
            k = (t - 1) % n
            if k == 0:
                past[active] = 0.0
                weighted[active] = 0.0
                eta_sum = 0.0
            eta = 0.85 if t >= nu else 0.9
            m1, v1 = conditional_parameters(pre, past[active], k)
            m2, v2 = conditional_parameters(post, past[active], k)
            true_mean, true_var = mixed_next_parameters(alpha, q, eta, eta_sum, weighted[active])
            x = rng.normal(true_mean, np.sqrt(true_var), len(active))
            ell = 0.5 * (np.log(v1 / v2) + (x - m1) ** 2 / v1 - (x - m2) ** 2 / v2)
            bw[active] = np.maximum(0.0, bw[active] + ell)
            past[active] += x
            weighted[active] += np.sqrt(eta) * (x - alpha * np.sqrt(eta))
            eta_sum += eta
            hit = bw[active] > threshold
            bs[active[hit]] = t
            be[active[hit]] = True
            active = active[~hit]
        runs[nu] = Runs(bs, be, horizon, nu)
    return runs, {
        "seed": seed,
        "burnin_seed": seeds[0],
        "continuation_seeds": {str(burnin + o): seeds[i + 1] for i, o in enumerate(offsets)},
        "shared_prefix": True,
        "burnin": burnin,
        "burnin_first_alarms": int(event.sum()),
    }
