"""Deterministic scalar Page-CUSUM delay conditional on surviving to nu.

Independent of Monte Carlo sampling and importance weights. Uses the same
exact Gaussian LLR CDF and finite-volume state convention as arl_solver.
Only iid scalar models apply; correlated conditional increments do not.
"""

import numpy as np
from scipy.linalg import solve
from .arl_solver import increment_cdf


def transition_matrix(pre, post, threshold, cells, sampling):
    """Return the substochastic Page transition matrix under a fixed sampling law.

    State zero has an explicit reflection atom; other states are cell
    midpoints. Missing row mass is absorption at a strict upper crossing.
    """
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be finite and nonnegative")
    if isinstance(cells, (bool, np.bool_)) or not isinstance(cells, (int, np.integer)) or cells < 8:
        raise ValueError("cells must be an integer >=8")
    if post.kl(pre) <= 0:
        raise ValueError("distinct pre/post distributions required")
    if threshold == 0:
        return np.array([[float(increment_cdf(0.0, pre, post, sampling))]])
    edges = np.linspace(0.0, threshold, cells + 1)
    states = np.r_[0.0, (edges[:-1] + edges[1:]) / 2]
    f = increment_cdf(edges[None, :] - states[:, None], pre, post, sampling)
    q = np.column_stack([f[:, 0], np.maximum(0.0, np.diff(f, axis=1))])
    if np.any(q.sum(axis=1) > 1 + 1e-12):
        raise ArithmeticError("invalid substochastic transition matrix")
    return q


def conditional_delay_numerical(
    pre, post, threshold, change_at=1005, cells=512, sampling_pre=None, sampling_post=None
):
    """Estimate E[T−nu | T>=nu] and P(T<nu) by first-step equations.

    Propagate the surviving pre-change distribution exactly nu−1 discrete
    steps, renormalizing each step and retaining log survival. Solve
    (I−Q_post)g=1 for expected post-change samples to the first alarm. Then
    delay = survivor_distribution @ g − 1. Subtracting one gives delay zero
    when the first changed pulse alarms. This is not worst-case Lorden delay.
    Optional sampling laws change the physical observations while pre/post
    remain the detector's assumed likelihood models (parameter mismatch).
    """
    if (
        isinstance(change_at, (bool, np.bool_))
        or not isinstance(change_at, (int, np.integer))
        or change_at < 1
    ):
        raise ValueError("change_at must be a positive integer")
    qpre = transition_matrix(
        pre, post, threshold, cells, pre if sampling_pre is None else sampling_pre
    )
    qpost = transition_matrix(
        pre, post, threshold, cells, post if sampling_post is None else sampling_post
    )
    distribution = np.zeros(qpre.shape[0])
    distribution[0] = 1.0
    log_survival = 0.0
    for _ in range(change_at - 1):
        following = distribution @ qpre
        mass = float(following.sum())
        if not 0 < mass <= 1 + 1e-12:
            raise ArithmeticError("invalid surviving mass")
        distribution = following / mass
        log_survival += np.log(min(mass, 1.0))
    matrix = np.eye(qpost.shape[0]) - qpost
    g = solve(matrix, np.ones(qpost.shape[0]), assume_a="gen", check_finite=False)
    if not np.all(np.isfinite(g)) or np.any(g < 1 - 1e-10):
        raise ArithmeticError("invalid post-change stopping-time solution")
    return {
        "threshold": float(threshold),
        "change_at": int(change_at),
        "cells": int(cells),
        "conditional_delay": float(distribution @ g - 1),
        "false_alarm_probability": float(-np.expm1(log_survival)),
        "log_survival_to_change": float(log_survival),
        "post_linear_residual": float(np.max(abs(matrix @ g - 1))),
        "method": "finite-volume absorbing pre-change propagation and post-change first-step solution",
    }
