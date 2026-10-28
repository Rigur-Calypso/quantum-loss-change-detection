"""Finite-volume ARL equations checked against Gaussian tails and direct sampling."""

import numpy as np
import pytest
from scipy.special import ndtr
from scipy.stats import ncx2
from lossqcd.models import Gaussian, coherent_pair, squeezed_pair
from lossqcd.arl_solver import (
    increment_cdf,
    increment_coefficients,
    mean_stopping_time,
    calibrate_numerical,
)
from lossqcd.simulation import simulate_stopping_times, summarize_arl, summarize_delay


@pytest.mark.parametrize(
    "pre,post",
    [
        (Gaussian(0, 1), Gaussian(1, 1)),
        (Gaussian(-1, 0.5), Gaussian(0.2, 2)),
        (Gaussian(-1, 2), Gaussian(0.2, 0.5)),
    ],
)
@pytest.mark.parametrize("sample_post", [False, True])
def test_llr_cdf_against_independent_normal_or_noncentral_chisquare(pre, post, sample_post):
    sampling = post if sample_post else pre
    a, b, c = increment_coefficients(pre, post, sampling)
    u = np.linspace(-10, 10, 101)
    if a == 0:
        expected = ndtr((u - c) / abs(b))
    else:
        vertex = c - b * b / (4 * a)
        q = (u - vertex) / a
        if a > 0:
            expected = ncx2.cdf(q, df=1, nc=(b / (2 * a)) ** 2)
        else:
            expected = ncx2.sf(q, df=1, nc=(b / (2 * a)) ** 2)
    np.testing.assert_allclose(increment_cdf(u, pre, post, sampling), expected, atol=1e-12)
    z = np.array([-2.0, -0.1, 0.0, 1.0, 3.0])
    x = sampling.mean + np.sqrt(sampling.variance) * z
    np.testing.assert_allclose(a * z * z + b * z + c, post.llr(x, pre), atol=1e-12)


def test_zero_threshold_geometric_first_positive_llr():
    pre, post = coherent_pair()
    delta = post.mean - pre.mean
    expected = 1 / ndtr(-abs(delta) / (2 * np.sqrt(pre.variance)))
    assert mean_stopping_time(pre, post, 0)["mean_stopping_time"] == pytest.approx(expected)


@pytest.mark.parametrize("pair", [coherent_pair(Na=0.1), squeezed_pair(Na=0.1)])
def test_integral_solver_against_seeded_arl_and_zero_state_delay(pair):
    pre, post = pair
    for sampling, nu in [(pre, None), (post, 1)]:
        runs = simulate_stopping_times(
            pre, post, 2.0, trials=50000, horizon=10000, change_at=nu, seed=802
        )
        numeric = mean_stopping_time(pre, post, 2.0, cells=256, sampling=sampling)
        if nu is None:
            stats = summarize_arl(runs)
            estimate, se = stats["arl"], stats["arl_se"]
        else:
            stats = summarize_delay(runs)
            estimate, se = stats["conditional_delay"] + 1, stats["conditional_delay_se"]
        assert abs(numeric["mean_stopping_time"] - estimate) < 5 * se
        assert numeric["relative_linear_residual"] < 1e-12


def test_calibration_and_mesh_guard():
    pre, post = coherent_pair()
    result = calibrate_numerical(
        pre, post, target_arl=1000, cells=128, check_cells=256, mesh_tolerance=0.01
    )
    assert result["mesh_check_passed"]
    assert result["grid_checks"][1]["mean_stopping_time"] == pytest.approx(1000, rel=1e-6)
    with pytest.raises(RuntimeError, match="mesh check failed"):
        calibrate_numerical(
            pre, post, target_arl=1000, cells=16, check_cells=32, mesh_tolerance=1e-6
        )


@pytest.mark.parametrize(
    "kwargs", [{"threshold": -1}, {"threshold": np.inf}, {"threshold": 1, "cells": 7}]
)
def test_invalid_solver_inputs(kwargs):
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        mean_stopping_time(pre, post, **kwargs)


def test_identical_model_and_cdf_endpoints():
    pre, post = coherent_pair()
    np.testing.assert_array_equal(increment_cdf([-np.inf, np.inf], pre, post), [0, 1])
    np.testing.assert_array_equal(increment_cdf([-1, 0, 1], pre, pre), [0, 1, 1])
    with pytest.raises(ValueError):
        mean_stopping_time(pre, pre, 1)
