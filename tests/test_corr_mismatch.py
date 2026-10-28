"""Independent dense-covariance oracle for the streaming mismatch detector."""

import numpy as np
import pytest
from lossqcd.corr_cusum import simulate_conditional
from lossqcd.entangled import mixed_block_moments, entangled_pair
from lossqcd.models import coherent_pair, squeezed_pair
from lossqcd.simulation import simulate_stopping_times


def dense_oracle(n, actual_eta2, change_at, seed, trials=700, horizon=80, actual_eta1=0.9):
    N, Na, eta1, eta2, threshold = 30.0, 0.3, 0.9, 0.85, 1.5
    pre, post = entangled_pair(N, Na, eta1, eta2, n)
    rng = np.random.default_rng(seed)
    stop = np.full(trials, horizon, dtype=int)
    event = np.zeros(trials, dtype=bool)
    active = np.arange(trials)
    history = np.zeros((trials, n))
    page = np.zeros(trials)

    def conditional(mean, covariance, k):
        if k == 0:
            return mean[0], covariance[0, 0]
        coeff = np.linalg.solve(covariance[:k, :k], covariance[:k, k])
        mu = mean[k] + (history[active, :k] - mean[:k]) @ coeff
        return mu, covariance[k, k] - covariance[k, :k] @ coeff

    for t in range(1, horizon + 1):
        k = (t - 1) % n
        if k == 0:
            history.fill(0)
            times = np.arange(t, t + n)
            etas = np.full(n, actual_eta1)
            if change_at is not None:
                etas[times >= change_at] = actual_eta2
            true_mean, true_cov = mixed_block_moments(N, Na, etas)
        mu, var = conditional(true_mean, true_cov, k)
        x = rng.normal(mu, np.sqrt(var), len(active))
        m1, v1 = conditional(pre.mean, pre.covariance, k)
        m2, v2 = conditional(post.mean, post.covariance, k)
        page[active] = np.maximum(
            0, page[active] + 0.5 * (np.log(v1 / v2) + (x - m1) ** 2 / v1 - (x - m2) ** 2 / v2)
        )
        history[active, k] = x
        crossed = active[page[active] > threshold]
        stop[crossed] = t
        event[crossed] = True
        active = active[page[active] <= threshold]
        if not len(active):
            break
    return stop, event


@pytest.mark.parametrize("n", [2, 4, 8])
@pytest.mark.parametrize("change_at", [1, 11, None])
def test_streaming_matches_dense_physical_covariance(n, change_at):
    expected_stop, expected_event = dense_oracle(n, 0.839, change_at, 714)
    actual = simulate_conditional(
        n=n,
        N=30,
        Na=0.3,
        threshold=1.5,
        trials=700,
        horizon=80,
        change_at=change_at,
        seed=714,
        sampling_eta2=0.839,
    )
    np.testing.assert_array_equal(actual.stop, expected_stop)
    np.testing.assert_array_equal(actual.event, expected_event)


@pytest.mark.parametrize("n", [1, 2, 4, 8])
def test_explicit_assumed_loss_preserves_seeded_result(n):
    args = dict(n=n, threshold=2, trials=500, horizon=400, change_at=15, seed=715)
    old = simulate_conditional(**args)
    new = simulate_conditional(**args, sampling_eta2=0.85)
    np.testing.assert_array_equal(old.stop, new.stop)
    np.testing.assert_array_equal(old.event, new.event)


@pytest.mark.parametrize("n", [1, 2, 4, 8])
def test_true_post_loss_does_not_change_no_change_runs(n):
    args = dict(n=n, threshold=2, trials=500, horizon=400, seed=716)
    a = simulate_conditional(**args)
    b = simulate_conditional(**args, sampling_eta2=0.83)
    np.testing.assert_array_equal(a.stop, b.stop)
    np.testing.assert_array_equal(a.event, b.event)


@pytest.mark.parametrize("kind,pair", [("coherent", coherent_pair), ("entangled", squeezed_pair)])
def test_scalar_mismatch_delegates_to_physical_scalar_law(kind, pair):
    pre, post = pair(Na=0.2)
    _, actual_post = pair(Na=0.2, eta2=0.839)
    a = simulate_stopping_times(
        pre, post, 2, trials=300, horizon=200, change_at=11, seed=717, sampling_post=actual_post
    )
    b = simulate_conditional(
        kind,
        1,
        Na=0.2,
        threshold=2,
        trials=300,
        horizon=200,
        change_at=11,
        seed=717,
        sampling_eta2=0.839,
    )
    np.testing.assert_array_equal(a.stop, b.stop)
    np.testing.assert_array_equal(a.event, b.event)


@pytest.mark.parametrize("eta", [-0.1, 1.1, np.nan, np.inf])
def test_invalid_physical_loss(eta):
    with pytest.raises(ValueError, match="sampling_eta2"):
        simulate_conditional(sampling_eta2=eta, trials=1)


@pytest.mark.parametrize("n", [2, 4, 8])
@pytest.mark.parametrize("change_at", [11, None])
def test_initial_channel_mismatch_matches_dense_mixed_covariance(n, change_at):
    expected_stop, expected_event = dense_oracle(n, 0.85, change_at, 719, actual_eta1=0.895)
    actual = simulate_conditional(
        n=n,
        N=30,
        Na=0.3,
        threshold=1.5,
        trials=700,
        horizon=80,
        change_at=change_at,
        seed=719,
        sampling_eta1=0.895,
    )
    np.testing.assert_array_equal(actual.stop, expected_stop)
    np.testing.assert_array_equal(actual.event, expected_event)


@pytest.mark.parametrize("eta", [-0.1, 1.1, np.nan, np.inf])
def test_invalid_initial_physical_loss(eta):
    with pytest.raises(ValueError, match="sampling_eta1"):
        simulate_conditional(sampling_eta1=eta, trials=1)
