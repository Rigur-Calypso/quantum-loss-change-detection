import numpy as np
import pytest
from scipy.special import ndtr
import lossqcd.rare_event as rare
from lossqcd.models import Gaussian, coherent_pair, squeezed_pair
from lossqcd.simulation import simulate_stopping_times, summarize_arl


def test_cycle_stopping_boundaries_keep_unreflected_sum(monkeypatch):
    class FakeRng:
        def __init__(self):
            self.x = iter([np.array([-0.5, 1.5]), np.array([1.5]), np.array([-1.5])])

        def normal(self, mean, sd, size):
            values = next(self.x)
            assert len(values) == size
            return values

    monkeypatch.setattr(rare.np.random, "default_rng", lambda seed: FakeRng())
    batch = rare.sample_cycles(Gaussian(0, 1), Gaussian(1, 1), 2.0, cycles=2)
    np.testing.assert_array_equal(batch.length, [1, 3])
    np.testing.assert_array_equal(batch.terminal_sum, [-1, 0])
    assert batch.complete.all()
    assert not batch.upper.any()


def test_weights_include_overshoot_and_zero_failed_cycles():
    pre, post = coherent_pair()
    p = rare.Cycles(
        np.array([1, 1, 3, 2]),
        np.array([-1.0, -0.2, 2.5, -0.3]),
        np.array([False, False, True, False]),
        np.ones(4, dtype=bool),
        pre,
        post,
        2.0,
        "pre",
        10,
    )
    q = rare.Cycles(
        np.array([2, 1, 3, 2]),
        np.array([2.5, -0.5, 3.0, -1.0]),
        np.array([True, False, True, False]),
        np.ones(4, dtype=bool),
        pre,
        post,
        2.0,
        "post",
        10,
    )
    weights = rare.success_weights(q)
    np.testing.assert_allclose(weights, [np.exp(-2.5), 0, np.exp(-3), 0])
    summary = rare.summarize_cycles(p, q)
    assert summary["cycle_alarm_probability"] == pytest.approx((np.exp(-2.5) + np.exp(-3)) / 4)
    assert summary["arl"] == pytest.approx(np.mean(p.length) / np.mean(weights))
    assert summary["largest_success_weight"] < np.exp(-2)


def test_censored_cycles_are_not_silently_treated_as_failures(monkeypatch):
    class FakeRng:
        def normal(self, mean, sd, size):
            return np.full(size, 0.6)

    monkeypatch.setattr(rare.np.random, "default_rng", lambda seed: FakeRng())
    pre, post = Gaussian(0, 1), Gaussian(1, 1)
    p = rare.sample_cycles(pre, post, 2, cycles=2, max_steps=2, proposal="pre")
    q = rare.sample_cycles(pre, post, 2, cycles=2, max_steps=2, proposal="post")
    assert not p.complete.any()
    with pytest.raises(RuntimeError, match="unfinished cycles"):
        rare.summarize_cycles(p, q)
    with pytest.raises(RuntimeError, match="unfinished cycles"):
        rare.success_weights(q)


@pytest.mark.parametrize("pair", [coherent_pair(Na=0.1), squeezed_pair(Na=0.1)])
def test_renewal_estimator_against_direct_no_change_simulation(pair):
    pre, post = pair
    direct = summarize_arl(
        simulate_stopping_times(pre, post, 2.0, trials=30000, horizon=10000, seed=900)
    )
    estimate, p, q = rare.estimate_regenerative_arl(pre, post, 2.0, cycles=100000, seed=901)
    combined = np.hypot(direct["arl_se"], estimate["arl_se"])
    assert abs(estimate["arl"] - direct["arl"]) < 4 * combined
    assert estimate["unfinished_pre_cycles"] == estimate["unfinished_post_cycles"] == 0
    assert estimate["stream_seeds"]["pre"] != estimate["stream_seeds"]["post"]
    assert estimate["success_weight_ess"] > 1000
    assert np.max(rare.success_weights(q)) <= np.exp(-2)


def test_zero_threshold_has_known_geometric_arl():
    pre, post = coherent_pair()
    expected = 1 / ndtr(-abs(post.mean - pre.mean) / (2 * np.sqrt(pre.variance)))
    estimate, _, _ = rare.estimate_regenerative_arl(pre, post, 0.0, cycles=100000, seed=902)
    assert abs(estimate["arl"] - expected) < 4 * estimate["arl_se"]
    assert estimate["mean_pre_cycle_length"] == 1
    assert estimate["mean_pre_cycle_length_se"] == 0


def test_scaled_moments_do_not_underflow_tiny_weights():
    mean, se = rare._mean_se(np.array([1e-218, 2e-218, 0.0, 0.0]))
    assert mean == pytest.approx(0.75e-218, rel=1e-12, abs=0)
    assert np.isfinite(se) and se > 0


def test_reproducible_seed():
    pre, post = coherent_pair()
    a, p, q = rare.estimate_regenerative_arl(pre, post, 4.0, cycles=100, seed=44)
    b, r, s = rare.estimate_regenerative_arl(pre, post, 4.0, cycles=100, seed=44)
    assert a == b
    for key, values in rare.cycle_arrays(p, q).items():
        np.testing.assert_array_equal(values, rare.cycle_arrays(r, s)[key])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"threshold": -1},
        {"threshold": 501},
        {"threshold": np.inf},
        {"threshold": 1, "cycles": 1},
        {"threshold": 1, "max_steps": 1},
        {"threshold": 1, "proposal": "wrong"},
    ],
)
def test_invalid_cycle_inputs(kwargs):
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        rare.sample_cycles(pre, post, **kwargs)


def test_identical_models_and_wrong_proposal_rejected():
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        rare.sample_cycles(pre, pre, 1, cycles=10)
    cycles = rare.sample_cycles(pre, post, 1, cycles=10)
    with pytest.raises(ValueError, match="post-proposal"):
        rare.success_weights(cycles)
