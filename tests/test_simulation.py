import numpy as np
import pytest
from lossqcd import coherent_pair, simulate_stopping_times, summarize_arl, summarize_delay
from lossqcd.cusum import stopping_time
from lossqcd.simulation import Runs


def test_known_false_alarm_delay_and_censoring_accounting():
    runs = Runs(np.array([2, 5, 8, 10]), np.array([True, True, True, False]), 10, 5)
    stats = summarize_delay(runs)
    assert stats['false_alarms'] == 1
    assert stats['false_alarm_probability'] == .25
    assert stats['survivors'] == 3
    assert stats['post_change_censored'] == 1
    assert stats['conditional_delay'] is None
    assert stats['restricted_conditional_delay'] == pytest.approx((0+3+5)/3)
    complete = Runs(np.array([2, 5, 8, 10]), np.ones(4, dtype=bool), 10, 5)
    assert summarize_delay(complete)['conditional_delay'] == pytest.approx(8/3)


def test_arl_never_drops_or_mislabels_censored_runs():
    runs = Runs(np.array([2, 8, 10]), np.array([True, True, False]), 10, None)
    assert summarize_arl(runs)['arl'] is None
    assert summarize_arl(runs)['restricted_mean_min_T_H'] == pytest.approx(20/3)
    complete = Runs(np.array([2, 8, 10]), np.ones(3, dtype=bool), 10, None)
    assert summarize_arl(complete)['arl'] == pytest.approx(20/3)


def test_no_survivors_and_all_censored():
    runs = Runs(np.array([1, 2]), np.ones(2, dtype=bool), 10, 5)
    assert summarize_delay(runs)['conditional_delay'] is None
    assert summarize_delay(runs)['survivors'] == 0
    pre, post = coherent_pair(eta2=.9)  # zero LLR identically
    runs = simulate_stopping_times(pre, post, 1, trials=8, horizon=10)
    assert summarize_arl(runs)['censored'] == 8
    assert summarize_arl(runs)['restricted_mean_min_T_H'] == 10


@pytest.mark.parametrize("nu", [None, 1, 7])
def test_single_run_matches_replayed_gaussian_observations(nu):
    pre, post = coherent_pair()
    rng = np.random.default_rng(39)
    x = []
    for t in range(1, 51):
        dist = post if nu is not None and t >= nu else pre
        x.append(rng.normal(dist.mean, np.sqrt(dist.variance)))
    expected = stopping_time(post.llr(x, pre), 2)
    result = simulate_stopping_times(pre, post, 2, trials=1, horizon=50, change_at=nu, seed=39)
    assert result.event[0] == (expected is not None)
    assert result.stop[0] == (50 if expected is None else expected)


def test_seed_reproducibility():
    pre, post = coherent_pair()
    a = simulate_stopping_times(pre, post, 2, trials=100, horizon=500, change_at=30, seed=99)
    b = simulate_stopping_times(pre, post, 2, trials=100, horizon=500, change_at=30, seed=99)
    np.testing.assert_array_equal(a.stop, b.stop)
    np.testing.assert_array_equal(a.event, b.event)


@pytest.mark.parametrize("kwargs", [{"trials":0}, {"horizon":0}, {"change_at":0},
                                      {"change_at":11, "horizon":10}, {"trials":2.5}])
def test_invalid_simulation_inputs(kwargs):
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        simulate_stopping_times(pre, post, 2, **kwargs)
