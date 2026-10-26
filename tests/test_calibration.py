import numpy as np
import pytest
import lossqcd.calibration as calibration
from lossqcd.models import coherent_pair, squeezed_pair
from lossqcd.simulation import simulate_stopping_times
from lossqcd.cusum import stopping_time


def test_bisection_on_known_monotone_arl(monkeypatch):
    seen = []
    def fake_simulate(pre, post, threshold, trials, horizon, seed):
        seen.append(seed)
        return threshold
    monkeypatch.setattr(calibration, 'simulate_stopping_times', fake_simulate)
    monkeypatch.setattr(calibration, 'summarize_arl',
                        lambda threshold: {'arl': float(np.exp(threshold)), 'arl_se':0., 'censored':0})
    pre, post = coherent_pair()
    result = calibration.calibrate_threshold(pre, post, target_arl=100, trials=2, iterations=12)
    assert result['threshold'] == pytest.approx(np.log(100), abs=.002)
    assert len({tuple(seed.spawn_key) for seed in seen}) == len(seen)
    assert result['threshold_bracket'][0] <= np.log(100) <= result['threshold_bracket'][1]


def test_censored_calibration_refuses_biased_arl(monkeypatch):
    monkeypatch.setattr(calibration, 'simulate_stopping_times', lambda *args, **kwargs: None)
    monkeypatch.setattr(calibration, 'summarize_arl', lambda runs: {'arl':None})
    pre, post = coherent_pair()
    with pytest.raises(RuntimeError, match='censored'):
        calibration.calibrate_threshold(pre, post)


def test_zero_signal_rejected():
    pre, post = coherent_pair(eta2=.9)
    with pytest.raises(ValueError, match='identical'):
        calibration.calibrate_threshold(pre, post)


@pytest.mark.parametrize('target', [0, 1, np.nan, np.inf])
def test_invalid_target_rejected(target):
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        calibration.calibrate_threshold(pre, post, target_arl=target)


def test_independent_ratio_uncertainty():
    result = calibration.independent_ratio(20, 1, 10, .5)
    assert result['ratio'] == 2
    assert result['se'] == pytest.approx(np.sqrt(.02))
    assert result['approx_95_interval'] == pytest.approx([2-1.96*np.sqrt(.02),2+1.96*np.sqrt(.02)])
    with pytest.raises(ValueError):
        calibration.independent_ratio(1, .1, 0, .1)


@pytest.mark.parametrize('nu', [1, 20])
def test_unequal_variance_squeezed_simulation_replays_observations(nu):
    pre, post = squeezed_pair(Na=.1)
    rng = np.random.default_rng(555)
    x = []
    for t in range(1, 300):
        dist = post if t >= nu else pre
        x.append(rng.normal(dist.mean, np.sqrt(dist.variance)))
    expected = stopping_time(post.llr(x, pre), 4.5)
    result = simulate_stopping_times(pre, post, 4.5, trials=1, horizon=299,
                                     change_at=nu, seed=555)
    assert result.event[0] == (expected is not None)
    assert result.stop[0] == (299 if expected is None else expected)
