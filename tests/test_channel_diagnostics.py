import numpy as np
import pytest
from lossqcd.arl_solver import increment_cdf, mean_stopping_time
from lossqcd.channel_diagnostics import fit_actual_eta2, laws
from lossqcd.delay_solver import conditional_delay_numerical
from lossqcd.models import Gaussian, coherent_pair
from lossqcd.simulation import simulate_stopping_times, summarize_delay


def test_explicit_matched_sampling_preserves_seeded_runs():
    pre, post = coherent_pair()
    args = dict(threshold=3, change_at=40, horizon=2000, trials=1000, seed=8401)
    a = simulate_stopping_times(pre, post, **args)
    b = simulate_stopping_times(pre, post, sampling_pre=pre, sampling_post=post, **args)
    np.testing.assert_array_equal(a.stop, b.stop)
    np.testing.assert_array_equal(a.event, b.event)


def test_mismatched_zero_threshold_has_geometric_solution():
    pre, post = coherent_pair()
    actual_pre, actual_post = Gaussian(pre.mean + 0.02, 0.3), Gaussian(post.mean - 0.03, 0.2)
    stats = conditional_delay_numerical(
        pre, post, 0, change_at=7, cells=32, sampling_pre=actual_pre, sampling_post=actual_post
    )
    p = float(increment_cdf(0, pre, post, actual_pre))
    q = float(increment_cdf(0, pre, post, actual_post))
    assert stats["false_alarm_probability"] == pytest.approx(1 - p**6)
    assert stats["conditional_delay"] == pytest.approx(q / (1 - q))


@pytest.mark.parametrize("kind", ["coherent", "entangled"])
def test_mismatch_simulation_matches_first_step_equation(kind):
    pre, post, actual = laws(kind, 100, 0.2, 0.9, 0.85, 0.84)
    stats = conditional_delay_numerical(pre, post, 3, change_at=30, cells=256, sampling_post=actual)
    runs = simulate_stopping_times(
        pre, post, 3, change_at=30, trials=40000, horizon=2000, seed=9720, sampling_post=actual
    )
    mc = summarize_delay(runs)
    assert mc["post_change_censored"] == 0
    assert (
        abs(stats["conditional_delay"] - mc["conditional_delay"]) < 4 * mc["conditional_delay_se"]
    )
    # Identical pre-change distribution and stopping rule => identical false alarms.
    matched = conditional_delay_numerical(pre, post, 3, change_at=30, cells=256)
    assert matched["false_alarm_probability"] == stats["false_alarm_probability"]


def test_true_pre_sampling_is_used_without_change():
    pre, post = Gaussian(0, 1), Gaussian(1, 1)
    actual_pre = Gaussian(0.3, 1)
    result = simulate_stopping_times(
        pre, post, 0, trials=60000, horizon=100, seed=954, sampling_pre=actual_pre
    )
    expected = mean_stopping_time(pre, post, 0, cells=32, sampling=actual_pre)["mean_stopping_time"]
    assert result.event.all()
    assert abs(result.stop.mean() - expected) < 4 * result.stop.std(ddof=1) / np.sqrt(
        len(result.stop)
    )


def test_loss_fit_recovers_known_physical_parameter_without_changing_score():
    pre, post, actual = laws("entangled", 100, 0.2, 0.9, 0.85, 0.842)
    factory = lambda eta: laws("entangled", 100, 0.2, 0.9, 0.85, eta)[2]
    target = conditional_delay_numerical(
        pre, post, 4, change_at=30, cells=64, sampling_post=actual
    )["conditional_delay"]
    fit = fit_actual_eta2(pre, post, 4, target, factory, change_at=30, cells=64)
    assert fit["actual_eta2"] == pytest.approx(0.842, abs=1e-8)
    assert factory(0.842).variance != post.variance


@pytest.mark.parametrize(
    "bracket,target", [((0.9, 0.8), 10), ((-0.1, 0.8), 10), ((0.8, 0.85), np.nan)]
)
def test_fit_rejects_invalid_inputs(bracket, target):
    pre, post = coherent_pair()
    with pytest.raises(ValueError):
        fit_actual_eta2(pre, post, 4, target, lambda eta: post, bracket=bracket)


def test_laws_reject_unknown_detector():
    with pytest.raises(ValueError):
        laws("typo", 100, 0.1, 0.9, 0.85, 0.84)
