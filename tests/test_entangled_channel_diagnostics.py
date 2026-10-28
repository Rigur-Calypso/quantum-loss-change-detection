import numpy as np
import pytest
from lossqcd.entangled_channel_diagnostics import pool_runs, repeat_z
from lossqcd.simulation import Runs, summarize_delay


def test_pooling_uses_all_surviving_first_alarms_not_mean_of_means():
    a = Runs(np.array([2, 5, 6, 10]), np.ones(4, dtype=bool), 20, 5)
    b = Runs(np.array([2, 3, 4, 15]), np.ones(4, dtype=bool), 20, 5)
    stats = summarize_delay(pool_runs([a, b]))
    assert stats["conditional_delay"] == 4  # [0,1,5,10], not average(2,10)
    assert stats["false_alarms"] == 4
    assert stats["survivors"] == 4
    assert stats["conditional_delay_se"] == pytest.approx(np.std([0, 1, 5, 10], ddof=1) / 2)


def test_pooling_keeps_censoring_visible():
    a = Runs(np.array([10, 20]), np.array([True, False]), 20, 5)
    stats = summarize_delay(pool_runs([a, a]))
    assert stats["conditional_delay"] is None
    assert stats["post_change_censored"] == 2
    assert repeat_z(stats, stats) is None


def test_pooling_rejects_inconsistent_change_points():
    a = Runs(np.array([10]), np.array([True]), 20, 5)
    b = Runs(np.array([10]), np.array([True]), 20, 6)
    with pytest.raises(ValueError):
        pool_runs([a, b])


def test_repeat_diagnostic_uses_independent_standard_errors():
    a = {"conditional_delay": 20.0, "conditional_delay_se": 0.3}
    b = {"conditional_delay": 21.0, "conditional_delay_se": 0.4}
    assert repeat_z(a, b) == pytest.approx(2.0)
