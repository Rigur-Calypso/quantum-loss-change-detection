import numpy as np
import pytest
from lossqcd.phase_simulation import simulate_phases
from lossqcd.entangled import entangled_pair, conditional_llrs
from lossqcd.conditional_detection import simulate_conditional
from lossqcd.simulation import summarize_delay
from lossqcd.cusum import stopping_time


def test_branch_replays_full_path_without_resetting_page(monkeypatch):
    import lossqcd.phase_simulation as module

    prefix = [9.5, 9.7, 9.6, 9.5]
    continuations = [[9.1, 8.8, 9.0, 9.0], [9.5, 9.1, 8.8, 9.0]]
    paths = iter([prefix, *continuations])

    class Rng:
        def __init__(self, values):
            self.values = iter(values)

        def normal(self, mean, sd, size):
            return np.full(size, next(self.values))

    monkeypatch.setattr(module.np.random, "default_rng", lambda seed: Rng(next(paths)))
    result, _ = simulate_phases(n=4, threshold=2, trials=2, burnin=4, offsets=[1, 3], horizon=8)
    a, b = entangled_pair(n=4, Na=0.1)
    for (nu, runs), continuation in zip(result.items(), continuations):
        increments = conditional_llrs(np.array(prefix + continuation).reshape(2, 4), a, b).ravel()
        t = stopping_time(increments, 2)
        assert np.all(runs.stop == (t if t is not None else 8))
        assert np.all(runs.event == (t is not None))


def test_prefix_alarms_are_preserved_and_never_restarted(monkeypatch):
    import lossqcd.phase_simulation as module

    class Rng:
        def normal(self, mean, sd, size):
            return np.full(size, 7.0)

    monkeypatch.setattr(module.np.random, "default_rng", lambda seed: Rng())
    result, meta = simulate_phases(n=4, threshold=1, trials=3, burnin=4, offsets=[1, 3], horizon=10)
    assert meta["burnin_first_alarms"] == 3
    for runs in result.values():
        assert np.all(runs.stop == 1)
        assert summarize_delay(runs)["survivors"] == 0


@pytest.mark.parametrize("n", [1, 4, 8])
def test_branched_marginals_match_full_physical_simulation(n):
    branches, _ = simulate_phases(
        n=n, threshold=3, trials=15000, burnin=16, offsets=[1, 3], horizon=500, seed=801
    )
    for nu, runs in branches.items():
        a = summarize_delay(runs)
        b = summarize_delay(
            simulate_conditional(
                n=n, threshold=3, trials=15000, horizon=500, change_at=nu, seed=802 + nu
            )
        )
        assert abs(a["conditional_delay"] - b["conditional_delay"]) < 4 * np.hypot(
            a["conditional_delay_se"], b["conditional_delay_se"]
        )
        assert abs(a["false_alarm_probability"] - b["false_alarm_probability"]) < 4 * np.hypot(
            a["false_alarm_probability_se"], b["false_alarm_probability_se"]
        )


def test_repeatability_and_censoring():
    a, seeds = simulate_phases(
        trials=100, burnin=8, offsets=[1, 3], horizon=11, threshold=50, seed=803
    )
    b, same = simulate_phases(
        trials=100, burnin=8, offsets=[1, 3], horizon=11, threshold=50, seed=803
    )
    assert seeds == same
    assert seeds["continuation_seeds"]["9"] != seeds["continuation_seeds"]["11"]
    for nu in a:
        np.testing.assert_array_equal(a[nu].stop, b[nu].stop)
        assert summarize_delay(a[nu])["conditional_delay"] is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"burnin": 3},
        {"trials": 0},
        {"offsets": []},
        {"offsets": [1, 1]},
        {"offsets": [0]},
        {"horizon": 1000},
        {"burnin": True},
    ],
)
def test_invalid_phase_inputs(kwargs):
    with pytest.raises(ValueError):
        simulate_phases(**kwargs)
