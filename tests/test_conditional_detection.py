import numpy as np
import pytest
from scipy.stats import norm
from lossqcd.entangled import entangled_pair, conditional_llrs, mixed_collective
from lossqcd.block_detection import simulate_entangled_blocks, minimum_arl_feasibility
from lossqcd.conditional_detection import simulate_conditional
from lossqcd.models import coherent_pair, squeezed_pair
from lossqcd.simulation import simulate_stopping_times
from lossqcd.cusum import stopping_time


@pytest.mark.parametrize(
    "kind,pair", [("coherent", coherent_pair(Na=0.1)), ("entangled", squeezed_pair(Na=0.1))]
)
def test_n1_conditional_matches_existing_scalar_detector(kind, pair):
    a, b = pair
    expected = simulate_stopping_times(a, b, 3, trials=100, horizon=1000, change_at=20, seed=70)
    actual = simulate_conditional(
        kind, n=1, threshold=3, trials=100, horizon=1000, change_at=20, seed=70
    )
    np.testing.assert_array_equal(actual.stop, expected.stop)
    np.testing.assert_array_equal(actual.event, expected.event)


def test_conditional_history_resets_at_block_boundaries_but_statistic_persists(monkeypatch):
    import lossqcd.conditional_detection as module

    pre, post = entangled_pair(n=4, Na=0.1)
    x = np.array([9.0, 9.8, 9.2, 9.1, 9.5, 9.0, 8.9, 9.3])
    llr = conditional_llrs(x.reshape(2, 4), pre, post).ravel()
    expected = stopping_time(llr, 2.0)

    class ScriptedRng:
        def __init__(self):
            self.iterator = iter(x)

        def normal(self, mean, sd, size):
            assert size == 1
            return np.array([next(self.iterator)])

    monkeypatch.setattr(module.np.random, "default_rng", lambda seed: ScriptedRng())
    runs = simulate_conditional(n=4, threshold=2, trials=1, horizon=8, seed=99)
    assert runs.event[0] == (expected is not None)
    assert runs.stop[0] == (8 if expected is None else expected)


@pytest.mark.parametrize("nu", [1, 3, 5])
def test_collective_block_detector_replays_mixed_block_pulse_timing(nu):
    pre, post = entangled_pair(n=4, Na=0.1)
    a, b = pre.collective, post.collective
    rng = np.random.default_rng(71)
    values = []
    for end in range(4, 65, 4):
        start = end - 3
        dist = (
            a if end < nu else (b if start >= nu else mixed_collective(100, 0.1, 4, end - nu + 1))
        )
        values.append(rng.normal(dist.mean, np.sqrt(dist.variance)))
    expected = stopping_time(b.llr(values, a), 2.0)
    runs = simulate_entangled_blocks(
        n=4, threshold=2, trials=1, horizon_pulses=64, change_at=nu, seed=71
    )
    assert runs.event[0] == (expected is not None)
    assert runs.stop[0] == (64 if expected is None else 4 * expected)
    assert runs.stop[0] % 4 == 0


def test_large_full_block_minimum_arl_uses_stable_small_tail():
    a, b = entangled_pair(n=1, Na=0.1)
    check = minimum_arl_feasibility(a.collective, b.collective, 1, 2e6)
    assert check["attainable"]
    a, b = entangled_pair(n=256, Na=0.1)
    check = minimum_arl_feasibility(a.collective, b.collective, 256, 2e6)
    assert not check["attainable"]
    assert check["minimum_arl_pulses"] > 1e30
    a, b = coherent_pair()
    d = b.mean - a.mean
    expected = 1 / norm.cdf(-abs(d) / (2 * np.sqrt(a.variance)))
    assert minimum_arl_feasibility(a, b, 1, 1000)["minimum_arl_pulses"] == pytest.approx(expected)


@pytest.mark.parametrize("kwargs", [{"n": 0}, {"horizon": 0}, {"change_at": 0}, {"threshold": -1}])
def test_invalid_conditional_inputs(kwargs):
    with pytest.raises(ValueError):
        simulate_conditional(**kwargs)


def test_zero_entanglement_detector_matches_ordinary_light():
    a = simulate_conditional(
        "coherent", 1, Na=0, threshold=3, trials=100, horizon=1000, change_at=20, seed=78
    )
    b = simulate_conditional(
        "entangled", 8, Na=0, threshold=3, trials=100, horizon=1000, change_at=20, seed=78
    )
    np.testing.assert_array_equal(a.stop, b.stop)
    np.testing.assert_array_equal(a.event, b.event)


def test_conditional_calibration_uses_distinct_direct_streams(monkeypatch):
    import lossqcd.conditional_detection as module

    calls = []

    def fake_simulate(*args, **kwargs):
        calls.append((args, kwargs))
        return kwargs["threshold"]

    monkeypatch.setattr(module, "simulate_conditional", fake_simulate)
    monkeypatch.setattr(
        module, "summarize_arl", lambda h: {"arl": float(np.exp(h)), "arl_se": 0.0, "censored": 0}
    )
    result = module.calibrate_conditional(
        "entangled", 8, target_arl=100, trials=20, seed=79, Na=0.5
    )
    assert result["threshold"] == pytest.approx(np.log(100), abs=0.003)
    assert all(args[0] == "entangled" and args[1] == 8 and args[3] == 0.5 for args, kw in calls)
    assert len({kw["seed"] for args, kw in calls}) == len(calls)


def test_conditional_calibration_rejects_censoring(monkeypatch):
    import lossqcd.conditional_detection as module

    monkeypatch.setattr(module, "simulate_conditional", lambda *args, **kw: None)
    monkeypatch.setattr(module, "summarize_arl", lambda runs: {"arl": None})
    with pytest.raises(RuntimeError, match="censored"):
        module.calibrate_conditional("entangled", 4, trials=20)
