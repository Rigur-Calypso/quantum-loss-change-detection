import numpy as np
import pytest
from scipy.stats import norm
from lossqcd.robust_thresholds import (
    GRID,
    calibrate,
    threshold_update,
    approximate_grid_lower_bound,
)


def test_threshold_update_raises_threshold_when_false_alarm_arl_is_too_short():
    assert threshold_update(10, 10, 1e6, 2e6) == pytest.approx(10 + np.log(2))
    assert threshold_update(10, 10, 4e6, 2e6) == 10


@pytest.mark.parametrize("values", [(0, 0, 0, 2e6), (-1, 2, 1e6, 2e6), (0, np.inf, 1e6, 2e6)])
def test_invalid_calibration_inputs(values):
    with pytest.raises(ValueError):
        threshold_update(*values)


def test_multiplicity_adjustment_is_more_conservative_than_single_point():
    a = approximate_grid_lower_bound(2.15e6, 20000, points=25)
    b = approximate_grid_lower_bound(2.15e6, 20000, points=1)
    assert a["lower"] < b["lower"] < 2.15e6
    assert norm.sf(a["z"]) * 25 == pytest.approx(0.05)
    assert a["lower"] == pytest.approx(2.15e6 * np.exp(-a["z"] * 20000 / 2.15e6))


@pytest.mark.parametrize("kwargs", [{"arl": 0}, {"se": -1}, {"points": True}, {"alpha": 1}])
def test_invalid_lower_bound_inputs(kwargs):
    args = {"arl": 2e6, "se": 10000} | kwargs
    with pytest.raises(ValueError):
        approximate_grid_lower_bound(**args)


def test_fixed_round_calibration_uses_all_grid_points_and_separate_seeds(monkeypatch):
    import lossqcd.robust_thresholds as module

    calls = []

    def fake_estimate(**kwargs):
        calls.append(kwargs)
        return (
            {
                "arl": float(np.exp(kwargs["threshold"]) * kwargs["sampling_eta1"]),
                "diagnostics_passed": True,
            },
            None,
            None,
        )

    monkeypatch.setattr(module, "estimate_conditional_arl", fake_estimate)
    model = {"kind": "entangled", "n": 4, "label": "Test", "calibration": {"threshold": 1.0}}
    result = calibrate(model, list(range(20)), cycles=20, design_target=1000)
    assert len(calls) == 20
    assert len({c["seed"] for c in calls}) == 20
    assert [c["sampling_eta1"] for c in calls] == list(GRID) * 4
    assert result["threshold"] == pytest.approx(np.log(1000 / min(GRID)))
    assert result["threshold"] == result["history"][-1]["threshold"]
    assert all(c["eta1"] == 0.9 and c["eta2"] == 0.85 and c["n"] == 4 for c in calls)


@pytest.mark.parametrize("seeds", [[1], [1] * 20])
def test_calibration_rejects_reused_or_missing_seed_slots(seeds):
    model = {"kind": "coherent", "n": 1, "label": "Test", "calibration": {"threshold": 1.0}}
    with pytest.raises(ValueError, match="independent seed"):
        calibrate(model, seeds, cycles=20)
