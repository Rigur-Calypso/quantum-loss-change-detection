import numpy as np
import pytest
from lossqcd.conditional_large_arl import calibrate_rare, assess_match, main


def test_calibration_has_fixed_budget_and_distinct_streams(monkeypatch):
    import lossqcd.conditional_large_arl as module

    calls = []

    def estimate(kind, n, h, **kwargs):
        calls.append((kind, n, h, kwargs))
        return {"arl": float(np.exp(h))}, None, None

    monkeypatch.setattr(module, "estimate_conditional_arl", estimate)
    result = calibrate_rare("entangled", 4, target=np.exp(10), cycles=100, seed=45, Na=0.1)
    assert result["threshold"] == pytest.approx(10)
    assert len(calls) == 4
    assert len({kw["seed"] for *_, kw in calls}) == 4
    assert all(kw["cycles"] == 100 for *_, kw in calls)


def test_matching_requires_entire_target_interval_and_weight_diagnostics():
    base = {
        "arl": 2e6,
        "arl_se": 1000,
        "diagnostics_passed": True,
        "approx_95_interval": [1.998e6, 2.002e6],
    }
    assert assess_match(base, base, 2e6)["passed"]
    assert not assess_match(base | {"diagnostics_passed": False}, base, 2e6)["passed"]
    assert not assess_match(base | {"approx_95_interval": [1.89e6, 2.11e6]}, base, 2e6)["passed"]


@pytest.mark.parametrize(
    "args",
    [
        ["--Na", "-1"],
        ["--validation-cycles", "1"],
        ["--target-arl", "nan"],
        ["--change-at", "0"],
        ["--change-at", "1", "1"],
        ["--tolerance", "2"],
    ],
)
def test_bad_cli_does_not_create_output(args, tmp_path):
    out = tmp_path / "output"
    with pytest.raises(SystemExit):
        main(["--output", str(out)] + args)
    assert not out.exists()
