import pytest
from lossqcd.rare_validation import compare_estimate, main


def test_agreement_includes_both_independent_sampling_errors():
    estimate = {"arl": 1040.0, "arl_se": 2.0}
    assert not compare_estimate(estimate, 1000.0, reference_se=1.0)["agreement_passed"]
    result = compare_estimate(estimate, 1000.0, reference_se=20.0)
    assert result["agreement_passed"]
    assert result["combined_se"] == pytest.approx((2**2 + 20**2) ** 0.5)


def test_missing_historical_sources_refused_before_output_creation(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(SystemExit):
        main(["--output", str(out), "--moderate-source", str(tmp_path / "missing.json")])
    assert not out.exists()
