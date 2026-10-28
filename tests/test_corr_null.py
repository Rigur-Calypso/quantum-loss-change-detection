import numpy as np
import pytest
from scipy.stats import norm
from lossqcd.corr_renewal import (
    sample_conditional_cycles,
    estimate_conditional_arl,
    summarize_conditional_cycles,
)
from lossqcd.corr_cusum import simulate_conditional
from lossqcd.entangled import entangled_pair, conditional_llrs
from lossqcd.simulation import summarize_arl
from lossqcd.null_channel import physical_laws
from lossqcd.arl_solver import mean_stopping_time


@pytest.mark.parametrize("proposal", ["pre", "post", "adaptive"])
def test_explicit_matched_initial_channel_preserves_arrays(proposal):
    args = dict(n=4, threshold=3, cycles=1000, seed=740, proposal=proposal)
    a = sample_conditional_cycles(**args)
    b = sample_conditional_cycles(**args, sampling_eta1=0.9)
    for field in ("length", "raw_sum", "page", "upper", "complete", "log_weight"):
        np.testing.assert_array_equal(getattr(a, field), getattr(b, field))
    assert a.parameters == b.parameters


@pytest.mark.parametrize("proposal", ["post", "adaptive"])
def test_density_weight_uses_actual_conditional_law_and_whole_mixture(monkeypatch, proposal):
    import lossqcd.corr_renewal as module

    values = np.array([9.8, 9.0, 9.8, 9.8])

    class Scripted:
        def __init__(self):
            self.index = 0

        def random(self, size):
            return np.zeros(size)

        def normal(self, mean, sd, size):
            result = np.full(size, values[self.index])
            self.index += 1
            return result

    monkeypatch.setattr(module.np.random, "default_rng", lambda seed: Scripted())
    batch = sample_conditional_cycles(
        n=4, threshold=100, cycles=2, max_steps=4, proposal=proposal, sampling_eta1=0.895
    )
    pre, post = entangled_pair(Na=0.1, n=4)
    actual, _ = entangled_pair(Na=0.1, n=4, eta1=0.895)
    increments = conditional_llrs(values, pre, post)
    page = 0.0
    expected = 0.0

    def density(block, k):
        cov = block.covariance
        mu = block.mean
        if k == 0:
            return norm.logpdf(values[k], mu[k], np.sqrt(cov[k, k]))
        coeff = np.linalg.solve(cov[:k, :k], cov[:k, k])
        mean = mu[k] + (values[:k] - mu[:k]) @ coeff
        variance = cov[k, k] - cov[k, :k] @ coeff
        return norm.logpdf(values[k], mean, np.sqrt(variance))

    for k in range(4):
        log0, log2 = density(actual, k), density(post, k)
        if proposal == "post":
            expected += log0 - log2
        else:
            beta = 1 / (1 + np.exp(2 - page))
            expected += log0 - np.logaddexp(np.log1p(-beta) + log0, np.log(beta) + log2)
        page = max(0, page + increments[k])
    np.testing.assert_allclose(batch.log_weight, expected, atol=1e-11)
    assert not np.allclose(batch.log_weight, -batch.raw_sum)


@pytest.mark.parametrize("proposal", ["post", "adaptive"])
def test_scalar_actual_null_estimate_matches_independent_solver(proposal):
    stats, _, _ = estimate_conditional_arl(
        n=1, threshold=3, cycles=80000, seed=741, sampling_eta1=0.899, proposal=proposal
    )
    pre, post, actual = physical_laws("entangled", 0.899)
    expected = mean_stopping_time(pre, post, 3, cells=512, sampling=actual)["mean_stopping_time"]
    assert abs(stats["arl"] - expected) < 4 * stats["arl_se"]


@pytest.mark.parametrize("n", [2, 4, 8])
def test_correlated_null_estimate_matches_direct_simulation(n):
    stats, _, _ = estimate_conditional_arl(
        n=n, threshold=2, cycles=80000, seed=742, sampling_eta1=0.899, proposal="adaptive"
    )
    direct = summarize_arl(
        simulate_conditional(
            n=n, threshold=2, trials=10000, horizon=30000, seed=743, sampling_eta1=0.899
        )
    )
    assert direct["arl"] is not None
    assert abs(stats["arl"] - direct["arl"]) < 4 * np.hypot(stats["arl_se"], direct["arl_se"])


def test_different_true_null_models_cannot_be_combined():
    a = sample_conditional_cycles(n=2, threshold=2, cycles=100, proposal="pre", sampling_eta1=0.899)
    b = sample_conditional_cycles(
        n=2, threshold=2, cycles=100, proposal="adaptive", sampling_eta1=0.901
    )
    with pytest.raises(ValueError, match="same detector"):
        summarize_conditional_cycles(a, b)
