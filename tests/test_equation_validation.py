"""Independent receiver mixture identities, rather than stored plot values."""

import numpy as np
import pytest
from lossqcd.equation_validation import integrated_counts, count_moments, sample_counts


@pytest.mark.parametrize("eta", [0.9, 0.5, 0.2])
def test_integrated_photon_distribution_has_analytic_mixed_poisson_moments(eta):
    p, tail = integrated_counts(101, eta, eta)
    k = np.arange(len(p))
    mean, var = count_moments(101, eta, eta)
    assert p.sum() == pytest.approx(1, abs=1e-12)
    assert tail < 1e-12
    assert p @ k == pytest.approx(mean, rel=1e-10)
    assert p @ (k - mean) ** 2 == pytest.approx(var, rel=1e-10)


def test_mixed_poisson_extra_variance_disappears_at_zero_phase_noise():
    mean, var = count_moments(101, 0.85, 0.9, sigma=0)
    assert mean == pytest.approx(101 * (np.sqrt(0.85) - 1.01 * np.sqrt(0.9)) ** 2)
    assert var == mean


def test_complex_amplitude_sampler_matches_marginal_distribution():
    p, _ = integrated_counts(101, 0.9, 1.0)
    counts = sample_counts(101, 0.9, 1.0, np.random.default_rng(20261024), 150000)
    # Check identifiable event probabilities, avoiding meaningless tiny tail bins.
    for k in (0, 1, 2, 3):
        error = abs(np.mean(counts == k) - p[k])
        se = np.sqrt(p[k] * (1 - p[k]) / len(counts))
        assert error < 4 * se
