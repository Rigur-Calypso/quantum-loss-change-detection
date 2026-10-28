"""Independent physical/numerical checks at the PDF legend settings."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import norm, poisson
from lossqcd.receivers import PhaseKennedy
from lossqcd.verified_figures import fixed_seed, figure12


@pytest.mark.parametrize("k", [0, 1, 3, 10])
def test_published_phase_variance_poisson_mixture_against_adaptive_integral(k):
    receiver = PhaseKennedy(101, 0.9, 1, 0.01, 0.1, 160)
    expected = quad(
        lambda z: (
            norm.pdf(z) * poisson.pmf(k, 101 * abs(np.sqrt(0.9) - 1.01 * np.exp(0.1j * z)) ** 2)
        ),
        -12,
        12,
        epsabs=1e-13,
        epsrel=1e-11,
    )[0]
    assert np.exp(receiver.logpmf(k)) == pytest.approx(expected, rel=1e-10, abs=1e-13)
    default = PhaseKennedy(101, 0.9, 1, 0.01, 0.1)
    assert np.exp(default.logpmf(k)) == pytest.approx(expected, rel=1e-8, abs=1e-12)


@pytest.mark.parametrize("db", [0, 6, 12, 18])
def test_fixed_seed_dense_gaussian_kl_and_ordering(db):
    from lossqcd.entangled import entangled_pair

    values = []
    for n in (1, 2, 4, 8):
        row = fixed_seed(db, n)
        pre, post = entangled_pair(Na=row["Na"], n=n)
        delta = post.mean - pre.mean
        dense = (
            0.5
            * (
                np.trace(np.linalg.solve(pre.covariance, post.covariance))
                - n
                + delta @ np.linalg.solve(pre.covariance, delta)
                + np.linalg.slogdet(pre.covariance)[1]
                - np.linalg.slogdet(post.covariance)[1]
            )
            / n
        )
        assert row["kl_per_pulse"] == pytest.approx(dense, rel=1e-12)
        assert row["kl_per_pulse"] == pytest.approx(row["independent_A_over_n_plus_B"], rel=1e-12)
        values.append(dense)
    assert np.all(np.diff(values) <= 1e-12)


def test_weak_probe_data_processing_and_zero_loss_boundary():
    from lossqcd.paper_models import Bernoulli

    assert np.isinf(Bernoulli(0.9).kl(Bernoulli(1.0)))
    for loss in (0.05, 1, 3, 5):
        row = figure12(loss)
        assert row["photon_counting"] >= row["photon_homodyne"]
        for residual in ("1e4", "1e5"):
            assert row[f"kennedy_pnr_{residual}"] >= row[f"kennedy_onoff_{residual}"]
