import numpy as np
import pytest
from scipy.integrate import quad
from lossqcd.models import Gaussian, coherent_pair, squeezed_pair, coherent_kl, squeezed_kl


@pytest.mark.parametrize("pre,post", [
    (Gaussian(9.4, .25), Gaussian(9.1, .25)),
    (Gaussian(-1, .7), Gaussian(2, 1.9)),
    (Gaussian(0, 3), Gaussian(0, .2)),
])
def test_llr_and_kl_against_density_integration(pre, post):
    x = np.linspace(-7, 12, 79)
    np.testing.assert_allclose(post.llr(x, pre), post.logpdf(x)-pre.logpdf(x), atol=1e-12)
    # Independent numerical integration in standardized coordinates.
    phi = lambda z: np.exp(-z*z/2) / np.sqrt(2*np.pi)
    expectation = quad(lambda z: phi(z) * (post.logpdf(post.mean+np.sqrt(post.variance)*z)
                       - pre.logpdf(post.mean+np.sqrt(post.variance)*z)), -12, 12)[0]
    assert expectation == pytest.approx(post.kl(pre), abs=1e-10)
    reverse_expectation = quad(lambda z: phi(z) * post.llr(
        pre.mean+np.sqrt(pre.variance)*z, pre), -12, 12)[0]
    assert reverse_expectation == pytest.approx(-pre.kl(post), abs=1e-10)


@pytest.mark.parametrize("Na", [0, .01, .1, 1, 5])
def test_paper_equations_and_energy_budget(Na):
    pre, post = coherent_pair(Na=Na)
    assert pre.mean == pytest.approx(np.sqrt(.9)*np.sqrt(100+Na))
    assert post.mean == pytest.approx(np.sqrt(.85)*np.sqrt(100+Na))
    assert pre.variance == post.variance == .25
    assert post.kl(pre) == pytest.approx(coherent_kl(100, Na, .9, .85))
    pre, post = squeezed_pair(Na=Na)
    assert pre.mean == pytest.approx(np.sqrt(.9)*10)
    assert post.variance == pytest.approx((.85*np.exp(-2*np.arcsinh(np.sqrt(Na)))+.15)/4)
    assert post.kl(pre) == pytest.approx(squeezed_kl(100, Na, .9, .85))


def test_zero_squeezing_and_paper_ratio():
    assert squeezed_pair(Na=0) == coherent_pair(Na=0)
    assert squeezed_kl(100, .1, .9, .85)/coherent_kl(100, .1, .9, .85) == pytest.approx(1.72, rel=.01)
    assert Gaussian(0, 1).kl(Gaussian(0, 1)) == 0


@pytest.mark.parametrize("mean,variance", [(0, 0), (0, -1), (float('nan'), 1), (0, float('inf'))])
def test_invalid_gaussian(mean, variance):
    with pytest.raises(ValueError):
        Gaussian(mean, variance)


@pytest.mark.parametrize("kwargs", [{"N":-1}, {"Na":-1}, {"eta1":1.1}, {"eta2":-.1}, {"Na":np.nan}])
def test_invalid_optical_parameters(kwargs):
    with pytest.raises(ValueError):
        coherent_pair(**kwargs)
    with pytest.raises(ValueError):
        squeezed_pair(**kwargs)
