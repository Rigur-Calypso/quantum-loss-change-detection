import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import kstest
from lossqcd.models import Gaussian, coherent_pair, squeezed_pair
from lossqcd.probes import (
    BPSK,
    bpsk_pair,
    bpsk_capacity,
    DisplacedPhoton,
    Bernoulli,
    Poisson,
    poisson_kl,
    kennedy_residual,
)
from lossqcd.receivers import PhaseHomodyne, PhaseKennedy, phase_quadrature
from lossqcd.iid_cusum import simulate_iid, estimate_iid_arl
from lossqcd.simulation import simulate_stopping_times, summarize_arl


@pytest.mark.parametrize("mean,var", [(0, 1), (0.4, 0.7), (9, 0.25)])
def test_bpsk_density_normalized_and_symmetric(mean, var):
    law = BPSK(mean, var)
    assert quad(
        lambda x: np.exp(law.logpdf(x)), -mean - 12 * np.sqrt(var), mean + 12 * np.sqrt(var)
    )[0] == pytest.approx(1, abs=1e-10)
    x = np.linspace(-12, 12, 41)
    np.testing.assert_allclose(law.logpdf(x), law.logpdf(-x))


def test_low_snr_bpsk_is_not_oracle_sign_recovery():
    pre, post = Gaussian(0.4, 0.4), Gaussian(0.2, 0.6)
    a, b = bpsk_pair(pre, post)
    assert 0 < b.kl(a) < post.kl(pre)
    assert b.llr(0.1, a) != pytest.approx(post.llr(0.1, pre))


def test_high_snr_bpsk_kl_and_score_match_unmodulated_law():
    pre, post = squeezed_pair(Na=0.1)
    a, b = bpsk_pair(pre, post)
    assert b.kl(a) == pytest.approx(post.kl(pre), abs=1e-12)
    np.testing.assert_allclose(
        b.llr(np.array([-9.0, 9.0, -10.0]), a), post.llr(np.array([9.0, 9.0, 10.0]), pre)
    )


def test_bpsk_capacity_limits_and_monotonicity():
    assert bpsk_capacity(0, 1) == pytest.approx(0, abs=1e-12)
    values = [bpsk_capacity(x, 1) for x in (0.1, 0.5, 1, 2, 10)]
    assert np.all(np.diff(values) > 0)
    assert values[-1] == pytest.approx(1, abs=1e-12)


@pytest.mark.parametrize("eta", [0, 0.3, 0.9, 1])
def test_photon_density_cdf_and_exact_sampler(eta):
    law = DisplacedPhoton(2, eta)
    assert quad(lambda x: np.exp(law.logpdf(x)), law.mean - 10, law.mean + 10, points=[law.mean])[
        0
    ] == pytest.approx(1, abs=1e-10)
    x = law.sample(np.random.default_rng(811), 60000)
    assert abs(x.mean() - law.mean) < 4 * np.sqrt(law.variance / len(x))
    assert np.var(x) == pytest.approx(law.variance, rel=0.025)
    assert kstest(x, law.cdf).pvalue > 0.001
    y = law.mean + 0.23
    h = 1e-5
    assert (law.cdf(y + h) - law.cdf(y - h)) / (2 * h) == pytest.approx(
        np.exp(law.logpdf(y)), rel=1e-7
    )


def test_dv_paper_entropy_benchmark_and_pure_vacuum_limit():
    a, b = DisplacedPhoton(10, 0.9), DisplacedPhoton(10, 0.85)
    assert b.kl(a) == pytest.approx(0.2063, abs=0.0001)
    assert DisplacedPhoton(0, 0).kl(DisplacedPhoton(0, 0)) == pytest.approx(0, abs=1e-12)


def test_discrete_relative_entropy_boundary_cases():
    assert Bernoulli(0.85).kl(Bernoulli(0.9)) > 0
    assert np.isinf(Bernoulli(0.1).kl(Bernoulli(0)))
    assert Bernoulli(0).kl(Bernoulli(0)) == 0
    assert poisson_kl(0, 0) == 0
    assert poisson_kl(0, 2) == 2
    assert np.isinf(kennedy_residual(100, 1, 0.9, 0.85, 0))


def test_generic_gaussian_detector_preserves_existing_seeded_arrays():
    a, b = coherent_pair()
    x = simulate_iid(a, b, 3, trials=500, horizon=2000, change_at=40, seed=812)
    y = simulate_stopping_times(a, b, 3, trials=500, horizon=2000, change_at=40, seed=812)
    np.testing.assert_array_equal(x.stop, y.stop)
    np.testing.assert_array_equal(x.event, y.event)


@pytest.mark.parametrize(
    "pair", [(BPSK(0.8, 0.25), BPSK(0.5, 0.3)), (DisplacedPhoton(2, 0.9), DisplacedPhoton(2, 0.7))]
)
def test_nongaussian_renewal_agrees_with_direct_false_alarms(pair):
    a, b = pair
    stats, _, _ = estimate_iid_arl(a, b, 2, cycles=80000, seed=813)
    direct = summarize_arl(simulate_iid(a, b, 2, trials=10000, horizon=100000, seed=814))
    assert direct["arl"] is not None
    assert abs(stats["arl"] - direct["arl"]) < 4 * np.hypot(stats["arl_se"], direct["arl_se"])


def test_phase_quadrature_normalization_and_truncated_second_moment():
    theta, w = phase_quadrature(0.02, 96)
    assert w.sum() == pytest.approx(1)
    assert np.sum(w * theta**2) == pytest.approx(0.02**2, rel=1e-10)


def test_zero_phase_homodyne_reduces_to_gaussian_model():
    pre, post = squeezed_pair(Na=0.2)
    a, b = PhaseHomodyne(100, 0.2, 0.9), PhaseHomodyne(100, 0.2, 0.85)
    np.testing.assert_allclose(a.logpdf([9, 10]), pre.logpdf([9, 10]))
    assert b.kl(a) == pytest.approx(post.kl(pre), abs=1e-10)


def test_zero_phase_kennedy_reduces_to_poisson_and_ideal_null_is_singular():
    a = PhaseKennedy(101, 0.9, 0.9, 0.01)
    b = PhaseKennedy(101, 0.85, 0.9, 0.01)
    assert b.kl(a) == pytest.approx(poisson_kl(b.components[0][0], a.components[0][0]), rel=1e-10)
    assert np.isinf(PhaseKennedy(101, 0.85, 0.9).kl(PhaseKennedy(101, 0.9, 0.9)))


def test_phase_mixture_quadrature_refinement():
    a, b = PhaseHomodyne(100, 0.2, 0.9, 0.03), PhaseHomodyne(100, 0.2, 0.85, 0.03)
    c, d = (
        PhaseHomodyne(100, 0.2, 0.9, 0.03, order=160),
        PhaseHomodyne(100, 0.2, 0.85, 0.03, order=160),
    )
    assert b.kl(a) == pytest.approx(d.kl(c), rel=1e-6)
    e, f = PhaseKennedy(101, 0.9, 0.9, 0.01, 0.03), PhaseKennedy(101, 0.85, 0.9, 0.01, 0.03)
    g, h = (
        PhaseKennedy(101, 0.9, 0.9, 0.01, 0.03, 160),
        PhaseKennedy(101, 0.85, 0.9, 0.01, 0.03, 160),
    )
    assert f.kl(e) == pytest.approx(h.kl(g), rel=1e-6)


@pytest.mark.parametrize("n", [1, 2, 4, 8, 32])
def test_hadamard_network_covariance_matches_analytic_entanglement(n):
    from lossqcd.probes import green_machine_covariance
    from lossqcd.entangled import entangled_pair

    pre, _ = entangled_pair(Na=0.2, n=n)
    np.testing.assert_allclose(green_machine_covariance(n, 0.2, 0.9), pre.covariance, atol=1e-15)


def test_shared_path_thresholds_keep_overshoot_and_strict_equality(monkeypatch):
    import lossqcd.iid_cusum as module

    values = iter((2.5, 0.7))  # scores 2.0, then 0.2
    monkeypatch.setattr(module, "draw", lambda law, rng, size: np.full(size, next(values)))
    runs = module.simulate_threshold_grid(
        Gaussian(0, 1), Gaussian(1, 1), [0, 0.5, 2, 2.1], trials=1, horizon=2
    )
    assert [r.stop[0] for r in runs] == [1, 1, 2, 2]
    assert all(r.event[0] for r in runs)


def test_shared_path_grid_retains_censoring_and_nested_alarm_times():
    from lossqcd.iid_cusum import simulate_threshold_grid

    a, b = coherent_pair()
    runs = simulate_threshold_grid(a, b, [1, 2, 3], trials=300, horizon=100, seed=815)
    assert np.all(np.diff(np.array([r.stop for r in runs]), axis=0) >= 0)
    assert np.all(runs[-1].event <= runs[0].event)
    assert summarize_arl(runs[-1])["arl"] is None


def test_main_text_change_point_argmin_includes_zero_prefix_and_earliest_tie():
    from lossqcd.cusum import detect_with_changepoint

    assert detect_with_changepoint([2, 2], 3) == {
        "stop": 2,
        "estimated_change_at": 1,
        "statistic": 4.0,
    }
    result = detect_with_changepoint([-2, 1, -1, 2, 3], 3)
    assert result == {"stop": 5, "estimated_change_at": 2, "statistic": 5.0}
    assert detect_with_changepoint([-1, -2], 3)["estimated_change_at"] is None


@pytest.mark.parametrize("Na", [0.001, 0.1, 2])
def test_squeezing_derivative_matches_finite_difference_of_eq3(Na):
    from lossqcd.probes import squeezed_energy_derivative
    from lossqcd.models import squeezed_kl

    step = Na * 1e-4
    expected = (squeezed_kl(100, Na + step, 0.9, 0.85) - squeezed_kl(100, Na - step, 0.9, 0.85)) / (
        2 * step
    )
    assert squeezed_energy_derivative(100, Na, 0.9, 0.85) == pytest.approx(expected, rel=1e-7)


def test_augmentation_crossover_retains_small_variance_kl_term():
    from lossqcd.probes import augmentation_crossover
    from lossqcd.models import squeezed_kl, coherent_kl

    result = augmentation_crossover()
    a = result["exact_Na"]
    assert result["approximate_Na"] == pytest.approx(900)
    assert squeezed_kl(100, a, 0.9, 0.85) == pytest.approx(coherent_kl(100, a, 0.9, 0.85))
    assert squeezed_kl(100, a * 1.1, 0.9, 0.85) < coherent_kl(100, a * 1.1, 0.9, 0.85)


@pytest.mark.parametrize("pair", [(Poisson(0.3), Poisson(0.8)), (Bernoulli(0.9), Bernoulli(0.81))])
def test_discrete_receiver_renewal_matches_direct_page_detection(pair):
    pre, post = pair
    stats, _, _ = estimate_iid_arl(pre, post, 2, cycles=100000, seed=816)
    direct = summarize_arl(simulate_iid(pre, post, 2, trials=8000, horizon=100000, seed=817))
    assert direct["arl"] is not None
    assert abs(stats["arl"] - direct["arl"]) < 4 * np.hypot(stats["arl_se"], direct["arl_se"])


def test_ideal_kennedy_null_never_false_alarms_but_waits_for_post_change_click():
    from lossqcd.simulation import summarize_delay

    pre, post = Poisson(0), Poisson(0.2)
    null = simulate_iid(pre, post, 2, trials=100, horizon=1000, seed=818)
    assert not null.event.any()
    changed = summarize_delay(
        simulate_iid(pre, post, 2, trials=40000, horizon=1000, change_at=10, seed=819)
    )
    assert changed["post_change_censored"] == 0 and changed["false_alarms"] == 0
    assert (
        abs(changed["conditional_delay"] - 1 / np.expm1(0.2)) < 4 * changed["conditional_delay_se"]
    )


def test_kennedy_binary_and_photon_resolved_observations_remain_distinct():
    a = PhaseKennedy(101, 0.9, 0.9, 0.01, 0.02)
    b = PhaseKennedy(101, 0.85, 0.9, 0.01, 0.02)
    assert 0 < b.onoff_kl(a) < b.kl(a)
    assert kennedy_residual(100, 1, 0.9, 0.85, 0.01, "onoff") < kennedy_residual(
        100, 1, 0.9, 0.85, 0.01, "pnr"
    )
