"""Paper Eqs. (2)–(3); all logarithms are natural, all widths are variances."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Gaussian:
    """Scalar homodyne law with finite mean and strictly positive variance.

    Methods use natural logarithms; post.llr(x, pre) scores a change, and
    post.kl(pre) measures expected post-change evidence per observation.
    """

    mean: float
    variance: float

    def __post_init__(self):
        if not np.isfinite(self.mean):
            raise ValueError("mean must be finite")
        if not np.isfinite(self.variance) or self.variance <= 0:
            raise ValueError("variance must be finite and positive")

    def logpdf(self, x):
        """Evaluate the Gaussian log density elementwise, preserving the input shape."""
        x = np.asarray(x, dtype=float)
        return -0.5 * (np.log(2 * np.pi * self.variance) + (x - self.mean) ** 2 / self.variance)

    def llr(self, x, pre):
        """log p_self(x)/p_pre(x); self is the post-change distribution.

        Centering at pre.mean avoids subtracting large squared means. Unlike
        the equal-variance linear formula this also includes the quadratic term.
        """
        y = np.asarray(x, dtype=float) - pre.mean
        delta = self.mean - pre.mean
        return (
            0.5 * np.log(pre.variance / self.variance)
            + 0.5 * y**2 * (1 / pre.variance - 1 / self.variance)
            + y * delta / self.variance
            - delta**2 / (2 * self.variance)
        )

    def kl(self, pre):
        """D(self || pre), in nats per observation."""
        ratio = self.variance / pre.variance
        return float(0.5 * (ratio - 1 - np.log(ratio) + (self.mean - pre.mean) ** 2 / pre.variance))


def _validate(N, Na, eta1, eta2):
    """Reject nonfinite energies, negative photon counts and transmission outside [0, 1]."""
    if not all(np.isfinite(x) for x in (N, Na, eta1, eta2)):
        raise ValueError("parameters must be finite")
    if N < 0 or Na < 0 or not (0 <= eta1 <= 1 and 0 <= eta2 <= 1):
        raise ValueError("energies must be nonnegative and transmissivities in [0,1]")


def coherent_pair(N=100.0, Na=0.0, eta1=0.9, eta2=0.85):
    """Equal total energy N+Na: alpha=sqrt(N+Na), vacuum variance=1/4."""
    _validate(N, Na, eta1, eta2)
    alpha = np.sqrt(N + Na)
    return (Gaussian(np.sqrt(eta1) * alpha, 0.25), Gaussian(np.sqrt(eta2) * alpha, 0.25))


def squeezed_pair(N=100.0, Na=0.0, eta1=0.9, eta2=0.85):
    """Independent displaced squeezed pulses: alpha=sqrt(N), sinh(r)^2=Na."""
    _validate(N, Na, eta1, eta2)
    r = np.arcsinh(np.sqrt(Na))

    def distribution(eta):
        return Gaussian(np.sqrt(eta * N), (eta * np.exp(-2 * r) + (1 - eta)) / 4)

    return distribution(eta1), distribution(eta2)


def coherent_kl(N, Na, eta1, eta2):
    """Paper Eq. (2), independent of the generic Gaussian KL implementation."""
    _validate(N, Na, eta1, eta2)
    return float(2 * (N + Na) * (np.sqrt(eta2) - np.sqrt(eta1)) ** 2)


def squeezed_kl(N, Na, eta1, eta2):
    """Paper Eq. (3); denominator v1 is the PRE-change variance."""
    pre, post = squeezed_pair(N, Na, eta1, eta2)
    ratio = post.variance / pre.variance
    return float(
        (ratio - np.log(ratio) - 1) / 2
        + N * (np.sqrt(eta2) - np.sqrt(eta1)) ** 2 / (2 * pre.variance)
    )
