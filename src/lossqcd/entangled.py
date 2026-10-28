"""Entanglement-augmented real-homodyne statistics, paper Eqs. (4)–(6).

The seed squeezing is s=asinh(sqrt(n*Na)), not asinh(sqrt(Na)). The paper
switches between s and r in this discussion; Appendix C.1 fixes seed energy.
K_eta = I/4 + (v_eta-1/4) uu^T, u=ones(n)/sqrt(n),
v_eta=[1-eta+eta*exp(-2s)]/4. All n-1 orthogonal modes have variance 1/4,
zero mean, and identical pre/post laws. Y=u^T X is therefore sufficient for
the homogeneous pre/post BLOCK likelihood ratio. KL/n is per pulse.
"""

from dataclasses import dataclass
import numpy as np
from .models import Gaussian, _validate


def validate_block_size(n):
    """Require a positive integer pulse count; booleans are not block sizes."""
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("block size must be a positive integer")


@dataclass(frozen=True)
class SymmetricGaussianBlock:
    """Uniform-mean Gaussian block with one squeezed collective mode.

    The other n-1 modes have vacuum variance 1/4. Homogeneous block LLRs
    therefore depend only on sum(x)/sqrt(n), while sampling retains all modes.
    """

    n: int
    pulse_mean: float
    collective_variance: float

    def __post_init__(self):
        validate_block_size(self.n)
        if not np.isfinite(self.pulse_mean):
            raise ValueError("pulse mean must be finite")
        if not np.isfinite(self.collective_variance) or self.collective_variance <= 0:
            raise ValueError("collective variance must be positive and finite")

    @property
    def collective(self):
        """Return the sufficient scalar law for the normalized sum of all n pulses."""
        return Gaussian(float(np.sqrt(self.n) * self.pulse_mean), self.collective_variance)

    @property
    def mean(self):
        return np.full(self.n, self.pulse_mean)

    @property
    def covariance(self):
        """Return the full rank-one covariance in pulse coordinates."""
        return np.eye(self.n) * 0.25 + np.full(
            (self.n, self.n), (self.collective_variance - 0.25) / self.n
        )

    def _observations(self, x):
        x = np.asarray(x, dtype=float)
        if x.ndim < 1 or x.shape[-1] != self.n or not np.all(np.isfinite(x)):
            raise ValueError("observations must be finite with last dimension n")
        return x

    def logpdf(self, x):
        """Evaluate the joint density by separating collective and vacuum modes."""
        x = self._observations(x)
        centered = x - self.pulse_mean
        collective_center = centered.sum(axis=-1) / np.sqrt(self.n)
        perpendicular = np.maximum(0.0, np.sum(centered**2, axis=-1) - collective_center**2)
        return -0.5 * (
            np.log(2 * np.pi * self.collective_variance)
            + collective_center**2 / self.collective_variance
        ) - 0.5 * ((self.n - 1) * np.log(2 * np.pi * 0.25) + perpendicular / 0.25)

    def llr(self, x, pre):
        """Score homogeneous block hypotheses through their sufficient collective mode."""
        if self.n != pre.n:
            raise ValueError("block dimensions must match")
        x = self._observations(x)
        return self.collective.llr(x.sum(axis=-1) / np.sqrt(self.n), pre.collective)

    def kl(self, pre):
        """Return joint-block KL in nats; divide by n for per-pulse information."""
        if self.n != pre.n:
            raise ValueError("block dimensions must match")
        return self.collective.kl(pre.collective)

    def sample(self, rng, size=1):
        """Draw full correlated blocks, returning an array with shape (size, n)."""
        if not isinstance(size, (int, np.integer)) or size < 1:
            raise ValueError("sample size must be a positive integer")
        z = rng.normal(size=(size, self.n))
        mean = z.mean(axis=-1, keepdims=True)
        return self.pulse_mean + 0.5 * (z - mean) + np.sqrt(self.collective_variance) * mean


def entangled_pair(N=100.0, Na=0.0, eta1=0.9, eta2=0.85, n=2):
    """Construct pre/post block laws with seed energy n*Na and displacement energy N."""
    _validate(N, Na, eta1, eta2)
    validate_block_size(n)
    if not np.isfinite(n * Na):
        raise ValueError("total seed energy n*Na must be finite")
    s = np.arcsinh(np.sqrt(n * Na))

    def state(eta):
        return SymmetricGaussianBlock(
            n, float(np.sqrt(eta * N)), float((eta * np.exp(-2 * s) + (1 - eta)) / 4)
        )

    return state(eta1), state(eta2)


def entangled_kl_per_pulse(N, Na, eta1, eta2, n):
    """Independent closed form, derived by diagonalizing Eqs. (4)–(6)."""
    _validate(N, Na, eta1, eta2)
    validate_block_size(n)
    s = np.arcsinh(np.sqrt(n * Na))
    v1 = (eta1 * np.exp(-2 * s) + (1 - eta1)) / 4
    v2 = (eta2 * np.exp(-2 * s) + (1 - eta2)) / 4
    if v1 <= 0 or v2 <= 0:
        raise ValueError("covariances must be positive definite")
    ratio = v2 / v1
    return float(
        (ratio - 1 - np.log(ratio)) / (2 * n) + N * (np.sqrt(eta2) - np.sqrt(eta1)) ** 2 / (2 * v1)
    )


def mixed_block_moments(N, Na, etas):
    """Output moments when transmission changes WITHIN an entangled block.

    Independent vacuum loss noise restores the I/4 term; the rank-one
    direction becomes sqrt(etas), not ones(n). No pre/post independence within
    a straddling block is assumed.
    """
    etas = np.asarray(etas, dtype=float)
    if (
        etas.ndim != 1
        or etas.size < 1
        or np.any(~np.isfinite(etas))
        or np.any((etas < 0) | (etas > 1))
    ):
        raise ValueError("etas must be a nonempty finite vector in [0,1]")
    _validate(N, Na, float(etas.min()), float(etas.max()))
    n = len(etas)
    s = np.arcsinh(np.sqrt(n * Na))
    direction = np.sqrt(etas)
    mean = np.sqrt(N) * direction
    covariance = np.eye(n) * 0.25 + np.expm1(-2 * s) / (4 * n) * np.outer(direction, direction)
    return mean, covariance


def mixed_collective(N, Na, n, post_pulses, eta1=0.9, eta2=0.85):
    """Exact law of Y=sum(X)/sqrt(n) in a straddling block.

    The detector still scores the homogeneous pre/post block hypotheses;
    this is not a within-block generalized likelihood-ratio detector.
    """
    _validate(N, Na, eta1, eta2)
    validate_block_size(n)
    if not isinstance(post_pulses, (int, np.integer)) or not 0 <= post_pulses <= n:
        raise ValueError("post_pulses must be an integer between 0 and n")
    root_average = ((n - post_pulses) * np.sqrt(eta1) + post_pulses * np.sqrt(eta2)) / n
    s = np.arcsinh(np.sqrt(n * Na))
    return Gaussian(
        float(np.sqrt(n * N) * root_average), float((1 + np.expm1(-2 * s) * root_average**2) / 4)
    )


def conditional_parameters(block, past_sum, seen):
    """Appendix B.4 (16)–(18), specialized to uniform rank-one covariance."""
    if not isinstance(seen, (int, np.integer)) or not 0 <= seen < block.n:
        raise ValueError("seen must be an integer in [0,n-1]")
    c = (block.collective_variance - 0.25) / block.n
    denominator = 0.25 + seen * c
    mean = block.pulse_mean + c / denominator * (np.asarray(past_sum) - seen * block.pulse_mean)
    variance = 0.25 * (0.25 + (seen + 1) * c) / denominator
    return mean, float(variance)


def conditional_llrs(x, pre, post):
    """Per-pulse conditional LLRs; their unreflected sum is the joint block LLR."""
    if pre.n != post.n:
        raise ValueError("block dimensions must match")
    x = pre._observations(x)
    increments = np.empty_like(x)
    past = np.zeros(x.shape[:-1])
    for k in range(pre.n):
        m1, v1 = conditional_parameters(pre, past, k)
        m2, v2 = conditional_parameters(post, past, k)
        z = x[..., k]
        increments[..., k] = 0.5 * (np.log(v1 / v2) + (z - m1) ** 2 / v1 - (z - m2) ** 2 / v2)
        past = past + z
    return increments


def mixed_next_parameters(alpha, q, eta, previous_eta_sum, weighted_residual_sum):
    """True next-pulse conditional law for K=I/4+q sqrt(eta)sqrt(eta)^T."""
    denominator = 0.25 + q * previous_eta_sum
    mean = alpha * np.sqrt(eta) + q * np.sqrt(eta) / denominator * weighted_residual_sum
    variance = 0.25 * (0.25 + q * (previous_eta_sum + eta)) / denominator
    return mean, float(variance)


def sample_mixed_blocks(N, Na, etas, rng, size=1):
    """Generate physical mixed-transmission blocks sequentially for validation."""
    mean, _ = mixed_block_moments(N, Na, etas)
    etas = np.asarray(etas, dtype=float)
    n = len(etas)
    if not isinstance(size, (int, np.integer)) or size < 1:
        raise ValueError("size must be a positive integer")
    q = np.expm1(-2 * np.arcsinh(np.sqrt(n * Na))) / (4 * n)
    weighted = np.zeros(size)
    eta_sum = 0.0
    out = np.empty((size, n))
    for k, eta in enumerate(etas):
        mu, var = mixed_next_parameters(np.sqrt(N), q, eta, eta_sum, weighted)
        x = rng.normal(mu, np.sqrt(var), size)
        out[:, k] = x
        weighted += np.sqrt(eta) * (x - mean[k])
        eta_sum += eta
    return out
