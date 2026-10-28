"""Scalar renewal ARL with separate assumed score and actual null law.

Cycles end when the unreflected assumed log-likelihood sum is <=0 or >h.
ARL=E_actual_null[length]/P_actual_null[upper]. Numerator cycles are drawn
from the actual null. Independent proposal cycles use a chosen Gaussian law.
Their weight is exp(sum(log p_actual_null - log p_proposal)), NOT in general
exp(-sum(assumed score)). All failures contribute zero to the denominator.
Weights are ordinary, not self-normalized; no clipping or bound is assumed.
Only iid scalar observations apply. No correlated-block claim is made.
"""

from dataclasses import dataclass
import numpy as np
from .cusum import validate_threshold
from .renewal import _mean_se


@dataclass(frozen=True)
class MismatchCycles:
    """Stopped assumed-score paths with a distinct actual-null/proposal log weight."""

    length: np.ndarray
    score: np.ndarray
    log_weight: np.ndarray
    upper: np.ndarray
    complete: np.ndarray


def sample_mismatch_cycles(
    pre, post, actual_null, proposal, threshold, cycles=200000, seed=2026102801, max_steps=100000
):
    """Sample Gaussian cycles using assumed scores and actual-null likelihood weights.

    The same observation updates both quantities, but they coincide only in
    the matched post-proposal case. Do not replace log_weight with -score.
    """
    validate_threshold(threshold)
    if post.kl(pre) <= 0:
        raise ValueError("distinct score laws required")
    for value in (cycles, max_steps):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or value < 2
        ):
            raise ValueError("cycles and max_steps must be integers >=2")
    rng = np.random.default_rng(seed)
    score = np.zeros(cycles)
    weight = np.zeros(cycles)
    length = np.full(cycles, max_steps, dtype=np.int64)
    upper = np.zeros(cycles, dtype=bool)
    complete = upper.copy()
    active = np.arange(cycles)
    for t in range(1, max_steps + 1):
        x = rng.normal(proposal.mean, np.sqrt(proposal.variance), len(active))
        score[active] += post.llr(x, pre)
        # Gaussian.llr evaluates this density difference stably.
        weight[active] += actual_null.llr(x, proposal)
        hit = score[active] > threshold
        finished = hit | (score[active] <= 0)
        upper[active[hit]] = True
        complete[active[finished]] = True
        length[active[finished]] = t
        active = active[~finished]
        if not len(active):
            break
    return MismatchCycles(length, score, weight, upper, complete)


def summarize_mismatch_cycles(numerator, denominator):
    """Form the renewal ratio from complete null/proposal ensembles and report weight diagnostics."""
    for batch in (numerator, denominator):
        size = len(batch.length)
        if size < 2 or any(
            a.shape != (size,) for a in (batch.score, batch.log_weight, batch.upper, batch.complete)
        ):
            raise ValueError("invalid cycle arrays")
        if not batch.complete.all():
            raise RuntimeError("unfinished cycles: unrestricted ARL unavailable")
        if not np.all(np.isfinite(batch.score)) or not np.all(np.isfinite(batch.log_weight)):
            raise ValueError("nonfinite score or likelihood weight")
    if not np.all(numerator.log_weight == 0):
        raise ValueError("numerator must be sampled under the actual null law")
    weights = np.zeros(len(denominator.length))
    with np.errstate(over="raise", under="ignore", invalid="raise"):
        weights[denominator.upper] = np.exp(denominator.log_weight[denominator.upper])
    mean, se = _mean_se(numerator.length)
    prob, pse = _mean_se(weights)
    if prob <= 0:
        raise RuntimeError("no positive weighted successes")
    arl = mean / prob
    rse = float(np.hypot(se / mean, pse / prob))
    scaled = weights / weights.max()
    ess = float(scaled.sum() ** 2 / (scaled @ scaled))
    fraction = float(scaled.max() / scaled.sum())
    return {
        "arl": float(arl),
        "arl_se": float(arl * rse),
        "relative_se": rse,
        "approx_95_interval": [float(arl * np.exp(-1.96 * rse)), float(arl * np.exp(1.96 * rse))],
        "mean_null_cycle_length": mean,
        "mean_null_cycle_length_se": se,
        "cycle_alarm_probability": prob,
        "cycle_alarm_probability_se": pse,
        "success_weight_ess": ess,
        "largest_weight_fraction": fraction,
        "diagnostics_passed": ess >= 1000 and fraction <= 0.01 and rse <= 0.02,
        "numerator_cycles": len(numerator.length),
        "proposal_cycles": len(denominator.length),
        "proposal_successes": int(denominator.upper.sum()),
        "uncertainty": "independent-means delta SE and approximate log-scale interval; diagnostics are not a finite-variance proof",
    }


def estimate_mismatch_arl(pre, post, actual_null, threshold, cycles=200000, seed=2026102801):
    """Use independent streams for actual-null cycle lengths and post-proposal alarm weights."""
    seeds = [
        int(s.generate_state(1, dtype=np.uint64)[0]) for s in np.random.SeedSequence(seed).spawn(2)
    ]
    numerator = sample_mismatch_cycles(
        pre, post, actual_null, actual_null, threshold, cycles, seeds[0]
    )
    denominator = sample_mismatch_cycles(pre, post, actual_null, post, threshold, cycles, seeds[1])
    result = summarize_mismatch_cycles(numerator, denominator)
    result["stream_seeds"] = seeds
    return result, numerator, denominator


def cycle_arrays(numerator, denominator):
    """Serialize numerator and proposal arrays with their distinct score and weight fields."""
    return {
        f"{prefix}_{field}": getattr(batch, field)
        for prefix, batch in (("null", numerator), ("proposal", denominator))
        for field in ("length", "score", "log_weight", "upper", "complete")
    }


def read_cycles(data, prefix):
    """Reconstruct a cycle record from an NPZ prefix without dropping event flags."""
    return MismatchCycles(
        *(
            data[f"{prefix}_{field}"]
            for field in ("length", "score", "log_weight", "upper", "complete")
        )
    )
