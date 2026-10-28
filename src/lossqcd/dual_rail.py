"""Appendix E.2 ideal dual-rail qubits, heralded loss and pair-clock detection.

Basis: input (early photon, late photon); output (vacuum, early, late).
Only total photon survival is heralded, not the occupied rail. This models
the ideal logical erasure channel, not the Duan-Kimble hardware dynamics or
the security of a QKD protocol. Readout is assumed available after BOTH bins.
"""

from dataclasses import dataclass
import numpy as np
from .paper_models import Bernoulli
from .cusum import validate_threshold
from .simulation import Runs, summarize_arl, summarize_delay


def validate_eta(eta):
    """Reject nonreal or nonfinite transmission outside the physical interval [0, 1]."""
    if not np.isscalar(eta) or not np.isreal(eta) or not np.isfinite(eta) or not 0 <= eta <= 1:
        raise ValueError("transmission must lie in [0,1]")


def density(state):
    """Validate and return a normalized positive Hermitian 2x2 qubit state."""
    rho = np.asarray(state, dtype=complex)
    if (
        rho.shape != (2, 2)
        or not np.all(np.isfinite(rho))
        or not np.allclose(rho, rho.conj().T, rtol=0, atol=1e-12)
        or not np.isclose(np.trace(rho), 1, rtol=0, atol=1e-12)
        or np.linalg.eigvalsh(rho).min() < -1e-12
    ):
        raise ValueError("a normalized positive Hermitian 2x2 density matrix is required")
    return rho


def pure_state(amplitudes):
    """Form a density matrix from two normalized complex rail amplitudes."""
    v = np.asarray(amplitudes, dtype=complex)
    if (
        v.shape != (2,)
        or not np.all(np.isfinite(v))
        or not np.isclose(np.vdot(v, v), 1, rtol=0, atol=1e-12)
    ):
        raise ValueError("two normalized finite qubit amplitudes required")
    return np.outer(v, v.conj())


def loss_kraus(eta_early, eta_late=None):
    """Return the survival and two loss operators in (vacuum, early, late) output order."""
    eta_late = eta_early if eta_late is None else eta_late
    validate_eta(eta_early)
    validate_eta(eta_late)
    return (
        np.array([[0, 0], [np.sqrt(eta_early), 0], [0, np.sqrt(eta_late)]], dtype=complex),
        np.array([[np.sqrt(1 - eta_early), 0], [0, 0], [0, 0]], dtype=complex),
        np.array([[0, np.sqrt(1 - eta_late)], [0, 0], [0, 0]], dtype=complex),
    )


def transmit(state, eta_early, eta_late=None):
    """Apply the trace-preserving pure-loss map, including the vacuum outcome."""
    rho = density(state)
    return sum(k @ rho @ k.conj().T for k in loss_kraus(eta_early, eta_late))


def heralded_state(state, eta_early, eta_late=None):
    """Return survival probability and surviving qubit state (None at zero)."""
    output = transmit(state, eta_early, eta_late)
    block = output[1:, 1:]
    probability = float(np.trace(block).real)
    return probability, None if probability == 0 else block / probability


@dataclass(frozen=True)
class DualRailRuns:
    """Pair-clock stopping records with a separate physical-bin change time."""

    pairs: Runs
    change_pulse: int | None

    @property
    def readout_pulses(self):
        """Convert pair alarms to the physical pulse at the end of each qubit."""
        return 2 * self.pairs.stop


def simulate_heralds(
    eta1,
    eta2,
    threshold,
    trials=30000,
    horizon_pairs=10000,
    change_pulse=None,
    encoding="haar",
    seed=2026102503,
):
    """Score binary survival once per qubit, at the end of its two time bins.

    Each rail undergoes its actual channel loss. A mid-pair change has survival
    p_early*eta_early+(1-p_early)*eta_late. The fixed detector uses Bernoulli
    eta1/eta2 likelihoods; it is not supplied the true change pulse or qubit.
    Haar-state early-bin population is uniform in [0,1]; phase does not affect
    total survival. Herald and encoding random streams are independent.
    """
    validate_eta(eta1)
    validate_eta(eta2)
    validate_threshold(threshold)
    if not 0 < eta1 < 1 or not 0 < eta2 < 1 or eta1 == eta2:
        raise ValueError("detector requires distinct interior transmissions")
    if any(
        isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x < 1
        for x in (trials, horizon_pairs)
    ):
        raise ValueError("positive integer trial/pair counts required")
    if change_pulse is not None and (
        isinstance(change_pulse, bool)
        or not isinstance(change_pulse, (int, np.integer))
        or not 1 <= change_pulse <= 2 * horizon_pairs
    ):
        raise ValueError("change pulse must lie in the physical observation window")
    if encoding not in ("haar", "early", "late", "plus"):
        raise ValueError("unknown qubit encoding")
    child = np.random.SeedSequence(seed).spawn(2)
    herald_rng, population_rng = [np.random.default_rng(s) for s in child]
    pre, post = Bernoulli(eta1), Bernoulli(eta2)
    active = np.arange(trials)
    page = np.zeros(trials)
    stop = np.full(trials, horizon_pairs, dtype=np.int64)
    event = np.zeros(trials, dtype=bool)
    for q in range(1, horizon_pairs + 1):
        early = eta2 if change_pulse is not None and 2 * q - 1 >= change_pulse else eta1
        late = eta2 if change_pulse is not None and 2 * q >= change_pulse else eta1
        if early == late:
            probability = early
        else:
            population = (
                population_rng.random(len(active))
                if encoding == "haar"
                else {"early": 1.0, "late": 0.0, "plus": 0.5}[encoding]
            )
            probability = population * early + (1 - population) * late
        observed = (herald_rng.random(len(active)) < probability).astype(int)
        page[active] = np.maximum(0, page[active] + post.llr(observed, pre))
        hit = page[active] > threshold
        stop[active[hit]] = q
        event[active[hit]] = True
        active = active[~hit]
        if not len(active):
            break
    nu = None if change_pulse is None else (change_pulse + 1) // 2
    return DualRailRuns(Runs(stop, event, horizon_pairs, nu), change_pulse)


def summarize_dual_rail(runs):
    """Keep pair-count statistics and physical-bin statistics explicitly separate."""
    if runs.change_pulse is None:
        pair_stats = summarize_arl(runs.pairs)
        physical = {
            key: (None if pair_stats[key] is None else 2 * pair_stats[key])
            for key in ("arl", "arl_se", "restricted_mean_min_T_H", "restricted_mean_se")
        }
    else:
        pair_stats = summarize_delay(runs.pairs)
        # An odd/even physical change has a different end-of-pair delay offset.
        offset = 2 * runs.pairs.change_at - runs.change_pulse
        physical = {
            key: (None if pair_stats[key] is None else 2 * pair_stats[key] + offset)
            for key in ("conditional_delay", "restricted_conditional_delay")
        }
        physical.update(
            {
                key: (None if pair_stats[key] is None else 2 * pair_stats[key])
                for key in ("conditional_delay_se", "restricted_conditional_delay_se")
            }
        )
        physical["false_alarms"] = pair_stats["false_alarms"]
        physical["post_change_censored"] = pair_stats["post_change_censored"]
    return {
        "pair_statistics": pair_stats,
        "physical_time_bin_statistics": physical,
        "change_pulse": runs.change_pulse,
        "readout_convention": "herald available at end of each two-bin qubit",
    }
