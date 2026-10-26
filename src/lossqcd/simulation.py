"""Single-shot CUSUM Monte Carlo with explicit false alarms and censoring.

The detector is never reset after an alarm. Each replicate stops at its FIRST
crossing. change_at=nu means observations 1,...,nu-1 are pre-change and nu onward
are post-change. None means no change. No sentinel can be mistaken for an event.
"""
from dataclasses import dataclass
import numpy as np
from .cusum import validate_threshold


@dataclass(frozen=True)
class Runs:
    stop: np.ndarray  # observed T, or horizon when censored
    event: np.ndarray  # True only when a threshold crossing was actually seen
    horizon: int
    change_at: int | None

    def __post_init__(self):
        stop = np.asarray(self.stop)
        event = np.asarray(self.event)
        if (not isinstance(self.horizon, (int, np.integer)) or self.horizon < 1
            or stop.ndim != 1 or stop.size < 1 or event.shape != stop.shape
            or event.dtype != bool or not np.issubdtype(stop.dtype, np.integer)
            or np.any(stop < 1) or np.any(stop > self.horizon)
            or np.any(stop[~event] != self.horizon)):
            raise ValueError("invalid stopping-time/event arrays or horizon")
        if self.change_at is not None and (not isinstance(self.change_at, (int, np.integer))
                or not 1 <= self.change_at <= self.horizon):
            raise ValueError("change_at must be within the observation window")


def simulate_stopping_times(pre, post, threshold, trials=2000, horizon=100000,
                            change_at=None, seed=20261020,
                            sampling_pre=None, sampling_post=None):
    """Sample Gaussian observations and score with the assumed pre/post laws.

    Optional sampling laws allow physical parameter mismatch. Defaults retain
    the original matched model and exactly the same seeded random draws.
    """
    validate_threshold(threshold)
    for name, value in (("trials", trials), ("horizon", horizon)):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if change_at is not None and (not isinstance(change_at, (int, np.integer))
                                  or not 1 <= change_at <= horizon):
        raise ValueError("change_at must be within the observation window")
    rng = np.random.default_rng(seed)
    actual_pre = pre if sampling_pre is None else sampling_pre
    actual_post = post if sampling_post is None else sampling_post
    stop = np.full(trials, horizon, dtype=np.int64)
    event = np.zeros(trials, dtype=bool)
    active = np.arange(trials)
    w = np.zeros(trials)
    for t in range(1, horizon + 1):
        dist = actual_post if change_at is not None and t >= change_at else actual_pre
        x = rng.normal(dist.mean, np.sqrt(dist.variance), active.size)
        w[active] = np.maximum(0.0, w[active] + post.llr(x, pre))
        crossed = w[active] > threshold
        stopped = active[crossed]
        stop[stopped], event[stopped] = t, True
        active = active[~crossed]
        if not active.size:
            break
    return Runs(stop, event, horizon, change_at)


def _mean_se(values):
    if not len(values):
        return None, None
    return (float(np.mean(values)),
            float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else None)


def summarize_arl(runs):
    """Uncensored sample mean estimates E_infinity[T].

    With censoring we expose only E[min(T,H)] and its Monte Carlo SE, a lower
    bound on ARL, never the mean of just the observed failures. Null means the
    requested unrestricted estimate is unavailable. SE is not a confidence bound.
    """
    if runs.change_at is not None:
        raise ValueError("ARL requires no-change runs")
    mean, se = _mean_se(runs.stop)
    censored = int(np.sum(~runs.event))
    return {"trials": len(runs.stop), "horizon": runs.horizon, "censored": censored,
            "arl": mean if censored == 0 else None,
            "arl_se": se if censored == 0 else None,
            "restricted_mean_min_T_H": mean, "restricted_mean_se": se}


def summarize_delay(runs):
    """E[T-nu | T>=nu], with zero delay for detection on observation nu.

    Reports P(T<nu) separately. Conditioning on survival is NOT Lorden's
    worst-case delay. With right censoring, unrestricted delay is unavailable;
    report E[min(T,H)-nu | T>=nu] as a restricted lower-bound estimate instead.
    """
    nu = runs.change_at
    if nu is None:
        raise ValueError("delay requires a finite change point")
    false = runs.event & (runs.stop < nu)
    survivors = ~false
    censored = int(np.sum(survivors & ~runs.event))
    mean, se = _mean_se(runs.stop[survivors] - nu)
    probability = float(np.mean(false))
    return {"trials": len(runs.stop), "change_at": nu, "horizon": runs.horizon,
            "false_alarms": int(false.sum()), "false_alarm_probability": probability,
            "false_alarm_probability_se": float(np.sqrt(probability*(1-probability)/len(false))),
            "survivors": int(survivors.sum()), "post_change_censored": censored,
            "conditional_delay": mean if censored == 0 else None,
            "conditional_delay_se": se if censored == 0 else None,
            "restricted_conditional_delay": mean, "restricted_conditional_delay_se": se}
