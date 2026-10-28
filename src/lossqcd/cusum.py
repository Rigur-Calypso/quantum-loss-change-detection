"""Page recursion W_0=0, W_t=max(0,W_{t-1}+ell_t), T=inf{t:W_t>h}."""

import numpy as np


def validate_threshold(threshold):
    """Require a finite, nonnegative Page threshold; alarms use a strict crossing."""
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be finite and nonnegative")


def page_cusum(llrs):
    """Return W_1,...,W_n. Empty suffix is allowed (the reflection at zero).

    For h>=0, positive threshold crossings match the paper's max over nonempty
    suffixes. A raw cumulative sum alone is not a Page CUSUM detector.
    """
    llrs = np.asarray(llrs, dtype=float)
    if llrs.ndim != 1 or not np.all(np.isfinite(llrs)):
        raise ValueError("llrs must be a finite one-dimensional sequence")
    out = np.empty(llrs.size)
    w = 0.0
    for i, increment in enumerate(llrs):
        w = max(0.0, w + increment)
        out[i] = w
    return out


def stopping_time(llrs, threshold):
    """One-based first STRICT crossing; None if not observed in the sequence."""
    validate_threshold(threshold)
    crossings = np.flatnonzero(page_cusum(llrs) > threshold)
    return int(crossings[0] + 1) if crossings.size else None


def detect_with_changepoint(llrs, threshold):
    """First alarm and ML change-point estimate from the main-text argmin rule.

    Use the earliest minimum prefix if tied, including S[0]=0. A change point
    is reported only at an observed alarm; it is an estimate, not the true nu.
    """
    validate_threshold(threshold)
    values = np.asarray(llrs, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise ValueError("llrs must be a finite one-dimensional sequence")
    total = 0.0
    minimum = 0.0
    minimum_index = 0
    for t, increment in enumerate(values, 1):
        total += increment
        w = max(0.0, total - minimum)
        if w > threshold:
            return {"stop": t, "estimated_change_at": minimum_index + 1, "statistic": float(w)}
        if total < minimum:
            minimum = total
            minimum_index = t
    return {"stop": None, "estimated_change_at": None, "statistic": None}
