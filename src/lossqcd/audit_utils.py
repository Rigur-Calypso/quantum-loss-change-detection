"""Shared source-file hashes and saved-cycle comparison tolerances."""

import hashlib
from pathlib import Path


def sha256(path):
    """Hash a file in bounded-memory chunks; accept a path object or string."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare_cycle_summaries(actual, expected):
    """Assert agreement of saved ARL, SE and importance-weight diagnostics.

    These tolerances allow floating-point rounding in pooled moments, not
    Monte Carlo disagreement between independently generated ensembles.
    """
    import numpy as np

    for key in ("arl", "arl_se", "success_weight_ess", "largest_weight_fraction"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-11, atol=1e-14, err_msg=key)
