"""Strict Page crossings and earliest-prefix change-point estimates."""

import numpy as np
import pytest
from lossqcd.cusum import page_cusum, stopping_time


def test_reflection_and_strict_one_based_stopping():
    llr = [-3, 1, 1, -5, 2, 1]
    np.testing.assert_array_equal(page_cusum(llr), [0, 1, 2, 0, 2, 3])
    assert stopping_time(llr, 2) == 6  # equality at samples 3 and 5 is not a crossing
    assert stopping_time(llr, 3) is None
    assert stopping_time([0.1], 0) == 1
    assert stopping_time([], 2) is None


def test_recursion_matches_exhaustive_suffix_search():
    rng = np.random.default_rng(44)
    llr = rng.normal(-0.1, 1, 100)
    brute = [max([0] + [llr[k:t].sum() for k in range(t)]) for t in range(1, 101)]
    np.testing.assert_allclose(page_cusum(llr), brute, atol=1e-13)


@pytest.mark.parametrize("h", [-1, np.inf, np.nan])
def test_invalid_threshold(h):
    with pytest.raises(ValueError):
        stopping_time([1], h)
