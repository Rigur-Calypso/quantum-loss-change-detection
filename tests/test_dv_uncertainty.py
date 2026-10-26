"""The grouped bootstrap must equal the ordinary empirical resampling law."""
import itertools
import numpy as np
import pytest
from scripts.compare_dv_uncertainty import bootstrap_mean


@pytest.mark.parametrize('data',[[-2,1,1,4],[-844,-129,17,34,50]])
def test_bootstrap_moments_against_exhaustive_individual_resampling(data):
    data=np.array(data)
    exact=np.array([np.mean(x) for x in itertools.product(data,repeat=len(data))])
    sample=bootstrap_mean(data,200000,20261025)
    mean_se=np.sqrt(exact.var()/len(sample))
    assert abs(sample.mean()-exact.mean())<4*mean_se
    assert sample.var()==pytest.approx(exact.var(),rel=.02)
    assert sample.min()>=exact.min() and sample.max()<=exact.max()
    # Early/negative observations must remain eligible in the signed bootstrap.
    expected_negative=np.mean(exact<0)
    se=np.sqrt(expected_negative*(1-expected_negative)/len(sample))
    assert abs(np.mean(sample<0)-expected_negative)<4*se


def test_bootstrap_constant_data_has_zero_uncertainty():
    assert np.array_equal(bootstrap_mean(np.repeat(7,500),20000,20261025),np.repeat(7.,20000))
