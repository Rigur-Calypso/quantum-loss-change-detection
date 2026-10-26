import numpy as np
import pytest
from lossqcd.models import Gaussian,coherent_pair,squeezed_pair
from lossqcd.mismatch_renewal import sample_mismatch_cycles,estimate_mismatch_arl,summarize_mismatch_cycles
from lossqcd.arl_solver import mean_stopping_time
from lossqcd.rare_event import sample_cycles


def test_matched_case_preserves_original_cycles_and_likelihood_identity():
    pre,post=squeezed_pair(Na=.2)
    a=sample_cycles(pre,post,3,cycles=5000,proposal='post',seed=721)
    b=sample_mismatch_cycles(pre,post,pre,post,3,cycles=5000,seed=721)
    np.testing.assert_array_equal(a.length,b.length)
    np.testing.assert_array_equal(a.terminal_sum,b.score)
    np.testing.assert_array_equal(a.upper,b.upper)
    np.testing.assert_allclose(b.log_weight,-b.score,atol=1e-12)


def test_mismatch_weight_is_actual_density_ratio_not_assumed_score(monkeypatch):
    import lossqcd.mismatch_renewal as module
    pre,post,actual=Gaussian(0,1),Gaussian(1,1),Gaussian(.2,.7)
    values=iter((.8,1.2))
    class Scripted:
        def normal(self,*args):
            return np.full(args[-1],next(values))
    monkeypatch.setattr(module.np.random,'default_rng',lambda seed:Scripted())
    batch=sample_mismatch_cycles(pre,post,actual,post,.5,cycles=2)
    assert np.all(batch.length==2)
    expected=sum(actual.logpdf(x)-post.logpdf(x) for x in (.8,1.2))
    np.testing.assert_allclose(batch.log_weight,expected)
    assert not np.allclose(batch.log_weight,-batch.score)


@pytest.mark.parametrize('pair',[coherent_pair,squeezed_pair])
def test_mismatched_renewal_matches_independent_first_step_solver(pair):
    pre,post=pair(Na=.1)
    actual,_=pair(Na=.1,eta1=.895)
    stats,_,_=estimate_mismatch_arl(pre,post,actual,3,cycles=100000,seed=722)
    expected=mean_stopping_time(pre,post,3,cells=512,sampling=actual)['mean_stopping_time']
    assert abs(expected-stats['arl'])<4*stats['arl_se']


def test_unfinished_cycles_are_rejected():
    pre,post=coherent_pair()
    batch=sample_mismatch_cycles(pre,post,pre,post,100,cycles=100,max_steps=2)
    with pytest.raises(RuntimeError,match='unfinished'):
        summarize_mismatch_cycles(batch,batch)


@pytest.mark.parametrize('cycles',[True,1,2.5])
def test_invalid_counts(cycles):
    pre,post=coherent_pair()
    with pytest.raises(ValueError):
        sample_mismatch_cycles(pre,post,pre,post,3,cycles=cycles)
