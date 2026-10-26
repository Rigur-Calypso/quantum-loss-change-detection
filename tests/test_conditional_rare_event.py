from dataclasses import replace
import numpy as np
import pytest
from lossqcd.conditional_rare_event import (
    sample_conditional_cycles,estimate_conditional_arl,summarize_conditional_cycles)
from lossqcd.rare_event import sample_cycles
from lossqcd.models import squeezed_pair
from lossqcd.entangled import entangled_pair,conditional_llrs
from lossqcd.arl_solver import mean_stopping_time


@pytest.mark.parametrize('proposal',['pre','post'])
def test_n1_reduces_to_existing_scalar_renewal(proposal):
    pre,post=squeezed_pair(Na=.1)
    expected=sample_cycles(pre,post,5,cycles=1000,proposal=proposal,seed=91)
    actual=sample_conditional_cycles(n=1,threshold=5,cycles=1000,proposal=proposal,seed=91)
    np.testing.assert_array_equal(actual.length,expected.length)
    np.testing.assert_array_equal(actual.upper,expected.upper)
    np.testing.assert_allclose(actual.raw_sum,expected.terminal_sum,atol=1e-12)


def scripted(monkeypatch,values):
    import lossqcd.conditional_rare_event as module
    class Rng:
        def __init__(self):
            self.values=iter(values)
        def normal(self,mean,sd,size):
            return np.full(size,next(self.values))
    monkeypatch.setattr(module.np.random,'default_rng',lambda seed:Rng())


def test_zero_inside_block_does_not_end_cycle(monkeypatch):
    # Strong positive residuals give negative evidence at every pulse.
    scripted(monkeypatch,[10.]*4)
    batch=sample_conditional_cycles(n=4,threshold=5,cycles=2,max_steps=4)
    assert np.all(batch.length==4)
    assert np.all(batch.complete & ~batch.upper)
    assert np.all(batch.page==0)
    assert np.all(batch.raw_sum<0)


def test_partial_block_alarm_preserves_prefix_and_raw_lr(monkeypatch):
    # First observation reflects to zero, then a low observation alarms.
    x=np.array([10.,8.,9.,9.])
    pre,post=entangled_pair(n=4,Na=.1)
    ell=conditional_llrs(x[None,:],pre,post)[0]
    assert ell[0]<0 and ell[1]>1
    scripted(monkeypatch,x)
    batch=sample_conditional_cycles(n=4,threshold=1,cycles=2,max_steps=4,proposal='post')
    assert np.all(batch.length==2)
    assert np.all(batch.upper & batch.complete)
    np.testing.assert_allclose(batch.raw_sum,ell[:2].sum(),atol=1e-12)
    np.testing.assert_allclose(batch.page,ell[1],atol=1e-12)
    # This distinction catches the tempting but incorrect exp(-Page) weight.
    assert np.all(batch.raw_sum<batch.page)
    prefix=x[:2];k=pre.covariance[:2,:2];l=post.covariance[:2,:2]
    def dense(mu,cov):
        delta=prefix-mu
        return -.5*(np.linalg.slogdet(cov)[1]+delta@np.linalg.solve(cov,delta))
    expected=dense(post.pulse_mean,l)-dense(pre.pulse_mean,k)
    np.testing.assert_allclose(batch.raw_sum,expected,atol=1e-12)


def test_unfinished_cycles_are_refused():
    s,pre,post=estimate_conditional_arl(n=2,threshold=3,cycles=100,seed=97)
    complete=post.complete.copy();complete[0]=False
    with pytest.raises(RuntimeError,match='unfinished'):
        summarize_conditional_cycles(pre,replace(post,complete=complete))
    with pytest.raises(ValueError,match='same detector'):
        summarize_conditional_cycles(pre,replace(post,n=4))


def test_independent_streams_and_scalar_integral_crosscheck():
    s,pre,post=estimate_conditional_arl(n=1,threshold=5,cycles=40000,seed=98)
    a,b=squeezed_pair(Na=.1)
    exact=mean_stopping_time(a,b,5,cells=512)['mean_stopping_time']
    assert abs(s['arl']-exact)<5*s['arl_se']+.002*exact
    assert s['stream_seeds'][0]!=s['stream_seeds'][1]
    assert np.all(pre.length[~pre.upper]%pre.n==0)


@pytest.mark.parametrize('kwargs',[{'cycles':1},{'cycles':True},{'max_steps':1},
    {'n':0},{'threshold':-1},{'proposal':'bad'},{'eta1':.9,'eta2':.9}])
def test_invalid_inputs(kwargs):
    with pytest.raises(ValueError):
        sample_conditional_cycles(**kwargs)


def test_denominator_uses_all_cycles_and_raw_path_weights():
    from lossqcd.conditional_rare_event import ConditionalCycles
    pre=ConditionalCycles(np.array([4,8]),np.array([-1.,-2.]),np.zeros(2),
        np.zeros(2,dtype=bool),np.ones(2,dtype=bool),'entangled',4,1.,'pre',(100.,.1,.9,.85))
    post=replace(pre,proposal='post',length=np.array([2,4]),raw_sum=np.array([.5,-1.]),
        page=np.array([2.,0.]),upper=np.array([True,False]))
    stats=summarize_conditional_cycles(pre,post)
    assert stats['cycle_alarm_probability']==pytest.approx(np.exp(-.5)/2)
    assert stats['arl']==pytest.approx(6/(np.exp(-.5)/2))
    assert not stats['diagnostics_passed']  # Two cycles cannot establish precision.


@pytest.mark.parametrize('n',[2,8])
def test_correlated_renewal_matches_direct_first_alarm(n):
    from lossqcd.conditional_detection import simulate_conditional
    from lossqcd.simulation import summarize_arl
    direct=summarize_arl(simulate_conditional(n=n,threshold=2,trials=20000,horizon=5000,seed=301))
    estimated,_,_=estimate_conditional_arl(n=n,threshold=2,cycles=200000,seed=302)
    assert abs(estimated['arl']-direct['arl'])<5*np.hypot(estimated['arl_se'],direct['arl_se'])
