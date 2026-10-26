from dataclasses import replace
import numpy as np
import pytest
from scipy.special import expit
from scipy.stats import norm
from lossqcd.conditional_detection import detector_models,simulate_conditional
from lossqcd.conditional_rare_event import (adaptive_log_weight,sample_conditional_cycles,
    estimate_conditional_arl,summarize_conditional_cycles,conditional_cycle_arrays)
from lossqcd.batched_rare_event import pool_cycle_summaries
from lossqcd.simulation import summarize_arl


def dense_conditional(block,previous):
    k=len(previous);cov=block.covariance
    if k==0:
        return block.pulse_mean,cov[0,0]
    coefficient=np.linalg.solve(cov[:k,:k],cov[:k,k])
    mean=block.pulse_mean+coefficient@(np.asarray(previous)-block.pulse_mean)
    var=cov[k,k]-cov[k,:k]@coefficient
    return mean,var


@pytest.mark.parametrize('page',[0.,2.,12.])
@pytest.mark.parametrize('seen',[0,3,7])
def test_step_weight_matches_independent_dense_gaussian_mixture(page,seen):
    (pre,post),_,_=detector_models('entangled',8,Na=.1)
    previous=np.linspace(9.2,9.8,seen)
    m1,v1=dense_conditional(pre,previous);m2,v2=dense_conditional(post,previous)
    x=np.linspace(7,11,31);b=2.;beta=expit(page-b)
    l1=norm.logpdf(x,m1,np.sqrt(v1));l2=norm.logpdf(x,m2,np.sqrt(v2))
    mixture=np.logaddexp(np.log1p(-beta)+l1,np.log(beta)+l2)
    np.testing.assert_allclose(adaptive_log_weight(page,l2-l1,b),l1-mixture,atol=2e-13)


def test_partial_block_weight_uses_mixture_density_before_page_reflection(monkeypatch):
    import lossqcd.conditional_rare_event as module
    path=[10.,8.]
    class Rng:
        def __init__(self):
            self.values=iter(path)
        def random(self,size):
            return np.full(size,.5)
        def normal(self,mean,sd,size):
            return np.full(size,next(self.values))
    monkeypatch.setattr(module.np.random,'default_rng',lambda seed:Rng())
    batch=sample_conditional_cycles(n=4,threshold=1,cycles=2,proposal='adaptive',max_steps=4)
    assert np.all(batch.upper) and np.all(batch.length==2)
    (pre,post),_,_=detector_models('entangled',4,Na=.1)
    w=0.;log_weight=0.;raw=0.
    for k,x in enumerate(path):
        m1,v1=dense_conditional(pre,path[:k]);m2,v2=dense_conditional(post,path[:k])
        l1=norm.logpdf(x,m1,np.sqrt(v1));l2=norm.logpdf(x,m2,np.sqrt(v2))
        beta=expit(w-2)
        log_weight+=l1-np.logaddexp(np.log1p(-beta)+l1,np.log(beta)+l2)
        ell=l2-l1;raw+=ell;w=max(0,w+ell)
    np.testing.assert_allclose(batch.log_weight,log_weight,atol=1e-13)
    assert np.all(abs(batch.log_weight+batch.raw_sum)>.1)
    assert 'post_log_weight' in conditional_cycle_arrays(batch,batch)


@pytest.mark.parametrize('Na',[.1,5.])
def test_adaptive_stopped_path_mass_and_weighted_cycle_length(Na):
    stats,pre,proposal=estimate_conditional_arl(n=8,Na=Na,threshold=4.5,
        cycles=150000,seed=910,proposal='adaptive')
    weight=np.exp(proposal.log_weight)
    se=weight.std(ddof=1)/np.sqrt(len(weight))
    assert abs(weight.mean()-1)<5*se
    length=weight*proposal.length
    length_se=length.std(ddof=1)/np.sqrt(len(length))
    assert abs(length.mean()-pre.length.mean())<5*np.hypot(length_se,stats['mean_pre_cycle_length_se'])
    direct=summarize_arl(simulate_conditional(n=8,Na=Na,threshold=4.5,trials=20000,horizon=20000,seed=911))
    assert abs(stats['arl']-direct['arl'])<5*np.hypot(stats['arl_se'],direct['arl_se'])


def test_adaptive_pooling_equals_concatenating_full_paths_and_rejects_mixed_methods():
    summaries=[];pres=[];posts=[]
    for seed in [912,913]:
        stats,pre,post=estimate_conditional_arl(n=8,threshold=5,cycles=5000,seed=seed,proposal='adaptive')
        summaries.append(stats);pres.append(pre);posts.append(post)
    def combine(batches):
        fields=['length','raw_sum','page','upper','complete']
        values={field:np.concatenate([getattr(b,field) for b in batches]) for field in fields}
        if batches[0].log_weight is not None:
            values['log_weight']=np.concatenate([b.log_weight for b in batches])
        return replace(batches[0],**values)
    expected=summarize_conditional_cycles(combine(pres),combine(posts))
    actual=pool_cycle_summaries(summaries)
    for key in ['arl','arl_se','success_weight_ess','largest_weight_fraction']:
        assert actual[key]==pytest.approx(expected[key],rel=1e-12)
    with pytest.raises(ValueError,match='same detector'):
        pool_cycle_summaries([summaries[0],summaries[1]|{'proposal_offset':4.}])


def test_missing_path_weight_is_refused():
    _,pre,post=estimate_conditional_arl(n=8,threshold=3,cycles=100,seed=914,proposal='adaptive')
    with pytest.raises(ValueError,match='path log weights'):
        summarize_conditional_cycles(pre,replace(post,log_weight=None))


@pytest.mark.parametrize('kwargs',[{'proposal_offset':-1},{'proposal_offset':np.nan},
    {'proposal_offset':11},{'proposal':'invalid'}])
def test_invalid_proposal_inputs(kwargs):
    with pytest.raises(ValueError):
        sample_conditional_cycles(proposal='adaptive',**kwargs) if 'proposal' not in kwargs else sample_conditional_cycles(**kwargs)
