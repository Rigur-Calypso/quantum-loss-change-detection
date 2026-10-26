import numpy as np
import pytest
from scipy.stats import multivariate_normal
from lossqcd.models import coherent_pair,squeezed_pair
from lossqcd.entangled import (entangled_pair,entangled_kl_per_pulse,mixed_block_moments,
    mixed_collective,conditional_parameters,conditional_llrs,sample_mixed_blocks)


@pytest.mark.parametrize('n',[1,2,4,16,256])
@pytest.mark.parametrize('Na',[0.,.1,1.])
def test_covariance_energy_and_kl_against_dense_gaussian(n,Na):
    pre,post=entangled_pair(n=n,Na=Na)
    s=np.arcsinh(np.sqrt(n*Na))
    for state,eta in [(pre,.9),(post,.85)]:
        expected_diag=eta*(np.exp(-2*s)+n-1)/(4*n)+(1-eta)/4
        assert np.diag(state.covariance)==pytest.approx(np.full(n,expected_diag))
        if n>1:
            assert state.covariance[0,1]==pytest.approx(eta*(np.exp(-2*s)-1)/(4*n))
        assert np.linalg.eigvalsh(state.covariance).min()>0
    K1,K2=pre.covariance,post.covariance
    delta=post.mean-pre.mean
    dense=.5*(np.trace(np.linalg.solve(K1,K2))-n+2*np.log(np.diag(np.linalg.cholesky(K1))).sum()-2*np.log(np.diag(np.linalg.cholesky(K2))).sum()
                 +delta@np.linalg.solve(K1,delta))
    assert post.kl(pre)==pytest.approx(dense,rel=1e-10,abs=1e-11)
    assert entangled_kl_per_pulse(100,Na,.9,.85,n)==pytest.approx(dense/n,rel=1e-10,abs=1e-11)
    if Na==0:
        a,b=coherent_pair(Na=0)
        assert post.kl(pre)/n==pytest.approx(b.kl(a))


@pytest.mark.parametrize('n',[1,2,8,16])
def test_joint_and_conditional_llrs_against_independent_dense_densities(n):
    pre,post=entangled_pair(n=n,Na=.1)
    x=post.sample(np.random.default_rng(32),size=79)
    actual=post.llr(x,pre)
    dense=multivariate_normal(post.mean,post.covariance).logpdf(x)-multivariate_normal(pre.mean,pre.covariance).logpdf(x)
    np.testing.assert_allclose(actual,dense,atol=2e-11)
    np.testing.assert_allclose(post.logpdf(x),multivariate_normal(post.mean,post.covariance).logpdf(x),atol=2e-11)
    np.testing.assert_allclose(conditional_llrs(x,pre,post).sum(axis=-1),dense,atol=2e-11)
    for k in range(n):
        mean,var=conditional_parameters(pre,x[:,:k].sum(axis=-1),k)
        if k==0:
            expected_mean=np.full(len(x),pre.mean[k]);expected_var=pre.covariance[k,k]
        else:
            weights=np.linalg.solve(pre.covariance[:k,:k],pre.covariance[:k,k])
            expected_mean=pre.mean[k]+(x[:,:k]-pre.mean[:k])@weights
            expected_var=pre.covariance[k,k]-pre.covariance[k,:k]@weights
        np.testing.assert_allclose(mean,expected_mean,atol=1e-12)
        assert var==pytest.approx(expected_var,abs=1e-12)


def test_n1_is_independent_squeezing():
    a,b=entangled_pair(n=1,Na=.1)
    c,d=squeezed_pair(Na=.1)
    assert a.collective.mean==pytest.approx(c.mean)
    assert b.collective.variance==pytest.approx(d.variance)
    assert b.kl(a)==pytest.approx(d.kl(c))


def test_correlated_sampler_and_mixed_loss_generator():
    pre,post=entangled_pair(n=4,Na=.1)
    x=post.sample(np.random.default_rng(33),100000)
    np.testing.assert_allclose(x.mean(axis=0),post.mean,atol=.005)
    np.testing.assert_allclose(np.cov(x,rowvar=False),post.covariance,atol=.0025)
    llr=post.llr(x,pre)
    assert abs(llr.mean()-post.kl(pre))<4*llr.std(ddof=1)/np.sqrt(len(llr))
    etas=np.array([.9,.9,.85,.85])
    mean,K=mixed_block_moments(100,.1,etas)
    x=sample_mixed_blocks(100,.1,etas,np.random.default_rng(34),100000)
    np.testing.assert_allclose(x.mean(axis=0),mean,atol=.005)
    np.testing.assert_allclose(np.cov(x,rowvar=False),K,atol=.0025)
    u=np.ones(4)/2
    collective=mixed_collective(100,.1,4,2)
    assert collective.mean==pytest.approx(u@mean)
    assert collective.variance==pytest.approx(u@K@u)


@pytest.mark.parametrize('n',[0,-1,2.5,True])
def test_invalid_dimensions(n):
    with pytest.raises(ValueError):
        entangled_pair(n=n,Na=.1)


def test_pair_default_energy_matches_existing_models():
    a,b=entangled_pair(n=1)
    c,d=squeezed_pair()
    assert a.collective.mean==pytest.approx(c.mean)
    assert b.collective.variance==pytest.approx(d.variance)
    assert b.kl(a)==pytest.approx(d.kl(c))
