"""Appendix A receiver imperfections, Eqs.7,10,11.

Phase is a normalized Gaussian truncated to [-pi,pi]. Quadrature drops at
most the |Z|>12 tails when the truncation is wider (mass <4e-33), and normalizes
the retained weights. Kennedy uses photon-number resolution here; its binary
on/off idealization in the main text is a different observation channel.
"""
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.special import logsumexp,gammaln,xlogy
from scipy.integrate import quad
from scipy.stats import truncnorm
from .models import _validate


def phase_quadrature(sigma,order=96):
    if not np.isfinite(sigma) or sigma<0 or isinstance(order,bool) or not isinstance(order,int) or order<16:
        raise ValueError('nonnegative phase SD and integer order >=16 required')
    if sigma==0:return np.array([0.]),np.array([1.])
    bound=min(np.pi/sigma,12.)
    z,w=np.polynomial.legendre.leggauss(order);z=z*bound
    weights=w*bound*np.exp(-z*z/2)/np.sqrt(2*np.pi);weights/=weights.sum()
    return sigma*z,weights


@dataclass(frozen=True)
class PhaseHomodyne:
    N: float
    Na: float
    eta: float
    phase_sd: float=0.
    thermal: float=0.
    order: int=96

    def __post_init__(self):
        _validate(self.N,self.Na,self.eta,self.eta)
        if not np.isfinite(self.thermal) or self.thermal<0:raise ValueError('thermal occupation must be nonnegative')
        phase_quadrature(self.phase_sd,self.order)

    @cached_property
    def components(self):
        theta,w=phase_quadrature(self.phase_sd,self.order);r=np.arcsinh(np.sqrt(self.Na))
        low=((1+2*self.thermal)*(1-self.eta)+self.eta*np.exp(-2*r))/4
        high=((1+2*self.thermal)*(1-self.eta)+self.eta*np.exp(2*r))/4
        return np.sqrt(self.eta*self.N)*np.cos(theta),low*np.cos(theta)**2+high*np.sin(theta)**2,w

    def logpdf(self,x):
        mean,var,w=self.components;y=np.asarray(x)[...,None]
        return logsumexp(-.5*(np.log(2*np.pi*var)+(y-mean)**2/var)+np.log(w),axis=-1)

    def llr(self,x,pre):return self.logpdf(x)-pre.logpdf(x)

    def sample(self,rng,size):
        theta=np.zeros(size) if self.phase_sd==0 else truncnorm.rvs(-np.pi/self.phase_sd,np.pi/self.phase_sd,
            scale=self.phase_sd,size=size,random_state=rng)
        r=np.arcsinh(np.sqrt(self.Na));noise=(1+2*self.thermal)*(1-self.eta)
        var=(noise+self.eta*(np.exp(-2*r)*np.cos(theta)**2+np.exp(2*r)*np.sin(theta)**2))/4
        return rng.normal(np.sqrt(self.eta*self.N)*np.cos(theta),np.sqrt(var),size)

    def kl(self,pre):
        mean,var,_=self.components;sd=np.sqrt(var.max());lo=float(mean.min()-12*sd);hi=float(mean.max()+12*sd)
        return float(quad(lambda x:float(np.exp(self.logpdf(x))*self.llr(x,pre)),lo,hi,
            epsabs=1e-9,epsrel=1e-8,limit=300)[0])


@dataclass(frozen=True)
class PhaseKennedy:
    photons: float
    eta: float
    reference_eta: float
    amplitude_error: float=0.
    phase_sd: float=0.
    order: int=96

    def __post_init__(self):
        _validate(self.photons,0,self.eta,self.reference_eta)
        if not np.isfinite(self.amplitude_error) or self.amplitude_error<=-1:raise ValueError('fractional amplitude error must exceed -1')
        phase_quadrature(self.phase_sd,self.order)

    @cached_property
    def components(self):
        theta,w=phase_quadrature(self.phase_sd,self.order)
        ref=np.sqrt(self.reference_eta)*(1+self.amplitude_error);signal=np.sqrt(self.eta)
        # Avoid subtracting almost identical cosines near a perfect null.
        rates=self.photons*((signal-ref)**2+4*signal*ref*np.sin(theta/2)**2)
        return rates,w

    def logpmf(self,k):
        rates,w=self.components;k=np.asarray(k)[...,None]
        return logsumexp(xlogy(k,rates)-rates-gammaln(k+1)+np.log(w),axis=-1)

    def logpdf(self,k):return self.logpmf(k)

    def llr(self,k,pre):return self.logpmf(k)-pre.logpmf(k)

    def sample(self,rng,size):
        theta=np.zeros(size) if self.phase_sd==0 else truncnorm.rvs(-np.pi/self.phase_sd,np.pi/self.phase_sd,
            scale=self.phase_sd,size=size,random_state=rng)
        ref=np.sqrt(self.reference_eta)*(1+self.amplitude_error);signal=np.sqrt(self.eta)
        rate=self.photons*((signal-ref)**2+4*signal*ref*np.sin(theta/2)**2)
        return rng.poisson(rate)

    def kl_details(self,pre):
        rates,_=self.components;other,_=pre.components
        maximum=float(max(rates.max(),other.max()))
        upper=int(np.ceil(maximum+14*np.sqrt(maximum+1)+50));k=np.arange(upper+1)
        lp,lq=self.logpmf(k),pre.logpmf(k);keep=np.isfinite(lp)
        value=float(np.sum(np.exp(lp[keep])*(lp[keep]-lq[keep])))
        return {'kl':value,'count_cutoff':upper,
            'post_omitted_mass':float(max(0.,1-np.exp(lp).sum())),
            'pre_omitted_mass':float(max(0.,1-np.exp(lq).sum()))}

    def kl(self,pre):return self.kl_details(pre)['kl']

    def click_probability(self):
        rates,w=self.components
        return float(np.sum(w*(-np.expm1(-rates))))

    def onoff_kl(self,pre):
        from .paper_models import Bernoulli
        return Bernoulli(self.click_probability()).kl(Bernoulli(pre.click_probability()))
