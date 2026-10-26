"""Additional scalar observation models from arXiv:2503.12276v2.

BPSK uses raw homodyne units: equal mixtures of N(+sqrt(eta)*alpha,v)
and N(-sqrt(eta)*alpha,v). No unknown true channel is used to rescale data.
The displaced one-photon density is Eq.22, expressed in centered coordinates
as normal(x; sqrt(eta)*alpha,1/4)*(1-eta+4*eta*(x-mean)**2).
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr,xlogy,xlog1py
from .models import Gaussian


def draw(law,rng,size):
    return rng.normal(law.mean,np.sqrt(law.variance),size) if isinstance(law,Gaussian) else law.sample(rng,size)


@dataclass(frozen=True)
class BPSK:
    mean: float
    variance: float

    def __post_init__(self):
        Gaussian(self.mean,self.variance)
        if self.mean<0:raise ValueError('BPSK amplitude must be nonnegative')

    def logpdf(self,x):
        y=np.abs(np.asarray(x,dtype=float))
        return Gaussian(self.mean,self.variance).logpdf(y)-np.log(2)+np.logaddexp(0,-2*self.mean*y/self.variance)

    def llr(self,x,pre):
        if not isinstance(pre,BPSK):return self.logpdf(x)-pre.logpdf(x)
        y=np.abs(np.asarray(x,dtype=float))
        return (Gaussian(self.mean,self.variance).llr(y,Gaussian(pre.mean,pre.variance))
            +np.logaddexp(0,-2*self.mean*y/self.variance)-np.logaddexp(0,-2*pre.mean*y/pre.variance))

    def sample(self,rng,size):
        values=rng.normal(self.mean,np.sqrt(self.variance),size)
        return values*rng.choice(np.array([-1.,1.]),size=size)

    def cdf(self,x):
        x=np.asarray(x);sd=np.sqrt(self.variance)
        return .5*(ndtr((x-self.mean)/sd)+ndtr((x+self.mean)/sd))

    def kl(self,pre):
        # The score is even, so one positive-mean Gaussian component suffices.
        sd=np.sqrt(self.variance)
        return float(quad(lambda z:np.exp(-z*z/2)/np.sqrt(2*np.pi)*float(self.llr(self.mean+sd*z,pre)),
            -12,12,epsabs=1e-11,epsrel=1e-10,limit=200)[0])


def bpsk_pair(pre,post):
    return BPSK(abs(pre.mean),pre.variance),BPSK(abs(post.mean),post.variance)


def bpsk_capacity(mean,variance):
    """Uniform-input BPSK mutual information in bits/pulse (main text)."""
    Gaussian(mean,variance)
    mean=abs(mean);sd=np.sqrt(variance)
    penalty=quad(lambda z:np.exp(-z*z/2)/np.sqrt(2*np.pi)
        *np.logaddexp(0,-2*mean*(mean+sd*z)/variance)/np.log(2),-12,12,
        epsabs=1e-11,limit=200)[0]
    return float(np.clip(1-penalty,0,1))


@dataclass(frozen=True)
class DisplacedPhoton:
    alpha: float
    eta: float

    def __post_init__(self):
        if not np.isfinite(self.alpha) or not np.isfinite(self.eta) or not 0<=self.eta<=1:
            raise ValueError('finite amplitude and eta in [0,1] required')

    @property
    def mean(self):return self.alpha*np.sqrt(self.eta)

    @property
    def variance(self):return .25+.5*self.eta

    def logpdf(self,x):
        y=np.asarray(x,dtype=float)-self.mean
        with np.errstate(divide='ignore'):
            return Gaussian(self.mean,.25).logpdf(x)+np.log(1-self.eta+4*self.eta*y*y)

    def llr(self,x,pre):
        y=np.asarray(x,dtype=float)
        if isinstance(pre,DisplacedPhoton):
            with np.errstate(divide='ignore'):
                return (Gaussian(self.mean,.25).llr(y,Gaussian(pre.mean,.25))
                    +np.log(1-self.eta+4*self.eta*(y-self.mean)**2)
                    -np.log(1-pre.eta+4*pre.eta*(y-pre.mean)**2))
        return self.logpdf(y)-pre.logpdf(y)

    def cdf(self,x):
        y=np.asarray(x,dtype=float)-self.mean
        return np.clip(ndtr(2*y)-2*self.eta*y*np.exp(-2*y*y)/np.sqrt(2*np.pi),0,1)

    def sample(self,rng,size):
        # Vacuum with probability 1-eta; otherwise the Fock-1 quadrature.
        # Its signed squared standardized radius has chi-square(3) law.
        fock=rng.random(size)<self.eta
        vacuum=rng.normal(0,.5,size)
        photon=rng.choice(np.array([-1.,1.]),size=size)*np.sqrt(rng.chisquare(3,size))/2
        return self.mean+np.where(fock,photon,vacuum)

    def kl(self,pre):
        # Integrate in the post-centered coordinate, splitting both density nodes.
        points=sorted(set([-10.,0.,10.]+[float(np.clip(pre.mean-self.mean,-10,10))]))
        def integrand(y):
            logp=float(self.logpdf(self.mean+y))
            return 0. if not np.isfinite(logp) else float(np.exp(logp)*self.llr(self.mean+y,pre))
        return float(sum(quad(integrand,a,b,epsabs=1e-10,limit=300)[0] for a,b in zip(points[:-1],points[1:])))


@dataclass(frozen=True)
class Bernoulli:
    p: float  # probability of outcome 1

    def __post_init__(self):
        if not np.isfinite(self.p) or not 0<=self.p<=1:raise ValueError('p must be in [0,1]')

    def sample(self,rng,size):return rng.binomial(1,self.p,size)

    def logpdf(self,x):
        x=np.asarray(x)
        with np.errstate(divide='ignore'):
            return np.where(x==1,np.log(self.p),np.where(x==0,np.log1p(-self.p),-np.inf))

    def llr(self,x,pre):return self.logpdf(x)-pre.logpdf(x)

    def kl(self,pre):
        if self.p==pre.p:return 0.
        return float(xlogy(self.p,self.p)-xlogy(self.p,pre.p)
            +xlog1py(1-self.p,-self.p)-xlog1py(1-self.p,-pre.p))


def poisson_kl(post_mean,pre_mean):
    if not np.isfinite(post_mean) or not np.isfinite(pre_mean) or min(post_mean,pre_mean)<0:
        raise ValueError('finite nonnegative means required')
    return float(xlogy(post_mean,post_mean)-xlogy(post_mean,pre_mean)+pre_mean-post_mean)


@dataclass(frozen=True)
class Poisson:
    rate: float

    def __post_init__(self):
        if not np.isfinite(self.rate) or self.rate<0:raise ValueError('nonnegative finite rate required')

    def sample(self,rng,size):return rng.poisson(self.rate,size)

    def logpdf(self,x):
        from scipy.special import gammaln
        x=np.asarray(x)
        return np.where((x>=0)&(x==np.floor(x)),xlogy(x,self.rate)-self.rate-gammaln(x+1),-np.inf)

    def llr(self,x,pre):return xlogy(x,self.rate)-xlogy(x,pre.rate)+pre.rate-self.rate

    def kl(self,pre):return poisson_kl(self.rate,pre.rate)


def green_machine_covariance(n,Na,eta):
    """Appendix C: propagate one squeezed input through a normalized Hadamard."""
    from scipy.linalg import hadamard
    from .models import _validate
    from .entangled import validate_block_size
    validate_block_size(n);_validate(0,Na,eta,eta)
    if n&(n-1):raise ValueError('Hadamard construction requires a power-of-two block')
    seed=np.eye(n)/4;seed[0,0]=np.exp(-2*np.arcsinh(np.sqrt(n*Na)))/4
    unitary=hadamard(n)/np.sqrt(n)
    return eta*(unitary@seed@unitary.T)+(1-eta)*np.eye(n)/4


def squeezed_energy_derivative(N,Na,eta1,eta2):
    """Differentiate Eq.3 directly, rather than a high-energy approximation."""
    from .models import _validate
    _validate(N,Na,eta1,eta2)
    delta=(np.sqrt(eta2)-np.sqrt(eta1))**2
    if Na==0:
        return float('inf') if N*eta1*delta>0 else float((eta1-eta2)**2)
    r=np.arcsinh(np.sqrt(Na));e=np.exp(-2*r)
    v1=(1-eta1+eta1*e)/4;v2=(1-eta2+eta2*e)/4
    d1=-eta1*e/2;d2=-eta2*e/2;q=v2/v1;dq=(d2*v1-v2*d1)/v1**2
    return float((.5*(1-1/q)*dq-N*delta*d1/(2*v1*v1))/(2*np.sinh(r)*np.cosh(r)))


def augmentation_crossover(N=100.,eta1=.9,eta2=.85):
    """Appendix B.3: exact positive energy crossover and its asymptotic estimate."""
    from scipy.optimize import brentq
    from .models import coherent_kl,squeezed_kl,_validate
    _validate(N,0,eta1,eta2)
    if N<=0 or not 0<=eta2<eta1<1:raise ValueError('need positive N and a lossy drop below eta1<1')
    approximate=N*eta1/(1-eta1)
    objective=lambda a:squeezed_kl(N,a,eta1,eta2)-coherent_kl(N,a,eta1,eta2)
    high=max(1.,2*approximate)
    while objective(high)>0:high*=2
    root=brentq(objective,1e-8,high)
    return {'exact_Na':float(root),'approximate_Na':float(approximate)}


def kennedy_residual(N,Na,eta1,eta2,residual_photons,mode='pnr'):
    """Appendix A Eq.7, absolute residual-amplitude convention epsilon=sqrt(Nerr)."""
    from .models import _validate
    _validate(N,Na,eta1,eta2)
    if not np.isfinite(residual_photons) or residual_photons<0:raise ValueError('nonnegative residual photons required')
    shift=np.sqrt(N+Na)*(np.sqrt(eta2)-np.sqrt(eta1))
    post_rate=(shift-np.sqrt(residual_photons))**2
    if mode=='pnr':return poisson_kl(post_rate,residual_photons)
    if mode=='onoff':return Bernoulli(-np.expm1(-post_rate)).kl(Bernoulli(-np.expm1(-residual_photons)))
    raise ValueError('mode must be pnr or onoff')
