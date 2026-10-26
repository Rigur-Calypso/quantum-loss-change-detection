"""Finite-volume discretization of the Page-CUSUM first-step integral equation.

For increment CDF F and threshold h, the mean remaining stopping time obeys
 g(w) = 1 + F(-w) g(0) + integral_(0,h] g(y) dF(y-w).
The atom at zero is explicit. Positive states use midpoint values; transition
probabilities integrate the EXACT increment CDF across each cell. Solve
(I-Q)g=1, including absorption when the statistic exceeds h. Grid refinement is
an empirical discretization diagnostic, not a certified error bound.

Optional dependency: scipy (install lossqcd[numerics] or lossqcd[dev]).
"""
import numpy as np
from scipy.linalg import solve
from scipy.optimize import brentq
from scipy.special import ndtr


def increment_coefficients(pre, post, sampling=None):
    """LLR(X)=a Z^2+b Z+c, Z~N(0,1), X drawn from sampling (default pre)."""
    s = pre if sampling is None else sampling
    a = .5*s.variance*(1/pre.variance-1/post.variance)
    b = np.sqrt(s.variance)*((s.mean-pre.mean)/pre.variance-(s.mean-post.mean)/post.variance)
    c = .5*np.log(pre.variance/post.variance)+.5*((s.mean-pre.mean)**2/pre.variance
                                                  -(s.mean-post.mean)**2/post.variance)
    return float(a), float(b), float(c)


def increment_cdf(values, pre, post, sampling=None):
    """CDF of a possibly quadratic Gaussian LLR; includes unequal variances."""
    a,b,c = increment_coefficients(pre,post,sampling)
    u=np.asarray(values,dtype=float)
    if np.any(np.isnan(u)):
        raise ValueError('CDF arguments must not be NaN')
    if a==0:
        return ndtr((u-c)/abs(b)) if b!=0 else (u>=c).astype(float)
    center=-b/(2*a)
    vertex=c-b*b/(4*a)
    radius=np.sqrt(np.maximum((u-vertex)/a,0.))
    low,high=center-radius,center+radius
    if a>0:
        # sf(low)-sf(high) is more accurate than subtracting two CDFs near 1.
        between=np.where(low>0,ndtr(-low)-ndtr(-high),ndtr(high)-ndtr(low))
        out=np.where(u<vertex,0.,between)
    else:
        outside=ndtr(low)+ndtr(-high)
        out=np.where(u>=vertex,1.,outside)
    return np.clip(out,0.,1.)


def mean_stopping_time(pre, post, threshold, cells=512, sampling=None):
    """Expected stopping time starting at W=0 under a fixed sampling law.

    pre sampling gives the no-change ARL. post sampling gives E[T] for a change
    on pulse 1, whose T-1 delay is one less. Does not calculate worst-case delay.
    """
    if not np.isfinite(threshold) or threshold<0:
        raise ValueError('threshold must be finite and nonnegative')
    if not isinstance(cells,(int,np.integer)) or cells<8:
        raise ValueError('cells must be an integer >=8')
    if post.kl(pre)<=0:
        raise ValueError('distinct pre/post distributions are required')
    if threshold==0:
        probability=1-float(increment_cdf(0.,pre,post,sampling))
        return {'mean_stopping_time':float(1/probability),'threshold':0.,'cells':cells,
                'absolute_linear_residual':0.,'relative_linear_residual':0.}
    edges=np.linspace(0.,threshold,cells+1)
    states=np.r_[0.,(edges[:-1]+edges[1:])/2]
    f=increment_cdf(edges[None,:]-states[:,None],pre,post,sampling)
    q=np.column_stack([f[:,0],np.maximum(0.,np.diff(f,axis=1))])
    if np.any(q.sum(axis=1)>1+1e-12):
        raise ArithmeticError('invalid substochastic transition matrix')
    matrix=np.eye(cells+1)-q
    g=solve(matrix,np.ones(cells+1),assume_a='gen',check_finite=False)
    if not np.all(np.isfinite(g)) or np.any(g<=0):
        raise ArithmeticError('nonpositive/nonfinite stopping-time solution')
    residual=float(np.max(np.abs(matrix@g-1)))
    return {'mean_stopping_time':float(g[0]),'threshold':float(threshold),'cells':cells,
            'absolute_linear_residual':residual,
            'relative_linear_residual':residual/float(np.max(np.abs(g)))}


def calibrate_numerical(pre,post,target_arl=2e6,cells=1024,check_cells=2048,
                        mesh_tolerance=.001):
    """Calibrate numerically, then check the SAME threshold on two other grids.

    Refinement changes, and fine-grid deviation from target, must both meet
    tolerance before the result is accepted. This is not independent MC
    validation of a rare-event ARL and not a rigorous continuum error bound.
    """
    if not np.isfinite(target_arl) or target_arl<=1:
        raise ValueError('target ARL must be finite and >1')
    if cells<16 or check_cells<=cells:
        raise ValueError('need cells>=16 and check_cells>cells')
    if not np.isfinite(mesh_tolerance) or not 0<mesh_tolerance<1:
        raise ValueError('mesh tolerance must lie in (0,1)')
    history=[]
    def objective(h):
        r=mean_stopping_time(pre,post,h,cells=cells)
        history.append(r)
        return np.log(r['mean_stopping_time']/target_arl)
    low=0.
    if objective(low)>=0:
        raise ValueError('target below nonnegative-threshold range')
    high=5.
    for _ in range(20):
        if objective(high)>0:
            break
        high*=1.5
    else:
        raise RuntimeError('cannot bracket numerical ARL target')
    h=float(brentq(objective,low,high,xtol=1e-7))
    grids=[mean_stopping_time(pre,post,h,cells=n) for n in [cells//2,cells,check_cells]]
    coarse,base,fine=[r['mean_stopping_time'] for r in grids]
    change=abs(fine/base-1)
    error=abs(fine/target_arl-1)
    passed=change<=mesh_tolerance and error<=mesh_tolerance
    result={'threshold':h,'target_arl':target_arl,'cells':cells,'check_cells':check_cells,
            'mesh_tolerance':mesh_tolerance,'base_to_fine_relative_change':change,
            'fine_relative_target_error':error,'mesh_check_passed':passed,
            'grid_checks':grids,'root_history':history,
            'coarse_to_base_relative_change':abs(base/coarse-1)}
    if not passed:
        raise RuntimeError(f'mesh check failed; increase cells: {result}')
    return result
