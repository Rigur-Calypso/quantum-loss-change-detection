"""Sample-based Page CUSUM and bounded-weight renewal ARL for general iid laws.

Applies to BPSK mixtures and Eq.22 single-photon homodyne distributions without
approximating them as Gaussians. Renewals occur at each return to Page zero;
this is for iid pulses, not correlated blocks. With post-proposal sampling,
successful likelihood weights exp(-raw_score) are bounded by exp(-threshold).
"""
import numpy as np
from .paper_models import draw
from .cusum import validate_threshold
from .simulation import Runs
from .mismatch_renewal import MismatchCycles,summarize_mismatch_cycles


def simulate_iid(pre,post,threshold,trials=30000,horizon=5000,change_at=None,seed=2026103101):
    validate_threshold(threshold)
    if any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) or v<1 for v in (trials,horizon)):
        raise ValueError('positive integer counts required')
    if change_at is not None and (isinstance(change_at,bool) or not isinstance(change_at,(int,np.integer)) or not 1<=change_at<=horizon):
        raise ValueError('change point must lie within the horizon')
    rng=np.random.default_rng(seed);active=np.arange(trials);page=np.zeros(trials)
    stop=np.full(trials,horizon,dtype=np.int64);event=np.zeros(trials,dtype=bool)
    for t in range(1,horizon+1):
        law=post if change_at is not None and t>=change_at else pre
        x=draw(law,rng,len(active));page[active]=np.maximum(0,page[active]+post.llr(x,pre))
        hit=page[active]>threshold
        stop[active[hit]]=t;event[active[hit]]=True;active=active[~hit]
        if not len(active):break
    return Runs(stop,event,horizon,change_at)


def sample_iid_cycles(pre,post,threshold,cycles=200000,proposal='pre',seed=2026103101,max_steps=100000):
    validate_threshold(threshold)
    if proposal not in ('pre','post'):raise ValueError('proposal must be pre or post')
    if any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) or v<2 for v in (cycles,max_steps)):
        raise ValueError('integer cycle counts >=2 required')
    rng=np.random.default_rng(seed);active=np.arange(cycles);score=np.zeros(cycles)
    length=np.full(cycles,max_steps,dtype=np.int64);complete=np.zeros(cycles,dtype=bool);upper=complete.copy()
    law=pre if proposal=='pre' else post
    for t in range(1,max_steps+1):
        x=draw(law,rng,len(active));score[active]+=post.llr(x,pre)
        if not np.all(np.isfinite(score[active])):raise ValueError('renewal estimation requires finite likelihood ratios')
        hit=score[active]>threshold;done=hit|(score[active]<=0)
        upper[active[hit]]=True;complete[active[done]]=True;length[active[done]]=t;active=active[~done]
        if not len(active):break
    return MismatchCycles(length,score,-score if proposal=='post' else np.zeros(cycles),upper,complete)


def estimate_iid_arl(pre,post,threshold,cycles=200000,seed=2026103101):
    seeds=[int(s.generate_state(1,dtype=np.uint64)[0]) for s in np.random.SeedSequence(seed).spawn(2)]
    a=sample_iid_cycles(pre,post,threshold,cycles,'pre',seeds[0])
    b=sample_iid_cycles(pre,post,threshold,cycles,'post',seeds[1])
    stats=summarize_mismatch_cycles(a,b)
    stats.update({'seed':seed,'stream_seeds':seeds,'success_weight_bound':float(np.exp(-threshold))})
    if np.any(b.log_weight[b.upper]>=-threshold):raise ArithmeticError('iid successful weight must be below exp(-h)')
    return stats,a,b


def simulate_threshold_grid(pre,post,thresholds,trials=1000,horizon=100000,seed=2026103101):
    """Appendix F: all first crossings on each shared, unreset null path.

    Every threshold crossed by a single overshooting pulse receives that SAME
    pulse index. Lower alarms do not reset the path. Output ensembles across
    thresholds are correlated; do not use independent-difference error bars.
    """
    thresholds=np.asarray(thresholds,dtype=float)
    if thresholds.ndim!=1 or not len(thresholds) or np.any(~np.isfinite(thresholds)) or np.any(thresholds<0) or np.any(np.diff(thresholds)<=0):
        raise ValueError('finite nonnegative strictly increasing thresholds required')
    if any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) or v<1 for v in (trials,horizon)):
        raise ValueError('positive integer counts required')
    rng=np.random.default_rng(seed);active=np.arange(trials);page=np.zeros(trials);crossed_count=np.zeros(trials,dtype=int)
    stop=np.full((trials,len(thresholds)),horizon,dtype=np.int64);event=np.zeros_like(stop,dtype=bool)
    for t in range(1,horizon+1):
        x=draw(pre,rng,len(active));page[active]=np.maximum(0,page[active]+post.llr(x,pre))
        following=np.maximum(crossed_count[active],np.searchsorted(thresholds,page[active],side='left'))
        for j in range(len(thresholds)):
            hit=active[(crossed_count[active]<=j)&(following>j)]
            stop[hit,j]=t;event[hit,j]=True
        crossed_count[active]=following;active=active[following<len(thresholds)]
        if not len(active):break
    return [Runs(stop[:,j],event[:,j],horizon,None) for j in range(len(thresholds))]
