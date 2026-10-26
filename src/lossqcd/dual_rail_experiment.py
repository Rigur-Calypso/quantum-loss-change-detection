"""Ideal Appendix E.2 state/erasure model with explicit physical time bins."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .dual_rail import pure_state,heralded_state,simulate_heralds,summarize_dual_rail
from .paper_models import Bernoulli
from .iid_detection import estimate_iid_arl
from .mismatch_renewal import cycle_arrays


def state_checks(seed,count=1500):
    rng=np.random.default_rng(seed);population=rng.random(count);phase=rng.uniform(0,2*np.pi,count)
    vectors=np.column_stack([np.sqrt(population),np.sqrt(1-population)*np.exp(1j*phase)])
    probability=[];uniform_fidelity=[];mixed_fidelity=[]
    for v in vectors:
        rho=pure_state(v);p,uniform=heralded_state(rho,.9);_,mixed=heralded_state(rho,.9,.81)
        probability.append(p);uniform_fidelity.append(float(np.vdot(v,uniform@v).real))
        mixed_fidelity.append(float(np.vdot(v,mixed@v).real))
    return {'vectors':vectors,'population':population,'uniform_probability':np.array(probability),
        'uniform_fidelity':np.array(uniform_fidelity),'mixed_fidelity':np.array(mixed_fidelity)}


def run(output,cycles=250000,trials=30000,seed=2026102503):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    seeds=[int(c.generate_state(1,dtype=np.uint64)[0]) for c in np.random.SeedSequence(seed).spawn(19)]
    pre,post=Bernoulli(.9),Bernoulli(.81);target_pairs=1e6;h=8.;history=[]
    report={'source_url':'https://arxiv.org/html/2503.12276v2#A5.SS2',
        'scope':'Ideal logical dual-rail erasure/heralding model; no hardware dynamics, QKD security or optimal-tradeoff claim.',
        'settings':{'eta1':.9,'eta2':.81,'target_arl_pairs':target_pairs,'target_arl_physical_bins':2*target_pairs,
            'cycles':cycles,'trials':trials,'seed':seed,'horizon_pairs':10000,'state_count':1500},
        'calibration':history,'validation':[],'changed':[]}
    for i in range(6):
        stats,_,_=estimate_iid_arl(pre,post,h,cycles=cycles,seed=seeds[i])
        history.append({'threshold':h,'summary':stats});h=max(0.,h+np.log(target_pairs/stats['arl']))
    report['threshold']=float(h)
    for i in range(2):
        stats,a,b=estimate_iid_arl(pre,post,h,cycles=cycles,seed=seeds[6+i])
        np.savez_compressed(output/f'validation_{i}.npz',**cycle_arrays(a,b));report['validation'].append(stats)
    a,b=report['validation'];z=float(abs(a['arl']-b['arl'])/np.hypot(a['arl_se'],b['arl_se']))
    report['validation_repeat_z']=z
    report['arl_confirmed']=bool(z<=4 and all(x['diagnostics_passed'] and .95e6<=x['approx_95_interval'][0]
        and x['approx_95_interval'][1]<=1.05e6 for x in report['validation']))
    control,num,den=estimate_iid_arl(pre,post,2,cycles=cycles,seed=seeds[8])
    direct=simulate_heralds(.9,.81,2,trials=8000,horizon_pairs=200000,seed=seeds[9])
    direct_stats=summarize_dual_rail(direct)
    np.savez_compressed(output/'control.npz',stop=direct.pairs.stop,event=direct.pairs.event,**cycle_arrays(num,den))
    scalar=direct_stats['pair_statistics']
    z=None if scalar['arl'] is None else float(abs(scalar['arl']-control['arl'])/np.hypot(scalar['arl_se'],control['arl_se']))
    report['control']={'weighted':control,'direct':direct_stats,'seed':seeds[9],'z':z,
        'passed':bool(z is not None and z<=4 and control['diagnostics_passed'])}
    for cp in (2001,2002):
        for encoding in ('early','late','plus','haar'):
            i=len(report['changed']);run=simulate_heralds(.9,.81,h,trials=trials,horizon_pairs=10000,
                change_pulse=cp,encoding=encoding,seed=seeds[10+i])
            filename=f'changed_{i}.npz';np.savez_compressed(output/filename,stop=run.pairs.stop,event=run.pairs.event)
            summary=summarize_dual_rail(run)
            report['changed'].append({'encoding':encoding,'change_pulse':cp,'seed':seeds[10+i],
                'file':filename,'summary':summary,'passed':summary['pair_statistics']['post_change_censored']==0})
            print(f'{encoding}, physical change bin {cp}: {summary["physical_time_bin_statistics"]["conditional_delay"]:.3f} bins',flush=True)
    states=state_checks(seeds[18]);np.savez_compressed(output/'states.npz',**states)
    report['states']={'seed':seeds[18],'max_survival_error':float(np.max(abs(states['uniform_probability']-.9))),
        'max_uniform_infidelity':float(np.max(abs(1-states['uniform_fidelity']))),
        'max_mid_pair_infidelity':float(np.max(1-states['mixed_fidelity']))}
    report['passed']=bool(report['arl_confirmed'] and report['control']['passed'] and all(c['passed'] for c in report['changed'])
        and report['states']['max_survival_error']<1e-12 and report['states']['max_uniform_infidelity']<1e-12)
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    rows=[]
    for c in report['changed']:
        pair=c['summary']['pair_statistics'];physical=c['summary']['physical_time_bin_statistics']
        rows.append({'encoding':c['encoding'],'change_physical_bin':c['change_pulse'],
            'false_alarms':pair['false_alarms'],'survivors':pair['survivors'],
            'conditional_delay_pairs':pair['conditional_delay'],'conditional_delay_physical_bins':physical['conditional_delay'],
            'delay_physical_se':physical['conditional_delay_se']})
    with (output/'delays.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    plot(output,report,states)
    return report


def plot(output,report,states):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    order=np.argsort(states['population']);x=states['population'][order]
    axes[0].plot(x,1-states['mixed_fidelity'][order],label='Change between rails: eta=0.90 / 0.81')
    axes[0].axhline(0,color='black',ls=':',label='Equal rail loss: state preserved')
    axes[0].set(xlabel='Input probability of photon in early rail',ylabel='Conditional state infidelity 1−F',
        title='A mid-qubit change is distinct from uniform loss')
    for cp,offset in ((2001,-.12),(2002,.12)):
        rows=[c for c in report['changed'] if c['change_pulse']==cp]
        axes[1].errorbar(np.arange(4)+offset,[c['summary']['physical_time_bin_statistics']['conditional_delay'] for c in rows],
            yerr=[1.96*c['summary']['physical_time_bin_statistics']['conditional_delay_se'] for c in rows],
            fmt='o',capsize=3,label=f'Change on bin {cp}'+(' (early)' if cp%2 else ' (late)'))
    axes[1].set_xticks(range(4),['early','late','plus','Haar'])
    axes[1].set(xlabel='Qubit encoding',ylabel='Conditional first-alarm delay (physical bins)',
        title='Readout after two bins; target null ARL 2 million bins')
    for ax in axes:ax.legend(fontsize=7);ax.grid(alpha=.2)
    for ext in ('png','svg'):fig.savefig(output/f'dual_rail.{ext}',dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='results_dual_rail');a=p.parse_args()
    r=run(a.output);print('All checks passed:',r['passed']);raise SystemExit(0 if r['passed'] else 1)


if __name__=='__main__':main()
