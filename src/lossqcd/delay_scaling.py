"""Finite-ARL check of Eq.1 for iid coherent and independently squeezed probes.

Compare cold-start E[T-1] with log(ARL)/D. No correlated-block worst-case or
Lorden theorem is inferred from these finite-grid numerical experiments.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .models import coherent_pair,squeezed_pair
from .arl_solver import calibrate_numerical,mean_stopping_time
from .delay_solver import conditional_delay_numerical
from .rare_event import estimate_regenerative_arl,cycle_arrays
from .simulation import simulate_stopping_times,summarize_delay

TARGETS=(100.,1000.,10000.,100000.,2000000.,10000000.,100000000.)


def pair(kind):
    if kind not in ('coherent','squeezed'):raise ValueError('unknown probe family')
    return (coherent_pair if kind=='coherent' else squeezed_pair)(N=100,Na=.1,eta1=.9,eta2=.85)


def run(output,cycles=150000,trials=40000,seed=2026102402):
    if any(isinstance(n,bool) or not isinstance(n,int) or n<2 for n in (cycles,trials)):
        raise ValueError('integer counts >=2 required')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    seeds=[int(c.generate_state(1,dtype=np.uint64)[0]) for c in np.random.SeedSequence(seed).spawn(44)]
    report={'source_url':'https://arxiv.org/html/2503.12276v2#S0.E1',
        'settings':{'cycles':cycles,'trials':trials,'seed':seed,'targets':list(TARGETS),
            'N':100,'Na':.1,'eta1':.9,'eta2':.85,'cells':512,'check_cells':1024,'mesh_tolerance':.005,
            'changed_horizon':5000},'cases':[],
        'scope':'Eq.1 finite-ARL scalar-model check, not an additional published figure reproduction or worst-case theorem.'}
    for kind in ('coherent','squeezed'):
        pre,post=pair(kind)
        for target in TARGETS:
            i=len(report['cases']);folder=output/f'case_{i:02d}';folder.mkdir()
            calibration=calibrate_numerical(pre,post,target,cells=512,check_cells=1024,mesh_tolerance=.005)
            h=calibration['threshold']
            delays=[mean_stopping_time(pre,post,h,cells=n,sampling=post) for n in (512,1024)]
            numeric=delays[-1]['mean_stopping_time']-1
            arl,num,den=estimate_regenerative_arl(pre,post,h,cycles=cycles,seed=seeds[3*i])
            np.savez_compressed(folder/'cycles.npz',**cycle_arrays(num,den))
            runs=simulate_stopping_times(pre,post,h,trials=trials,horizon=5000,change_at=1,seed=seeds[3*i+1])
            np.savez_compressed(folder/'changed.npz',stop=runs.stop,event=runs.event)
            mc=summarize_delay(runs)
            fine_arl=calibration['grid_checks'][-1]['mean_stopping_time']
            arl_z=float(abs(arl['arl']-fine_arl)/arl['arl_se'])
            delay_z=None if mc['conditional_delay'] is None else float(abs(mc['conditional_delay']-numeric)/mc['conditional_delay_se'])
            refinement=abs(delays[-1]['mean_stopping_time']/delays[0]['mean_stopping_time']-1)
            prediction=float(np.log(fine_arl)/post.kl(pre))
            case={'kind':kind,'target_arl':target,'folder':folder.name,'calibration':calibration,
                'post_grid_checks':delays,'numeric_arl':fine_arl,'numeric_delay':numeric,
                'numeric_post_samples':numeric+1,'delay_grid_relative_change':float(refinement),
                'kl':post.kl(pre),'first_order_prediction':prediction,'delay_to_first_order_ratio':numeric/prediction,
                'arl_mc':arl,'delay_mc':mc,'arl_z':arl_z,'delay_z':delay_z,
                'seeds':{'arl':seeds[3*i],'cold':seeds[3*i+1]},
                'passed':bool(arl_z<=4 and delay_z is not None and delay_z<=4 and refinement<.001
                    and arl['success_weight_ess']>=1000 and arl['unfinished_pre_cycles']==0
                    and arl['unfinished_post_cycles']==0 and mc['false_alarms']==0
                    and mc['post_change_censored']==0)}
            if target==2000000.:
                warm_seed=seeds[42 if kind=='coherent' else 43]
                warm=simulate_stopping_times(pre,post,h,trials=trials,horizon=5000,change_at=1000,seed=warm_seed)
                np.savez_compressed(folder/'warm_changed.npz',stop=warm.stop,event=warm.event)
                stats=summarize_delay(warm)
                numerical=conditional_delay_numerical(pre,post,h,change_at=1000,cells=1024)
                z=None if stats['conditional_delay'] is None else float(abs(stats['conditional_delay']-numerical['conditional_delay'])/stats['conditional_delay_se'])
                case['warm']={'seed':warm_seed,'mc':stats,'numeric':numerical,'z':z}
                case['passed']=bool(case['passed'] and z is not None and z<=4 and stats['post_change_censored']==0)
            report['cases'].append(case)
            print(f'{kind}, ARL={target:g}: passed={case["passed"]}',flush=True)
    report['passed']=all(c['passed'] for c in report['cases'])
    report['rows']=[{k:c[k] for k in ('kind','target_arl','numeric_arl','numeric_delay','numeric_post_samples',
        'kl','first_order_prediction','delay_to_first_order_ratio','arl_z','delay_z','passed')} for c in report['cases']]
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (output/'scaling.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(report['rows'][0]));w.writeheader();w.writerows(report['rows'])
    plot(output,report)
    return report


def plot(output,report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for kind in ('coherent','squeezed'):
        rows=[c for c in report['cases'] if c['kind']==kind];x=[c['numeric_arl'] for c in rows]
        line,=axes[0].plot(x,[c['numeric_delay'] for c in rows],label=f'{kind}: first-step solution')
        color=line.get_color()
        axes[0].plot(x,[c['first_order_prediction'] for c in rows],'--',color=color,label=f'{kind}: log(ARL)/KL')
        axes[0].errorbar(x,[c['delay_mc']['conditional_delay'] for c in rows],
            yerr=[1.96*c['delay_mc']['conditional_delay_se'] for c in rows],fmt='o',color=color,capsize=3)
        axes[1].plot(x,[c['delay_to_first_order_ratio'] for c in rows],'-o',color=color,label=kind)
    axes[0].set(xscale='log',xlabel='No-change average run length (pulses)',ylabel='Cold-start delay E[T−1] (pulses)',title='Finite-ARL delay versus first-order prediction')
    axes[1].axhline(1,color='black',ls=':',label='First-order asymptotic reference')
    axes[1].set(xscale='log',xlabel='No-change average run length (pulses)',ylabel='Calculated delay / [log(ARL)/KL]',title='Approach toward the asymptotic scale')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    for ext in ('png','svg'):fig.savefig(output/f'delay_scaling.{ext}',dpi=180)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='results_delay_scaling');a=p.parse_args()
    r=run(a.output);raise SystemExit(0 if r['passed'] else 1)


if __name__=='__main__':main()
