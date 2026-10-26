"""Check the finite-ARL Eq.1 experiment, including first alarms and raw cycles."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .delay_scaling import pair,TARGETS
from .arl_solver import mean_stopping_time
from .delay_solver import conditional_delay_numerical
from .rare_event import Cycles,summarize_cycles,estimate_regenerative_arl,cycle_arrays
from .simulation import Runs,summarize_delay,simulate_stopping_times


def audit(root,regenerate=False):
    root=Path(root);r=json.loads((root/'report.json').read_text());s=r['settings'];checks=[]
    def check(name,value):checks.append({'name':name,'passed':bool(value)})
    seeds=[int(c.generate_state(1,dtype=np.uint64)[0]) for c in np.random.SeedSequence(s['seed']).spawn(44)]
    check('complete predeclared grid',[(c['kind'],c['target_arl']) for c in r['cases']]==[(k,t) for k in ('coherent','squeezed') for t in TARGETS])
    for i,c in enumerate(r['cases']):
        prefix=f'case {i}: ';pre,post=pair(c['kind']);h=c['calibration']['threshold'];folder=root/c['folder']
        check(prefix+'independent seeds',c['seeds']=={'arl':seeds[3*i],'cold':seeds[3*i+1]})
        check(prefix+'KL',c['kl']==post.kl(pre))
        grids=[mean_stopping_time(pre,post,h,cells=n) for n in (256,512,1024)]
        check(prefix+'numerical ARL grids',grids==c['calibration']['grid_checks'])
        check(prefix+'target calibrated',abs(grids[-1]['mean_stopping_time']/c['target_arl']-1)<=s['mesh_tolerance'])
        postgrids=[mean_stopping_time(pre,post,h,cells=n,sampling=post) for n in (512,1024)]
        check(prefix+'post first-step grids',postgrids==c['post_grid_checks'])
        delay=postgrids[-1]['mean_stopping_time']-1;arl=grids[-1]['mean_stopping_time']
        check(prefix+'sample-count convention',delay==c['numeric_delay'] and delay+1==c['numeric_post_samples'])
        check(prefix+'Eq1 prediction',np.log(arl)/post.kl(pre)==c['first_order_prediction'])
        check(prefix+'Eq1 ratio',delay/c['first_order_prediction']==c['delay_to_first_order_ratio'])
        with np.load(folder/'cycles.npz',allow_pickle=False) as raw:
            cycles=[Cycles(*(raw[f'{law}_{f}'] for f in ('length','terminal_sum','upper','complete')),
                pre,post,h,law,c['arl_mc']['max_steps']) for law in ('pre','post')]
            stats=summarize_cycles(*cycles)
            check(prefix+'cycle counts',all(len(x.length)==s['cycles'] for x in cycles))
            check(prefix+'raw weighted ARL',all(stats[k]==c['arl_mc'][k] for k in stats))
        z=abs(stats['arl']-arl)/stats['arl_se']
        check(prefix+'ARL independent agreement',z==c['arl_z'] and z<=4)
        with np.load(folder/'changed.npz',allow_pickle=False) as raw:
            run=Runs(raw['stop'],raw['event'],s['changed_horizon'],1);mc=summarize_delay(run)
            check(prefix+'cold raw first-alarm summary',mc==c['delay_mc'])
            check(prefix+'cold complete and no pre-change alarms',mc['false_alarms']==0 and mc['post_change_censored']==0 and mc['trials']==s['trials'])
            if regenerate:
                fresh=simulate_stopping_times(pre,post,h,trials=s['trials'],horizon=s['changed_horizon'],change_at=1,seed=seeds[3*i+1])
                check(prefix+'exact cold sample regeneration',np.array_equal(fresh.stop,run.stop) and np.array_equal(fresh.event,run.event))
        delay_z=abs(mc['conditional_delay']-delay)/mc['conditional_delay_se']
        check(prefix+'delay independent agreement',delay_z==c['delay_z'] and delay_z<=4)
        if regenerate:
            stats2,a,b=estimate_regenerative_arl(pre,post,h,cycles=s['cycles'],seed=seeds[3*i])
            check(prefix+'exact ARL regeneration',stats2==c['arl_mc'])
            with np.load(folder/'cycles.npz',allow_pickle=False) as raw:
                check(prefix+'exact cycle-array regeneration',all(np.array_equal(raw[k],v) for k,v in cycle_arrays(a,b).items()))
        if 'warm' in c:
            warm=c['warm'];seed=seeds[42 if c['kind']=='coherent' else 43]
            check(prefix+'warm seed',warm['seed']==seed)
            numerical=conditional_delay_numerical(pre,post,h,change_at=1000,cells=1024)
            check(prefix+'warm first-step solution',numerical==warm['numeric'])
            with np.load(folder/'warm_changed.npz',allow_pickle=False) as raw:
                run=Runs(raw['stop'],raw['event'],s['changed_horizon'],1000);mc=summarize_delay(run)
                check(prefix+'warm raw summary',mc==warm['mc'])
                check(prefix+'warm first alarms and censoring',mc['post_change_censored']==0 and mc['survivors']+mc['false_alarms']==s['trials'])
                if regenerate:
                    fresh=simulate_stopping_times(pre,post,h,trials=s['trials'],horizon=s['changed_horizon'],change_at=1000,seed=seed)
                    check(prefix+'exact warm samples',np.array_equal(fresh.stop,run.stop) and np.array_equal(fresh.event,run.event))
            z=abs(mc['conditional_delay']-numerical['conditional_delay'])/mc['conditional_delay_se']
            check(prefix+'warm independent agreement',z==warm['z'] and z<=4)
        check(prefix+'reported status',c['passed'])
    with (root/'scaling.csv').open() as f:rows=list(csv.DictReader(f))
    check('CSV rows',rows==[{k:str(v) for k,v in row.items()} for row in r['rows']])
    check('report passed',r['passed'] and all(c['passed'] for c in r['cases']))
    result={'passed':all(c['passed'] for c in checks),'count':len(checks),'regenerated':regenerate,'checks':checks}
    (root/'audit.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root');p.add_argument('--regenerate',action='store_true');a=p.parse_args()
    r=audit(a.root,a.regenerate);print(f"{sum(x['passed'] for x in r['checks'])}/{r['count']} checks passed")
    raise SystemExit(0 if r['passed'] else 1)


if __name__=='__main__':main()
