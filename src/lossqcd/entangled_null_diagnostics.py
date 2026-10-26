"""Correlated-detector robustness to true initial-channel mismatch.

Renewals occur at block boundaries with zero Page evidence, or at the first
alarm. Adaptive proposals mix the actual null and assumed post conditional
laws. Their likelihood correction is distinct from the assumed detector score.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .conditional_detection import simulate_conditional
from .conditional_rare_event import estimate_conditional_arl,conditional_cycle_arrays,ConditionalCycles
from .batched_rare_event import pool_cycle_summaries
from .simulation import summarize_delay,summarize_arl


def models_for_study(source):
    return [m for m in source['models'] if m['figure']==7 and m['n'] in (2,4,8) and np.isclose(m['Na'],.1)]


def parameters(n,eta):
    return dict(kind='entangled',n=n,N=100.,Na=.1,eta1=.9,eta2=.85,sampling_eta1=eta)


def read_cycle_pair(data,n,eta,h):
    params=(100.,.1,.9,.85)+((eta,) if eta!=.9 else ())
    def build(prefix,proposal):
        arrays=[data[f'{prefix}_{field}'] for field in ('length','raw_sum','page','upper','complete')]
        return ConditionalCycles(*arrays,'entangled',n,h,proposal,params,
            data['post_log_weight'] if prefix=='post' else None,2. if proposal=='adaptive' else None)
    return build('pre','pre'),build('post','adaptive')


def run(source_root,output,cycles=250000,trials=30000,seed=2026102901):
    source_root,output=Path(source_root).resolve(),Path(output).resolve()
    if any(isinstance(x,bool) or not isinstance(x,int) or x<2 for x in (cycles,trials)):
        raise ValueError('cycles and trials must be integers >=2')
    source=json.loads((source_root/'report.json').read_text());models=models_for_study(source)
    if len(models)!=3 or not all(m['confirmed'] for m in models):
        raise ValueError('expected three confirmed Na=0.1 correlated source detectors')
    if any(source['settings'][k]!=v for k,v in (('N',100),('eta1',.9),('eta2',.85))):
        raise ValueError('expected N=100 and assumed eta=.9/.85')
    hashes={str(source_root/'report.json'):sha256(source_root/'report.json')}
    for m in models:
        path=source_root/Path(m['folder']).name/'changed_runs.npz';hashes[str(path)]=sha256(path)
    output.mkdir(parents=True,exist_ok=False)
    report={'source_root':str(source_root),'source_hashes':hashes,'rows':[],'controls':[],
        'settings':{'cycles_per_law_per_repeat':cycles,'changed_trials':trials,'seed':seed,
            'true_eta1_grid':[.899,.9,.901],'proposal':'adaptive','proposal_offset':2.,
            'N':100.,'Na':.1,'assumed_eta1':.9,'eta2':.85,'horizon':5000,
            'control_cycles':100000,'control_trials':10000,'control_horizon':100000},
        'interpretation':'Robustness diagnostic, no author mean fitted. Original thresholds and likelihoods fixed. '
            'Fresh actual-null ARL estimates; no iid scalar solver applied to correlated data. '
            'Approximate Monte Carlo uncertainty, not certified bounds.'}
    children=np.random.SeedSequence(seed).spawn(39)
    make_seed=lambda i:int(children[i].generate_state(1,dtype=np.uint64)[0])
    for m in models:
        for eta in report['settings']['true_eta1_grid']:
            i=len(report['rows']);h=m['calibration']['threshold'];repeats=[];files=[];seeds=[]
            for repeat in range(2):
                run_seed=make_seed(3*i+repeat)
                stats,pre,post=estimate_conditional_arl(threshold=h,cycles=cycles,seed=run_seed,
                    proposal='adaptive',proposal_offset=2.,**parameters(m['n'],eta))
                filename=f'cycles_{i:02d}_{repeat}.npz'
                np.savez_compressed(output/filename,**conditional_cycle_arrays(pre,post))
                repeats.append(stats);files.append(filename);seeds.append(run_seed)
            pooled=pool_cycle_summaries(repeats)
            z=float(abs(repeats[0]['arl']-repeats[1]['arl'])/np.hypot(repeats[0]['arl_se'],repeats[1]['arl_se']))
            changed_seed=make_seed(3*i+2)
            runs=simulate_conditional(threshold=h,trials=trials,horizon=5000,
                change_at=m['change_at'],seed=changed_seed,**parameters(m['n'],eta))
            changed=summarize_delay(runs);filename=f'changed_{i:02d}.npz'
            np.savez_compressed(output/filename,stop=runs.stop,event=runs.event)
            passed=all(r['diagnostics_passed'] for r in repeats) and pooled['diagnostics_passed']
            passed=passed and z<=4 and changed['conditional_delay'] is not None
            report['rows'].append({'n':m['n'],'label':m['label'],'true_eta1':eta,'threshold':h,
                'change_at':m['change_at'],'repeats':repeats,'pooled':pooled,'repeat_z':z,
                'seeds':seeds,'cycle_files':files,'changed_seed':changed_seed,'changed_file':filename,
                'changed':changed,'checks_passed':bool(passed)})
            print(f"{i+1}/9 n={m['n']} eta1={eta}: ARL={pooled['arl']:.0f} ±{pooled['arl_se']:.0f}; repeat z={z:.2f}",flush=True)
    for m in models:
        for eta in (.899,.901):
            i=len(report['controls']);direct_seed=make_seed(27+2*i);estimate_seed=make_seed(28+2*i)
            runs=simulate_conditional(threshold=3,trials=10000,horizon=100000,seed=direct_seed,**parameters(m['n'],eta))
            direct=summarize_arl(runs)
            stats,pre,post=estimate_conditional_arl(threshold=3,cycles=100000,seed=estimate_seed,
                proposal='adaptive',proposal_offset=2.,**parameters(m['n'],eta))
            z=None if direct['arl'] is None else float(abs(direct['arl']-stats['arl'])/np.hypot(direct['arl_se'],stats['arl_se']))
            filename=f'control_{i}.npz';np.savez_compressed(output/filename,
                stop=runs.stop,event=runs.event,**conditional_cycle_arrays(pre,post))
            report['controls'].append({'n':m['n'],'true_eta1':eta,'threshold':3.,'direct_seed':direct_seed,
                'estimate_seed':estimate_seed,'file':filename,'direct':direct,'estimate':stats,'comparison_z':z,
                'checks_passed':bool(z is not None and z<=4 and stats['diagnostics_passed'])})
    report['source_unchanged']=all(sha256(p)==v for p,v in hashes.items())
    report['all_checks_passed']=report['source_unchanged'] and all(r['checks_passed'] for r in report['rows']+report['controls'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    write_outputs(report,output)
    return report


def exported_rows(report):
    result=[]
    for r in report['rows']:
        center=next(c for c in report['rows'] if c['n']==r['n'] and c['true_eta1']==.9)
        result.append({'n':r['n'],'true_eta1':r['true_eta1'],'threshold':r['threshold'],
            'arl':r['pooled']['arl'],'arl_se':r['pooled']['arl_se'],
            'arl_ratio_to_fresh_matched':r['pooled']['arl']/center['pooled']['arl'],
            'delay':r['changed']['conditional_delay'],'delay_se':r['changed']['conditional_delay_se'],
            'false_alarms':r['changed']['false_alarms'],'false_alarm_probability':r['changed']['false_alarm_probability'],
            'repeat_z':r['repeat_z'],'checks_passed':r['checks_passed']})
    return result


def write_outputs(report,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    exported=exported_rows(report)
    with (output/'sensitivity.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(exported[0]));writer.writeheader();writer.writerows(exported)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    for n in (2,4,8):
        rows=[r for r in report['rows'] if r['n']==n];x=[r['true_eta1'] for r in rows]
        axes[0].errorbar(x,[r['pooled']['arl']/1e6 for r in rows],
            yerr=[1.96*r['pooled']['arl_se']/1e6 for r in rows],fmt='o-',capsize=3,label=f'n={n}')
        axes[1].errorbar(x,[r['changed']['conditional_delay'] for r in rows],
            yerr=[1.96*r['changed']['conditional_delay_se'] for r in rows],fmt='o-',capsize=3,label=f'n={n}')
    axes[0].axhline(2,color='black',ls='--',label='Original target')
    axes[0].set(ylabel='No-change ARL (million pulses)',title='Fresh correlated renewal estimates')
    axes[1].set(ylabel='Conditional delay (pulses)',title='True post-change eta2 stays 0.85')
    for ax in axes:
        ax.set_xticks([.899,.9,.901],['0.899','0.900','0.901']);ax.set_xlabel('True pre-change transmission eta1')
        ax.axvline(.9,color='gray',ls=':');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Entangled initial-channel robustness at original thresholds\nAssumed eta1=0.9, eta2=0.85; N=100, Na=0.1; error bars: ±1.96 Monte Carlo SE',fontsize=11)
    for ext in ('png','svg'):
        fig.savefig(output/f'entangled_null_sensitivity.{ext}',dpi=180)
    plt.close(fig)
    table=[]
    for n in (2,4,8):
        rows=[r for r in exported if r['n']==n]
        table.append(f"| {n} | "+' | '.join(f"{r['arl']:,.0f} ± {r['arl_se']:,.0f}" for r in rows)+' |')
    table='\n'.join(table)
    ratios=exported_rows(report)
    low=[r['arl_ratio_to_fresh_matched'] for r in ratios if r['true_eta1']==.899]
    high=[r['arl_ratio_to_fresh_matched'] for r in ratios if r['true_eta1']==.901]
    delay_span=max(max(r['changed']['conditional_delay'] for r in report['rows'] if r['n']==n)
        -min(r['changed']['conditional_delay'] for r in report['rows'] if r['n']==n) for n in (2,4,8))
    (output/'RESULTS.md').write_text(f"""# Entangled initial-channel robustness

True initial transmission is 89.9%, 90% or 90.1%. The detector still assumes
90% before the change and 85% after it. True post-change transmission stays
85%; N=100, Na=0.1. Every threshold is frozen at the source calibration.
No author latency is fitted. Fresh estimates are required because the no-change
law changes. ARLs below are in pulses, ±1 Monte Carlo SE:

| Entangled block n | True eta1=0.899 | True eta1=0.900 | True eta1=0.901 |
| --- | --- | --- | --- |
{table}

Relative to each fresh matched-channel estimate, the lower-transmission cases
reduce ARL by {100*(1-max(low)):.1f}–{100*(1-min(low)):.1f}%, while higher-transmission
cases increase it by {100*(min(high)-1):.1f}–{100*(max(high)-1):.1f}%.
The largest observed delay span across the three initial channels for any block
size is {delay_span:.3f} pulse. Those are Monte Carlo point estimates, with SEs
saved separately; they are not deterministic bounds on true delay changes.

Nine settings each have two independent {report['settings']['cycles_per_law_per_repeat']:,}-cycle
ensembles per sampling law and a {report['settings']['changed_trials']:,}-trial changed
ensemble. Numerator and denominator moments are pooled before forming their
ratio; ARL ratios are never averaged. Six h=3 controls compare fresh renewal
estimates against independent 10,000-trial direct no-change simulations.
All diagnostics, repetitions, controls and censoring checks pass:
**{report['all_checks_passed']}**. Approximate intervals and diagnostic thresholds
are not finite-variance proofs or certified coverage guarantees.

Each renewal begins with an empty block and W=0, ending at a first alarm or a
block boundary with W=0. Reflection within a block retains its observation
history. The adaptive proposal mixes the actual null conditional law with the
assumed post conditional law. The exact path weight uses that whole mixture,
not its sampled component. Page evidence still uses the assumed pre/post
score. The weight and score therefore use different likelihood ratios under
mismatch and are saved separately. No clipping or self-normalization is used.

Changed observations retain the exact physical covariance for mixed-loss blocks.
The change is at pulse 1005; delay is E[T−1005 | T>=1005]. Early first alarms
are counted separately and never restarted; censored survivors make unrestricted
mean delay unavailable. Dense-covariance oracle tests check the physical sampler
and conditional likelihood weights; scalar reductions match numerical equations.
No iid solver is applied to correlated increments.

These are robustness results, not recovered author settings or a reconciliation
of the published latency gap. Original source evidence is checksum-protected.
Full cycle arrays, path weights, changed/control runs, seeds, CSV and plots are saved.

Next: use these measured sensitivities to design a detector that tolerates initial
channel uncertainty, with fresh false-alarm calibration and held-out validation.
""")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',default='results_author_comparison')
    p.add_argument('--output',default='results_entangled_null_diagnostics')
    p.add_argument('--cycles',type=int,default=250000)
    p.add_argument('--trials',type=int,default=30000)
    p.add_argument('--seed',type=int,default=2026102901)
    a=p.parse_args();run(a.source,a.output,a.cycles,a.trials,a.seed)


if __name__=='__main__':
    main()
