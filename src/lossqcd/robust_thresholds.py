"""Conservative threshold calibration and independent sampled-grid validation.

The detector's score models are fixed; only its threshold changes. Design
targets a 7.5% ARL buffer over two million at the least favorable sampled
initial transmission. Five fixed grid points are tested independently after
selection. These empirical checks do not prove a continuum minimax guarantee.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from .channel_diagnostics import sha256
from .conditional_detection import simulate_conditional
from .conditional_rare_event import estimate_conditional_arl,conditional_cycle_arrays,ConditionalCycles
from .batched_rare_event import pool_cycle_summaries
from .simulation import summarize_delay


GRID=(.899,.8995,.9,.9005,.901)


def detector_args(model,eta):
    return dict(kind=model['kind'],n=model['n'],N=100.,Na=.1,eta1=.9,eta2=.85,sampling_eta1=eta)


def threshold_update(original,current,worst_arl,design_target):
    if not all(np.isfinite(v) for v in (original,current,worst_arl,design_target)) or min(original,current)<0 or min(worst_arl,design_target)<=0:
        raise ValueError('finite nonnegative thresholds and positive ARLs required')
    return float(max(original,current+np.log(design_target/worst_arl)))


def approximate_grid_lower_bound(arl,se,points=25,alpha=.05):
    if (not np.isfinite(arl) or arl<=0 or not np.isfinite(se) or se<0
        or isinstance(points,bool) or not isinstance(points,int) or points<1 or not 0<alpha<1):
        raise ValueError('invalid ARL, SE or multiplicity settings')
    z=float(norm.isf(alpha/points))
    return {'lower':float(arl*np.exp(-z*se/arl)),'z':z,'points':points,'alpha':alpha,
        'method':'Bonferroni one-sided nominal 95% family coverage using asymptotic log-scale delta SE; not a certified bound'}


def calibrate(model,seeds,cycles=100000,design_target=2150000.,iterations=4):
    if len(seeds)!=len(GRID)*iterations or len(set(seeds))!=len(seeds):
        raise ValueError('one independent seed per calibration evaluation required')
    original=model['calibration']['threshold'];h=original;history=[]
    for iteration in range(iterations):
        evaluations=[]
        for j,eta in enumerate(GRID):
            stats,_,_=estimate_conditional_arl(threshold=h,cycles=cycles,
                seed=seeds[iteration*len(GRID)+j],proposal='adaptive',proposal_offset=2.,**detector_args(model,eta))
            evaluations.append({'true_eta1':eta,'summary':stats})
        worst=min(evaluations,key=lambda e:e['summary']['arl'])
        history.append({'threshold':h,'evaluations':evaluations,'worst_eta1':worst['true_eta1'],
            'minimum_arl':worst['summary']['arl']})
        print(f"Calibration {model['label']} {iteration+1}/{iterations}: h={h:.5f}, minimum sampled ARL={worst['summary']['arl']:.0f}",flush=True)
        if iteration<iterations-1:
            h=threshold_update(original,h,worst['summary']['arl'],design_target)
    return {'threshold':h,'original_threshold':original,'history':history,'seeds':seeds,
        'design_target':design_target,'iterations':iterations,'cycles_per_law':cycles,
        'method':'fixed four-round sampled-grid log-ARL updates; no holdout or author-delay fitting'}


def cycle_pair(data,summary):
    p=summary['parameters'];params=tuple(p[k] for k in ('N','Na','eta1','eta2'))
    if 'sampling_eta1' in p:params+=(p['sampling_eta1'],)
    def build(prefix,proposal):
        arrays=[data[f'{prefix}_{f}'] for f in ('length','raw_sum','page','upper','complete')]
        return ConditionalCycles(*arrays,summary['kind'],summary['n'],summary['threshold'],proposal,params,
            data['post_log_weight'] if prefix=='post' else None,2. if proposal=='adaptive' else None)
    return build('pre','pre'),build('post','adaptive')


def source_models(source):
    return [m for m in source['models'] if m['figure']==7 and np.isclose(m['Na'],.1)]


def run(source_root,output,calibration_cycles=100000,validation_cycles=250000,trials=20000,seed=2026103001):
    for value in (calibration_cycles,validation_cycles,trials):
        if isinstance(value,bool) or not isinstance(value,int) or value<2:raise ValueError('sample counts must be integers >=2')
    source_root,output=Path(source_root).resolve(),Path(output).resolve()
    source=json.loads((source_root/'report.json').read_text());models=source_models(source)
    if len(models)!=5 or not all(m['confirmed'] for m in models):raise ValueError('five confirmed source models required')
    if any(source['settings'][k]!=v for k,v in (('N',100),('eta1',.9),('eta2',.85))):raise ValueError('unsupported source parameters')
    hashes={str(source_root/'report.json'):sha256(source_root/'report.json')}
    for m in models:
        p=source_root/Path(m['folder']).name/'changed_runs.npz';hashes[str(p)]=sha256(p)
    output.mkdir(parents=True,exist_ok=False)
    report={'source_root':str(source_root),'source_hashes':hashes,'models':[],'rows':[],
        'settings':{'seed':seed,'calibration_cycles':calibration_cycles,'validation_cycles':validation_cycles,
            'changed_trials':trials,'horizon':5000,'grid':list(GRID),'target_arl':2000000.,
            'design_target':2150000.,'calibration_iterations':4,'validation_points':25,'family_alpha':.05},
        'scope':'Sampled-grid conservative thresholds at N=100, Na=0.1, assumed eta1=.9/eta2=.85, true eta2=.85. '
            'This is a robustness extension, not a paper-figure reproduction or a continuum guarantee.'}
    children=np.random.SeedSequence(seed).spawn(200)
    seed_at=lambda i:int(children[i].generate_state(1,dtype=np.uint64)[0])
    # Selection uses streams 0:100; validation/delay use 100:200 only.
    for index,model in enumerate(models):
        fitted=calibrate(model,[seed_at(20*index+j) for j in range(20)],calibration_cycles)
        report['models'].append({'kind':model['kind'],'n':model['n'],'label':model['label'],'calibration':fitted})
    for model_index,model in enumerate(models):
        h=report['models'][model_index]['calibration']['threshold'];baseline_h=model['calibration']['threshold']
        for eta in GRID:
            i=len(report['rows']);repeats=[];files=[];seeds=[]
            for repeat in range(2):
                run_seed=seed_at(100+4*i+repeat)
                stats,pre,post=estimate_conditional_arl(threshold=h,cycles=validation_cycles,
                    seed=run_seed,proposal='adaptive',proposal_offset=2.,**detector_args(model,eta))
                filename=f'cycles_{i:02d}_{repeat}.npz';np.savez_compressed(output/filename,**conditional_cycle_arrays(pre,post))
                repeats.append(stats);files.append(filename);seeds.append(run_seed)
            pooled=pool_cycle_summaries(repeats)
            z=float(abs(repeats[0]['arl']-repeats[1]['arl'])/np.hypot(repeats[0]['arl_se'],repeats[1]['arl_se']))
            lower=approximate_grid_lower_bound(pooled['arl'],pooled['arl_se'])
            delays={}
            for j,(label,threshold) in enumerate((('baseline',baseline_h),('conservative',h))):
                run_seed=seed_at(100+4*i+2+j)
                runs=simulate_conditional(threshold=threshold,trials=trials,horizon=5000,
                    change_at=model['change_at'],seed=run_seed,**detector_args(model,eta))
                filename=f'delay_{i:02d}_{label}.npz';np.savez_compressed(output/filename,stop=runs.stop,event=runs.event)
                delays[label]={'seed':run_seed,'file':filename,'summary':summarize_delay(runs)}
            a,b=delays['baseline']['summary'],delays['conservative']['summary']
            cost=None if a['conditional_delay'] is None or b['conditional_delay'] is None else b['conditional_delay']-a['conditional_delay']
            cost_se=None if cost is None else float(np.hypot(a['conditional_delay_se'],b['conditional_delay_se']))
            passed=all(r['diagnostics_passed'] for r in repeats) and pooled['diagnostics_passed'] and z<=4
            passed=passed and lower['lower']>=2000000 and cost is not None
            row={'label':model['label'],'kind':model['kind'],'n':model['n'],'true_eta1':eta,
                'change_at':model['change_at'],'threshold':h,'baseline_threshold':baseline_h,
                'seeds':seeds,'cycle_files':files,'repeats':repeats,'pooled':pooled,'repeat_z':z,
                'grid_lower_bound':lower,'delays':delays,'delay_cost':cost,'delay_cost_se':cost_se,'passed':bool(passed)}
            report['rows'].append(row)
            print(f"Holdout {i+1}/25 {model['label']} eta1={eta}: ARL={pooled['arl']:.0f}, adjusted lower={lower['lower']:.0f}, passed={passed}",flush=True)
    report['source_unchanged']=all(sha256(p)==v for p,v in hashes.items())
    report['all_passed']=report['source_unchanged'] and all(r['passed'] for r in report['rows'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    write_outputs(report,output)
    return report


def exported_rows(report):
    return [{'label':r['label'],'n':r['n'],'true_eta1':r['true_eta1'],'baseline_threshold':r['baseline_threshold'],
        'conservative_threshold':r['threshold'],'arl':r['pooled']['arl'],'arl_se':r['pooled']['arl_se'],
        'adjusted_lower_arl':r['grid_lower_bound']['lower'],'baseline_delay':r['delays']['baseline']['summary']['conditional_delay'],
        'conservative_delay':r['delays']['conservative']['summary']['conditional_delay'],
        'delay_cost':r['delay_cost'],'delay_cost_se':r['delay_cost_se'],'passed':r['passed']} for r in report['rows']]


def write_outputs(report,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=exported_rows(report)
    with (output/'validation.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for model in report['models']:
        subset=[r for r in report['rows'] if r['label']==model['label']];x=[r['true_eta1'] for r in subset]
        axes[0].plot(x,[r['grid_lower_bound']['lower']/1e6 for r in subset],'o-',label=model['label'])
        axes[1].errorbar(x,[r['delay_cost'] for r in subset],yerr=[1.96*r['delay_cost_se'] for r in subset],
            fmt='o-',capsize=3,label=model['label'])
    axes[0].axhline(2,color='black',ls='--',label='Required ARL floor')
    axes[0].set(ylabel='Adjusted approximate lower ARL (million pulses)',title='Independent holdouts at frozen thresholds')
    axes[1].set(ylabel='Additional conditional delay (pulses)',title='Cost relative to original source threshold')
    for ax in axes:
        ax.set_xticks(GRID,[f'{eta:.4f}' for eta in GRID]);ax.set_xlabel('True initial transmission eta1')
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Conservative sampled-grid calibration: false-alarm protection and delay cost\nN=100, Na=0.1; fixed assumed channels 0.9/0.85; true post-change channel 0.85',fontsize=11)
    for ext in ('png','svg'):fig.savefig(output/f'protection_and_cost.{ext}',dpi=180)
    plt.close(fig)
    table=[]
    for m in report['models']:
        subset=[r for r in report['rows'] if r['label']==m['label']];cal=m['calibration']
        table.append(f"| {m['label']} | {cal['original_threshold']:.5f} | {cal['threshold']:.5f} | "
            f"{min(r['grid_lower_bound']['lower'] for r in subset):,.0f} | "
            f"{min(r['delay_cost'] for r in subset):.2f}–{max(r['delay_cost'] for r in subset):.2f} |")
    table='\n'.join(table)
    (output/'RESULTS.md').write_text(f"""# Conservative thresholds for uncertain initial transmission

This extends the paper implementation with a robustness experiment. Score laws
remain at eta1=0.9/eta2=0.85, N=100, Na=0.1; true post-change eta2 stays 0.85.
Five initial transmissions, 0.899 through 0.901 in steps of 0.0005, are sampled.
Original results remain unchanged. Only the Page alarm threshold is raised.

Four fixed calibration rounds use independent {report['settings']['calibration_cycles']:,}-cycle ensembles per law
at every grid point, targeting a minimum sampled ARL of 2.15 million. These
calibration streams are separate from all holdout and changed-run streams.
Thresholds are then frozen. No holdout results or published delays are fitted.

Independent validation passes at **{sum(r['passed'] for r in report['rows'])}/25 sampled settings**;
overall status: **{report['all_passed']}**. Each setting has two independent
{report['settings']['validation_cycles']:,}-cycle ensembles per sampling law. Raw
cycle moments are pooled before forming the ARL ratio. Each adjusted approximate
lower ARL must exceed two million, with ESS, largest-weight, relative-SE and
repeatability checks also passing.

| Detector | Original h | Conservative h | Smallest adjusted lower ARL | Observed delay cost (pulses) |
| --- | --- | --- | --- | --- |
{table}

The lower bounds use one-sided Bonferroni allocation of nominal 5% error over
25 points and asymptotic log-scale delta SEs. This is approximate sampled-grid
evidence, not certified coverage or a guarantee for every point in the continuous
interval. No monotonicity in true transmission is assumed for correlated scores.
Renewal boundaries and exact actual-null conditional mixture corrections are
retained; no iid solver is substituted for correlated data.

Delay costs use independent fresh {report['settings']['changed_trials']:,}-trial
ensembles at the original and conservative thresholds, at each true initial
channel. The change is at pulse 1005; delay is E[T−1005 | T>=1005]. SEs for the
differences combine the two independent errors. Early first alarms stay separate,
and censored survivors make unrestricted means unavailable. The plotted delay
error bars are pointwise ±1.96 Monte Carlo SE, not simultaneous bounds.

The source evidence is checksum-protected. Calibration histories and seeds,
all holdout cycle arrays and path weights, baseline/conservative first-alarm runs,
CSV and plots are saved. Calibration arrays can be regenerated from recorded
seeds. This does not reconcile the publication's latency gap or identify author
parameters. The conservative thresholds have deliberately unequal ARLs away
from the least favorable sampled channel, so this is not exact ARL matching
across detectors at every channel.

Next: test channels between the calibration points with fresh samples, then
study estimating the initial channel from a separate training period.
""")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',default='results_author_comparison');p.add_argument('--output',default='results_robust_thresholds')
    p.add_argument('--calibration-cycles',type=int,default=100000);p.add_argument('--validation-cycles',type=int,default=250000)
    p.add_argument('--trials',type=int,default=20000);p.add_argument('--seed',type=int,default=2026103001)
    a=p.parse_args();run(a.source,a.output,a.calibration_cycles,a.validation_cycles,a.trials,a.seed)


if __name__=='__main__':main()
