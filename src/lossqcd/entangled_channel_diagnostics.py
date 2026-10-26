"""Held-out correlated-block predictions at the existing scalar-fitted loss.

No new fitting, threshold tuning or ARL calibration. Two independent direct
ensembles per setting quantify simulation repeatability, not author uncertainty.
"""
import argparse
import csv
import json
import platform
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .conditional_detection import simulate_conditional
from .simulation import Runs, summarize_delay


def correlated_models(source):
    return [m for m in source['models'] if m['figure']==7 and m['n']>1]


def simulation_arguments(model, settings, eta, trials, seed):
    return dict(kind=model['kind'], n=model['n'], N=settings['N'], Na=model['Na'],
        eta1=settings['eta1'], eta2=settings['eta2'], threshold=model['calibration']['threshold'],
        trials=trials, horizon=settings['horizon'], change_at=model['change_at'],
        seed=seed, sampling_eta2=eta)


def pool_runs(runs):
    if len({(r.horizon,r.change_at) for r in runs})!=1:
        raise ValueError('ensembles must share horizon and change point')
    return Runs(np.concatenate([r.stop for r in runs]),np.concatenate([r.event for r in runs]),
        runs[0].horizon,runs[0].change_at)


def repeat_z(a,b):
    if a['conditional_delay'] is None or b['conditional_delay'] is None:
        return None
    denominator=np.hypot(a['conditional_delay_se'],b['conditional_delay_se'])
    return float(abs(a['conditional_delay']-b['conditional_delay'])/denominator) if denominator>0 else None


def run(scalar_root, output, trials=30000, seed=2026102701):
    scalar_root,output=Path(scalar_root).resolve(),Path(output).resolve()
    if not isinstance(trials,int) or trials<2:
        raise ValueError('trials must be an integer >=2')
    scalar=json.loads((scalar_root/'report.json').read_text())
    if not scalar['all_checks_passed'] or len(scalar['rows'])!=18:
        raise ValueError('expected completed scalar diagnostic')
    source_root=Path(scalar['source_root'])
    source=json.loads((source_root/'report.json').read_text())
    settings=source['settings'];models=correlated_models(source)
    if len(models)!=27 or not all(m['confirmed'] for m in models):
        raise ValueError('expected 27 confirmed correlated source settings')
    hashes=dict(scalar['source_hashes'])
    paths=[scalar_root/'report.json',scalar_root/'predictions.csv']
    paths.extend(scalar_root/r['runs'] for r in scalar['rows'])
    paths.extend(source_root/Path(m['folder']).name/'changed_runs.npz' for m in models)
    hashes.update({str(p):sha256(p) for p in paths})
    if not all(sha256(p)==h for p,h in hashes.items()):
        raise ValueError('source evidence changed')
    output.mkdir(parents=True,exist_ok=False)
    report={'scalar_root':str(scalar_root),'source_root':str(source_root),'source_hashes':hashes,
        'settings':{'trials_per_repeat':trials,'seed':seed,'N':settings['N'],
            'eta1':settings['eta1'],'assumed_eta2':settings['eta2'],
            'actual_eta2':scalar['fit']['actual_eta2'],'horizon':settings['horizon']},
        'versions':{'python':platform.python_version(),'numpy':np.__version__},'rows':[],
        'interpretation':'No new fits. Conditional-LLR Page detector; true mixed-loss covariance. '
            'Pre-change law and threshold unchanged: original ARL evidence reused, not re-estimated. '
            'Author uncertainties and actual Figure 7 eta are unavailable. This is a parameter diagnostic.'}
    streams=np.random.SeedSequence(seed).spawn(2*len(models))
    for index,model in enumerate(models):
        runs=[];summaries=[];seeds=[];files=[]
        for repeat in range(2):
            run_seed=int(streams[2*index+repeat].generate_state(1,dtype=np.uint64)[0])
            batch=simulate_conditional(**simulation_arguments(model,settings,report['settings']['actual_eta2'],trials,run_seed))
            filename=f'runs_{index:02d}_{repeat}.npz'
            np.savez_compressed(output/filename,stop=batch.stop,event=batch.event)
            runs.append(batch);summaries.append(summarize_delay(batch));seeds.append(run_seed);files.append(filename)
        pooled=summarize_delay(pool_runs(runs));z=repeat_z(*summaries)
        published=next(r for r in source['figure7_rows'] if r['label']==model['label'] and np.isclose(r['Na'],model['Na']))
        relative=None if pooled['conditional_delay'] is None else pooled['conditional_delay']/published['author_delay']-1
        row={'Na':model['Na'],'n':model['n'],'label':model['label'],'kind':model['kind'],
            'change_at':model['change_at'],'threshold':model['calibration']['threshold'],
            'author_delay':published['author_delay'],'original_mc_delay':published['our_delay'],
            'source_arl':model['validation']['arl'],'source_arl_se':model['validation']['arl_se'],
            'seeds':seeds,'files':files,'repeats':summaries,'pooled':pooled,'repeat_z':z,
            'relative_difference':relative,'within_5_percent':relative is not None and abs(relative)<=.05,
            'checks_passed':z is not None and z<=4 and pooled['post_change_censored']==0}
        report['rows'].append(row)
        print(f"{index+1}/27 {model['label']} Na={model['Na']:g}: delay={pooled['conditional_delay']}, repeat z={z:.2f}",flush=True)
    report['entangled_within_5_percent']=sum(r['within_5_percent'] for r in report['rows'])
    report['combined_heldout_within_5_percent']=scalar['heldout_within_5_percent']+report['entangled_within_5_percent']
    report['source_unchanged']=all(sha256(p)==h for p,h in hashes.items())
    report['all_checks_passed']=report['source_unchanged'] and all(r['checks_passed'] for r in report['rows'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    write_outputs(report,scalar,output)
    return report


def combined_rows(report,scalar):
    rows=[]
    for r in scalar['rows']:
        rows.append({'Na':r['Na'],'label':r['label'],'n':1,'training':r['training'],
            'author_delay':r['author_delay'],'original_delay':r['original_mc_delay'],
            'prediction':r['numerical_1024']['conditional_delay'],
            'prediction_method':'scalar deterministic 1024 cells','mc_se':r['mc']['conditional_delay_se'],
            'relative_difference':r['relative_difference'],'within_5_percent':r['within_5_percent']})
    for r in report['rows']:
        rows.append({'Na':r['Na'],'label':r['label'],'n':r['n'],'training':False,
            'author_delay':r['author_delay'],'original_delay':r['original_mc_delay'],
            'prediction':r['pooled']['conditional_delay'],'prediction_method':'two pooled direct MC ensembles',
            'mc_se':r['pooled']['conditional_delay_se'],'relative_difference':r['relative_difference'],
            'within_5_percent':r['within_5_percent']})
    return rows


def write_outputs(report,scalar,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=combined_rows(report,scalar)
    with (output/'combined_predictions.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(1,3,figsize=(13,4.6),layout='constrained')
    for ax,n in zip(axes,(2,4,8)):
        subset=[r for r in rows if r['n']==n];x=[r['Na'] for r in subset]
        ax.plot(x,[r['author_delay'] for r in subset],'o-',label='Published means')
        ax.plot(x,[r['original_delay'] for r in subset],':',label='Original matched simulation')
        ax.errorbar(x,[r['prediction'] for r in subset],yerr=[1.96*r['mc_se'] for r in subset],
            fmt='s-',capsize=3,label='Fixed-loss diagnostic ±1.96 SE')
        ax.set(xlabel='Added energy Na',ylabel='Conditional delay (pulses)',title=f'Entangled block n={n}')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle(f"True eta2={report['settings']['actual_eta2']:.6f}; all detectors still assume 0.85\nNo entangled means fitted; original thresholds and pre-change laws fixed",fontsize=11)
    for ext in ('png','svg'):
        fig.savefig(output/f'entangled_channel_predictions.{ext}',dpi=180)
    plt.close(fig)
    group_summary='\n'.join(f"- n={n}: {sum(r['within_5_percent'] for r in report['rows'] if r['n']==n)}/9 within ±5%." for n in (2,4,8))
    gap_table='\n'.join(
        f"| {n} | {100*np.mean([abs(r['original_mc_delay']/r['author_delay']-1) for r in report['rows'] if r['n']==n]):.2f}% | "
        f"{100*np.mean([abs(r['relative_difference']) for r in report['rows'] if r['n']==n]):.2f}% |"
        for n in (2,4,8))
    (output/'RESULTS.md').write_text(f"""# Entangled channel-loss diagnostic

We reused **eta2={report['settings']['actual_eta2']:.8f}**, previously fitted using
only the coherent Na=0.1 mean. No additional author means or thresholds were
fitted. Detector likelihoods still assume eta1=0.9 and eta2=0.85; N=100.
The true pre-change eta1 remains 0.9. The exact mixed-loss covariance is retained
when the physical change at pulse 1005 falls inside a block. Conditional history
resets at block boundaries; Page evidence persists, including through reflections.

**{report['entangled_within_5_percent']}/27 entangled means** fall within a descriptive ±5%
of published values. Including the earlier scalar diagnostic gives
**{report['combined_heldout_within_5_percent']}/44 held-out means** within ±5%; the one
fitted classical point is excluded from this total.

{group_summary}

Mean absolute relative differences across each nine-point energy grid:

| Block size | Original matched simulation | Common-loss diagnostic |
| --- | --- | --- |
{gap_table}

Each correlated setting has two independent {report['settings']['trials_per_repeat']:,}-trial
ensembles, pooled at the raw stopping-time level. All repeatability/censoring
checks pass: **{report['all_checks_passed']}**. Reported uncertainty is our Monte Carlo
SE only; author uncertainties are unavailable. ±5% is not a significance test.
Delay means E[T−1005 | T>=1005]; early first alarms are counted separately and
never restarted. A censored survivor makes unrestricted mean delay unavailable.

The entire no-change stopping law is unchanged, so original near-two-million
pulse ARL evidence applies unchanged; it is not re-estimated here. Unit tests
compare the streaming implementation with an independent dense-covariance
conditional-Gaussian calculation, including within-block changes and first alarms.
No iid scalar delay solver is applied to correlated blocks.

The fitted loss improves agreement but is a diagnostic assumption, not an
estimate of the authors' actual channel. Figure 7 does not save eta separately.
Figure 3's specified channel is not refitted. Original source files and scalar
results are protected by checksums; new raw ensembles, seeds and summaries are saved.
See combined_predictions.csv and entangled_channel_predictions.png.

Next: consolidate the remaining discrepancies and test sensitivity to the true
pre-change channel, which can also change false-alarm performance. Such a study
needs new ARL checks; the current unchanged-null argument would no longer apply.
""")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar-source',default='results_channel_diagnostics')
    parser.add_argument('--output',default='results_entangled_channel_diagnostics')
    parser.add_argument('--trials',type=int,default=30000)
    parser.add_argument('--seed',type=int,default=2026102701)
    args=parser.parse_args()
    run(args.scalar_source,args.output,args.trials,args.seed)


if __name__=='__main__':
    main()
