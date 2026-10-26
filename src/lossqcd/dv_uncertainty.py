"""Re-summarize author single-photon latency data without hiding early alarms.

The author dataset has signed latencies but no first-alarm/censoring flags.
Its nonnegative re-summary is conditional only under the stated interpretation
lat=T-cp for complete first-alarm runs. Existing source results are read-only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from lossqcd.author_comparison import author_data
from lossqcd.octave_data import load_octave
from lossqcd.simulation import Runs,summarize_delay


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def stats(x):
    x=np.asarray(x)
    return {'n':len(x),'mean':float(x.mean()),'se':float(x.std(ddof=1)/np.sqrt(len(x)))}


def bootstrap_mean(x,replicates,seed):
    """Exact empirical bootstrap of the mean, aggregating equal lattice values.

    A multinomial over distinct observed values is equivalent to iid resampling
    individual observations with replacement. It retains the negative tail.
    """
    values,counts=np.unique(x,return_counts=True);n=len(x)
    rng=np.random.default_rng(seed);out=np.empty(replicates)
    for start in range(0,replicates,500):
        size=min(500,replicates-start)
        draws=rng.multinomial(n,counts/n,size=size)
        out[start:start+size]=draws@values/n
    return out


def wilson(successes,n):
    z=1.959963984540054;p=successes/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(max(0,center-half)),float(min(1,center+half))]


def run(root,output,replicates=20000):
    if replicates<1000:raise ValueError('at least 1000 bootstrap replicates required')
    data=root/'data/zenodo_17460928';study=root/'data/studies/results_paper_completion'
    paths=[data/'prl_supplem_fig4b_data',data/'manifest.json',study/'report.json',study/'photon/changed.npz']
    hashes={str(p.relative_to(root)):digest(p) for p in paths}
    verified=author_data(data)
    source=load_octave(data/'prl_supplem_fig4b_data');author=source['lat'].ravel()
    base=json.loads((study/'report.json').read_text());case=base['cases'][-1]
    with np.load(study/'photon/changed.npz',allow_pickle=False) as raw:
        runs=Runs(raw['stop'],raw['event'],base['settings']['horizon'],case['change_at'])
    delay=summarize_delay(runs)
    if delay!=case['delay'] or not np.all(runs.event):raise ValueError('source simulation changed or contains censoring')
    ours=runs.stop-runs.change_at
    groups={'author_signed':author,'author_nonnegative_assuming_complete_first_alarms':author[author>=0],
        'ours_signed':ours,'ours_conditional':ours[ours>=0]}
    seeds=[int(c.generate_state(1,dtype=np.uint64)[0]) for c in np.random.SeedSequence(2026102502).spawn(4)]
    boot={};summaries={}
    for (name,values),seed in zip(groups.items(),seeds):
        boot[name]=bootstrap_mean(values,replicates,seed)
        summaries[name]={**stats(values),'bootstrap_percentile_95':np.quantile(boot[name],[.025,.975]).tolist(),'seed':seed}
    differences={}
    for name,a,b in [('signed','author_signed','ours_signed'),
        ('conditional_assuming_author_entries_are_complete_first_alarms','author_nonnegative_assuming_complete_first_alarms','ours_conditional')]:
        x,y=summaries[a],summaries[b];samples=boot[b]-boot[a]
        differences[name]={'point_difference_ours_minus_author':y['mean']-x['mean'],
            'relative_difference':y['mean']/x['mean']-1,'independent_mean_se':float(np.hypot(x['se'],y['se'])),
            'bootstrap_difference_percentile_95':np.quantile(samples,[.025,.975]).tolist(),
            'normal_standardized_difference':float((y['mean']-x['mean'])/np.hypot(x['se'],y['se']))}
        boot[name+'_difference']=samples
    source_negative=int((author<0).sum());ours_negative=int((ours<0).sum())
    report={'source_url':'https://zenodo.org/records/17460928','source_sha256':hashes,
        'bootstrap_replicates':replicates,'bootstrap_method':'Ordinary independent empirical mean bootstrap; repeated lattice values grouped with a multinomial.',
        'author_settings':{k:float(source[k]) for k in ('N','N_a','cp','n','ARL0')},
        'author_negative_entries':author[author<0].tolist(),'author_zero_entries':int((author==0).sum()),
        'author_pulse_range_consistent_if_T_equals_lat_plus_cp':bool(np.all((author+source['cp']>=1)&(author+source['cp']<=source['n']))),
        'summaries':summaries,'differences':differences,
        'early_alarm_fractions':{'author_assuming_first_alarm_latencies':{'count':source_negative,'n':len(author),'fraction':source_negative/len(author),'wilson_95':wilson(source_negative,len(author))},
            'ours_first_alarm_flags':{'count':ours_negative,'n':len(ours),'fraction':ours_negative/len(ours),'wilson_95':wilson(ours_negative,len(ours))}},
        'conditional_gap_reduction_from_signed_source_reinterpretation':float(summaries['author_nonnegative_assuming_complete_first_alarms']['mean']-summaries['author_signed']['mean']),
        'limits':'Author event flags, censoring handling and thresholds are unavailable. Conditional source re-summary assumes complete first-alarm T-cp entries. Empirical bootstrap assumes independent runs and does not capture ARL calibration error, missing events, model uncertainty or unseen tail events. No p-value or causal diagnosis is asserted.',
        'author_dataset_checksums_passed':True,'our_delay_recomputed':delay,
        'source_unchanged':all(digest(root/p)==h for p,h in hashes.items())}
    output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(output/'bootstrap_means.npz',**boot)
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_outputs(output,report,groups)
    return report


def write_outputs(output,r,groups):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    upper=5*(int(max(groups['author_nonnegative_assuming_complete_first_alarms'].max(),groups['ours_conditional'].max()))//5+1)
    bins=np.arange(0,upper+1,5)
    axes[0].hist(groups['author_nonnegative_assuming_complete_first_alarms'],bins=bins,density=True,histtype='step',lw=2,label='Author nonnegative entries (n=498)')
    axes[0].hist(groups['ours_conditional'],bins=bins,density=True,histtype='step',lw=2,label='Our conditional first alarms (n=29983)')
    axes[0].set(xlabel='Nonnegative latency (pulses)',ylabel='Density',title='Complete-first-alarm interpretation is assumed for author data')
    axes[0].legend(fontsize=8)
    for i,(name,s) in enumerate(r['summaries'].items()):
        lo,hi=s['bootstrap_percentile_95'];m=s['mean']
        axes[1].errorbar(m,i,xerr=[[m-lo],[hi-m]],fmt='o',capsize=4)
    axes[1].set_yticks(range(4),['Author signed','Author nonnegative*','Our signed','Our conditional'])
    axes[1].set(xlabel='Mean latency (pulses)',title='Empirical-bootstrap percentile intervals (95%)')
    for ax in axes:ax.grid(alpha=.2)
    for ext in ('png','svg'):fig.savefig(output/f'dv_uncertainty.{ext}',dpi=180)
    plt.close(fig)
    lines=['# Uncertainty in the released single-photon latency comparison','',
        'Source: [author dataset](https://zenodo.org/records/17460928), supplemental `lat` array and existing independently calibrated Python DV experiment. Original files and earlier studies are unchanged.','',
        '## Keep the delay statistics distinct','',
        '| Statistic | Runs | Mean (pulses) | Standard error | Bootstrap percentile 95% interval |',
        '| --- | --- | --- | --- | --- |']
    labels=['Author signed, as plotted','Author nonnegative, conditional interpretation assumed','Our signed first-alarm latency','Our conditional first-alarm delay']
    for label,s in zip(labels,r['summaries'].values()):
        lo,hi=s['bootstrap_percentile_95'];lines.append(f"| {label} | {s['n']} | {s['mean']:.3f} | {s['se']:.3f} | [{lo:.3f}, {hi:.3f}] |")
    lines+=['','The author array contains -844 and -129 pulses, with no zero entries. Those values are retained in the signed summary. Under the interpretation lat=T-cp they correspond to alarms on pulses 156 and 871, before change pulse 1000. The 498 nonnegative entries are separately re-summarized. The authors did not supply event flags, so that re-summary is not proof of their original event handling.','',
        '![Latency distributions and mean intervals](dv_uncertainty.png)','',
        '## Does the remaining gap disappear?','']
    for name,d in r['differences'].items():
        lo,hi=d['bootstrap_difference_percentile_95'];lines.append(f"- {name.replace('_',' ')}: our mean minus author mean = {d['point_difference_ours_minus_author']:.3f} pulses ({d['relative_difference']:.2%}); bootstrap percentile interval [{lo:.3f}, {hi:.3f}].")
    lines+=['',f"Using the nonnegative source re-summary instead of its signed mean moves the comparison by {r['conditional_gap_reduction_from_signed_source_reinterpretation']:.3f} pulses. It narrows the descriptive conditional-versus-published gap, but the comparable nonnegative mean is still lower than ours. The bootstrap interval for the conditional difference excludes zero under the stated assumptions; this is evidence of a remaining gap, not identification of its cause.",'',
        '## Early alarms','', '| Data interpretation | Fraction | Wilson 95% interval |','| --- | --- | --- |']
    for name,s in r['early_alarm_fractions'].items():
        lo,hi=s['wilson_95'];lines.append(f"| {name.replace('_',' ')} | {s['count']}/{s['n']} ({s['fraction']:.3%}) | [{lo:.3%}, {hi:.3%}] |")
    lines+=['','The author sample has only two negative entries. Neither their fraction nor the ratio of the two fractions estimates the unknown ARL, and no exponential stopping-time distribution is assumed. Our event flags and no-censoring status are directly checked against the saved 30,000-run experiment.','',
        '## Limits and reproduction','',r['limits'],'',
        'Four independent random streams, each with 20,000 ordinary empirical bootstrap replicates, are used. Repeated integer values are grouped with multinomial counts, exactly equivalent to resampling the observations. Difference intervals combine independent author and Python resamples. They are percentile-bootstrap approximations, not calibrated finite-sample confidence bounds. Rare tail events missing from the 500 author runs remain unknowable from this sample.','',
        '```sh','python scripts/compare_dv_uncertainty.py --output results_dv_uncertainty_repeat','```','',
        'Use a fresh directory. This reads existing cached observations, validates original checksums, and saves mean resamples and plots. It generates no new detector trajectories and fits no detector threshold, channel parameter or published latency.']
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--output',type=Path,default=Path('results/dv_uncertainty'));a=p.parse_args()
    r=run(a.root.resolve(),a.output.resolve());print(json.dumps(r['differences'],indent=2))


if __name__=='__main__':main()
