"""Parameter sweep of equal-energy coherent and independent-squeezing probes."""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .models import coherent_pair, squeezed_pair
from .calibration import calibrate_threshold, independent_ratio
from .simulation import simulate_stopping_times, summarize_arl, summarize_delay


def run_condition(Na, target_arl, seed, calibration_trials=30000,
                  validation_trials=50000, delay_trials=30000, tolerance=.05):
    """One prespecified setting; no threshold tuning on the validation ensemble.

    Zero squeezing uses one common threshold because the two models are
    identical, with independent validation and delay streams for both labels.
    Censored ensembles cannot support an unrestricted comparison.
    """
    horizon = max(100000, int(100*target_arl))
    models = {'coherent': coherent_pair(Na=Na), 'squeezed': squeezed_pair(Na=Na)}
    streams = np.random.SeedSequence(seed).spawn(8)
    seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in streams]
    entries, raw = {}, {}
    for i, (name, (pre, post)) in enumerate(models.items()):
        if Na == 0 and name == 'squeezed':
            baseline = models['coherent']
            np.testing.assert_allclose([pre.mean, pre.variance, post.mean, post.variance],
                [baseline[0].mean, baseline[0].variance, baseline[1].mean, baseline[1].variance],
                rtol=0, atol=1e-14)
            # Use exactly the same floating-point model and threshold at the null.
            pre, post = baseline
            calibration = entries['coherent']['calibration'] | {'reused_from':'coherent'}
        else:
            calibration = calibrate_threshold(pre, post, target_arl, trials=calibration_trials,
                                               horizon=horizon, seed=seeds[i])
        h = calibration['threshold']
        runs = simulate_stopping_times(pre, post, h, trials=validation_trials,
                                        horizon=horizon, seed=seeds[2+i])
        validation = summarize_arl(runs)
        if validation['arl'] is None:
            raise RuntimeError('censored ARL validation; enlarge the horizon')
        raw[f'{name}_validation_stop'], raw[f'{name}_validation_event'] = runs.stop, runs.event
        delay = {}
        for j, nu in enumerate([1, 100]):
            runs = simulate_stopping_times(pre, post, h, trials=delay_trials,
                horizon=horizon, change_at=nu, seed=seeds[4+2*j+i])
            delay[str(nu)] = summarize_delay(runs)
            if delay[str(nu)]['conditional_delay'] is None:
                raise RuntimeError('censored/no-survivor delay ensemble')
            raw[f'{name}_nu{nu}_stop'], raw[f'{name}_nu{nu}_event'] = runs.stop, runs.event
        entries[name] = {'calibration':calibration, 'validation':validation, 'delay':delay,
                         'pre': {'mean':float(pre.mean), 'variance':float(pre.variance)},
                         'post': {'mean':float(post.mean), 'variance':float(post.variance)},
                         'kl':post.kl(pre)}
    a, b = entries['coherent'], entries['squeezed']
    arl_ratio = independent_ratio(a['validation']['arl'], a['validation']['arl_se'],
                                  b['validation']['arl'], b['validation']['arl_se'])
    lo, hi = arl_ratio['approx_95_interval']
    matched = (1-tolerance <= lo and hi <= 1+tolerance and
               all(abs(e['validation']['arl']/target_arl-1) <= tolerance for e in entries.values()))
    ratios = {}
    for nu in ['1','100']:
        x, y = a['delay'][nu], b['delay'][nu]
        r = independent_ratio(x['conditional_delay'],x['conditional_delay_se'],
                              y['conditional_delay'],y['conditional_delay_se'])
        ratios[nu] = r | {'delay_reduction_percent':100*(1-1/r['ratio'])}
    report = {'Na':Na, 'target_arl':target_arl, 'seed':seed, 'horizon':horizon,
              'stream_seeds':seeds, 'models':entries, 'arl_ratio':arl_ratio,
              'arl_match_tolerance':tolerance, 'arl_matched':matched, 'delay_ratios':ratios,
              'kl_ratio':b['kl']/a['kl']}
    return report, raw


def summary_rows(report):
    rows = []
    for nu in ['1','100']:
        a, b = report['models']['coherent'], report['models']['squeezed']
        x, y = a['delay'][nu], b['delay'][nu]
        r = report['delay_ratios'][nu]
        rows.append({'Na':report['Na'], 'target_arl':report['target_arl'], 'change_at':int(nu),
            'arl_matched':report['arl_matched'],
            'coherent_threshold':a['calibration']['threshold'], 'squeezed_threshold':b['calibration']['threshold'],
            'coherent_arl':a['validation']['arl'], 'coherent_arl_se':a['validation']['arl_se'],
            'squeezed_arl':b['validation']['arl'], 'squeezed_arl_se':b['validation']['arl_se'],
            'arl_ratio':report['arl_ratio']['ratio'], 'arl_ratio_se':report['arl_ratio']['se'],
            'coherent_delay':x['conditional_delay'], 'coherent_delay_se':x['conditional_delay_se'],
            'squeezed_delay':y['conditional_delay'], 'squeezed_delay_se':y['conditional_delay_se'],
            'delay_ratio':r['ratio'], 'delay_ratio_se':r['se'],
            'ratio_95_low':r['approx_95_interval'][0], 'ratio_95_high':r['approx_95_interval'][1],
            'delay_reduction_percent':r['delay_reduction_percent'], 'kl_ratio':report['kl_ratio'],
            'coherent_false_alarm_probability':x['false_alarm_probability'],
            'squeezed_false_alarm_probability':y['false_alarm_probability'],
            'coherent_survivors':x['survivors'], 'squeezed_survivors':y['survivors']})
    return rows


def write_summary(root, reports, settings):
    rows = [row for report in reports for row in summary_rows(report)]
    with (root/'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    overall = {'settings':settings, 'paper':'https://arxiv.org/html/2503.12276v2',
        'parameters':{'N':100, 'eta1':.9, 'eta2':.85, 'change_at':[1,100]},
        'python':sys.version, 'platform':platform.platform(),
        'versions':{k:importlib.metadata.version(k) for k in ['numpy','matplotlib']},
        'conditions':reports}
    (root/'report.json').write_text(json.dumps(overall, indent=2)+'\n')
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    nas = settings['Na']
    colors = ['#52677c','#c5413c','#258875','#8658a6']
    for k, target in enumerate(settings['target_arl']):
        data = [row for row in rows if row['target_arl']==target and row['change_at']==100]
        for i, row in enumerate(data):
            marker = 'o' if row['arl_matched'] else 'x'
            axes[0].errorbar(i, row['delay_ratio'], yerr=1.96*row['delay_ratio_se'],
                marker=marker, color=colors[k%len(colors)], capsize=3,
                label=f'ARL target {target:g}' if i==0 else None)
        axes[0].plot(range(len(data)),[row['delay_ratio'] for row in data], color=colors[k%len(colors)], alpha=.5)
        axes[1].errorbar(np.arange(len(data))+(k-(len(settings['target_arl'])-1)/2)*.12,
            [row['arl_ratio'] for row in data], yerr=[1.96*row['arl_ratio_se'] for row in data],
            fmt='o', capsize=3, color=colors[k%len(colors)])
    reference = [r for r in reports if r['target_arl']==settings['target_arl'][0]]
    axes[0].plot(range(len(nas)), [r['kl_ratio'] for r in reference], 'k--', label='KL ratio: large-ARL reference')
    axes[0].axhline(1, color='grey', ls=':')
    axes[0].set(ylabel='Ordinary delay / squeezed delay', title='Conditional detection gain: change at pulse 100')
    axes[0].legend(fontsize=9)
    tol=settings['match_tolerance']
    axes[1].axhspan(1-tol,1+tol,color='#258875',alpha=.10,label=f'±{tol*100:g}% equivalence band')
    axes[1].axhline(1, color='black', ls='--', lw=1)
    axes[1].set(ylabel='Ordinary ARL / squeezed ARL',title='Independent false-alarm validation')
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.set(xticks=range(len(nas)),xticklabels=[f'{x:g}' for x in nas], xlabel='Tested Na (photons/pulse; categorical spacing)')
        ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
    fig.suptitle('Equal energy per setting · pointwise approximate 95% Monte Carlo intervals')
    fig.tight_layout(); fig.savefig(root/'parameter_sweep.png',dpi=180)
    fig.savefig(root/'parameter_sweep.svg'); plt.close(fig)
    matched=sum(r['arl_matched'] for r in reports)
    lines=['# Parameter sweep: independent squeezing versus coherent light','',
        f"Grid: Na={nas}; target ARL={settings['target_arl']}. Fixed N=100, eta1=0.9, eta2=0.85. Same total energy N+Na within each comparison.",'',
        f"**{matched}/{len(reports)} settings passed the prespecified approximate ARL equivalence criterion:** the entire pointwise 95% ARL-ratio interval lies within ±{100*tol:g}%, and both validated ARL estimates lie within that tolerance of target. Unmatched settings, if any, remain visible and are exploratory only.",'',
        '| Na | Target ARL | Validated ordinary / squeezed ARL | Delay ordinary / squeezed (nu=100) | Delay gain ratio | Approx. 95% ratio interval | Reduction | Matched |',
        '|---|---:|---:|---:|---:|---:|---:|---|']
    for row in rows:
        if row['change_at']!=100:
            continue
        lines.append(f"| {row['Na']:g} | {row['target_arl']:g} | {row['coherent_arl']:.1f} / {row['squeezed_arl']:.1f} | {row['coherent_delay']:.2f} / {row['squeezed_delay']:.2f} | {row['delay_ratio']:.3f} | [{row['ratio_95_low']:.3f}, {row['ratio_95_high']:.3f}] | {row['delay_reduction_percent']:.1f}% | {row['arl_matched']} |")
    lines += ['', '## Control and interpretation','',
        'At Na=0, both physical models and thresholds are identical. Only the observations are drawn independently. The delay ratio should be near 1; the 95% intervals need not all contain 1 by chance. This control checks that the pipeline does not manufacture a squeezing advantage.', '',
        'Na=0.01, 0.1 and 1 test progressively stronger squeezing. Each condition includes changes at pulse 1 and pulse 100. The CSV and JSON report conditional delays, survivor counts and early-alarm probabilities for both. nu=1 has no pre-change alarms; nu=100 conditions on no earlier alarm. This is not a worst-case delay calculation.', '',
        'The plot uses equally spaced categories for the tested Na values. The dashed KL ratio is an asymptotic reference, not an exact finite-ARL delay prediction. A growing delay gain as ARL increases is suggestive of that limit, not a proof of convergence. Exact equality of true ARLs is not established by finite Monte Carlo samples.', '',
        f"Calibration uses {settings['calibration_trials']:,} trials/evaluation, separate validation {settings['validation_trials']:,} per method/setting, and {settings['delay_trials']:,} delay trials per method/change point. The broader sweep uses a 5% equivalence band by default; the earlier high-precision single-setting comparison retained its separate 2% band.", '',
        'All intervals are pointwise approximate Monte Carlo intervals; there is no multiple-comparison correction and no simultaneous confidence claim for the whole grid. Threshold calibration uncertainty is not included in the fixed-threshold delay intervals. Every condition was prespecified; no failing validation was discarded or used to retune its thresholds.', '',
        '## Saved evidence and scope','',
        'Each condition folder contains report.json with model parameters, stream seeds and calibration history, plus raw_runs.npz containing all validation/delay stopping times and event flags. The aggregate report.json, summary.csv and figures use those records. A fresh output folder is required. No relevant calibration/validation/delay ensemble was censored; the runner refuses unrestricted comparisons otherwise.', '',
        'This extends the moderate-ARL numerical validation. It does not reproduce Figure 3 at ARL=2,000,000, implement entangled blocks, or validate a real optical experiment. Next: efficient large-ARL calibration and the correlated entangled model.', '']
    (root/'RESULTS.md').write_text('\n'.join(lines))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('results_sweep'))
    p.add_argument('--seed',type=int,default=20261023)
    p.add_argument('--Na',type=float,nargs='+',default=[0.,.01,.1,1.])
    p.add_argument('--target-arl',type=float,nargs='+',default=[500.,1000.,5000.])
    p.add_argument('--calibration-trials',type=int,default=30000)
    p.add_argument('--validation-trials',type=int,default=50000)
    p.add_argument('--delay-trials',type=int,default=30000)
    p.add_argument('--match-tolerance',type=float,default=.05)
    args=p.parse_args(argv)
    if any(not np.isfinite(x) or x<0 for x in args.Na) or len(set(args.Na))!=len(args.Na):
        p.error('Na values must be finite, nonnegative and unique')
    if any(not np.isfinite(x) or x<=1 for x in args.target_arl) or len(set(args.target_arl))!=len(args.target_arl):
        p.error('target ARLs must be finite, >1 and unique')
    if min(args.calibration_trials,args.validation_trials,args.delay_trials)<2:
        p.error('ensemble sizes must be >=2')
    if not np.isfinite(args.match_tolerance) or not 0<args.match_tolerance<1:
        p.error('match tolerance must lie in (0,1)')
    args.output.mkdir(parents=True,exist_ok=False)
    settings=vars(args)|{'output':str(args.output)}
    (args.output/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    grid=[(a,g) for g in args.target_arl for a in args.Na]
    streams=np.random.SeedSequence(args.seed).spawn(len(grid))
    reports=[]
    for i,((a,g),stream) in enumerate(zip(grid,streams)):
        print(f'Setting {i+1}/{len(grid)}: Na={a:g}, target ARL={g:g}',flush=True)
        report,raw=run_condition(a,g,int(stream.generate_state(1,dtype=np.uint64)[0]),
            args.calibration_trials,args.validation_trials,args.delay_trials,args.match_tolerance)
        folder=args.output/f'condition_{i:02d}'
        folder.mkdir()
        (folder/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        np.savez_compressed(folder/'raw_runs.npz',**raw)
        reports.append(report)
        print(f"  ARL matched: {report['arl_matched']}; delay ratio at nu=100: {report['delay_ratios']['100']['ratio']:.3f}",flush=True)
    write_summary(args.output,reports,settings)
    print(f"Completed {len(reports)} settings; {sum(r['arl_matched'] for r in reports)} passed ARL equivalence.",flush=True)


if __name__=='__main__':
    main()
