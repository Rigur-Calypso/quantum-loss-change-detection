"""Equal-energy, approximately equal-ARL coherent versus squeezed comparison."""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import sys
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .models import coherent_pair, squeezed_pair
from .calibration import calibrate_threshold, independent_ratio
from .simulation import simulate_stopping_times, summarize_arl, summarize_delay


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=Path('results_equal_arl'))
    p.add_argument('--seed', type=int, default=20261022)
    p.add_argument('--target-arl', type=float, default=1000.)
    p.add_argument('--calibration-trials', type=int, default=100000)
    p.add_argument('--validation-trials', type=int, default=100000)
    p.add_argument('--delay-trials', type=int, default=50000)
    p.add_argument('--horizon', type=int, default=100000)
    p.add_argument('--change-at', type=int, nargs='+', default=[1, 100, 1000])
    p.add_argument('--Na', type=float, default=.1)
    p.add_argument('--match-tolerance', type=float, default=.02)
    args = p.parse_args(argv)
    if min(args.calibration_trials, args.validation_trials, args.delay_trials) < 2:
        p.error('each ensemble requires at least two trials')
    if any(not 1 <= nu <= args.horizon for nu in args.change_at):
        p.error('change points must lie in the observation window')
    if len(set(args.change_at)) != len(args.change_at):
        p.error('change points must be unique')
    if not np.isfinite(args.match_tolerance) or not 0 < args.match_tolerance < 1:
        p.error('match tolerance must be between 0 and 1')
    # A fresh output directory prevents accidental replacement of previous studies.
    args.output.mkdir(parents=True, exist_ok=False)
    streams = np.random.SeedSequence(args.seed).spawn(4+2*len(args.change_at))
    stream_seeds = [int(s.generate_state(1, dtype=np.uint64)[0]) for s in streams]
    models = {'coherent': coherent_pair(Na=args.Na), 'squeezed': squeezed_pair(Na=args.Na)}
    report = {'paper': 'https://arxiv.org/html/2503.12276v2',
              'parameters': {'N':100., 'Na':args.Na, 'eta1':.9, 'eta2':.85,
                             'total_transmitted_energy_per_pulse':100+args.Na},
              'settings': vars(args) | {'output': str(args.output)},
              'stream_seeds': stream_seeds, 'rng': 'default_rng/PCG64; SeedSequence.spawn',
              'python': sys.version, 'platform': platform.platform(),
              'versions': {k:importlib.metadata.version(k) for k in ['numpy','matplotlib']},
              'models': {}}
    raw, delay_rows = {}, []
    for i, (name, (pre, post)) in enumerate(models.items()):
        print(f'Calibrating {name} to ARL={args.target_arl:g} ...', flush=True)
        calibration = calibrate_threshold(pre, post, args.target_arl,
            trials=args.calibration_trials, horizon=args.horizon, seed=stream_seeds[i])
        h = calibration['threshold']
        print(f'{name}: fixed threshold {h:.6f}; independent validation ...', flush=True)
        validation = simulate_stopping_times(pre, post, h, args.validation_trials,
                                            args.horizon, seed=stream_seeds[2+i])
        stats = summarize_arl(validation)
        if stats['arl'] is None:
            raise RuntimeError('censored validation: increase horizon; no equal-ARL claim can be made')
        raw[f'{name}_validation_stop'] = validation.stop
        raw[f'{name}_validation_event'] = validation.event
        entry = {'pre': {'mean':float(pre.mean), 'variance':float(pre.variance)},
                 'post': {'mean':float(post.mean), 'variance':float(post.variance)},
                 'kl_nats_per_pulse': post.kl(pre), 'calibration':calibration,
                 'validation':stats, 'delay':{}}
        for j, nu in enumerate(args.change_at):
            runs = simulate_stopping_times(pre, post, h, args.delay_trials, args.horizon,
                                           change_at=nu, seed=stream_seeds[4+2*j+i])
            delay = summarize_delay(runs)
            if delay['conditional_delay'] is None:
                raise RuntimeError('censored/no-survivor delay ensemble: unrestricted comparison unavailable')
            entry['delay'][str(nu)] = delay
            delay_rows.append({'model':name, 'threshold':h, **delay})
            raw[f'{name}_nu{nu}_stop'] = runs.stop
            raw[f'{name}_nu{nu}_event'] = runs.event
        report['models'][name] = entry
        # Persist calibration history as soon as a model completes.
        (args.output/f'{name}_calibration.json').write_text(json.dumps(calibration, indent=2)+'\n')
        print(f"{name}: validated ARL={stats['arl']:.2f} +/- {stats['arl_se']:.2f} MC SE", flush=True)

    classical, squeezed = report['models']['coherent'], report['models']['squeezed']
    a, b = classical['validation'], squeezed['validation']
    arl_ratio = independent_ratio(a['arl'], a['arl_se'], b['arl'], b['arl_se'])
    tolerance = args.match_tolerance
    # Require the entire approximate 95% ratio interval inside the equivalence
    # band, plus each point estimate close to target. Do not retune on holdout.
    interval = arl_ratio['approx_95_interval']
    matched = (1-tolerance <= interval[0] and interval[1] <= 1+tolerance
               and abs(a['arl']/args.target_arl-1) <= tolerance
               and abs(b['arl']/args.target_arl-1) <= tolerance)
    comparisons = {}
    for nu in args.change_at:
        x, y = classical['delay'][str(nu)], squeezed['delay'][str(nu)]
        ratio = independent_ratio(x['conditional_delay'], x['conditional_delay_se'],
                                  y['conditional_delay'], y['conditional_delay_se'])
        comparisons[str(nu)] = ratio | {'delay_reduction_percent':100*(1-1/ratio['ratio'])}
    report['arl_comparison'] = arl_ratio
    report['arl_equivalence_tolerance'] = tolerance
    report['approximately_equal_arl_validated'] = matched
    report['delay_comparisons'] = comparisons
    report['kl_ratio_asymptotic_reference_only'] = squeezed['kl_nats_per_pulse']/classical['kl_nats_per_pulse']
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(args.output/'raw_runs.npz', **raw)
    with (args.output/'delay_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(delay_rows[0]))
        writer.writeheader(); writer.writerows(delay_rows)

    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5))
    colors = ['#52677c', '#c5413c']
    axs[0].bar(['Ordinary light','Squeezed light'], [a['arl'],b['arl']],
               yerr=[1.96*a['arl_se'],1.96*b['arl_se']], color=colors, capsize=5)
    axs[0].axhline(args.target_arl, color='black', ls='--', lw=1, label='Target')
    axs[0].set(ylabel='Pulses before a false alarm', title='Independent ARL validation')
    axs[0].legend()
    index = np.arange(len(args.change_at))
    for i, (name, label) in enumerate([('coherent','Ordinary light'),('squeezed','Squeezed light')]):
        data = report['models'][name]['delay']
        axs[1].bar(index+(i-.5)*.36, [data[str(nu)]['conditional_delay'] for nu in args.change_at],
                   width=.36, color=colors[i], label=label,
                   yerr=[1.96*data[str(nu)]['conditional_delay_se'] for nu in args.change_at], capsize=3)
    axs[1].set(xticks=index, xticklabels=args.change_at, xlabel='Change at pulse',
               ylabel='Conditional detection delay (pulses)', title='Delay given no earlier alarm')
    axs[1].legend()
    axs[2].errorbar(index, [comparisons[str(nu)]['ratio'] for nu in args.change_at],
                    yerr=[1.96*comparisons[str(nu)]['se'] for nu in args.change_at],
                    fmt='o', color=colors[1], capsize=5, label='Simulated delay ratio')
    axs[2].axhline(1, color='grey', ls=':')
    axs[2].axhline(report['kl_ratio_asymptotic_reference_only'], color='black', ls='--', label='KL ratio (asymptotic reference)')
    axs[2].set(xticks=index, xticklabels=args.change_at, xlabel='Change at pulse',
               ylabel='Ordinary delay / squeezed delay', title='Gain from independent squeezing')
    axs[2].legend(fontsize=8)
    for ax in axs:
        ax.grid(axis='y', alpha=.2); ax.set_axisbelow(True)
    fig.suptitle(f'Equal energy: N=100, Na={args.Na:g}; target ARL={args.target_arl:g}; bars: approx. 95% MC intervals')
    fig.tight_layout(); fig.savefig(args.output/'equal_arl_comparison.png', dpi=180)
    fig.savefig(args.output/'equal_arl_comparison.svg'); plt.close(fig)
    lines = ['# Coherent versus independently squeezed detection', '',
             'Same total energy per pulse; separately calibrated thresholds; independent ARL validation and delay ensembles.', '',
             f"Parameters: N=100, Na={args.Na:g}, eta1=0.9, eta2=0.85; seed {args.seed}. Target ARL {args.target_arl:g} pulses.", '',
             '| Method | Threshold | Validated ARL ± MC SE |', '|---|---:|---:|']
    for name, label in [('coherent','Ordinary light'),('squeezed','Independent squeezing')]:
        e = report['models'][name]
        lines.append(f"| {label} | {e['calibration']['threshold']:.6f} | {e['validation']['arl']:.2f} ± {e['validation']['arl_se']:.2f} |")
    lines += ['', f"Approximate ARL ratio 95% interval: [{interval[0]:.4f}, {interval[1]:.4f}]. Entire interval inside ±{100*tolerance:g}% equivalence band: **{matched}**. Each validated ARL point estimate is also required within that tolerance of target. Validation was not reused to retune thresholds.", '',
              '| Change at pulse | Ordinary delay ± SE | Squeezed delay ± SE | Delay ratio (approx. 95% interval) | Delay reduction | Early alarms: ordinary / squeezed |',
              '|---|---:|---:|---:|---:|---:|']
    for nu in args.change_at:
        x, y = classical['delay'][str(nu)], squeezed['delay'][str(nu)]
        c = comparisons[str(nu)]; lo, hi = c['approx_95_interval']
        lines.append(f"| {nu} | {x['conditional_delay']:.3f} ± {x['conditional_delay_se']:.3f} | {y['conditional_delay']:.3f} ± {y['conditional_delay_se']:.3f} | {c['ratio']:.3f} [{lo:.3f}, {hi:.3f}] | {c['delay_reduction_percent']:.1f}% | {100*x['false_alarm_probability']:.2f}% / {100*y['false_alarm_probability']:.2f}% |")
    lines += ['', '## Interpretation and limits', '',
              'Delay means E[T-nu | T>=nu], with zero delay for an alarm on the first changed pulse. Early alarms are excluded from this conditional mean and reported separately; this is not a worst-case delay. nu=1 is a fresh detector with no pre-change false alarms. Longer pre-change periods can produce substantial early-alarm rates even at the same ARL.', '',
              f"There were no censored runs in calibration, validation or the delay ensembles (the script refuses unrestricted comparisons if any are censored). Calibration uses {args.calibration_trials:,} trials/evaluation; validation {args.validation_trials:,} per method; delays {args.delay_trials:,} per method/change point. Log-likelihoods include unequal variances for squeezed measurements.", '',
              f"The analytical KL ratio is {report['kl_ratio_asymptotic_reference_only']:.6f}. It is a large-ARL reference, not an exact finite-threshold delay prediction. Uncertainty intervals describe Monte Carlo sampling at the fixed chosen thresholds; they do not include a hypothetical recalibration's variability or optical modeling errors. Exact equality of true ARLs is not established by finite simulations.", '',
              'This is a moderate-ARL, single-Na computational comparison. It does not reproduce Figure 3 at ARL=2,000,000, implement entangled probes, or constitute a hardware experiment. Next: validate a range of Na and ARL settings, then implement efficient large-ARL calibration.', '',
              '## Reproducibility', '', 'The JSON report records settings, stream seeds, calibration histories, distribution parameters, package versions and all summaries. Raw per-run stopping times and event flags are in raw_runs.npz. Streams 0/1 calibrate coherent/squeezed, 2/3 validate, and 4+2j/5+2j simulate the jth change point. All ensembles are independent. Each calibration stream spawns independent streams for its threshold evaluations.', '']
    (args.output/'RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps({'equal_arl_validated':matched, 'arl_ratio':arl_ratio, 'delay_comparisons':comparisons}, indent=2), flush=True)
    if not matched:
        print('ARL equivalence criterion failed: comparison is exploratory, not a validated equal-ARL result.', flush=True)


if __name__ == '__main__':
    main()
