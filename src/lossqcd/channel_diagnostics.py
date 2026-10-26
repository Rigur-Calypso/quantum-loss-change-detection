"""One-point channel-mismatch fit and held-out scalar predictions.

This is an inverse-model diagnostic, not a reproduction or an estimate of
the authors' actual channel. Only the true post-change transmissivity varies.
The assumed likelihoods, thresholds and true pre-change law stay fixed, so
the original no-change ARL evidence applies to exactly the same stopping law.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from .delay_solver import conditional_delay_numerical
from .models import coherent_pair, squeezed_pair
from .simulation import simulate_stopping_times, summarize_delay


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def laws(kind, N, Na, eta1, assumed_eta2, actual_eta2):
    if kind not in ('coherent', 'entangled'):
        raise ValueError('expected coherent or entangled (scalar squeezed)')
    pair = coherent_pair if kind == 'coherent' else squeezed_pair
    pre, post = pair(N=N, Na=Na, eta1=eta1, eta2=assumed_eta2)
    _, actual_post = pair(N=N, Na=Na, eta1=eta1, eta2=actual_eta2)
    return pre, post, actual_post


def fit_actual_eta2(pre, post, threshold, target_delay, post_factory,
                    bracket=(.80, .85), change_at=1005, cells=256):
    """Fit one physical parameter, holding the detector likelihoods fixed."""
    if not np.isfinite(target_delay) or target_delay <= 0:
        raise ValueError('target delay must be finite and positive')
    if len(bracket) != 2 or not 0 <= bracket[0] < bracket[1] <= 1:
        raise ValueError('bracket must be increasing within [0,1]')
    history = []
    def objective(eta):
        delay = conditional_delay_numerical(pre, post, threshold,
            change_at=change_at, cells=cells,
            sampling_post=post_factory(eta))['conditional_delay']
        history.append({'actual_eta2': float(eta), 'delay': delay})
        return delay-target_delay
    eta = brentq(objective, *bracket, xtol=1e-10)
    return {'actual_eta2': float(eta), 'target_delay': float(target_delay),
            'cells': cells, 'bracket': list(bracket), 'history': history}


def scalar_models(source):
    return [m for m in source['models'] if m['figure'] == 7 and m['n'] == 1]


def numerical_row(model, settings, actual_eta2):
    pre, post, actual = laws(model['kind'], settings['N'], model['Na'],
        settings['eta1'], settings['eta2'], actual_eta2)
    args = dict(change_at=model['change_at'])
    h = model['calibration']['threshold']
    coarse = conditional_delay_numerical(pre, post, h, cells=512,
        sampling_post=actual, **args)
    fine = conditional_delay_numerical(pre, post, h, cells=1024,
        sampling_post=actual, **args)
    matched = conditional_delay_numerical(pre, post, h, cells=512, **args)
    return pre, post, actual, coarse, fine, matched


def run(source_root, output, trials=30000, seed=2026102601):
    source_root, output = Path(source_root).resolve(), Path(output).resolve()
    if trials < 2:
        raise ValueError('at least two trials required')
    source_path = source_root/'report.json'
    source = json.loads(source_path.read_text())
    models = scalar_models(source)
    if len(models) != 18 or not all(m['confirmed'] for m in models):
        raise ValueError('expected 18 confirmed scalar Figure 7 source settings')
    settings = source['settings']
    training = next(m for m in models if m['kind'] == 'coherent' and np.isclose(m['Na'], .1))
    author = next(r for r in source['figure7_rows'] if r['label'] == training['label'] and np.isclose(r['Na'], .1))
    pre, post, _ = laws('coherent', settings['N'], .1, settings['eta1'], settings['eta2'], settings['eta2'])
    factory = lambda eta: coherent_pair(N=settings['N'], Na=.1,
        eta1=settings['eta1'], eta2=eta)[1]
    fit = fit_actual_eta2(pre, post, training['calibration']['threshold'],
        author['author_delay'], factory, change_at=training['change_at'])
    output.mkdir(parents=True, exist_ok=False)
    hashes = {str(source_path): sha256(source_path)}
    # Hash the original evidence and author values without altering any source.
    for m in models:
        path = source_root/Path(m['folder']).name/'changed_runs.npz'
        hashes[str(path)] = sha256(path)
    data = Path(settings['data']).resolve()
    for name in ('prl_fig7_data', 'plot_PRL_figs.m', 'manifest.json'):
        path = data/name
        hashes[str(path)] = sha256(path)
    report = {'source_root': str(source_root), 'source_hashes': hashes,
        'settings': {'trials': trials, 'seed': seed, 'horizon': settings['horizon'],
            'N': settings['N'], 'eta1': settings['eta1'], 'assumed_eta2': settings['eta2']},
        'fit': fit, 'training': {'label': training['label'], 'Na': .1}, 'rows': [],
        'scope': 'Scalar coherent and independent squeezed only; correlated blocks not tested.',
        'interpretation': 'Diagnostic one-point fit, not inferred author parameters. Author errors unknown. '
            'Figure 7 eta is not saved in the dataset. Original ARL evidence is reused, not newly estimated: '
            'the true pre-change law, score models and threshold are identical.'}
    streams = np.random.SeedSequence(seed).spawn(len(models))
    for index, (model, stream) in enumerate(zip(models, streams)):
        pre, post, actual, coarse, fine, matched = numerical_row(model, settings, fit['actual_eta2'])
        run_seed = int(stream.generate_state(1, dtype=np.uint64)[0])
        runs = simulate_stopping_times(pre, post, model['calibration']['threshold'],
            trials=trials, horizon=settings['horizon'], change_at=model['change_at'],
            seed=run_seed, sampling_post=actual)
        stats = summarize_delay(runs)
        filename = f'runs_{index:02d}.npz'
        np.savez_compressed(output/filename, stop=runs.stop, event=runs.event)
        published = next(r for r in source['figure7_rows'] if r['label'] == model['label'] and np.isclose(r['Na'], model['Na']))
        delta = abs(fine['conditional_delay']-coarse['conditional_delay'])
        row = {'Na': model['Na'], 'label': model['label'], 'kind': model['kind'],
            'training': model is training, 'threshold': model['calibration']['threshold'],
            'change_at': model['change_at'], 'seed': run_seed, 'runs': filename,
            'author_delay': published['author_delay'], 'original_mc_delay': published['our_delay'],
            'matched_numerical': matched, 'numerical_512': coarse, 'numerical_1024': fine,
            'mesh_difference': delta, 'mc': stats,
            'relative_difference': fine['conditional_delay']/published['author_delay']-1,
            'within_5_percent': abs(fine['conditional_delay']/published['author_delay']-1) <= .05,
            'source_arl': model['validation']['arl'], 'source_arl_se': model['validation']['arl_se'],
            'checks_passed': stats['conditional_delay'] is not None and delta < .05
                and abs(fine['conditional_delay']-stats['conditional_delay']) <= 4*stats['conditional_delay_se']+delta
                and matched['false_alarm_probability'] == coarse['false_alarm_probability']}
        report['rows'].append(row)
        print(f"{index+1}/18 {model['label']} Na={model['Na']:g}: predicted {fine['conditional_delay']:.3f}, author {published['author_delay']:.3f}", flush=True)
    report['source_unchanged'] = all(sha256(p) == h for p, h in hashes.items())
    report['all_checks_passed'] = report['source_unchanged'] and all(r['checks_passed'] for r in report['rows'])
    report['heldout_within_5_percent'] = sum(r['within_5_percent'] for r in report['rows'] if not r['training'])
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    write_outputs(report, output)
    return report


def write_outputs(report, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = report['rows']
    with (output/'predictions.csv').open('w', newline='') as handle:
        fields = ['Na', 'label', 'training', 'threshold', 'author_delay', 'original_mc_delay',
                  'prediction', 'mc_delay', 'mc_se', 'relative_difference', 'within_5_percent']
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for r in rows:
            writer.writerow({k: r[k] for k in fields if k in r} |
                {'prediction': r['numerical_1024']['conditional_delay'],
                 'mc_delay': r['mc']['conditional_delay'], 'mc_se': r['mc']['conditional_delay_se']})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained')
    for ax, kind, title in zip(axes, ('coherent', 'entangled'), ('Coherent light', 'Independent squeezed light')):
        subset = [r for r in rows if r['kind'] == kind]
        x = [r['Na'] for r in subset]
        ax.plot(x, [r['author_delay'] for r in subset], 'o-', label='Published means')
        ax.plot(x, [r['original_mc_delay'] for r in subset], ':', label='Original matched simulation')
        ax.plot(x, [r['numerical_1024']['conditional_delay'] for r in subset], '-', label='Common loss-change prediction')
        ax.errorbar(x, [r['mc']['conditional_delay'] for r in subset],
            yerr=[1.96*r['mc']['conditional_delay_se'] for r in subset], fmt='x', capsize=3, label='Fresh MC ±1.96 SE')
        ax.set(xlabel='Added energy Na', ylabel='Conditional delay (pulses)', title=title)
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle(f"Diagnostic: true eta2={report['fit']['actual_eta2']:.6f}; detector assumes 0.85\nFit uses only coherent Na=0.1; all thresholds and pre-change laws fixed", fontsize=11)
    for extension in ('png', 'svg'):
        fig.savefig(output/f'channel_predictions.{extension}', dpi=180)
    plt.close(fig)
    (output/'RESULTS.md').write_text(f"""# Channel-parameter diagnostic

A single true post-change transmissivity was fitted to the published coherent
mean at Na=0.1: **eta2={report['fit']['actual_eta2']:.8f}**, instead of the
detector's assumed 0.85. Pre-change eta1 stays 0.9, N stays 100. Every detector
keeps its original likelihoods and ARL-calibrated threshold. This changes the
physical post-change observations, including squeezed-state variance.

The fitted point is training data; **{report['heldout_within_5_percent']}/17 held-out scalar means**
are predicted within ±5%. This descriptive tolerance is not a statistical test:
author uncertainty estimates are unavailable. Correlated entangled blocks are
outside this experiment's scope. The fit does not establish the authors' loss.
Figure 7 source data do not separately save eta; Figure 3 explicitly specifies
0.9/0.85 and is not refitted here.

All 18 predictions use 512/1024-cell deterministic calculations and fresh
{report['settings']['trials']:,}-trial simulations, with independent seeds per setting.
All numerical/MC and mesh checks pass: **{report['all_checks_passed']}**.
Early first alarms are recorded separately; delay is E[T−1005 | T>=1005].
Censored delays are never reported as unrestricted means.

The no-change law is mathematically identical to the source study because
only the true post-change law changes. Existing ARL evidence (near two million
pulses) is reused unchanged; no new ARL estimate or calibration is claimed.
Original report, raw scalar runs and author files are checksum-protected.
See predictions.csv, report.json, raw runs and channel_predictions.png.

The common change reduces the gap but does not bring all held-out scalar means
within the chosen tolerance. Next: extend the physical-mismatch diagnostic to
correlated blocks, preserving their exact mixed-loss covariance and prefix
likelihoods. Requesting the authors' original
thresholds and simulator would resolve assumptions more directly; no message
has been sent.
""")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default='results_author_comparison')
    parser.add_argument('--output', default='results_channel_diagnostics')
    parser.add_argument('--trials', type=int, default=30000)
    parser.add_argument('--seed', type=int, default=2026102601)
    args = parser.parse_args()
    run(args.source, args.output, args.trials, args.seed)


if __name__ == '__main__':
    main()
