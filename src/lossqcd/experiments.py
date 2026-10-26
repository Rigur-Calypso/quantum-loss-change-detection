"""Reproduce Figure 2 classical/squeezed CRE and an initial classical MC study."""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .models import coherent_pair, coherent_kl, squeezed_kl
from .cusum import page_cusum
from .simulation import simulate_stopping_times, summarize_arl, summarize_delay


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--seed', type=int, default=20261020)
    parser.add_argument('--trials', type=int, default=4000)
    parser.add_argument('--horizon', type=int, default=100000)
    parser.add_argument('--change-at', type=int, default=100)
    args = parser.parse_args(argv)
    if args.trials < 2 or not 1 <= args.change_at <= args.horizon:
        parser.error('need trials>=2 and 1<=change-at<=horizon')
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    N, eta1, eta2 = 100., .9, .85
    # Fine near Na=0, where Eq. (3) has an infinite right derivative.
    Na = np.unique(np.r_[0, np.geomspace(1e-6, .1, 100), np.linspace(.1, 5, 491)])
    classical = np.array([coherent_kl(N, a, eta1, eta2) for a in Na])
    squeezed = np.array([squeezed_kl(N, a, eta1, eta2) for a in Na])
    np.savetxt(args.output/'figure2_data.csv', np.column_stack([Na, classical, squeezed]),
               delimiter=',', header='Na,classical_Eq2_nats,squeezed_Eq3_nats', comments='')
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Na, classical, 'r--', label='Classical coherent baseline · Eq. (2)')
    ax.plot(Na, squeezed, 'r-', label='Independent squeezing · Eq. (3)')
    ax.set(xlabel=r'Augmentation energy $N_a$ (photons/pulse)',
           ylabel=r'$D(P_2\Vert P_1)$ (nats/pulse)', xlim=(0, 5), ylim=(0, None),
           title=r'Figure 2: classical and independent squeezing curves')
    ax.text(.02, .96, r'$N=100,\ \eta_1=0.9,\ \eta_2=0.85$', transform=ax.transAxes, va='top')
    ax.legend(loc='center right'); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(args.output/'figure2_classical_squeezed.png', dpi=180)
    fig.savefig(args.output/'figure2_classical_squeezed.svg'); plt.close(fig)

    pre, post = coherent_pair(N=N, Na=0, eta1=eta1, eta2=eta2)
    thresholds = [2., 3., 4.]
    children = np.random.SeedSequence(args.seed).spawn(7)
    rows, raw = [], {}
    for i, h in enumerate(thresholds):
        arl = simulate_stopping_times(pre, post, h, args.trials, args.horizon,
                                     seed=children[2*i])
        delay = simulate_stopping_times(pre, post, h, args.trials, args.horizon,
                                       change_at=args.change_at, seed=children[2*i+1])
        row = {'threshold': h, **summarize_arl(arl), **summarize_delay(delay)}
        rows.append(row)
        for label, runs in [('arl', arl), ('delay', delay)]:
            raw[f'{label}_{i}_stop'] = runs.stop
            raw[f'{label}_{i}_event'] = runs.event
        print(json.dumps(row), flush=True)
    np.savez_compressed(args.output/'monte_carlo_runs.npz', **raw)
    with (args.output/'monte_carlo_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    specs = [('arl', 'arl_se', 'No-change average run length'),
             ('conditional_delay', 'conditional_delay_se', r'Conditional delay $E[T-\nu\mid T\geq\nu]$'),
             ('false_alarm_probability', 'false_alarm_probability_se', r'False alarms $P(T<\nu)$')]
    for ax, (key, se, label) in zip(axs, specs):
        values = [np.nan if row[key] is None else row[key] for row in rows]
        errors = [0 if row[se] is None else row[se] for row in rows]
        ax.errorbar(thresholds, values, yerr=errors, fmt='o-', capsize=4)
        ax.set(xlabel='CUSUM threshold h', ylabel=label); ax.grid(alpha=.2)
    fig.suptitle(f'Classical baseline · {args.trials:,} runs/setting · change at {args.change_at} · bars ±1 MC SE')
    fig.tight_layout(); fig.savefig(args.output/'classical_monte_carlo.png', dpi=180); plt.close(fig)

    rng = np.random.default_rng(children[6]); nu = args.change_at
    length = max(nu+300, 500)
    x = np.r_[rng.normal(pre.mean, np.sqrt(pre.variance), nu-1),
              rng.normal(post.mean, np.sqrt(post.variance), length-nu+1)]
    w = page_cusum(post.llr(x, pre))
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(np.arange(1, length+1), w, lw=1, label='Page statistic (continued for illustration)')
    ax.axvline(nu, color='k', ls=':', label='True change')
    ax.axhline(4, color='r', ls='--', label='h=4; detector stops at first crossing')
    ax.set(xlabel='Pulse index (one based)', ylabel='Page CUSUM statistic', title='One seeded classical trajectory')
    ax.legend(); fig.tight_layout(); fig.savefig(args.output/'classical_trace.png', dpi=180); plt.close(fig)
    report = {'paper': 'https://arxiv.org/html/2503.12276v2',
              'parameters': {'N':N, 'eta1':eta1, 'eta2':eta2, 'simulation_Na':0},
              'seed': args.seed, 'rng': 'NumPy default_rng/PCG64; SeedSequence.spawn(7)',
              'python': sys.version, 'platform': platform.platform(),
              'versions': {name: importlib.metadata.version(name) for name in ['numpy','matplotlib']},
              'figure2_endpoints': {'classical_Na0': float(classical[0]),
                                    'classical_Na5': float(classical[-1]),
                                    'squeezed_Na5': float(squeezed[-1]),
                                    'squeezed_to_classical_Na0_1': squeezed_kl(N,.1,eta1,eta2)/coherent_kl(N,.1,eta1,eta2)},
              'monte_carlo': rows}
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
