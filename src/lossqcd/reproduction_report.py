"""Assemble and check an offline report of the main-paper evidence.

Existing auditors run in temporary directories containing links to source
links, so their new audit.json files never replace historical audit records.
Research extensions are linked in the coverage document, not counted as paper
reproductions. A successful evidence audit does not imply published agreement.
"""
import argparse
import csv
import hashlib
import html
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
from urllib.parse import quote

import numpy as np
from .models import coherent_kl, squeezed_kl
from .entangled import entangled_kl_per_pulse

STUDIES = ('results_entangled', 'results_author_comparison',
           'results_paper_completion', 'results_verified_figures_v2',
           'results_equation_validation', 'results_delay_scaling', 'results_dual_rail')
AUDITORS = (
    ('results_author_comparison', 'lossqcd.audit_author_comparison'),
    ('results_paper_completion', 'lossqcd.audit_paper_completion'),
    ('results_verified_figures_v2', 'lossqcd.audit_verified_figures'),
    ('results_equation_validation', 'lossqcd.audit_equation_validation'),
    ('results_delay_scaling', 'lossqcd.audit_delay_scaling'),
    ('results_dual_rail', 'lossqcd.audit_dual_rail'),
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def inventory(root):
    """Relative paths are portable; only the declared core studies count."""
    root = Path(root)
    paths = []
    for folder in (*STUDIES, 'src/lossqcd', 'tests', 'data/zenodo_17460928',
                   'notebooks', 'walkthrough_run_v2', 'scripts', 'source_audit',
                   'data/arxiv_source_audit', 'results_dv_uncertainty_final'):
        paths.extend(p for p in (root / folder).rglob('*')
                     if p.is_file() and '__pycache__' not in p.parts)
    paths.extend(root / name for name in
                 ('pyproject.toml', 'requirements-lock.txt', 'PAPER_IMPLEMENTATION.md',
                  'README.md', 'RESULTS.md', 'REPRODUCING.md'))
    return {str(p.relative_to(root)): digest(p) for p in sorted(set(paths))}


def verify_manifest(root, report):
    root, report = Path(root).resolve(), Path(report).resolve()
    manifest = json.loads((report / 'manifest.json').read_text())
    problems = []
    for group, base in (('inputs', root), ('outputs', report)):
        for name, expected in manifest[group].items():
            path = (base / name).resolve()
            if not path.is_relative_to(base):
                problems.append(f'{group}: unsafe path {name}')
            elif not path.is_file():
                problems.append(f'{group}: missing {name}')
            elif digest(path) != expected:
                problems.append(f'{group}: changed {name}')
    return {'passed': not problems, 'files_checked': sum(len(manifest[g]) for g in ('inputs', 'outputs')),
            'problems': problems}


def run_auditor(root, folder, module, destination):
    """Only the temporary audit.json can be written by the existing auditor."""
    source = root / folder
    with tempfile.TemporaryDirectory(prefix='lossqcd-audit-') as temp:
        stage = Path(temp)
        for path in source.iterdir():
            if path.name != 'audit.json':
                (stage / path.name).symlink_to(path.resolve(), target_is_directory=path.is_dir())
        command = [sys.executable, '-m', module, str(stage)]
        completed = subprocess.run(command, cwd=root, capture_output=True, text=True,
                                   env={**os.environ, 'MPLCONFIGDIR': os.environ.get('MPLCONFIGDIR', '/private/tmp/tiju-mpl')})
        log = completed.stdout + completed.stderr
        (destination / f'{folder}.log').write_text(log)
        audit_path = stage / 'audit.json'
        if completed.returncode != 0 or not audit_path.exists():
            return {'study': folder, 'passed': False, 'count': 0,
                    'error': f'Auditor exited {completed.returncode}; inspect saved log.'}
        result = json.loads(audit_path.read_text())
        (destination / f'{folder}.json').write_text(json.dumps(result, indent=2) + '\n')
        passed = result.get('passed', result.get('data_audit_passed', False))
        checks = result.get('checks', [])
        return {'study': folder, 'passed': bool(passed), 'count': len(checks)}


def check_figure2(root):
    path = root / 'results_entangled/figure2_all_curves.csv'
    data = np.genfromtxt(path, delimiter=',', names=True)
    checks = []
    for name in data.dtype.names[1:]:
        if name == 'coherent':
            fresh = [coherent_kl(100, float(a), .9, .85) for a in data['Na']]
        elif name == 'squeezed_n1':
            fresh = [squeezed_kl(100, float(a), .9, .85) for a in data['Na']]
        else:
            n = int(name.removeprefix('entangled_n'))
            fresh = [entangled_kl_per_pulse(100, float(a), .9, .85, n) for a in data['Na']]
        checks.append({'name': name, 'passed': bool(np.allclose(fresh, data[name], rtol=1e-12, atol=1e-14)),
                       'points': len(fresh)})
    checks.append({'name': 'Na covers 0 through 5', 'passed': bool(data['Na'][0] == 0 and data['Na'][-1] == 5)})
    return {'study': 'results_entangled Figure 2 table', 'passed': all(c['passed'] for c in checks),
            'count': len(checks), 'checks': checks}


def coverage():
    # Deliberately no aggregate reproduction percentage: evidence types differ.
    return [
        ('Eq.1', 'Finite-ARL scaling checked', 'Fourteen scalar settings; rare-event ARL and direct delay validation.', 'Asymptotic trend, not an exact finite-delay formula or proof of a worst-case theorem.', 'results_delay_scaling/delay_scaling.png'),
        ('1', 'Observation models implemented', 'Coherent, squeezed and entangled transmitter laws.', 'Hardware schematic, not an optical implementation.', ''),
        ('2', 'Equations checked', 'All coherent, squeezed and n=2–256 per-pulse entropy curves.', 'No original numeric figure array for pointwise comparison.', 'results_entangled/figure2_entangled.png'),
        ('3', 'Author-data difference remains', 'Eight unmodulated and eight BPSK energy settings; ARL validated.', 'Delay ratios differ from published data.', 'results_author_comparison/figure3_overlay.png'),
        ('4', 'Published ordering unresolved', 'Fixed-seed equations, dense covariance and random checks agree.', 'Equations give n=1 ≥ n=2; published inset gives n=2 > n=1.', 'results_verified_figures_v2/fixed_seed_equation_check.png'),
        ('5', 'Published Kennedy curve unresolved', 'Legend settings verified; independent Eq.10 integration and sampling agree.', 'About 12-fold Kennedy gap; roughly 10-photon diagnostic is not the stated setting.', 'results_equation_validation/independent_checks.png'),
        ('6', 'Receiver laws implemented', 'Homodyne, photon-resolving and binary receiver statistics.', 'No optical or electronic hardware implementation.', ''),
        ('7', 'Author-data difference remains', '45 settings: five detector families at nine energies; ARL validated.', 'Our mean delays are approximately 13.6–34.4% higher by family.', 'results_author_comparison/figure7_overlay.png'),
        ('8–10', 'Statistical network checked', 'BPSK mixtures, displaced squeezed laws, Hadamard covariance propagation.', 'No capacity-achieving code or physical Green Machine.', ''),
        ('11(a)', 'Equation benchmark checked', 'Single-photon quadrature law, CDF, sampler and entropy.', '0.206340 agrees with quoted 0.2063; plot uses rounded tap.', 'results_paper_completion/photon_probe_comparisons.png'),
        ('11(b)', 'Latency statistic differs', 'DV detector independently calibrated; CV and BPSK comparisons saved.', 'Our conditional first-alarm delay differs from source signed latency mean.', 'results_paper_completion/figure11b_comparison.png'),
        ('12', 'Sparse pixel agreement', 'Seven curves within 2% at 1 dB; legend settings verified.', 'Sparse approximate pixels, not full numeric source data.', 'results_verified_figures_v2/verified_receiver_curves.png'),
        ('App.E.2', 'Ideal logical channel implemented', 'Arbitrary dual-rail states, heralded loss, Page detection and two-bin timing.', 'No Duan-Kimble hardware dynamics, QKD security or optimal-tradeoff solution.', 'results_dual_rail/dual_rail.png'),
    ]


def evidence(root):
    def read(name): return json.loads((root / name / 'report.json').read_text())
    author, paper, verified, independent = [read(name) for name in STUDIES[1:5]]
    scaling = read('results_delay_scaling')
    slopes = {}
    for kind in ('coherent', 'squeezed'):
        low, high = [c for c in scaling['cases'] if c['kind'] == kind and c['target_arl'] in (2e6, 1e8)]
        slopes[kind] = low['kl'] * (high['numeric_delay'] - low['numeric_delay']) / np.log(high['numeric_arl'] / low['numeric_arl'])
    bpsk = float(np.mean([abs(r['relative_difference']) for r in paper['bpsk_rows']]))
    return {
        'author_arl_confirmed': sum(m['confirmed'] for m in author['models']),
        'author_arl_settings': len(author['models']),
        'new_arl_confirmed': sum(c['arl_confirmed'] for c in paper['cases']),
        'new_arl_settings': len(paper['cases']),
        'figure3_unmodulated_mean_absolute_ratio_difference': author['figure3_metrics']['mean_absolute_relative_ratio_difference'],
        'figure3_bpsk_mean_absolute_ratio_difference': bpsk,
        'figure7_by_family': author['figure7_metrics'],
        'dv_entropy': paper['analytic']['benchmarks']['dv_kl'],
        'dv_delay': paper['cases'][-1]['delay'],
        'figure12_max_sparse_pixel_difference': max(abs(x['relative_difference']) for x in verified['figure12_digitization']),
        'counterfactual_kennedy_fit_photons': independent['fitted_energy'],
        'independent_draws': independent['trials'] * (len(independent['kennedy']) + len(independent['figure4_dense'])),
        'delay_scaling_confirmed': sum(c['passed'] for c in scaling['cases']),
        'delay_scaling_cases': len(scaling['cases']),
        'delay_scaling_high_arl_normalized_slopes': slopes,
    }


def write_report(root, output, report):
    rows = coverage(); e = report['evidence']
    def link(path): return quote(os.path.relpath(root / path, output), safe='/')
    lines = ['# Paper reproduction report', '',
             '**Main paper models implemented; complete numerical reproduction remains unresolved.**', '',
             '[Paper](https://arxiv.org/html/2503.12276v2) · [Coverage details](' + link('PAPER_IMPLEMENTATION.md') + ')', '',
             '## Evidence at a glance', '',
             f"- Original author-grid study: {e['author_arl_confirmed']}/{e['author_arl_settings']} settings internally confirmed at ARL near two million.",
             f"- New BPSK/DV study: {e['new_arl_confirmed']}/{e['new_arl_settings']} settings independently confirmed.",
             f"- Figure 3 mean absolute delay-ratio differences: {e['figure3_unmodulated_mean_absolute_ratio_difference']:.2%} unmodulated, {e['figure3_bpsk_mean_absolute_ratio_difference']:.2%} BPSK.",
             f"- DV entropy: {e['dv_entropy']:.6f} nats, versus the paper's rounded 0.2063.",
             f"- DV conditional delay: {e['dv_delay']['conditional_delay']:.3f} ± {e['dv_delay']['conditional_delay_se']:.3f} pulses (one SE).", '',
             'ARL confirmation and audit passes validate our internal calculations. They do not establish published latency agreement.', '',
             '## Figure-by-figure coverage', '', '| Figure | Evidence status | Implemented | Remaining limit |', '| --- | --- | --- | --- |']
    for number, status, implementation, limit, plot in rows:
        label = f'[{number}]({link(plot)})' if plot else number
        lines.append(f'| {label} | {status} | {implementation} | {limit} |')
    lines += ['', '## Interpretation and limits', '',
              'Delay means are E[T−nu | T>=nu], with early first alarms and censoring tracked separately. These are not Lorden worst-case delays. ARL intervals are approximate Monte Carlo evidence, not certified finite-sample bounds.', '',
              'The source single-photon latency array contains two negative values and its plotting script uses a signed mean. Our first-alarm conditional mean is a different statistic. Author uncertainty estimates and original Figure 3/7 thresholds are unavailable.', '',
              '[Additional DV uncertainty re-analysis](' + link('results_dv_uncertainty_final/RESULTS.md') + ') separately compares signed and conditional summaries using empirical bootstrap intervals. Its author conditional interpretation assumes complete first-alarm entries because the source event flags are missing. This ancillary analysis is included in the checksum inventory but is not counted as an additional core study or as part of the component-audit total below.', '',
              '[Public source/version audit and specific author questions](' + link('source_audit/REPORT.md') + ') document the bounded search for missing original code; no message was sent.', '',
              f"The {e['counterfactual_kennedy_fit_photons']:.4f}-photon Kennedy fit uses one pixel point; three others are held out. It is a non-unique diagnostic, not an inferred author parameter or a correction to N=100.", '',
              'Robust-threshold, channel-mismatch and uncertainty studies remain separate research extensions. Hardware, capacity-achieving codes and the paper’s open optimal tradeoff problem are not claimed as implemented.', '',
              '## Direct check of Eq.1', '',
              f"All {e['delay_scaling_confirmed']}/{e['delay_scaling_cases']} coherent/squeezed settings pass independent ARL and cold-start delay checks, from target ARL 100 to 100 million pulses. At the largest targets the change in delay per log(ARL) agrees with 1/KL: normalized slopes {e['delay_scaling_high_arl_normalized_slopes']['coherent']:.6f} coherent and {e['delay_scaling_high_arl_normalized_slopes']['squeezed']:.6f} squeezed.", '',
              'The absolute cold-start delay remains below log(ARL)/KL at finite targets. Delay counts E[T−1]; counting post-change samples instead adds one. This is supporting numerical scaling evidence, not a proof of a worst-case theorem or an explanation of published numerical latency gaps.', '',
              '[Scaling data, grid checks, first-alarm handling and exact-regeneration commands](' + link('results_delay_scaling/RESULTS.md') + ')', '',
              '## Combined verification', '',
              f"Current audit result: **{'PASS' if report['passed'] else 'FAIL'}**; {report['check_count']} component checks. This run reads saved ensembles and recalculates summary/formula checks; it does not rerun the historical random draws.", '',
              '| Component | Passed | Checks |', '| --- | --- | --- |']
    for check in report['audits']:
        lines.append(f"| {check['study']} | {check['passed']} | {check['count']} |")
    lines += ['', '## Reproduction', '',
              '[Start with the guided notebook](' + link('notebooks/01_paper_walkthrough.ipynb') + ') or [read its executed copy](' + link('walkthrough_run_v2/01_paper_walkthrough.ipynb') + '). This teaching demonstration uses ARL 1,000 and is separate from the seven core paper studies.', '',
              'Run from the existing project environment:', '', '```sh',
              'python -m lossqcd.reproduction_report --output paper_reproduction_repeat',
              'python -m lossqcd.reproduction_report --check paper_reproduction_repeat',
              'python -m pytest -q', '```', '',
              'Each build needs a fresh output directory. The build runs six original audit modules in temporary linked directories, adds the Figure 2 table check, and saves new audit records here. Original audit.json files remain unchanged. Manifest checking detects changes to the recorded source/code/data and generated report files; the manifest itself is not an external signature.', '',
              'For exact random regeneration use the separate study commands in their linked results documents. Cached study reports contain absolute source paths; moving the project may require regenerating those studies with updated source paths. No new paper PDF download is required.', '',
              '## Next source evidence needed', '',
              '1. Figure 4: original fixed-seed covariance and energy-conversion calculation.',
              '2. Figure 5: original probe amplitude/energy and numerical plotting array.',
              '3. Figures 3/7/11: thresholds, false-alarm calibration, first-alarm arrays and event conventions.', '',
              '[Detailed independent checks and source questions](' + link('results_equation_validation/RESULTS.md') + ')', '',
              '## Selected plots', '']
    for number in ('Eq.1', '2', '3', '7', '11(b)', '12', '5', 'App.E.2'):
        row = next(r for r in rows if r[0] == number)
        lines += [f'### Figure {number}', '', f'![Figure {number}: {row[1]}]({link(row[4])})', '']
    (output / 'REPORT.md').write_text('\n'.join(lines))
    # HTML uses the same evidence content with escaped text and local assets.
    body = ['<header><p>TRANSMISSION-LOSS CHANGE DETECTION</p><h1>Paper reproduction report</h1>',
            '<p>Main models implemented. Complete numerical reproduction remains unresolved.</p></header>',
            '<main><section><h2>Evidence at a glance</h2><div class="cards">',
            f'<article><strong>{e["author_arl_confirmed"]}/{e["author_arl_settings"]}</strong><span>Author-grid ARL checks</span></article>',
            f'<article><strong>{e["new_arl_confirmed"]}/{e["new_arl_settings"]}</strong><span>New BPSK/DV ARL checks</span></article>',
            f'<article><strong>{e["dv_entropy"]:.6f}</strong><span>Single-photon KL, nats</span></article></div>',
            '<p>Audit passes support the implementation. They do not prove agreement with the published curves.</p></section>',
            '<section><h2>Figure-by-figure evidence</h2><div class="table"><table><thead><tr><th>Figure</th><th>Status</th><th>Implemented</th><th>Remaining limit</th></tr></thead><tbody>']
    for number, status, implementation, limit, plot in rows:
        label = f'<a href="{link(plot)}">{html.escape(number)}</a>' if plot else html.escape(number)
        body.append('<tr><td>' + label + '</td>' + ''.join('<td>' + html.escape(t) + '</td>' for t in (status, implementation, limit)) + '</tr>')
    body += ['</tbody></table></div></section><section><h2>How to read these results</h2>',
             '<p>Detection delay is conditional on surviving until the change. Early alarms and censored runs are tracked separately. Reported Monte Carlo intervals are approximate. Figure 12 has sparse pixel agreement; Figures 4, 5 and the absolute latency comparisons remain unresolved.</p>',
             '<p>The roughly ten-photon Kennedy comparison is a diagnostic hypothesis. It does not replace the stated bright-pulse model. Robustness extensions, hardware and open research problems are outside the numerical reproduction claim.</p>',
             f'<p>Combined saved-evidence audit: <strong>{"PASS" if report["passed"] else "FAIL"}</strong>, {report["check_count"]} component checks. Historical random draws were not rerun.</p>',
             '<p><a href="REPORT.md">Full report, uncertainty and reproduction commands</a> · <a href="summary.json">Machine-readable results</a> · <a href="manifest.json">SHA-256 inventory</a></p></section>',
             '<section><h2>Selected comparisons</h2>']
    for number in ('Eq.1', '2', '3', '7', '11(b)', '12', '5', 'App.E.2'):
        row = next(r for r in rows if r[0] == number)
        body.append(f'<figure><img src="{link(row[4])}" alt="Figure {html.escape(number)} comparison"><figcaption>Figure {html.escape(number)} — {html.escape(row[1])}. {html.escape(row[3])}</figcaption></figure>')
    body += ['</section><footer>Guha, John, Gong and Basu · <a href="https://arxiv.org/html/2503.12276v2">arXiv:2503.12276v2</a>. Generated from local saved evidence; external sources are not fetched.</footer></main>']
    css = 'body{margin:0;background:#f4f6f8;color:#172433;font:16px/1.6 system-ui,sans-serif}header{background:#162b42;color:white;padding:3rem max(5vw,1rem)}header p:first-child{letter-spacing:.14em;font-size:.8rem}h1{font-size:clamp(2rem,4vw,3rem);margin:.4rem 0}main{max-width:1180px;margin:auto;padding:1.5rem}section{background:white;border:1px solid #dae1e8;border-radius:12px;padding:1.5rem;margin-bottom:1.5rem}a{color:#165f95}h2{font-size:1.4rem}.cards{display:flex;gap:1rem;flex-wrap:wrap}.cards article{padding:1rem;background:#edf4f8;flex:1;min-width:180px;border-radius:8px}.cards strong{display:block;font-size:2rem}.cards span{display:block}.table{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{text-align:left;vertical-align:top;padding:.8rem;border-bottom:1px solid #dce3e9}th{background:#edf4f8}td:first-child{white-space:nowrap}figure{margin:1.5rem 0;border:1px solid #dce3e9;padding:.6rem;border-radius:8px}img{width:100%;height:auto}figcaption{font-size:.9rem;padding:.8rem}footer{font-size:.9rem;padding:1rem;color:#4b5c6b}'
    (output / 'report.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Paper reproduction report</title><style>' + css + '</style></head><body>' + ''.join(body) + '</body></html>')
    with (output / 'figure_coverage.csv').open('w', newline='') as handle:
        writer = csv.writer(handle); writer.writerow(('figure', 'status', 'implemented', 'remaining_limit', 'plot'))
        writer.writerows(rows)


def build(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if any(output.is_relative_to(root / study) for study in STUDIES):
        raise ValueError('Report output must be outside the existing study directories.')
    if not sys.flags.optimize == 0:
        raise RuntimeError('Run without Python -O: some historical auditors use assertions.')
    before = inventory(root)
    output.mkdir(parents=True, exist_ok=False)
    destination = output / 'audits'; destination.mkdir()
    audits = [check_figure2(root)]
    for folder, module in AUDITORS:
        print(f'Checking {folder}...', flush=True)
        audits.append(run_auditor(root, folder, module, destination))
    unchanged = inventory(root) == before
    result = {'paper_url': 'https://arxiv.org/html/2503.12276v2',
              'passed': unchanged and all(a['passed'] for a in audits),
              'published_reproduction_complete': False, 'inputs_unchanged': unchanged,
              'check_count': sum(a['count'] for a in audits), 'audits': audits,
              'evidence': evidence(root), 'python': platform.python_version(),
              'versions': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'matplotlib')},
              'scope': 'Seven core paper studies; research extensions excluded. Saved-data recalculation, not exact regeneration. Dual-rail addition is an ideal Appendix E.2 model, not a new published numerical figure.'}
    (output / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    write_report(root, output, result)
    outputs = {str(p.relative_to(output)): digest(p) for p in sorted(output.rglob('*')) if p.is_file()}
    (output / 'manifest.json').write_text(json.dumps({'algorithm': 'SHA-256', 'inputs': before, 'outputs': outputs}, indent=2) + '\n')
    verification = verify_manifest(root, output)
    if not verification['passed']: raise RuntimeError(verification['problems'])
    return result, verification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--output', type=Path, default=Path('paper_reproduction'))
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        result = verify_manifest(args.root, args.check)
        print(json.dumps(result, indent=2))
        raise SystemExit(0 if result['passed'] else 1)
    result, verification = build(args.root, args.output)
    print(f"{result['check_count']} component checks; passed={result['passed']}; {verification['files_checked']} file hashes verified.")
    raise SystemExit(0 if result['passed'] else 1)


if __name__ == '__main__': main()
