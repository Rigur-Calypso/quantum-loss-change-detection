"""Audit channel diagnostic source integrity, predictions and raw simulations."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256, scalar_models, numerical_row
from .simulation import Runs, simulate_stopping_times, summarize_delay


def audit(root, regenerate=False):
    root = Path(root)
    report = json.loads((root/'report.json').read_text())
    checks = []
    def check(name, passed):
        checks.append({'name': name, 'passed': bool(passed)})
    for path, expected in report['source_hashes'].items():
        check('source hash: '+path, sha256(path) == expected)
    source = json.loads((Path(report['source_root'])/'report.json').read_text())
    models = scalar_models(source)
    check('scalar coverage', len(models) == len(report['rows']) == 18)
    check('one training point', sum(r['training'] for r in report['rows']) == 1)
    csv_rows = list(csv.DictReader((root/'predictions.csv').open()))
    check('CSV coverage', len(csv_rows) == 18)
    settings = source['settings']
    check('physical settings', all(report['settings'][k] == settings[k] for k in ('N', 'eta1', 'horizon')))
    check('assumed post loss', report['settings']['assumed_eta2'] == settings['eta2'])
    streams = np.random.SeedSequence(report['settings']['seed']).spawn(18)
    for index, (model, row, exported, stream) in enumerate(zip(models, report['rows'], csv_rows, streams)):
        prefix = f'row {index}: '
        check(prefix+'model and threshold', all(row[k] == model[k] for k in ('Na', 'label', 'kind', 'change_at'))
              and row['threshold'] == model['calibration']['threshold'] and model['confirmed'])
        check(prefix+'seed', row['seed'] == int(stream.generate_state(1, dtype=np.uint64)[0]))
        published = next(r for r in source['figure7_rows'] if r['label'] == row['label'] and np.isclose(r['Na'], row['Na']))
        check(prefix+'author mapping', row['author_delay'] == published['author_delay'])
        check(prefix+'original ARL evidence', row['source_arl'] == model['validation']['arl']
              and row['source_arl_se'] == model['validation']['arl_se'])
        pre, post, actual, coarse, fine, matched = numerical_row(model, settings, report['fit']['actual_eta2'])
        check(prefix+'recomputed numerical predictions', all(
            np.isclose(saved[key], rebuilt[key], rtol=1e-10, atol=1e-10)
            for saved, rebuilt in ((row['numerical_512'],coarse),(row['numerical_1024'],fine),(row['matched_numerical'],matched))
            for key in ('conditional_delay','false_alarm_probability','log_survival_to_change')))
        with np.load(root/row['runs']) as data:
            runs = Runs(data['stop'], data['event'], report['settings']['horizon'], row['change_at'])
        stats = summarize_delay(runs)
        check(prefix+'raw delay statistics', stats == row['mc'])
        mesh = abs(fine['conditional_delay']-coarse['conditional_delay'])
        check(prefix+'independent numerical/MC agreement', stats['conditional_delay'] is not None and mesh < .05
              and abs(fine['conditional_delay']-stats['conditional_delay']) <= 4*stats['conditional_delay_se']+mesh)
        check(prefix+'unchanged pre-change first-alarm law', coarse['false_alarm_probability'] == matched['false_alarm_probability'])
        relative = fine['conditional_delay']/published['author_delay']-1
        check(prefix+'descriptive comparison', np.isclose(relative,row['relative_difference'])
              and (abs(relative)<=.05) == row['within_5_percent'])
        expected_training = row['kind'] == 'coherent' and np.isclose(row['Na'], .1)
        check(prefix+'training allocation', row['training'] == expected_training)
        if expected_training:
            fit = report['fit']
            check('training fit', abs(fine['conditional_delay']-fit['target_delay']) < .05
                  and fit['target_delay'] == row['author_delay'] and .8 <= fit['actual_eta2'] <= .85)
        check(prefix+'CSV values', exported['label'] == row['label'] and all(
            float(exported[k]) == value for k,value in (
                ('Na',row['Na']),('threshold',row['threshold']),('author_delay',row['author_delay']),
                ('original_mc_delay',row['original_mc_delay']),('prediction',fine['conditional_delay']),
                ('mc_delay',stats['conditional_delay']),('mc_se',stats['conditional_delay_se']),
                ('relative_difference',relative)))
              and exported['training'] == str(row['training'])
              and exported['within_5_percent'] == str(row['within_5_percent']))
        if regenerate:
            fresh = simulate_stopping_times(pre, post, row['threshold'], trials=report['settings']['trials'],
                horizon=report['settings']['horizon'], change_at=row['change_at'], seed=row['seed'], sampling_post=actual)
            check(prefix+'seeded raw regeneration', np.array_equal(fresh.stop,runs.stop) and np.array_equal(fresh.event,runs.event))
    check('held-out count', report['heldout_within_5_percent'] == sum(
        r['within_5_percent'] for r in report['rows'] if not r['training']))
    check('report status', report['all_checks_passed'] and report['source_unchanged'] and all(r['checks_passed'] for r in report['rows']))
    result = {'passed': all(c['passed'] for c in checks), 'checks': checks, 'regenerated': regenerate}
    (root/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    if not result['passed']:
        raise AssertionError([c['name'] for c in checks if not c['passed']])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', nargs='?', default='results_channel_diagnostics')
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    audit(args.root, args.regenerate)


if __name__ == '__main__':
    main()
