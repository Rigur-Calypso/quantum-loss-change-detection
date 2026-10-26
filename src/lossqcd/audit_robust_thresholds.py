"""Audit conservative threshold selection, independent validation and delay costs."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .robust_thresholds import GRID,source_models,detector_args,threshold_update,approximate_grid_lower_bound,cycle_pair,exported_rows
from .conditional_rare_event import summarize_conditional_cycles,estimate_conditional_arl,conditional_cycle_arrays
from .conditional_detection import simulate_conditional
from .batched_rare_event import pool_cycle_summaries
from .simulation import Runs,summarize_delay


def audit(root,regenerate=False,regenerate_calibration=False):
    root=Path(root);report=json.loads((root/'report.json').read_text());s=report['settings']
    source=json.loads((Path(report['source_root'])/'report.json').read_text());models=source_models(source);checks=[]
    def check(name,passed):checks.append({'name':name,'passed':bool(passed)})
    for path,digest in report['source_hashes'].items():check('source '+path,sha256(path)==digest)
    check('scope',len(models)==len(report['models'])==5 and len(report['rows'])==25 and s['grid']==list(GRID)
        and s['validation_points']==25 and s['target_arl']==2000000 and s['design_target']==2150000 and s['family_alpha']==.05)
    children=np.random.SeedSequence(s['seed']).spawn(200)
    seed=lambda i:int(children[i].generate_state(1,dtype=np.uint64)[0])
    calibration_seeds=[];heldout_seeds=[]
    for mi,(model,saved) in enumerate(zip(models,report['models'])):
        prefix=f'model {mi}: ';fit=saved['calibration'];original=model['calibration']['threshold']
        expected=[seed(20*mi+j) for j in range(20)];calibration_seeds.extend(expected)
        check(prefix+'mapping and seed allocation',all(saved[k]==model[k] for k in ('kind','n','label')) and fit['seeds']==expected)
        check(prefix+'four fixed selection rounds',len(fit['history'])==fit['iterations']==4 and fit['cycles_per_law']==s['calibration_cycles']
            and fit['design_target']==s['design_target'] and fit['original_threshold']==original)
        h=original
        for iteration,step in enumerate(fit['history']):
            check(prefix+f'round {iteration} threshold',step['threshold']==h and len(step['evaluations'])==5)
            for j,evaluation in enumerate(step['evaluations']):
                stats=evaluation['summary'];eta=GRID[j]
                expected_params={k:v for k,v in detector_args(model,eta).items() if k not in ('kind','n') and (k!='sampling_eta1' or eta!=.9)}
                check(prefix+f'round {iteration} grid {j}',evaluation['true_eta1']==eta and stats['seed']==expected[5*iteration+j]
                    and stats['threshold']==h and stats['parameters']==expected_params
                    and stats['kind']==model['kind'] and stats['n']==model['n']
                    and stats['pre_cycles']==stats['post_cycles']==s['calibration_cycles'])
                if regenerate_calibration:
                    fresh,_,_=estimate_conditional_arl(threshold=h,cycles=s['calibration_cycles'],seed=stats['seed'],
                        proposal='adaptive',proposal_offset=2.,**detector_args(model,eta))
                    check(prefix+f'exact calibration regeneration {iteration}/{j}',fresh==stats)
            worst=min(step['evaluations'],key=lambda e:e['summary']['arl'])
            check(prefix+f'round {iteration} worst grid point',step['minimum_arl']==worst['summary']['arl'] and step['worst_eta1']==worst['true_eta1'])
            if iteration<3:h=threshold_update(original,h,step['minimum_arl'],s['design_target'])
        check(prefix+'frozen conservative threshold',fit['threshold']==h and h>=original)
    for i,row in enumerate(report['rows']):
        model=models[i//5];eta=GRID[i%5];fit=report['models'][i//5]['calibration'];prefix=f'holdout {i}: '
        expected=[seed(100+4*i+j) for j in range(4)];heldout_seeds.extend(expected)
        check(prefix+'mapping and thresholds',all(row[k]==model[k] for k in ('kind','n','label','change_at')) and row['true_eta1']==eta
            and row['threshold']==fit['threshold'] and row['baseline_threshold']==model['calibration']['threshold'])
        check(prefix+'independent seeds',row['seeds']==expected[:2])
        for repeat in range(2):
            saved=row['repeats'][repeat]
            with np.load(root/row['cycle_files'][repeat]) as data:
                pre,post=cycle_pair(data,saved);stats=summarize_conditional_cycles(pre,post)
                check(prefix+f'raw moments and weights {repeat}',all(saved[k]==v for k,v in stats.items()))
                expected_params={k:v for k,v in detector_args(model,eta).items() if k not in ('kind','n') and (k!='sampling_eta1' or eta!=.9)}
                check(prefix+f'count and metadata {repeat}',len(pre.length)==len(post.length)==s['validation_cycles']
                    and saved['parameters']==expected_params and saved['kind']==model['kind'] and saved['n']==model['n']
                    and saved['threshold']==row['threshold'] and saved['seed']==expected[repeat])
                streams=[int(c.generate_state(1,dtype=np.uint64)[0]) for c in np.random.SeedSequence(expected[repeat]).spawn(2)]
                check(prefix+f'law streams and diagnostics {repeat}',saved['stream_seeds']==streams and streams[0]!=streams[1] and stats['diagnostics_passed'])
                if regenerate:
                    _,a,b=estimate_conditional_arl(threshold=row['threshold'],cycles=s['validation_cycles'],seed=expected[repeat],
                        proposal='adaptive',proposal_offset=2.,**detector_args(model,eta))
                    check(prefix+f'exact cycle regeneration {repeat}',all(np.array_equal(data[k],v) for k,v in conditional_cycle_arrays(a,b).items()))
        pooled=pool_cycle_summaries(row['repeats'])
        check(prefix+'pooled sufficient moments',pooled==row['pooled'])
        lower=approximate_grid_lower_bound(pooled['arl'],pooled['arl_se'],s['validation_points'],s['family_alpha'])
        check(prefix+'adjusted lower ARL',lower==row['grid_lower_bound'] and lower['lower']>=s['target_arl'])
        z=float(abs(row['repeats'][0]['arl']-row['repeats'][1]['arl'])/np.hypot(row['repeats'][0]['arl_se'],row['repeats'][1]['arl_se']))
        check(prefix+'independent repetitions',z==row['repeat_z'] and z<=4 and pooled['diagnostics_passed'])
        summaries=[]
        for j,(label,h) in enumerate((('baseline',row['baseline_threshold']),('conservative',row['threshold']))):
            saved=row['delays'][label]
            with np.load(root/saved['file']) as data:runs=Runs(data['stop'],data['event'],s['horizon'],row['change_at'])
            stats=summarize_delay(runs);summaries.append(stats)
            check(prefix+label+' raw first alarms',len(runs.stop)==s['changed_trials'] and stats==saved['summary']
                and saved['seed']==expected[2+j] and stats['conditional_delay'] is not None)
            if regenerate:
                fresh=simulate_conditional(threshold=h,trials=s['changed_trials'],horizon=s['horizon'],change_at=row['change_at'],
                    seed=saved['seed'],**detector_args(model,eta))
                check(prefix+label+' exact regeneration',np.array_equal(runs.stop,fresh.stop) and np.array_equal(runs.event,fresh.event))
        a,b=summaries
        check(prefix+'independent delay difference',row['delay_cost']==b['conditional_delay']-a['conditional_delay']
            and row['delay_cost_se']==float(np.hypot(a['conditional_delay_se'],b['conditional_delay_se'])))
        check(prefix+'classification',row['passed'])
    check('selection independent of holdouts',len(set(calibration_seeds+heldout_seeds))==200 and set(calibration_seeds).isdisjoint(heldout_seeds))
    with (root/'validation.csv').open() as handle:exported=list(csv.DictReader(handle))
    check('CSV',exported==[{k:str(v) for k,v in row.items()} for row in exported_rows(report)])
    check('status',report['source_unchanged'] and report['all_passed'])
    result={'passed':all(c['passed'] for c in checks),'regenerated':regenerate,
        'regenerated_calibration':regenerate_calibration,'checks':checks}
    (root/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    if not result['passed']:raise AssertionError([c['name'] for c in checks if not c['passed']])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',nargs='?',default='results_robust_thresholds');p.add_argument('--regenerate',action='store_true')
    p.add_argument('--regenerate-calibration',action='store_true')
    a=p.parse_args();audit(a.root,a.regenerate,a.regenerate_calibration)


if __name__=='__main__':main()
