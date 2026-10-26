"""Recompute saved sweep summaries and audit batches, without rerunning the sweep."""
import argparse
import json
from pathlib import Path
import numpy as np
from .conditional_rare_event import ConditionalCycles,summarize_conditional_cycles
from .batched_rare_event import pool_cycle_summaries
from .conditional_sweep import make_rows
from .simulation import Runs,summarize_delay


def same_numbers(actual,expected,keys):
    for key in keys:
        np.testing.assert_allclose(actual[key],expected[key],rtol=1e-11,atol=1e-14,err_msg=key)


def audit(root):
    root=Path(root);report=json.loads((root/'report.json').read_text());checks=[]
    for index,m in enumerate(report['models']):
        folder=root/f'energy_{index//5:02d}_model_{index%5:02d}'
        if not m['shared_arl']:
            repetitions=[]
            for j in range(2):
                saved=json.loads((folder/f'validation_{j}.json').read_text())
                pooled=pool_cycle_summaries(saved['batches'])
                same_numbers(pooled,saved['summary'],['arl','arl_se','success_weight_ess','largest_weight_fraction'])
                first=saved['batches'][0];arrays=[]
                with np.load(folder/f'audit_cycles_{j}.npz') as raw:
                    for law in ['pre','post']:
                        batch=ConditionalCycles(*(raw[f'{law}_{field}'] for field in
                            ['length','raw_sum','page','upper','complete']),m['kind'],m['n'],
                            m['calibration']['threshold'],law,(100.,m['Na'],.9,.85))
                        assert np.all(batch.complete)
                        assert np.all(batch.page[batch.upper]>batch.threshold)
                        assert np.all(batch.page[~batch.upper]==0.)
                        assert np.all(batch.length[~batch.upper]%m['n']==0)
                        assert np.all(batch.raw_sum<=batch.page+1e-9)
                        arrays.append(batch)
                recalculated=summarize_conditional_cycles(*arrays)
                same_numbers(recalculated,first,['arl','arl_se','success_weight_ess','largest_weight_fraction'])
                repetitions.append(pooled)
                checks.append(f'{index}: repetition {j}, pooled moments and raw audit batch')
            pooled=pool_cycle_summaries(repetitions)
            same_numbers(pooled,m['validation'],['arl','arl_se','success_weight_ess','largest_weight_fraction'])
            weighted=json.loads((folder/'moderate_batches.json').read_text())
            same_numbers(pool_cycle_summaries(weighted['batches']),m['moderate_check']['weighted'],['arl','arl_se'])
            checks.append(f'{index}: combined holdouts and moderate moments')
        else:
            base=report['models'][0]
            assert m['Na']==0 and m['calibration']==base['calibration'] and m['validation']==base['validation']
            checks.append(f'{index}: explicitly shared zero-energy evidence')
        with np.load(folder/'changed_runs.npz') as raw:
            prefix_stops=[]
            for nu in m['delay']:
                run=Runs(raw[f'nu{nu}_stop'],raw[f'nu{nu}_event'],report['settings']['horizon'],int(nu))
                assert summarize_delay(run)==m['delay'][nu]
                prefix_stops.append(np.where(run.stop<=report['settings']['burnin'],run.stop,-1))
                checks.append(f'{index}: changed-run summary on pulse {nu}')
            assert all(np.array_equal(prefix_stops[0],p) for p in prefix_stops[1:])
    assert make_rows(report['models'])==report['rows']
    result={'data_audit_passed':True,'checks':checks,'check_count':len(checks),
        'confirmed_settings':report['confirmed_settings'],'total_settings':report['total_settings'],
        'unconfirmed':[{'Na':m['Na'],'label':m['label'],'arl_matched':m['arl_match']['passed'],
            'repeat_passed':m['repeat_passed'],'moderate_passed':m['moderate_check']['passed']}
            for m in report['models'] if not m['confirmed']],
        'zero_control_passed':report['zero_control']['passed']}
    (root/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',type=Path)
    args=p.parse_args(argv)
    result=audit(args.results)
    print(f"{result['check_count']} data checks passed; {result['confirmed_settings']}/{result['total_settings']} settings confirmed.")


if __name__=='__main__':
    main()
