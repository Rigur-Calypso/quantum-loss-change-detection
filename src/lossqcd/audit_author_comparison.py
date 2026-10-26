"""Audit source checksums, frozen ARL calibration and author comparison arrays."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .author_comparison import author_data
from .batched_rare_event import pool_cycle_summaries
from .conditional_rare_event import ConditionalCycles,summarize_conditional_cycles
from .conditional_large_arl import assess_match
from .conditional_sweep import target_match
from .simulation import Runs,summarize_delay
from .calibration import independent_ratio


def compare(actual,expected):
    for key in ['arl','arl_se','success_weight_ess','largest_weight_fraction']:
        np.testing.assert_allclose(actual[key],expected[key],rtol=1e-11,atol=1e-14,err_msg=key)


def audit(root):
    root=Path(root);report=json.loads((root/'report.json').read_text());checks=[]
    authors=author_data(report['settings']['data'])
    assert authors==report['author_data'];checks.append('author checksums, metadata and analytical KL coordinates')
    for i,m in enumerate(report['models']):
        folder=root/Path(m['folder']).name;h=12.
        assert len(m['calibration']['history'])==5
        for entry in m['calibration']['history']:
            assert entry['threshold']==h
            h=max(0.,h+float(np.clip(np.log(2e6/entry['arl']),-2.,2.)))
        assert h==m['calibration']['threshold']
        repetitions=[]
        for j in range(2):
            batch=json.loads((folder/f'validation_{j}.json').read_text())
            summary=pool_cycle_summaries(batch['batches']);compare(summary,batch['summary'])
            assert summary['threshold']==h
            arrays=[]
            with np.load(folder/f'audit_{j}.npz') as raw:
                for law in ['pre','post']:
                    c=ConditionalCycles(*(raw[f'{law}_{field}'] for field in
                        ['length','raw_sum','page','upper','complete']),m['kind'],m['n'],h,
                        'pre' if law=='pre' else 'adaptive',(100.,m['Na'],.9,.85),
                        None if law=='pre' else raw['post_log_weight'],None if law=='pre' else 2.)
                    assert np.all(c.complete)
                    assert np.all(c.page[c.upper]>h)
                    assert np.all(c.page[~c.upper]==0)
                    assert np.all(c.length[~c.upper]%m['n']==0)
                    arrays.append(c)
            compare(summarize_conditional_cycles(*arrays),batch['batches'][0])
            repetitions.append(summary)
            checks.append(f'{i}: repetition {j}, frozen threshold, raw audit and batch moments')
        combined=pool_cycle_summaries(repetitions);compare(combined,m['validation'])
        z=float(abs(repetitions[0]['arl']-repetitions[1]['arl'])/np.hypot(repetitions[0]['arl_se'],repetitions[1]['arl_se']))
        np.testing.assert_allclose(z,m['repeat_z'],rtol=1e-11)
        assert bool(z<=4 and all(s['diagnostics_passed'] for s in repetitions))==m['repeat_passed']
        with np.load(folder/'changed_runs.npz') as raw:
            assert summarize_delay(Runs(raw['stop'],raw['event'],5000,m['change_at']))==m['delay']
        base=next(b for b in report['models'] if b['figure']==m['figure'] and b['Na']==m['Na'] and b['kind']=='coherent')
        match=({'passed':target_match(combined,2e6,.05),'ratio':1.} if m is base else assess_match(combined,base['validation'],2e6,.05))
        assert match==m['arl_match']
        expected=bool(match['passed'] and m['repeat_passed'] and (m is base or base['confirmed']))
        assert expected==m['confirmed']
        checks.append(f'{i}: pooled ARL, saved delay and confirmation classification')
    for name,rows in [('figure7_comparison.csv',report['figure7_rows']),('figure3_comparison.csv',report['figure3_rows'])]:
        with (root/name).open(newline='') as f:
            assert list(csv.DictReader(f))==[{k:str(v) for k,v in row.items()} for row in rows]
    for r in report['figure7_rows']:
        ei=authors['figure7']['energies'].index(r['Na'])
        index=['Coherent','Squeezed n=1','Entangled n=2','Entangled n=4','Entangled n=8'].index(r['label'])
        assert r['author_delay']==authors['figure7']['means'][index][ei]
        assert r['relative_difference']==r['our_delay']/r['author_delay']-1
    for r in report['figure3_rows']:
        index=authors['figure3']['energies'].index(r['Na'])
        assert r['author_ratio']==authors['figure3']['unmodulated_ratios'][index]
        assert r['kl_ratio']==authors['figure3']['kl_ratios'][index]
        pair=[next(m for m in report['models'] if m['figure']==3 and m['Na']==r['Na'] and m['kind']==kind)
            for kind in ['coherent','entangled']]
        c,q=[m['delay'] for m in pair]
        ratio=independent_ratio(c['conditional_delay'],c['conditional_delay_se'],q['conditional_delay'],q['conditional_delay_se'])
        assert ratio['ratio']==r['our_ratio'] and ratio['se']==r['our_ratio_se']
        assert r['relative_difference']==r['our_ratio']/r['author_ratio']-1
    checks.append('all author-to-column mappings, numerical differences and CSV rows')
    result={'data_audit_passed':True,'checks':checks,'check_count':len(checks),
        'simulation_settings':len(report['models']),
        'confirmed_arl_settings':sum(m['confirmed'] for m in report['models']),
        'author_kl_max_error':authors['figure3']['max_absolute_kl_error'],
        'latency_reproduction_claimed':False}
    (root/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    args=p.parse_args(argv);result=audit(args.results)
    print(f"{result['check_count']} data checks passed; {result['confirmed_arl_settings']}/{result['simulation_settings']} ARL settings confirmed. Latency gaps retained.")


if __name__=='__main__':
    main()
