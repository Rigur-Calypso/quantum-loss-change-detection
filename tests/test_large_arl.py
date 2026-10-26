import json
import pytest
from lossqcd.models import coherent_pair
from lossqcd.arl_solver import mean_stopping_time
from lossqcd.large_arl import crosscheck_existing_monte_carlo,main


def test_historical_crosscheck_reports_agreement_and_disagreement(tmp_path):
    pre,post=coherent_pair()
    value=mean_stopping_time(pre,post,2.,cells=1024)['mean_stopping_time']
    model={'pre':{'mean':pre.mean,'variance':pre.variance},
           'post':{'mean':post.mean,'variance':post.variance},
           'calibration':{'threshold':2.},'validation':{'arl':value,'arl_se':.1,'censored':0}}
    data={'conditions':[{'Na':0,'target_arl':value,'models':{'coherent':model}}]}
    path=tmp_path/'report.json';path.write_text(json.dumps(data))
    result=crosscheck_existing_monte_carlo(path)
    assert result['all_passed']
    assert len(result['source_sha256'])==64
    model['validation']['arl']=value*1.2
    path.write_text(json.dumps(data))
    assert not crosscheck_existing_monte_carlo(path)['all_passed']
    model['validation']['censored']=1
    path.write_text(json.dumps(data))
    with pytest.raises(RuntimeError,match='censored'):
        crosscheck_existing_monte_carlo(path)


def test_missing_benchmark_refused_before_output_created(tmp_path):
    output=tmp_path/'output'
    with pytest.raises(SystemExit):
        main(['--output',str(output),'--validation-source',str(tmp_path/'missing.json')])
    assert not output.exists()
