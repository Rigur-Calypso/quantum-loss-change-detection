import numpy as np
import pytest
from lossqcd.gap_diagnostics import delay_conventions,equivalent_threshold,main
from lossqcd.simulation import Runs
from lossqcd.models import coherent_pair
from lossqcd.delay_solver import conditional_delay_numerical


def test_conventions_preserve_distinct_estimands_and_early_alarms():
    runs=Runs(np.array([3,10,12]),np.ones(3,dtype=bool),20,10)
    x=delay_conventions(runs)
    assert x['conditional']==1
    assert x['signed_all_first_alarms']==pytest.approx(-5/3)
    assert x['clipped_all_first_alarms']==pytest.approx(2/3)
    assert x['conditional_inclusive']==2 and x['conditional_minus_one']==0
    assert x['early_alarm_fraction']==pytest.approx(1/3)


@pytest.mark.parametrize('runs',[
    Runs(np.array([20,12]),np.array([False,True]),20,10),
    Runs(np.array([1,2]),np.ones(2,dtype=bool),20,10),
    Runs(np.array([1,2]),np.ones(2,dtype=bool),20,None)])
def test_unavailable_metrics_are_refused(runs):
    with pytest.raises(ValueError):
        delay_conventions(runs)


def test_inverse_threshold_recovers_a_known_value_without_latency_mc():
    pre,post=coherent_pair(Na=.1)
    target=conditional_delay_numerical(pre,post,4.,change_at=20,cells=64)['conditional_delay']
    result=equivalent_threshold(pre,post,target,6.,change_at=20,cells=64)
    assert result['threshold']==pytest.approx(4.,abs=1e-7)
    assert 'counterfactual' in result['purpose']
    with pytest.raises(ValueError,match='not bracketed'):
        equivalent_threshold(pre,post,1000,6.,change_at=20,cells=64)


@pytest.mark.parametrize('args',[['--trials','1'],['--cycles','123'],['--batch-size','1']])
def test_invalid_diagnostic_budgets_do_not_create_outputs(args,tmp_path):
    out=tmp_path/'out'
    with pytest.raises(SystemExit):
        main(['--output',str(out)]+args)
    assert not out.exists()
