import numpy as np
import pytest
from lossqcd.models import coherent_pair,squeezed_pair
from lossqcd.arl_solver import mean_stopping_time,increment_cdf
from lossqcd.delay_solver import conditional_delay_numerical
from lossqcd.simulation import simulate_stopping_times,summarize_delay


@pytest.mark.parametrize('pair',[coherent_pair(Na=.1),squeezed_pair(Na=.1)])
def test_nu1_is_post_first_step_time_minus_one(pair):
    pre,post=pair
    x=conditional_delay_numerical(pre,post,3,change_at=1,cells=128)
    g=mean_stopping_time(pre,post,3,cells=128,sampling=post)
    assert x['conditional_delay']==pytest.approx(g['mean_stopping_time']-1,rel=1e-12)
    assert x['false_alarm_probability']==0


def test_zero_threshold_geometric_survival_and_zero_based_delay():
    pre,post=coherent_pair()
    x=conditional_delay_numerical(pre,post,0,change_at=20,cells=64)
    fpre=float(increment_cdf(0,pre,post));fpost=float(increment_cdf(0,pre,post,post))
    assert x['false_alarm_probability']==pytest.approx(1-fpre**19)
    assert x['conditional_delay']==pytest.approx(fpost/(1-fpost))


@pytest.mark.parametrize('pair',[coherent_pair(Na=.1),squeezed_pair(Na=.1)])
def test_conditioned_pre_history_matches_direct_simulation(pair):
    pre,post=pair
    x=conditional_delay_numerical(pre,post,3,change_at=40,cells=256)
    stats=summarize_delay(simulate_stopping_times(pre,post,3,trials=60000,horizon=3000,change_at=40,seed=934))
    assert abs(x['conditional_delay']-stats['conditional_delay'])<4*stats['conditional_delay_se']
    assert abs(x['false_alarm_probability']-stats['false_alarm_probability'])<4*stats['false_alarm_probability_se']


@pytest.mark.parametrize('kwargs',[{'change_at':0},{'change_at':True},{'cells':7},{'threshold':-1}])
def test_invalid_delay_inputs(kwargs):
    pre,post=coherent_pair()
    with pytest.raises(ValueError):
        conditional_delay_numerical(pre,post,**({'threshold':3}|kwargs))
