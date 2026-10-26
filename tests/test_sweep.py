import numpy as np
import pytest
import lossqcd.sweep as sweep
from lossqcd.simulation import Runs


@pytest.mark.parametrize('Na,calibration_count', [(0.,1),(.1,2)])
def test_threshold_reuse_only_for_identical_models_and_independent_ensembles(monkeypatch,Na,calibration_count):
    calibrations, simulation_seeds = [], []
    def fake_calibration(pre,post,target,**kwargs):
        calibrations.append(kwargs['seed'])
        return {'threshold':1.5,'history':[]}
    real_simulate = sweep.simulate_stopping_times
    def spy_simulate(*args,**kwargs):
        simulation_seeds.append(kwargs['seed'])
        return real_simulate(*args,**kwargs)
    monkeypatch.setattr(sweep,'calibrate_threshold',fake_calibration)
    monkeypatch.setattr(sweep,'simulate_stopping_times',spy_simulate)
    report,raw=sweep.run_condition(Na,500,seed=44,calibration_trials=20,
                                   validation_trials=100,delay_trials=100)
    assert len(calibrations)==calibration_count
    assert len(set(calibrations+simulation_seeds))==calibration_count+6
    if Na==0:
        assert report['models']['squeezed']['calibration']['reused_from']=='coherent'
        assert report['kl_ratio']==pytest.approx(1)
        assert report['models']['coherent']['pre']==report['models']['squeezed']['pre']
        assert report['models']['coherent']['post']==report['models']['squeezed']['post']
    else:
        assert 'reused_from' not in report['models']['squeezed']['calibration']
    assert len(raw)==12
    rows=sweep.summary_rows(report)
    assert [row['change_at'] for row in rows]==[1,100]
    assert rows[0]['coherent_false_alarm_probability']==0
    assert rows[0]['squeezed_false_alarm_probability']==0
    assert rows[1]['coherent_delay']==report['models']['coherent']['delay']['100']['conditional_delay']


def test_censored_sweep_refuses_unrestricted_comparison(monkeypatch):
    monkeypatch.setattr(sweep,'calibrate_threshold',lambda *a,**kw:{'threshold':1.5})
    def censored(*args,**kw):
        n,h=kw['trials'],kw['horizon']
        return Runs(np.full(n,h),np.zeros(n,dtype=bool),h,None)
    monkeypatch.setattr(sweep,'simulate_stopping_times',censored)
    with pytest.raises(RuntimeError,match='censored ARL'):
        sweep.run_condition(.1,1000,seed=1,calibration_trials=2,validation_trials=2,delay_trials=2)


@pytest.mark.parametrize('argv', [['--Na','-1'],['--Na','0','0'],['--target-arl','1'],
                                  ['--validation-trials','1'],['--match-tolerance','0']])
def test_invalid_sweep_settings_rejected_before_creating_outputs(argv,tmp_path):
    output=tmp_path/'no_output'
    with pytest.raises(SystemExit):
        sweep.main(['--output',str(output)]+argv)
    assert not output.exists()
