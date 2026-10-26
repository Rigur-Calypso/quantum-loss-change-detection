"""State-channel identities and physical clocks must agree independently."""
import numpy as np
import pytest
from lossqcd.dual_rail import pure_state,loss_kraus,transmit,heralded_state,simulate_heralds,summarize_dual_rail,DualRailRuns
from lossqcd.simulation import Runs


@pytest.mark.parametrize('etas',[(0,0),(1,1),(.9,.9),(.9,.2),(0,1)])
def test_kraus_channel_is_trace_preserving_and_positive(etas):
    ks=loss_kraus(*etas)
    np.testing.assert_allclose(sum(k.conj().T@k for k in ks),np.eye(2),atol=1e-14)
    rho=np.array([[.4,.1+.2j],[.1-.2j,.6]])
    out=transmit(rho,*etas)
    assert np.trace(out)==pytest.approx(1)
    assert np.linalg.eigvalsh(out).min()>=-1e-14


@pytest.mark.parametrize('eta',[.01,.2,.9,1.])
def test_uniform_loss_is_encoding_independent_erasure_and_preserves_coherence(eta):
    rho=pure_state([np.sqrt(.3),1j*np.sqrt(.7)])
    probability,conditional=heralded_state(rho,eta)
    assert probability==pytest.approx(eta)
    np.testing.assert_allclose(conditional,rho,atol=1e-14)
    out=transmit(rho,eta)
    np.testing.assert_allclose(out[1:,1:],eta*rho,atol=1e-14)
    assert out[0,0]==pytest.approx(1-eta)


def test_zero_survival_does_not_fabricate_a_received_qubit():
    probability,conditional=heralded_state(np.eye(2)/2,0)
    assert probability==0 and conditional is None


def test_mid_pair_differential_loss_is_not_state_independent():
    a,b=1/np.sqrt(2),1j/np.sqrt(2)
    probability,rho=heralded_state(pure_state([a,b]),.9,.3)
    expected=np.array([np.sqrt(.9)*a,np.sqrt(.3)*b])/np.sqrt(.6)
    assert probability==pytest.approx(.6)
    np.testing.assert_allclose(rho,np.outer(expected,expected.conj()),atol=1e-14)
    assert rho[0,0].real==pytest.approx(.75)


@pytest.mark.parametrize('change_pulse',[201,202])
def test_physical_delay_uses_two_bins_and_change_phase(change_pulse):
    runs=DualRailRuns(Runs(np.array([50,101,104]),np.ones(3,dtype=bool),200,101),change_pulse)
    stats=summarize_dual_rail(runs)
    assert stats['pair_statistics']['false_alarms']==1
    expected=np.mean(np.array([202,208])-change_pulse)
    assert stats['physical_time_bin_statistics']['conditional_delay']==expected


def test_censoring_keeps_unrestricted_physical_arl_unavailable():
    runs=DualRailRuns(Runs(np.array([2,10]),np.array([True,False]),10,None),None)
    result=summarize_dual_rail(runs)
    assert result['physical_time_bin_statistics']['arl'] is None
    assert result['physical_time_bin_statistics']['restricted_mean_min_T_H']==12


def test_uniform_loss_heralds_are_identical_for_all_encodings_with_shared_coin_stream():
    examples=[simulate_heralds(.9,.7,2,trials=100,horizon_pairs=1000,change_pulse=101,encoding=e,seed=11)
        for e in ('haar','early','late','plus')]
    for other in examples[1:]:
        np.testing.assert_array_equal(other.pairs.stop,examples[0].pairs.stop)
        np.testing.assert_array_equal(other.pairs.event,examples[0].pairs.event)


@pytest.mark.parametrize('state',[np.eye(3)/3,np.array([[1,1],[0,0]]),np.diag([1.1,-.1]),np.eye(2)])
def test_invalid_qubit_density_rejected(state):
    with pytest.raises(ValueError):transmit(state,.9)


@pytest.mark.parametrize('kwargs',[{'eta1':-1},{'eta2':1},{'change_pulse':0},{'encoding':'unknown'},{'trials':True}])
def test_invalid_detector_inputs_rejected(kwargs):
    with pytest.raises(ValueError):simulate_heralds(**({'eta1':.9,'eta2':.8,'threshold':2}|kwargs))
