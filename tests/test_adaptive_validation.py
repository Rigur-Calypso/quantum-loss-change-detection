import numpy as np
import pytest
from lossqcd.adaptive_validation import audit_mass,main


@pytest.mark.parametrize('args',[['--cycles','1'],['--batch-size','1'],
    ['--cycles','123'],['--direct-trials','1'],['--moderate-cycles','-1']])
def test_invalid_followup_budget_does_not_create_output(args,tmp_path):
    out=tmp_path/'results'
    with pytest.raises(SystemExit):
        main(['--output',str(out)]+args)
    assert not out.exists()


def test_mass_and_cycle_length_identities_and_zero_uncertainty(tmp_path):
    path=tmp_path/'audit.npz'
    def save(weights):
        np.savez(path,post_log_weight=np.log(weights),post_length=np.array([2,2]),
            pre_length=np.array([2,2]),pre_complete=np.ones(2,dtype=bool),post_complete=np.ones(2,dtype=bool))
    save(np.array([.5,1.5]))
    assert audit_mass(path)['passed']
    save(np.array([1.,1.]))
    assert audit_mass(path)['passed']
    save(np.array([.1,.2]))
    assert not audit_mass(path)['passed']
