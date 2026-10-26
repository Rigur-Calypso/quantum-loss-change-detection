import pytest

from lossqcd.experiments import main


def test_baseline_cli_refuses_to_overwrite_existing_evidence(tmp_path):
    saved = tmp_path / 'report.json'
    saved.write_bytes(b'original saved evidence')
    with pytest.raises(FileExistsError):
        main(['--output', str(tmp_path)])
    assert saved.read_bytes() == b'original saved evidence'
    assert list(tmp_path.iterdir()) == [saved]
