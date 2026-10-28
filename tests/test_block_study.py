import pytest
from lossqcd.block_study import main


@pytest.mark.parametrize("args", [["--Na", "-1"], ["--block-sizes", "1"], ["--change-at", "0"]])
def test_invalid_experiment_refused_before_output_creation(args, tmp_path):
    out = tmp_path / "output"
    with pytest.raises(SystemExit):
        main(["--output", str(out)] + args)
    assert not out.exists()
