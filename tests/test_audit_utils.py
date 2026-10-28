"""Golden source hashes and rejection of inconsistent saved diagnostics."""

import pytest

from lossqcd.audit_utils import compare_cycle_summaries, sha256


@pytest.mark.parametrize(
    "content, expected",
    [
        (b"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        (b"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"),
    ],
)
def test_sha256_matches_known_vectors_for_both_path_types(tmp_path, content, expected):
    path = tmp_path / "source"
    path.write_bytes(content)
    assert sha256(path) == sha256(str(path)) == expected


@pytest.mark.parametrize(
    "field", ["arl", "arl_se", "success_weight_ess", "largest_weight_fraction"]
)
def test_cycle_comparison_rejects_each_inconsistent_quantity(field):
    saved = {
        "arl": 2e6,
        "arl_se": 1000.0,
        "success_weight_ess": 2000.0,
        "largest_weight_fraction": 0.001,
    }
    compare_cycle_summaries(saved, dict(saved))
    changed = {**saved, field: saved[field] * 1.01}
    with pytest.raises(AssertionError, match=field):
        compare_cycle_summaries(changed, saved)
