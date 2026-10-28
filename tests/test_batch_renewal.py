from dataclasses import replace
import numpy as np
import pytest
from lossqcd.batch_renewal import pool_moments, pool_cycle_summaries, batched_estimate
from lossqcd.corr_renewal import estimate_conditional_arl, summarize_conditional_cycles


def test_pooled_moments_match_concatenated_sample():
    arrays = [np.array([1.0, 2.0, 3.0]), np.array([10.0, 11.0, 12.0, 13.0])]
    ensembles = [(len(x), x.mean(), x.std(ddof=1) / np.sqrt(len(x))) for x in arrays]
    n, mean, se = pool_moments(ensembles)
    full = np.concatenate(arrays)
    assert n == len(full)
    assert mean == pytest.approx(full.mean())
    assert se == pytest.approx(full.std(ddof=1) / np.sqrt(len(full)))


def test_pooling_cycle_moments_equals_full_raw_array_estimate():
    batches = []
    pres = []
    posts = []
    for seed in [91, 92, 93]:
        stats, pre, post = estimate_conditional_arl(n=4, threshold=5, cycles=2000, seed=seed)
        batches.append(stats)
        pres.append(pre)
        posts.append(post)

    def combine(data):
        return replace(
            data[0],
            **{
                k: np.concatenate([getattr(b, k) for b in data])
                for k in ["length", "raw_sum", "page", "upper", "complete"]
            },
        )

    expected = summarize_conditional_cycles(combine(pres), combine(posts))
    actual = pool_cycle_summaries(batches)
    for key in [
        "arl",
        "arl_se",
        "relative_se",
        "success_weight_ess",
        "largest_weight_fraction",
        "mean_pre_cycle_length_se",
        "cycle_alarm_probability_se",
    ]:
        assert actual[key] == pytest.approx(expected[key], rel=1e-12)
    with pytest.raises(ValueError, match="same detector"):
        pool_cycle_summaries([batches[0], batches[1] | {"threshold": 6}])
    with pytest.raises(RuntimeError, match="unfinished"):
        pool_cycle_summaries([batches[0] | {"unfinished_pre_cycles": 1}])


def test_batched_audit_and_repeatability(tmp_path):
    target = tmp_path / "audit.npz"
    a = batched_estimate(
        "entangled", 2, 3, 0.1, cycles=2000, batch_size=1000, seed=83, audit_path=target
    )
    b = batched_estimate("entangled", 2, 3, 0.1, cycles=2000, batch_size=1000, seed=83)
    assert a == b
    assert a["batches"][0]["seed"] != a["batches"][1]["seed"]
    with np.load(target) as raw:
        assert len(raw["pre_length"]) == 1000
        assert np.all(raw["post_complete"])


@pytest.mark.parametrize("ensembles", [[], [(1, 1.0, 0.0)], [(10, np.nan, 0.0)], [(10, 1.0, -1.0)]])
def test_invalid_moments(ensembles):
    with pytest.raises(ValueError):
        pool_moments(ensembles)
