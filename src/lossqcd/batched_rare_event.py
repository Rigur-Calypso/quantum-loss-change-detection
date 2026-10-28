"""Pool independent renewal ensembles using sufficient first/second moments.

No self-normalization of importance weights and no averaging of ARL ratios.
Each batch keeps its seed and full summary; optional audit batches keep arrays.
"""

import numpy as np
from .conditional_rare_event import estimate_conditional_arl, conditional_cycle_arrays


def pool_moments(ensembles):
    """(count, mean, SE) -> exact pooled sample mean and SE, up to rounding."""
    if not ensembles or any(
        n < 2 or not np.isfinite(mu) or not np.isfinite(se) or se < 0 for n, mu, se in ensembles
    ):
        raise ValueError("need finite means/SEs and counts >=2")
    total = sum(n for n, _, _ in ensembles)
    mean = sum(n * mu for n, mu, _ in ensembles) / total
    centered = sum((n - 1) * n * se**2 + n * (mu - mean) ** 2 for n, mu, se in ensembles)
    return total, float(mean), float(np.sqrt(centered / ((total - 1) * total)))


def pool_cycle_summaries(batches):
    """Pool raw cycle moments before forming the ARL ratio and weight diagnostics.

    All batches must describe the same fixed detector and proposal. Averaging
    individual ARL ratios would give a different estimator.
    """
    if not batches:
        raise ValueError("need batches")
    for b in batches[1:]:
        if any(
            b.get(k) != batches[0].get(k)
            for k in [
                "kind",
                "n",
                "threshold",
                "parameters",
                "max_steps",
                "proposal",
                "proposal_offset",
            ]
        ):
            raise ValueError("batches must describe the same detector")
    nc, mean, se = pool_moments(
        [
            (b["pre_cycles"], b["mean_pre_cycle_length"], b["mean_pre_cycle_length_se"])
            for b in batches
        ]
    )
    nw, p, pse = pool_moments(
        [
            (b["post_cycles"], b["cycle_alarm_probability"], b["cycle_alarm_probability_se"])
            for b in batches
        ]
    )
    if p <= 0:
        raise RuntimeError("no weighted successes")
    arl = mean / p
    rse = float(np.hypot(se / mean, pse / p))
    ess = float(nw * p * p / ((nw - 1) * pse * pse + p * p))
    largest_log = max(b["max_log_success_weight"] for b in batches)
    fraction = float(np.exp(largest_log - np.log(nw) - np.log(p)))
    out = {
        "arl": float(arl),
        "arl_se": float(arl * rse),
        "relative_se": rse,
        "approx_95_interval": [float(arl * np.exp(-1.96 * rse)), float(arl * np.exp(1.96 * rse))],
        "mean_pre_cycle_length": mean,
        "mean_pre_cycle_length_se": se,
        "cycle_alarm_probability": p,
        "cycle_alarm_probability_se": pse,
        "success_weight_ess": ess,
        "max_log_success_weight": largest_log,
        "largest_weight_fraction": fraction,
        "diagnostics_passed": bool(ess >= 1000 and fraction <= 0.01 and rse <= 0.02),
        "pre_cycles": nc,
        "post_cycles": nw,
        "longest_pre_cycle": max(b["longest_pre_cycle"] for b in batches),
        "longest_post_cycle": max(b["longest_post_cycle"] for b in batches),
        "uncertainty_method": "pooled independent cycle moments; delta SE; log-scale approximate 95% interval",
    }
    for key in [
        "pre_observed_alarm_cycles",
        "post_proposal_successes",
        "pre_pulses_simulated",
        "post_pulses_simulated",
        "unfinished_pre_cycles",
        "unfinished_post_cycles",
    ]:
        out[key] = sum(b[key] for b in batches)
    if out["unfinished_pre_cycles"] or out["unfinished_post_cycles"]:
        raise RuntimeError("unfinished cycles cannot be pooled")
    for key in ["kind", "n", "threshold", "parameters", "max_steps", "proposal", "proposal_offset"]:
        out[key] = batches[0].get(key)
    return out


def batched_estimate(
    kind,
    n,
    threshold,
    Na,
    cycles=1000000,
    batch_size=100000,
    seed=20261022,
    audit_path=None,
    proposal="post",
    proposal_offset=2.0,
):
    """Run independent child-seeded batches and optionally retain the first raw audit batch."""
    if cycles < 2 or batch_size < 2 or cycles % batch_size:
        raise ValueError("cycles must be a positive multiple of batch_size, both >=2")
    batches = []
    for i, child in enumerate(np.random.SeedSequence(seed).spawn(cycles // batch_size)):
        child_seed = int(child.generate_state(1, dtype=np.uint64)[0])
        stats, pre, post = estimate_conditional_arl(
            kind,
            n,
            threshold,
            cycles=batch_size,
            seed=child_seed,
            Na=Na,
            proposal=proposal,
            proposal_offset=proposal_offset,
        )
        if i == 0 and audit_path is not None:
            np.savez_compressed(audit_path, **conditional_cycle_arrays(pre, post))
        batches.append(stats)
    return {
        "seed": seed,
        "batch_size": batch_size,
        "batches": batches,
        "summary": pool_cycle_summaries(batches),
    }
