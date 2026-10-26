"""Monte Carlo threshold calibration; holdout validation belongs to the caller."""
import numpy as np
from .simulation import simulate_stopping_times, summarize_arl


def calibrate_threshold(pre, post, target_arl=1000., trials=20000,
                        horizon=100000, seed=20261021, iterations=10):
    """Noisy bisection in expected ARL, with an independent stream per evaluation.

    The search is approximate: a narrow threshold bracket does not imply a
    precise ARL. The selected threshold must be checked on a fresh ensemble.
    Any censored calibration run raises an error rather than biasing the ARL.
    """
    if not np.isfinite(target_arl) or target_arl <= 1:
        raise ValueError('target_arl must be finite and greater than 1')
    for name, value in [('trials', trials), ('horizon', horizon), ('iterations', iterations)]:
        if not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError(f'{name} must be an integer >=2')
    if post.kl(pre) <= 0:
        raise ValueError('identical distributions cannot be calibrated')
    streams = np.random.SeedSequence(seed).spawn(iterations + 24)
    history = []

    def evaluate(h):
        runs = simulate_stopping_times(pre, post, h, trials, horizon,
                                       seed=streams[len(history)])
        stats = summarize_arl(runs)
        if stats['arl'] is None:
            raise RuntimeError('censored calibration runs: increase horizon')
        history.append({'threshold': h, **stats})
        return stats['arl']

    low, high = 0., 5.
    if evaluate(low) >= target_arl:
        raise ValueError('target ARL is below the attainable nonnegative-threshold range')
    for _ in range(20):
        if evaluate(high) >= target_arl:
            break
        high += 1.
    else:
        raise RuntimeError('could not bracket target ARL')
    for _ in range(iterations):
        mid = (low + high) / 2
        if evaluate(mid) < target_arl:
            low = mid
        else:
            high = mid
    return {'threshold': (low + high) / 2, 'target_arl': target_arl,
            'threshold_bracket': [low, high], 'history': history,
            'method': 'noisy Monte Carlo bisection; independent stream per evaluation',
            'calibration_trials_per_evaluation': trials,
            'calibration_seed': seed}


def independent_ratio(numerator, numerator_se, denominator, denominator_se):
    """Ratio and first-order MC SE for two independent positive mean estimates."""
    if any(not np.isfinite(x) for x in (numerator, numerator_se, denominator, denominator_se)):
        raise ValueError('ratio inputs must be finite')
    if numerator <= 0 or denominator <= 0 or numerator_se < 0 or denominator_se < 0:
        raise ValueError('means must be positive and SEs nonnegative')
    ratio = numerator / denominator
    se = ratio * np.sqrt((numerator_se/numerator)**2 + (denominator_se/denominator)**2)
    return {'ratio': float(ratio), 'se': float(se),
            'approx_95_interval': [float(ratio-1.96*se), float(ratio+1.96*se)]}
