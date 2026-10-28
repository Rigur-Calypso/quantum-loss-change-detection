"""Calibrate and independently validate Appendix B.4 ARL in pulse units.

The saved default experiment is a Na=0.1 slice, not the full Figure 7 sweep.
Thresholds freeze after calibration. Validation never feeds back into them.
"""

import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .corr_renewal import (
    estimate_conditional_arl,
    conditional_cycle_arrays,
    summarize_conditional_cycles,
)
from .corr_cusum import simulate_conditional
from .simulation import summarize_delay
from .calibration import independent_ratio


def calibrate_rare(kind, n, target=2e6, cycles=500000, seed=20261021, iterations=4, Na=0.1):
    """Fixed-budget log-ARL correction; a fresh holdout is mandatory.

    ARL is locally approximately exponential in h in this rare-event regime.
    The unit-slope update is a search heuristic, not an analytical ARL claim.
    """
    if not np.isfinite(target) or target <= 1:
        raise ValueError("target must be finite and >1")
    if not isinstance(iterations, int) or iterations < 1:
        raise ValueError("iterations must be positive")
    h = 12.0
    history = []
    for child in np.random.SeedSequence(seed).spawn(iterations):
        child_seed = int(child.generate_state(1, dtype=np.uint64)[0])
        stats, _, _ = estimate_conditional_arl(kind, n, h, cycles=cycles, seed=child_seed, Na=Na)
        history.append(stats)
        h = max(0.0, h + float(np.clip(np.log(target / stats["arl"]), -2.0, 2.0)))
    return {
        "threshold": h,
        "history": history,
        "cycles_per_ensemble": cycles,
        "seed": seed,
        "target_arl_pulses": target,
        "method": "fixed-budget log-ARL corrections; threshold frozen before holdout",
    }


def summarize_validation(pre_batches, post_batches):
    """Concatenate independent renewal arrays before recomputing the ARL ratio."""
    from dataclasses import replace

    def combine(batches):
        return replace(
            batches[0],
            **{
                field: np.concatenate([getattr(b, field) for b in batches])
                for field in ["length", "raw_sum", "page", "upper", "complete"]
            },
        )

    return summarize_conditional_cycles(combine(pre_batches), combine(post_batches))


def assess_match(stats, baseline, target, tolerance=0.05):
    """Require target and relative-ARL intervals inside tolerance, with weight diagnostics."""
    ratio = independent_ratio(stats["arl"], stats["arl_se"], baseline["arl"], baseline["arl_se"])
    low, high = ratio["approx_95_interval"]
    # Require each target interval as well as the inter-detector interval inside band.
    on_target = all(
        (1 - tolerance) * target <= s["approx_95_interval"][0]
        and s["approx_95_interval"][1] <= (1 + tolerance) * target
        for s in (stats, baseline)
    )
    return ratio | {
        "passed": bool(
            on_target
            and 1 - tolerance <= low
            and high <= 1 + tolerance
            and stats["diagnostics_passed"]
            and baseline["diagnostics_passed"]
        )
    }


def write_outputs(root, report):
    """Write the pulse-ARL holdouts, moderate controls and conditional-delay comparisons."""
    models = report["models"]
    rows = []
    baseline = models[0]
    for m in models:
        for nu, d in m["delay"].items():
            b = baseline["delay"][nu]
            gain = (
                independent_ratio(
                    b["conditional_delay"],
                    b["conditional_delay_se"],
                    d["conditional_delay"],
                    d["conditional_delay_se"],
                )
                if m is not baseline
                else {"ratio": 1.0, "se": 0.0}
            )
            rows.append(
                {
                    "label": m["label"],
                    "n": m["n"],
                    "threshold": m["calibration"]["threshold"],
                    "arl": m["validation"]["arl"],
                    "arl_se": m["validation"]["arl_se"],
                    "match_passed": m["arl_match"]["passed"],
                    "change_at": int(nu),
                    "change_phase": (int(nu) - 1) % m["n"] + 1,
                    "delay": d["conditional_delay"],
                    "delay_se": d["conditional_delay_se"],
                    "early_alarm_probability": d["false_alarm_probability"],
                    "survivors": d["survivors"],
                    "coherent_delay_ratio": gain["ratio"],
                    "ratio_se": gain["se"],
                }
            )
    report["summary"] = rows
    report["all_checks_passed"] = all(
        m["arl_match"]["passed"] and m["replicates_agree"] and m["moderate_check"]["passed"]
        for m in models
    )
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    with (root / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    x = np.arange(len(models))
    labels = [m["label"] for m in models]
    axes[0].errorbar(
        x,
        [m["validation"]["arl"] / 1e6 for m in models],
        yerr=[1.96 * m["validation"]["arl_se"] / 1e6 for m in models],
        fmt="o",
        capsize=4,
    )
    axes[0].axhline(report["settings"]["target_arl"] / 1e6, color="black", ls="--")
    axes[0].set(ylabel="No-change ARL (million pulses)", title="Independent weighted Monte Carlo")
    for nu in report["settings"]["change_at"]:
        selected = [r for r in rows if r["change_at"] == nu]
        axes[1].errorbar(
            x,
            [r["delay"] for r in selected],
            yerr=[1.96 * r["delay_se"] for r in selected],
            fmt="o-",
            capsize=3,
            label=f"Change on pulse {nu}",
        )
    axes[1].set(
        ylabel="Conditional detection delay (pulses)", title="Pulse-by-pulse conditional detector"
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x, labels, rotation=25)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(f"N=100, Na={report['settings']['Na']}; bars: approximate 95% MC intervals")
    fig.savefig(root / "corr_arl.png", dpi=180)
    fig.savefig(root / "corr_arl.svg")
    plt.close(fig)
    lines = [
        "# Correlated conditional detector at large false-alarm ARL",
        "",
        "Source model and detector: [paper Appendix B.4](https://arxiv.org/html/2503.12276v2#A2.SS4).",
        "The regenerative extension and its validation here are our implementation work; they are not a method reported in the optical paper. General background: [Calvin, Glynn & Nakayama (2001)](https://web.stanford.edu/~glynn/papers/2001/CalvinGNakayama01.pdf).",
        "",
        f"N=100, Na={report['settings']['Na']}, eta1=0.9, eta2=0.85. Target ARL={report['settings']['target_arl']:,.0f} pulses. All checks passed: **{report['all_checks_passed']}**.",
        "",
        "| Detector | Threshold | ARL ± SE (pulses) | Matched | Repeat z | Moderate check z |",
        "|---|---:|---:|---|---:|---:|",
    ]
    for m in models:
        lines.append(
            f"| {m['label']} | {m['calibration']['threshold']:.6f} | {m['validation']['arl']:,.0f} ± {m['validation']['arl_se']:,.0f} | {m['arl_match']['passed']} | {m['replicate_z']:.2f} | {m['moderate_check']['z']:.2f} |"
        )
    lines += [
        "",
        "## Detection delays",
        "",
        "Each run stops at its first alarm. Delay is E[T−nu | T>=nu]; early alarms are reported separately and never counted as detections. Pulse nu is the first changed pulse; block phase is (nu−1) mod n + 1. No changed-run delay in this report is censored.",
        "",
        "| Detector | Change pulse | Phase | Delay ± SE | Early alarms | Coherent/detector delay |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['label']} | {r['change_at']} | {r['change_phase']} | {r['delay']:.3f} ± {r['delay_se']:.3f} | {100 * r['early_alarm_probability']:.3f}% | {r['coherent_delay_ratio']:.3f} |"
        )
    lines += [
        "",
        "## Why the rare-event calculation is different",
        "",
        "A cycle starts with zero Page evidence and an empty block prefix. It ends either at a strict alarm on any pulse, including inside a block, or at a block end with Page evidence exactly zero. Only that block-end zero renews the full state. Reflections inside a block retain the measured prefix. Blocks are independent under the homogeneous no-change law, so renewal gives ARL = mean pre-change cycle length / pre-change cycle alarm probability.",
        "",
        "The numerator uses pre-change cycles. An independent post-change proposal estimates the denominator with ordinary weights exp(−sum of ALL raw conditional LLRs), multiplied by the cycle-alarm indicator. The sum includes increments before Page reflections and stops at the actual alarm pulse. The weight is not exp(−Page), not a complete-block likelihood at a partial-block alarm, and not normalized by the sum of weights. The ratio estimator is consistent, not exactly unbiased for finite ensembles.",
        "",
        "Post-proposal weights can be heavy-tailed; no iid exp(−h) bound is assumed. We report ESS, largest weight fraction, relative SE and two independent repetitions. Screening requires ESS>=1000, largest contribution<=1%, and relative SE<=2%. Repeats must agree within 4 combined SE; the new method must also agree within 4 combined SE with an independent direct simulation at the earlier approximately 1000-pulse setting. These diagnostics do not prove finite variance or rule out unobserved tail contributions.",
        "",
        f"Thresholds use {report['settings']['calibration_cycles']:,} cycles per ensemble per search evaluation (three times as many for n>=8), with four fixed log-ARL corrections. They freeze before two independent holdout ensembles, each using {report['settings']['validation_cycles']:,} cycles per law (three times as many for n>=8). Validation never retunes thresholds. Target and inter-detector approximate 95% intervals must lie fully inside ±{100 * report['settings']['tolerance']:g}%. All intervals are pointwise Monte Carlo approximations conditional on thresholds, not rigorous bounds or model uncertainty.",
        "",
        "Raw renewal-cycle lengths, raw sums, final Page scores and event/completion flags are saved separately for both validation repetitions and both laws. Changed-run stopping times and event flags, calibration histories, random seeds, model settings, versions, CSV and plots are saved alongside. Output folders are never overwritten.",
        "",
        "This is a single Na=0.1 slice with the five detector types in Figure 7, not the full figure or a pixel-level reproduction. The physical model remains ideal homodyne with known losses and no excess receiver noise; the conditional detector is not claimed optimal. One-based indexing explicitly fixes the phase; pulse 1005 is phase 1 for n=2 and n=4 and phase 5 for n=8.",
        "",
        "Next: sweep Na and change phase using these matched-ARL checks. Larger blocks (such as n=16) need improved proposal design or substantially more samples: exploratory post-proposal pilots had unstable weights, so they are excluded from this confirmed comparison.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("results_conditional_large_arl"))
    p.add_argument("--reference", type=Path, default=Path("results_entangled"))
    p.add_argument("--seed", type=int, default=2026102160)
    p.add_argument("--Na", type=float, default=0.1)
    p.add_argument("--target-arl", type=float, default=2e6)
    p.add_argument("--calibration-cycles", type=int, default=500000)
    p.add_argument("--validation-cycles", type=int, default=1000000)
    p.add_argument("--delay-trials", type=int, default=50000)
    p.add_argument("--change-at", type=int, nargs="+", default=[1, 1005])
    p.add_argument("--horizon", type=int, default=100000)
    p.add_argument("--tolerance", type=float, default=0.05)
    args = p.parse_args(argv)
    if min(args.calibration_cycles, args.validation_cycles, args.delay_trials) < 2:
        p.error("ensemble sizes must be >=2")
    if (
        not np.isfinite(args.target_arl)
        or args.target_arl <= 1
        or not np.isfinite(args.Na)
        or args.Na < 0
    ):
        p.error("target ARL must be finite and >1; Na finite and nonnegative")
    if (
        not 0 < args.tolerance < 1
        or any(not 1 <= nu <= args.horizon for nu in args.change_at)
        or len(set(args.change_at)) != len(args.change_at)
    ):
        p.error("invalid tolerance or change points")
    # Validate source settings before creating outputs; never silently compare mismatched Na.
    references = [
        json.loads((args.reference / f"model_{i:02d}" / "report.json").read_text())
        for i in range(5)
    ]
    choices = [("coherent", 1, "Coherent"), ("entangled", 1, "Squeezed n=1")] + [
        ("entangled", n, f"Entangled n={n}") for n in [2, 4, 8]
    ]
    for m, (kind, n, _) in zip(references, choices):
        if (m["kind"], m["n"]) != (kind, n) or m["calibration"]["parameters"] != {
            "N": 100.0,
            "Na": args.Na,
            "eta1": 0.9,
            "eta2": 0.85,
        }:
            p.error("direct reference detector/model parameters differ")
        if m["validation"]["arl"] is None or m["validation"]["censored"]:
            p.error("direct reference has censored ARL")
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "settings": vars(args) | {"output": str(args.output), "reference": str(args.reference)},
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {key: importlib.metadata.version(key) for key in ["numpy", "matplotlib"]},
        "models": [],
    }
    streams = np.random.SeedSequence(args.seed).spawn(5 * (4 + len(args.change_at)))
    for i, (kind, n, label) in enumerate(choices):
        count = 4 + len(args.change_at)
        seeds = [
            int(child.generate_state(1, dtype=np.uint64)[0])
            for child in streams[i * count : (i + 1) * count]
        ]
        folder = args.output / f"model_{i:02d}"
        folder.mkdir()
        factor = 3 if n >= 8 else 1
        print(f"{label}: independent moderate-ARL crosscheck", flush=True)
        reference = references[i]
        moderate, _, _ = estimate_conditional_arl(
            kind,
            n,
            reference["calibration"]["threshold"],
            cycles=args.validation_cycles * factor,
            seed=seeds[0],
            Na=args.Na,
        )
        direct = reference["validation"]
        z = abs(moderate["arl"] - direct["arl"]) / np.hypot(moderate["arl_se"], direct["arl_se"])
        check = {
            "estimate": moderate,
            "direct": direct,
            "z": float(z),
            "passed": bool(z <= 4 and moderate["diagnostics_passed"]),
        }
        print(f"{label}: calibrating, moderate z={z:.2f}", flush=True)
        calibration = calibrate_rare(
            kind, n, args.target_arl, args.calibration_cycles * factor, seeds[1], Na=args.Na
        )
        h = calibration["threshold"]
        repeats = []
        pre_batches = []
        post_batches = []
        for j in range(2):
            stats, pre, post = estimate_conditional_arl(
                kind, n, h, cycles=args.validation_cycles * factor, seed=seeds[2 + j], Na=args.Na
            )
            np.savez_compressed(
                folder / f"validation_{j}.npz", **conditional_cycle_arrays(pre, post)
            )
            repeats.append(stats)
            pre_batches.append(pre)
            post_batches.append(post)
            print(
                f"  holdout {j + 1}: ARL {stats['arl']:,.0f} +/- {stats['arl_se']:,.0f}, ESS {stats['success_weight_ess']:,.0f}",
                flush=True,
            )
        combined = summarize_validation(pre_batches, post_batches)
        repeat_z = abs(repeats[0]["arl"] - repeats[1]["arl"]) / np.hypot(
            repeats[0]["arl_se"], repeats[1]["arl_se"]
        )
        del pre_batches, post_batches, pre, post
        delay = {}
        raw = {}
        for j, nu in enumerate(args.change_at):
            runs = simulate_conditional(
                kind,
                n,
                Na=args.Na,
                threshold=h,
                trials=args.delay_trials,
                horizon=args.horizon,
                change_at=nu,
                seed=seeds[4 + j],
            )
            delay[str(nu)] = summarize_delay(runs)
            if delay[str(nu)]["conditional_delay"] is None:
                raise RuntimeError(
                    "censored/no-survivor changed runs; unrestricted delay unavailable"
                )
            raw[f"nu{nu}_stop"], raw[f"nu{nu}_event"] = runs.stop, runs.event
        np.savez_compressed(folder / "delay_runs.npz", **raw)
        match = (
            assess_match(
                combined, report["models"][0]["validation"], args.target_arl, args.tolerance
            )
            if i
            else {
                "ratio": 1.0,
                "approx_95_interval": [1.0, 1.0],
                "passed": bool(
                    combined["diagnostics_passed"]
                    and (1 - args.tolerance) * args.target_arl <= combined["approx_95_interval"][0]
                    and combined["approx_95_interval"][1] <= (1 + args.tolerance) * args.target_arl
                ),
            }
        )
        model = {
            "kind": kind,
            "n": n,
            "label": label,
            "seeds": seeds,
            "moderate_check": check,
            "calibration": calibration,
            "validation_repeats": repeats,
            "validation": combined,
            "replicate_z": float(repeat_z),
            "replicates_agree": bool(
                repeat_z <= 4 and all(s["diagnostics_passed"] for s in repeats)
            ),
            "arl_match": match,
            "delay": delay,
        }
        report["models"].append(model)
        (folder / "report.json").write_text(json.dumps(model, indent=2) + "\n")
    write_outputs(args.output, report)
    print(f"Completed. All predeclared checks passed: {report['all_checks_passed']}", flush=True)


if __name__ == "__main__":
    main()
