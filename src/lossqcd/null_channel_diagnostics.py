"""Scalar sensitivity to true pre-change loss at frozen source thresholds.

No published latency is fitted. Both physical post-change and assumed detector
eta2 stay 0.85. Only true eta1 varies. Fresh numerical and renewal ARL estimates
are mandatory; matched-null ARL evidence is not reused for mismatched cases.
"""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .models import coherent_pair, squeezed_pair
from .arl_solver import mean_stopping_time
from .delay_solver import conditional_delay_numerical
from .mismatch_renewal import estimate_mismatch_arl, cycle_arrays
from .simulation import simulate_stopping_times, summarize_delay, summarize_arl


def physical_laws(kind, eta):
    """Return assumed score laws and the true initial-channel law at fixed post-change loss."""
    pair = coherent_pair if kind == "coherent" else squeezed_pair if kind == "entangled" else None
    if pair is None:
        raise ValueError("unknown scalar model")
    pre, post = pair(N=100, Na=0.1, eta1=0.9, eta2=0.85)
    actual, _ = pair(N=100, Na=0.1, eta1=eta, eta2=0.85)
    return pre, post, actual


def selected_models(source):
    """Select coherent and independently squeezed Na=0.1 source controls."""
    return [
        m for m in source["models"] if m["figure"] == 7 and m["n"] == 1 and np.isclose(m["Na"], 0.1)
    ]


def numerical_checks(pre, post, actual, threshold, change_at, arl_cells=(1024, 2048)):
    """Refine actual-null ARL and conditional-delay grids without changing the threshold."""
    arl = [mean_stopping_time(pre, post, threshold, cells=c, sampling=actual) for c in arl_cells]
    delay = [
        conditional_delay_numerical(
            pre, post, threshold, change_at=change_at, cells=c, sampling_pre=actual
        )
        for c in (512, 1024)
    ]
    return arl, delay


def run(source_root, output, cycles=200000, trials=30000, seed=2026102801, refine_from=None):
    """Re-estimate false alarms and changed delays over a prespecified true-initial-channel grid.

    Optional refinement reuses exactly the earlier raw ensembles; it changes
    only the numerical mesh checks, not the sampled evidence.
    """
    source_root, output = Path(source_root).resolve(), Path(output).resolve()
    if not isinstance(cycles, int) or cycles < 2 or not isinstance(trials, int) or trials < 2:
        raise ValueError("cycles and trials must be integers >=2")
    source = json.loads((source_root / "report.json").read_text())
    models = selected_models(source)
    s = source["settings"]
    if (
        len(models) != 2
        or (s["N"], s["eta1"], s["eta2"]) != (100, 0.9, 0.85)
        or not all(m["confirmed"] for m in models)
    ):
        raise ValueError(
            "expected the two confirmed Na=0.1 scalar source models at N=100, eta=.9/.85"
        )
    hashes = {str(source_root / "report.json"): sha256(source_root / "report.json")}
    for model in models:
        p = source_root / Path(model["folder"]).name / "changed_runs.npz"
        hashes[str(p)] = sha256(p)
    prior = None
    if refine_from is not None:
        prior = Path(refine_from).resolve()
        previous = json.loads((prior / "report.json").read_text())
        ps = previous["settings"]
        if (ps["seed"], ps["cycles_per_law_per_repeat"], ps["trials"]) != (seed, cycles, trials):
            raise ValueError("refinement must retain the original seeds and sample sizes")
        hashes[str(prior / "report.json")] = sha256(prior / "report.json")
        for path in sorted(prior.glob("*.npz")):
            hashes[str(path)] = sha256(path)
    output.mkdir(parents=True, exist_ok=False)
    etas = [0.899, 0.8995, 0.9, 0.9005, 0.901]
    report = {
        "source_root": str(source_root),
        "source_hashes": hashes,
        "settings": {
            "N": 100,
            "Na": 0.1,
            "assumed_eta1": 0.9,
            "eta2": 0.85,
            "true_eta1_grid": etas,
            "cycles_per_law_per_repeat": cycles,
            "trials": trials,
            "seed": seed,
            "horizon": 5000,
        },
        "rows": [],
        "direct_checks": [],
        "scope": "Scalar coherent and independent squeezed only. No author mean fitted. "
        "No-change law changes with actual eta1, requiring fresh ARL validation. "
        "This is robustness analysis, not an explanation inferred for author data.",
    }
    report["refinement"] = {
        "prior_root": str(prior) if prior is not None else None,
        "arl_cells": [1024, 2048],
        "delay_cells": [512, 1024],
        "reason": "Initial 512/1024 ARL grid changed by 0.124–0.135% for coherent cases, above the fixed 0.1% tolerance. "
        "Refine the numerical grid while retaining all seeds, raw samples, physical parameters and thresholds."
        if prior is not None
        else None,
    }
    streams = np.random.SeedSequence(seed).spawn(40)
    make_seed = lambda i: int(streams[i].generate_state(1, dtype=np.uint64)[0])
    index = 0
    for model in models:
        for eta in etas:
            pre, post, actual = physical_laws(model["kind"], eta)
            h = model["calibration"]["threshold"]
            nu = model["change_at"]
            arl, delay = numerical_checks(pre, post, actual, h, nu)
            repeats = []
            files = []
            seeds = []
            for repeat in range(2):
                run_seed = make_seed(3 * index + repeat)
                stats, num, den = estimate_mismatch_arl(pre, post, actual, h, cycles, run_seed)
                filename = f"cycles_{index:02d}_{repeat}.npz"
                np.savez_compressed(output / filename, **cycle_arrays(num, den))
                repeats.append(stats)
                files.append(filename)
                seeds.append(run_seed)
            changed_seed = make_seed(3 * index + 2)
            changed = simulate_stopping_times(
                pre,
                post,
                h,
                trials=trials,
                horizon=5000,
                change_at=nu,
                seed=changed_seed,
                sampling_pre=actual,
            )
            filename = f"changed_{index:02d}.npz"
            np.savez_compressed(output / filename, stop=changed.stop, event=changed.event)
            stats = summarize_delay(changed)
            z = abs(repeats[0]["arl"] - repeats[1]["arl"]) / np.hypot(
                repeats[0]["arl_se"], repeats[1]["arl_se"]
            )
            arl_mesh = abs(arl[1]["mean_stopping_time"] / arl[0]["mean_stopping_time"] - 1)
            delay_mesh = abs(delay[1]["conditional_delay"] - delay[0]["conditional_delay"])
            passed = (
                arl_mesh < 0.001
                and delay_mesh < 0.05
                and z <= 4
                and all(
                    r["diagnostics_passed"]
                    and abs(r["arl"] - arl[1]["mean_stopping_time"])
                    <= 4 * r["arl_se"]
                    + abs(arl[1]["mean_stopping_time"] - arl[0]["mean_stopping_time"])
                    for r in repeats
                )
            )
            passed = (
                passed
                and stats["conditional_delay"] is not None
                and abs(stats["conditional_delay"] - delay[1]["conditional_delay"])
                <= 4 * stats["conditional_delay_se"] + delay_mesh
            )
            report["rows"].append(
                {
                    "kind": model["kind"],
                    "label": model["label"],
                    "true_eta1": eta,
                    "threshold": h,
                    "change_at": nu,
                    "arl_numerical": arl,
                    "delay_numerical": delay,
                    "arl_mesh_relative_change": arl_mesh,
                    "delay_mesh_change": delay_mesh,
                    "repeat_z": float(z),
                    "repeats": repeats,
                    "cycle_files": files,
                    "seeds": seeds,
                    "changed_file": filename,
                    "changed_seed": changed_seed,
                    "changed": stats,
                    "checks_passed": bool(passed),
                }
            )
            print(
                f"{index + 1}/10 {model['label']} eta1={eta:g}: ARL={arl[1]['mean_stopping_time']:.0f}, delay={stats['conditional_delay']:.3f}",
                flush=True,
            )
            index += 1
    for model in models:
        for eta in (etas[0], etas[-1]):
            pre, post, actual = physical_laws(model["kind"], eta)
            i = len(report["direct_checks"])
            run_seed = make_seed(30 + 2 * i)
            direct = simulate_stopping_times(
                pre, post, 3, trials=10000, horizon=100000, seed=run_seed, sampling_pre=actual
            )
            direct_stats = summarize_arl(direct)
            estimate_seed = make_seed(31 + 2 * i)
            est, num, den = estimate_mismatch_arl(pre, post, actual, 3, cycles, estimate_seed)
            num_arl = mean_stopping_time(pre, post, 3, cells=512, sampling=actual)[
                "mean_stopping_time"
            ]
            filename = f"control_{i}.npz"
            np.savez_compressed(
                output / filename, stop=direct.stop, event=direct.event, **cycle_arrays(num, den)
            )
            passed = (
                direct_stats["arl"] is not None
                and abs(direct_stats["arl"] - num_arl) <= 4 * direct_stats["arl_se"]
            )
            passed = (
                passed
                and est["diagnostics_passed"]
                and abs(est["arl"] - num_arl) <= 4 * est["arl_se"]
            )
            report["direct_checks"].append(
                {
                    "kind": model["kind"],
                    "true_eta1": eta,
                    "threshold": 3.0,
                    "direct_seed": run_seed,
                    "estimate_seed": estimate_seed,
                    "file": filename,
                    "direct": direct_stats,
                    "estimate": est,
                    "numerical_arl": num_arl,
                    "checks_passed": bool(passed),
                }
            )
    report["source_unchanged"] = all(sha256(p) == h for p, h in hashes.items())
    if prior is not None:
        report["refinement"]["all_raw_samples_identical"] = all(
            sha256(path) == sha256(output / path.name) for path in prior.glob("*.npz")
        )
        report["source_unchanged"] = (
            report["source_unchanged"] and report["refinement"]["all_raw_samples_identical"]
        )
    report["all_checks_passed"] = report["source_unchanged"] and all(
        r["checks_passed"] for r in report["rows"] + report["direct_checks"]
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    write_outputs(report, output)
    return report


def write_outputs(report, output):
    """Save actual-null sensitivity, renewal checks and direct moderate-ARL controls."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = report["rows"]
    exported = [
        {
            "label": r["label"],
            "true_eta1": r["true_eta1"],
            "threshold": r["threshold"],
            "numerical_arl": r["arl_numerical"][1]["mean_stopping_time"],
            "arl_repeat0": r["repeats"][0]["arl"],
            "arl_se0": r["repeats"][0]["arl_se"],
            "arl_repeat1": r["repeats"][1]["arl"],
            "arl_se1": r["repeats"][1]["arl_se"],
            "delay": r["changed"]["conditional_delay"],
            "delay_se": r["changed"]["conditional_delay_se"],
            "early_alarm_probability": r["delay_numerical"][1]["false_alarm_probability"],
        }
        for r in rows
    ]
    with (output / "sensitivity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(exported[0]))
        writer.writeheader()
        writer.writerows(exported)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    for kind in ("coherent", "entangled"):
        subset = [r for r in rows if r["kind"] == kind]
        x = [r["true_eta1"] for r in subset]
        color = "C0" if kind == "coherent" else "C1"
        label = subset[0]["label"]
        axes[0].plot(
            x,
            [r["arl_numerical"][1]["mean_stopping_time"] / 1e6 for r in subset],
            color=color,
            label=label,
        )
        for repeat, marker in ((0, "o"), (1, "x")):
            axes[0].errorbar(
                x,
                [r["repeats"][repeat]["arl"] / 1e6 for r in subset],
                yerr=[1.96 * r["repeats"][repeat]["arl_se"] / 1e6 for r in subset],
                fmt=marker,
                color=color,
                capsize=3,
            )
        axes[1].plot(
            x,
            [r["delay_numerical"][1]["conditional_delay"] for r in subset],
            color=color,
            label=label,
        )
        axes[1].errorbar(
            x,
            [r["changed"]["conditional_delay"] for r in subset],
            yerr=[1.96 * r["changed"]["conditional_delay_se"] for r in subset],
            fmt="o",
            color=color,
            capsize=3,
        )
    axes[0].axhline(2, color="black", ls="--", label="Original ARL target")
    axes[0].set(
        ylabel="No-change ARL (million pulses)", title="Fresh numerical and renewal estimates"
    )
    axes[1].set(ylabel="Conditional delay (pulses)", title="True post-change eta2 stays 0.85")
    for ax in axes:
        ax.axvline(0.9, color="gray", ls=":")
        ax.set_xlabel("True pre-change transmissivity eta1")
        ax.set_xticks(
            report["settings"]["true_eta1_grid"],
            [f"{eta:.4f}" for eta in report["settings"]["true_eta1_grid"]],
        )
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    fig.suptitle(
        "Initial-channel sensitivity: fixed detector assumptions and original thresholds\nN=100, Na=0.1; lines: numerical; points: Monte Carlo ±1.96 SE",
        fontsize=11,
    )
    for ext in ("png", "svg"):
        fig.savefig(output / f"null_channel_sensitivity.{ext}", dpi=180)
    plt.close(fig)
    lines = []
    for kind in ("coherent", "entangled"):
        subset = [r for r in rows if r["kind"] == kind]
        low, center, high = subset[0], subset[2], subset[-1]
        lines.append(
            f"| {low['label']} | {low['arl_numerical'][1]['mean_stopping_time']:,.0f} | "
            f"{center['arl_numerical'][1]['mean_stopping_time']:,.0f} | {high['arl_numerical'][1]['mean_stopping_time']:,.0f} |"
        )
    table = "\n".join(lines)
    refinement = report["refinement"]
    provenance = (
        f"The initial failed 512/1024-grid run remains preserved in `{refinement['prior_root']}`. "
        f"All raw samples are identical: **{refinement['all_raw_samples_identical']}**. "
        "Only the ARL grid changed to 1024/2048 cells, retaining the same 0.1% tolerance."
        if refinement["prior_root"] is not None
        else "No prior run was supplied; the ARL grid uses 1024/2048 cells."
    )
    (output / "RESULTS.md").write_text(f"""# Initial-channel sensitivity

Detector assumptions stay eta1=0.9, eta2=0.85; N=100, Na=0.1. True eta1 varies
from 0.899 to 0.901 (89.9% to 90.1% transmission). The true post-change eta2
stays 0.85. Original source thresholds stay fixed. No author mean is fitted,
and the earlier post-change diagnostic fit is not used in this separate study.

Fresh numerical no-change ARLs, in pulses:

| Detector | True eta1=0.899 | True eta1=0.900 | True eta1=0.901 |
| --- | --- | --- | --- |
{table}

Each of 10 scalar settings uses 1024/2048-cell ARL solutions and 512/1024-cell delay
solutions, two independent
{report["settings"]["cycles_per_law_per_repeat"]:,}-cycle ensembles per sampling law,
and a fresh {report["settings"]["trials"]:,}-trial changed ensemble. Four moderate-threshold
controls also compare direct no-change simulation against renewal estimates
and the numerical solution. All checks pass: **{report["all_checks_passed"]}**.

The renewal numerator samples the actual null law. The independent denominator
uses the detector's assumed post law as proposal, with full path log weight
sum(log p_actual_null − log p_proposal). The assumed detector score is stored
separately. Substituting exp(−sum(assumed score)) would be wrong under mismatch.
No likelihood weights are clipped or self-normalized; unfinished cycles are
rejected. ESS, largest-weight fraction and relative-SE diagnostics are recorded.
Intervals and mesh refinements are approximate diagnostics, not rigorous bounds.

Changed delay is E[T−1005 | T>=1005], with early first alarms reported separately.
Censoring is explicit. Both initial-state survival and post-change timing are
checked numerically. The earlier matched-null ARL evidence is not substituted
for fresh estimates at the mismatched initial channel.

This scalar robustness study does not determine the authors' actual channel or
resolve the publication's latency discrepancy. Correlated entangled blocks are
outside its scope; they need distinct block-boundary renewals and corrected
conditional path weights under their actual null law. Original report and raw
source runs are checksum-protected. Seeds, full renewal arrays, changed runs,
control runs, CSV and plots are saved here.

{provenance}

Next: extend actual-null weighting to the correlated block detector and validate
it independently before drawing conclusions about entangled robustness.
""")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="results_author_comparison")
    parser.add_argument("--output", default="results_null_channel_diagnostics")
    parser.add_argument("--cycles", type=int, default=200000)
    parser.add_argument("--trials", type=int, default=30000)
    parser.add_argument("--seed", type=int, default=2026102801)
    parser.add_argument("--refine-from", default=None)
    args = parser.parse_args()
    run(args.source, args.output, args.cycles, args.trials, args.seed, args.refine_from)


if __name__ == "__main__":
    main()
