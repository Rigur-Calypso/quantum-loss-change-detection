"""Controlled diagnostics of the latency gap; production results stay frozen.

Counterfactual thresholds are fitted only to ask what false-alarm ARL would
accompany an author mean in OUR model. They are not author thresholds, fair
ARL calibrations, or replacements for validated results.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from .models import coherent_pair, squeezed_pair
from .delay_solver import conditional_delay_numerical
from .arl_solver import mean_stopping_time
from .conditional_detection import simulate_conditional
from .phase_simulation import simulate_phases
from .simulation import Runs, summarize_delay
from .batched_rare_event import batched_estimate, pool_cycle_summaries


def delay_conventions(runs):
    """Compare zero-based, inclusive and conditional summaries of the same first alarms."""
    if runs.change_at is None or not np.all(runs.event):
        raise ValueError("convention comparison requires uncensored changed runs")
    d = runs.stop - runs.change_at
    survive = d >= 0
    if not np.any(survive):
        raise ValueError("conditional delay unavailable without survivors")
    conditional = float(d[survive].mean())
    return {
        "conditional": conditional,
        "conditional_minus_one": conditional - 1,
        "conditional_inclusive": conditional + 1,
        "signed_all_first_alarms": float(d.mean()),
        "clipped_all_first_alarms": float(np.maximum(d, 0).mean()),
        "early_alarm_fraction": float(np.mean(~survive)),
        "signed_shift_from_conditional": float(d.mean() - conditional),
    }


def equivalent_threshold(pre, post, target_delay, upper, change_at=1005, cells=256):
    """Fit a counterfactual threshold to one latency, exposing its implied false-alarm cost."""
    if not np.isfinite(target_delay) or target_delay < 0:
        raise ValueError("target delay must be finite and nonnegative")
    history = []

    def objective(h):
        x = conditional_delay_numerical(pre, post, h, change_at, cells)
        history.append(x)
        return x["conditional_delay"] - target_delay

    if objective(0) >= 0 or objective(upper) <= 0:
        raise ValueError("author mean is not bracketed by diagnostic thresholds")
    h = float(brentq(objective, 0, upper, xtol=1e-8))
    return {
        "threshold": h,
        "author_delay_target": target_delay,
        "cells": cells,
        "history": history,
        "purpose": "counterfactual inverse delay in our fixed model; not an estimate of the author threshold",
    }


def find_model(report, Na, kind, n=1):
    """Locate one saved Figure 7 detector by augmentation energy, family and block size."""
    return next(
        m
        for m in report["models"]
        if m["figure"] == 7
        and m["kind"] == kind
        and m["n"] == n
        and np.isclose(m["Na"], Na, rtol=0, atol=1e-12)
    )


def write_outputs(root, report):
    """Save alternative delay conventions and fitted-threshold diagnostics as separate evidence."""
    for name, key in [
        ("conventions.csv", "conventions"),
        ("deterministic_checks.csv", "deterministic_checks"),
    ]:
        rows = report[key]
        with (root / name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    rows = report["conventions"]
    report["max_signed_convention_shift"] = max(
        abs(r["signed_shift_from_conditional"]) for r in rows
    )
    report["convention_closes_any_gap"] = any(
        abs(r["closest_difference"]) <= 0.05 * r["author_delay"] for r in rows
    )
    report["all_deterministic_checks_passed"] = all(
        r["passed"] for r in report["deterministic_checks"]
    )
    report["counterfactual_validation_passed"] = all(
        m["checks_passed"] for m in report["counterfactuals"]
    )
    report["all_counterfactual_arl_intervals_below_target"] = all(
        m["validation"]["approx_95_interval"][1] < 1.9e6 for m in report["counterfactuals"]
    )
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    selected = [r for r in rows if np.isclose(r["Na"], 0.1)]
    x = np.arange(len(selected))
    for offset, key, label in [
        (-0.3, "author_delay", "Authors’ mean"),
        (-0.1, "conditional", "Our conditional delay"),
        (0.1, "signed_all_first_alarms", "Include signed early alarms"),
        (0.3, "conditional_minus_one", "Subtract one pulse"),
    ]:
        ax.bar(x + offset, [r[key] for r in selected], width=0.19, label=label)
    ax.set(
        xticks=x,
        xticklabels=[r["label"] for r in selected],
        ylabel="Delay / signed timing statistic (pulses)",
        title="Same saved runs and thresholds: convention sensitivity at Na=0.1",
    )
    ax.tick_params(axis="x", rotation=15)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(root / "convention_sensitivity.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    cases = report["counterfactuals"]
    x = np.arange(len(cases))
    labels = [f"{m['label']}\nNa={m['Na']:g}" for m in cases]
    axes[0].plot(
        x, [m["original_threshold"] for m in cases], "o-", label="Original ARL-calibrated h"
    )
    axes[0].plot(
        x, [m["fit"]["threshold"] for m in cases], "s-", label="Diagnostic latency-matching h"
    )
    axes[0].set(
        ylabel="Page threshold h", title="Counterfactual fit; original detector remains unchanged"
    )
    axes[0].legend(fontsize=8)
    axes[1].errorbar(
        x,
        [m["validation"]["arl"] / 1e6 for m in cases],
        yerr=[1.96 * m["validation"]["arl_se"] / 1e6 for m in cases],
        fmt="o",
        capsize=4,
    )
    axes[1].axhline(2.0, color="black", ls="--", label="Required original ARL: 2 million")
    axes[1].set(
        ylabel="No-change ARL (million pulses)",
        title="False-alarm cost at the diagnostic thresholds",
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x, labels, rotation=20)
        ax.grid(axis="y", alpha=0.2)
    fig.savefig(root / "counterfactual_thresholds.png", dpi=180)
    fig.savefig(root / "counterfactual_thresholds.svg")
    plt.close(fig)
    phase = report["phase_check"]
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    ax.errorbar(
        range(1, 9),
        [r["conditional_delay"] for r in phase["delay"].values()],
        yerr=[1.96 * r["conditional_delay_se"] for r in phase["delay"].values()],
        fmt="o-",
        capsize=4,
    )
    ax.axhline(phase["author_delay"], color="black", ls="--", label="Author n=8 mean at Na=0.1")
    ax.set(
        xticks=range(1, 9),
        xlabel="Position of first changed pulse in n=8 block",
        ylabel="Conditional delay (pulses)",
        title="Block-phase probe at the frozen author-comparison threshold",
    )
    ax.legend()
    ax.grid(alpha=0.2)
    fig.savefig(root / "phase_sensitivity.png", dpi=180)
    plt.close(fig)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    low = min(m["validation"]["arl"] for m in cases)
    high = max(m["validation"]["arl"] for m in cases)
    lines = [
        "# Controlled diagnostics of the author-latency gap",
        "",
        "Reference: [authors’ supplemental data](https://doi.org/10.5281/zenodo.17460928) and [paper Appendix F](https://arxiv.org/html/2503.12276v2#A6). The source comparison report and raw runs are preserved; this folder contains separate diagnostic experiments.",
        "",
        "## Findings",
        "",
        f"- Across all 45 Figure 7 settings, including signed early alarms changes our saved conditional mean by at most **{report['max_signed_convention_shift']:.3f} pulses**. A one-pulse index shift or other listed convention brings any point within ±5% of its author mean: **{report['convention_closes_any_gap']}**. These checks concern our observed first alarms; the authors’ event processing remains unknown.",
        f"- An independent deterministic calculation agrees with the scalar coherent/squeezed Monte Carlo delays at both meshes for all 18 settings: **{report['all_deterministic_checks_passed']}**. The method does not use importance sampling or Monte Carlo observations. It shares the Gaussian LLR CDF and physical assumptions with the implementation, so it is not an independent check of every optical assumption.",
        f"- For six representative scalar cases, fitting a LOWER diagnostic threshold to the published latency mean yields independently estimated ARLs of **{low:,.0f}–{high:,.0f} pulses**, rather than two million. All diagnostic ARL intervals lie below 95% of the original target: **{report['all_counterfactual_arl_intervals_below_target']}**. All fresh diagnostic validation checks pass: **{report['counterfactual_validation_passed']}**.",
        "- The n=8 probe varies the first changed pulse across all eight block positions at one unchanged threshold. Its phase ranges and the author reference are shown separately below. Phase estimates share a burn-in; their range is descriptive, not a worst-case Lorden bound.",
        "",
        "These results rule out the tested convention changes as a sufficient explanation within our fixed model. A lower threshold can imitate a published mean, but sacrifices the validated false-alarm interval. This does not establish that the authors used a lower threshold or that their implementation is incorrect. Their thresholds, event arrays, calibration tables and uncertainty estimates are not supplied.",
        "",
        "## Deterministic scalar method",
        "",
        "Finite-volume transition matrices include the exact Gaussian LLR CDF cell probabilities, the Page atom at zero and first-alarm absorption. Starting at W=0, propagate the pre-change distribution for nu−1 pulses, renormalizing survivors and accumulating their log probability. Solve (I−Q_post)g=1 for the expected number of post-change samples to the first alarm. The conditional delay is survivor_distribution @ g − 1, giving zero delay for detection on the first changed pulse. No restart after an early alarm is allowed. This method applies to scalar iid coherent/squeezed laws; it is not applied to correlated n>1 increments.",
        "",
        "The saved source thresholds are checked with 512 and 1024 cells by default. Agreement requires |MC−fine| <= 4 MC SE + |fine−coarse| and mesh change below 0.05 pulse. Grid differences are empirical refinement checks, not rigorous discretization bounds. Source pre-change false-alarm fractions are also compared against deterministic survival probabilities in deterministic_checks.csv.",
        "",
        "## Counterfactual thresholds: diagnostics only",
        "",
        "Each inverse calculation fixes N=100, eta1=0.9, eta2=0.85, the probe model and change pulse 1005. A root finder chooses h to match the author point mean on a 256-cell delay grid; 512- and 1024-cell evaluations then inspect that SAME h. The author mean has unknown uncertainty, so this inverse fit has no claimed confidence interval for h. These fitted thresholds are deliberately not calibrated to two-million-pulse ARL and are never written into production/source results.",
        "",
        f"Each fitted threshold is tested on {report['settings']['trials']:,} fresh changed runs and two independent adaptive-importance-sampling ensembles of {report['settings']['cycles']:,} cycles per law. Raw audit cycles, ordinary importance-weight moments, raw changed stopping times and seeds are saved. Agreement of numerical and MC delay is checked using four MC SE plus observed mesh differences; ARL is checked against a scalar numerical first-step solution and between independent repetitions. Weight diagnostics retain ESS>=1000, largest contribution<=1%, relative SE<=2%. All intervals are pointwise approximate Monte Carlo intervals conditional on fixed fitted h, not simultaneous bounds or uncertainty in author means.",
        "",
        "| Probe | Na | Original h | Diagnostic h | Author mean | New MC delay ± SE | Diagnostic ARL ± SE |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for m in cases:
        lines.append(
            f"| {m['label']} | {m['Na']:g} | {m['original_threshold']:.6f} | {m['fit']['threshold']:.6f} | {m['author_delay']:.4f} | {m['delay']['conditional_delay']:.4f} ± {m['delay']['conditional_delay_se']:.4f} | {m['validation']['arl']:,.0f} ± {m['validation']['arl_se']:,.0f} |"
        )
    values = list(phase["delay"].values())
    means = [d["conditional_delay"] for d in values]
    lines += [
        "",
        "## Phase probe",
        "",
        f"Na=0.1, n=8, frozen h={phase['threshold']:.9f}. Across change pulses 1001..1008, mean conditional delay ranges from {min(means):.3f} to {max(means):.3f} pulses. The author mean at pulse 1005 is {phase['author_delay']:.4f}. This probe changes the physical change time and mixed block accordingly, while retaining the same ARL-calibrated threshold. It does not assume the paper used the same one-based indexing convention.",
        "",
        "The convention table compares different estimands: conditional delay, inclusive (+1), shifted (−1), signed all-first-alarm times, and clipping early alarms to zero. Signed/clipped means are diagnostics and are not replacements for conditional delay. Event censoring and no-survivor cases are refused.",
        "",
        "Next: obtain the original simulator or per-setting threshold/ARL tables, or test explicitly documented alternative optical/receiver assumptions. No messages to authors are sent by this workflow.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path("results_author_comparison"))
    p.add_argument("--output", type=Path, default=Path("results_gap_diagnostics"))
    p.add_argument("--seed", type=int, default=2026102501)
    p.add_argument("--trials", type=int, default=50000)
    p.add_argument("--cycles", type=int, default=500000)
    p.add_argument("--batch-size", type=int, default=100000)
    args = p.parse_args(argv)
    if args.trials < 2 or args.batch_size < 2 or args.cycles < 2 or args.cycles % args.batch_size:
        p.error("invalid budgets; cycles must be a positive multiple of batch size")
    source_path = args.source / "report.json"
    source_bytes = source_path.read_bytes()
    source = json.loads(source_bytes)
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "settings": vars(args) | {"source": str(args.source), "output": str(args.output)},
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "raw_source_hashes": {},
        "python": sys.version,
        "versions": {k: importlib.metadata.version(k) for k in ["numpy", "scipy", "matplotlib"]},
        "conventions": [],
        "deterministic_checks": [],
        "counterfactuals": [],
    }
    for m in source["models"]:
        if m["figure"] != 7:
            continue
        raw_path = args.source / Path(m["folder"]).name / "changed_runs.npz"
        report["raw_source_hashes"][str(raw_path)] = hashlib.sha256(
            raw_path.read_bytes()
        ).hexdigest()
        with np.load(raw_path) as raw:
            runs = Runs(raw["stop"], raw["event"], 5000, 1005)
            stats = delay_conventions(runs)
        author = next(
            r["author_delay"]
            for r in source["figure7_rows"]
            if r["Na"] == m["Na"] and r["label"] == m["label"]
        )
        alternatives = [
            stats[k]
            for k in [
                "conditional",
                "conditional_minus_one",
                "conditional_inclusive",
                "signed_all_first_alarms",
                "clipped_all_first_alarms",
            ]
        ]
        report["conventions"].append(
            {
                "label": m["label"],
                "Na": m["Na"],
                "n": m["n"],
                "author_delay": author,
                **stats,
                "closest_difference": float(
                    min(alternatives, key=lambda x: abs(x - author)) - author
                ),
            }
        )
        if m["n"] == 1:
            pre, post = (coherent_pair if m["kind"] == "coherent" else squeezed_pair)(Na=m["Na"])
            h = m["calibration"]["threshold"]
            coarse = conditional_delay_numerical(pre, post, h, 1005, 512)
            fine = conditional_delay_numerical(pre, post, h, 1005, 1024)
            mesh = abs(fine["conditional_delay"] - coarse["conditional_delay"])
            error = abs(fine["conditional_delay"] - m["delay"]["conditional_delay"])
            report["deterministic_checks"].append(
                {
                    "label": m["label"],
                    "Na": m["Na"],
                    "threshold": h,
                    "mc_delay": m["delay"]["conditional_delay"],
                    "mc_se": m["delay"]["conditional_delay_se"],
                    "coarse_delay": coarse["conditional_delay"],
                    "fine_delay": fine["conditional_delay"],
                    "mesh_change": mesh,
                    "numeric_false_alarm_probability": fine["false_alarm_probability"],
                    "mc_false_alarm_probability": m["delay"]["false_alarm_probability"],
                    "passed": bool(
                        mesh < 0.05 and error <= 4 * m["delay"]["conditional_delay_se"] + mesh
                    ),
                }
            )
    children = np.random.SeedSequence(args.seed).spawn(19)
    counter = 0
    for energy in [0.1, 0.3, 0.5]:
        for kind, label, pair in [
            ("coherent", "Coherent", coherent_pair),
            ("entangled", "Squeezed", squeezed_pair),
        ]:
            original = find_model(source, energy, kind)
            a = original["Na"]
            author = next(
                r["author_delay"]
                for r in source["figure7_rows"]
                if r["Na"] == a and r["label"] == original["label"]
            )
            pre, post = pair(Na=a)
            original_h = original["calibration"]["threshold"]
            print(
                f"Counterfactual {label}, Na={a:g}: invert published mean {author:.4f}", flush=True
            )
            fit = equivalent_threshold(pre, post, author, original_h)
            h = fit["threshold"]
            folder = args.output / f"case_{counter:02d}"
            folder.mkdir()
            seeds = [
                int(s.generate_state(1, dtype=np.uint64)[0])
                for s in children[counter * 3 : (counter + 1) * 3]
            ]
            repetitions = []
            for j in range(2):
                batch = batched_estimate(
                    kind,
                    1,
                    h,
                    a,
                    cycles=args.cycles,
                    batch_size=args.batch_size,
                    seed=seeds[j],
                    proposal="adaptive",
                    proposal_offset=2.0,
                    audit_path=folder / f"audit_{j}.npz",
                )
                (folder / f"validation_{j}.json").write_text(json.dumps(batch, indent=2) + "\n")
                repetitions.append(batch["summary"])
            validation = pool_cycle_summaries(repetitions)
            runs = simulate_conditional(
                kind,
                1,
                Na=a,
                threshold=h,
                trials=args.trials,
                horizon=5000,
                change_at=1005,
                seed=seeds[2],
            )
            delay = summarize_delay(runs)
            if delay["conditional_delay"] is None:
                raise RuntimeError("censored/no-survivor diagnostic delay")
            np.savez_compressed(folder / "changed_runs.npz", stop=runs.stop, event=runs.event)
            meshes = [
                conditional_delay_numerical(pre, post, h, 1005, cells) for cells in [256, 512, 1024]
            ]
            arls = [
                mean_stopping_time(pre, post, h, cells=cells)["mean_stopping_time"]
                for cells in [512, 1024]
            ]
            delay_mesh = abs(meshes[-1]["conditional_delay"] - meshes[-2]["conditional_delay"])
            arl_mesh = abs(arls[-1] - arls[-2])
            repeat_z = abs(repetitions[0]["arl"] - repetitions[1]["arl"]) / np.hypot(
                repetitions[0]["arl_se"], repetitions[1]["arl_se"]
            )
            passed = bool(
                abs(delay["conditional_delay"] - meshes[-1]["conditional_delay"])
                <= 4 * delay["conditional_delay_se"] + delay_mesh
                and abs(validation["arl"] - arls[-1]) <= 4 * validation["arl_se"] + arl_mesh
                and repeat_z <= 4
                and all(s["diagnostics_passed"] for s in repetitions)
                and validation["diagnostics_passed"]
                and delay_mesh < 0.05
                and arl_mesh / arls[-1] < 0.001
            )
            case = {
                "kind": kind,
                "label": label,
                "Na": a,
                "original_threshold": original_h,
                "author_delay": author,
                "fit": fit,
                "seeds": seeds,
                "validation": validation,
                "repeat_z": float(repeat_z),
                "delay": delay,
                "delay_grid_checks": meshes,
                "numerical_arls": arls,
                "checks_passed": passed,
            }
            report["counterfactuals"].append(case)
            (folder / "report.json").write_text(json.dumps(case, indent=2) + "\n")
            print(
                f"  diagnostic h={h:.4f}; ARL={validation['arl']:,.0f}; validation passed={passed}",
                flush=True,
            )
            counter += 1
    model = find_model(source, 0.1, "entangled", 8)
    phase_seed = int(children[-1].generate_state(1, dtype=np.uint64)[0])
    runs, phase_seeds = simulate_phases(
        n=8,
        Na=0.1,
        threshold=model["calibration"]["threshold"],
        trials=args.trials,
        seed=phase_seed,
    )
    delay = {str(nu): summarize_delay(r) for nu, r in runs.items()}
    if any(s["conditional_delay"] is None for s in delay.values()):
        raise RuntimeError("unavailable phase delay")
    np.savez_compressed(
        args.output / "phase_runs.npz",
        **{
            f"nu{nu}_{field}": getattr(r, field)
            for nu, r in runs.items()
            for field in ["stop", "event"]
        },
    )
    report["phase_check"] = {
        "threshold": model["calibration"]["threshold"],
        "Na": 0.1,
        "n": 8,
        "seeds": phase_seeds,
        "delay": delay,
        "author_delay": next(
            r["author_delay"]
            for r in source["figure7_rows"]
            if r["n"] == 8 and np.isclose(r["Na"], 0.1)
        ),
    }
    if hashlib.sha256(source_path.read_bytes()).hexdigest() != report["source_sha256"]:
        raise RuntimeError("source report changed")
    write_outputs(args.output, report)
    print("Diagnostic study complete; original thresholds and source runs preserved.", flush=True)


if __name__ == "__main__":
    main()
