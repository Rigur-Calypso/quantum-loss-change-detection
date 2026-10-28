"""Sampled Figure 7-style energy sweep and all eight block-change phases.

Every nonzero-energy detector is calibrated separately. Holdouts never change
thresholds. At zero energy all detectors are exactly the same scalar law, so
they share one threshold/ARL validation but use independent changed ensembles.
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
from .conditional_large_arl import calibrate_rare, assess_match
from .conditional_rare_event import estimate_conditional_arl
from .batched_rare_event import batched_estimate, pool_cycle_summaries
from .conditional_detection import simulate_conditional
from .phase_simulation import simulate_phases
from .simulation import summarize_arl, summarize_delay
from .calibration import independent_ratio

CHOICES = [
    ("coherent", 1, "Coherent"),
    ("entangled", 1, "Squeezed n=1"),
    ("entangled", 2, "Entangled n=2"),
    ("entangled", 4, "Entangled n=4"),
    ("entangled", 8, "Entangled n=8"),
]


def target_match(s, target, tol):
    """Require the entire approximate ARL interval inside the target band, plus weight checks."""
    return bool(
        s["diagnostics_passed"]
        and (1 - tol) * target <= s["approx_95_interval"][0]
        and s["approx_95_interval"][1] <= (1 + tol) * target
    )


def make_rows(models):
    """Flatten detector and phase results, retaining failed confirmations and the shared zero control."""
    rows = []
    for m in models:
        base = next(b for b in models if b["Na"] == m["Na"] and b["kind"] == "coherent")
        for nu, d in m["delay"].items():
            bd = base["delay"][nu]
            ratio = (
                {"ratio": 1.0, "se": 0.0, "approx_95_interval": [1.0, 1.0]}
                if m is base
                else independent_ratio(
                    bd["conditional_delay"],
                    bd["conditional_delay_se"],
                    d["conditional_delay"],
                    d["conditional_delay_se"],
                )
            )
            rows.append(
                {
                    "Na": m["Na"],
                    "label": m["label"],
                    "n": m["n"],
                    "kind": m["kind"],
                    "change_at": int(nu),
                    "phase": (int(nu) - 1) % m["n"] + 1,
                    "threshold": m["calibration"]["threshold"],
                    "arl": m["validation"]["arl"],
                    "arl_se": m["validation"]["arl_se"],
                    "arl_matched": m["arl_match"]["passed"],
                    "confirmed": m["confirmed"],
                    "delay": d["conditional_delay"],
                    "delay_se": d["conditional_delay_se"],
                    "early_alarms": d["false_alarms"],
                    "early_alarm_probability": d["false_alarm_probability"],
                    "survivors": d["survivors"],
                    "gain": ratio["ratio"],
                    "gain_se": ratio["se"],
                    "gain_95_low": ratio["approx_95_interval"][0],
                    "gain_95_high": ratio["approx_95_interval"][1],
                }
            )
    return rows


def write_outputs(root, report):
    """Write sweep tables and plots with unsupported settings visibly marked."""
    rows = make_rows(report["models"])
    report["rows"] = rows
    report["confirmed_settings"] = sum(m["confirmed"] for m in report["models"])
    report["total_settings"] = len(report["models"])
    # A zero-energy control must be equivalent, not merely nonsignificant.
    controls = [r for r in rows if r["Na"] == 0 and r["kind"] != "coherent"]
    report["zero_control"] = {
        "tolerance": 0.03,
        "passed": bool(
            all(0.97 <= r["gain_95_low"] and r["gain_95_high"] <= 1.03 for r in controls)
        )
        if controls
        else None,
        "gain_range": [min(r["gain"] for r in controls), max(r["gain"] for r in controls)]
        if controls
        else None,
        "criterion": "each pointwise approximate 95% delay-ratio interval wholly inside [0.97,1.03]",
    }
    with (root / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    fig, ax = plt.subplots(figsize=(9, 5.2), layout="constrained")
    for i, (_, _, label) in enumerate(CHOICES):
        data = [r for r in rows if r["label"] == label and r["change_at"] == 1005]
        x = [r["Na"] for r in data]
        y = [r["delay"] if r["confirmed"] else np.nan for r in data]
        ax.errorbar(
            x,
            y,
            yerr=[1.96 * r["delay_se"] for r in data],
            fmt="o--" if i == 0 else "o-",
            capsize=3,
            label=label,
            color=plt.cm.tab10(i),
        )
        bad = [r for r in data if not r["confirmed"]]
        if bad:
            ax.scatter(
                [r["Na"] for r in bad],
                [r["delay"] for r in bad],
                marker="x",
                s=65,
                color=plt.cm.tab10(i),
            )
    if any(not m["confirmed"] for m in report["models"]):
        ax.scatter([], [], marker="x", color="black", label="Unconfirmed ARL / weight check")
    ax.set_xscale("symlog", linthresh=0.01)
    nas_ticks = report["settings"]["energies"]
    ax.set_xticks(nas_ticks, [f"{a:g}" for a in nas_ticks])
    ax.set(
        xlabel="Augmentation energy Na (photons/pulse)",
        ylabel="Conditional detection delay (pulses)",
        title="Sampled Figure 7-style comparison: change on pulse 1005\nN=100; target ARL=2 million pulses; bars: approximate 95% MC intervals",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.savefig(root / "energy_sweep.png", dpi=180)
    fig.savefig(root / "energy_sweep.svg")
    plt.close(fig)
    nas = report["settings"]["energies"]
    phases = report["settings"]["offsets"]
    matrix = np.array(
        [
            [
                next(
                    r["delay"]
                    for r in rows
                    if r["Na"] == a
                    and r["n"] == 8
                    and r["change_at"] == report["settings"]["burnin"] + p
                )
                for p in phases
            ]
            for a in nas
        ]
    )
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    im = ax.imshow(matrix, aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(phases)), phases)
    ax.set_yticks(range(len(nas)), [f"{a:g}" for a in nas])
    ax.set(
        xlabel="Position of first changed pulse within n=8 block",
        ylabel="Augmentation energy Na",
        title="n=8 phase study: conditional delay in pulses\n* marks an unconfirmed ARL / weight check",
    )
    for i in range(len(nas)):
        confirmed = next(
            m["confirmed"] for m in report["models"] if m["Na"] == nas[i] and m["n"] == 8
        )
        for j in range(len(phases)):
            red, green, blue, _ = im.cmap(im.norm(matrix[i, j]))
            luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
            ax.text(
                j,
                i,
                f"{matrix[i, j]:.1f}" + (" *" if not confirmed else ""),
                ha="center",
                va="center",
                color="black" if luminance > 0.5 else "white",
                fontsize=9,
            )
    fig.colorbar(im, ax=ax, label="Conditional detection delay (pulses)")
    fig.savefig(root / "phase_sweep.png", dpi=180)
    fig.savefig(root / "phase_sweep.svg")
    plt.close(fig)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Energy and block-phase sweep",
        "",
        "Model and pulse-by-pulse detector: [paper Figure 7 and Appendix B.4](https://arxiv.org/html/2503.12276v2#A2.SS4).",
        "This is an independent, sampled Figure 7-style study with explicit ARL matching and one-based indexing. It is not a digitized or supplemental-data error comparison.",
        "",
        f"**{report['confirmed_settings']}/{report['total_settings']} settings confirmed.** Zero-energy delay-equivalence control: **{report['zero_control']['passed']}**.",
        f"N=100, eta1=0.9, eta2=0.85; target ARL={report['settings']['target_arl']:,.0f} pulses. Na values: {nas}. Each detector uses offsets {phases} after {report['settings']['burnin']} pre-change pulses.",
        "",
        "## Change on pulse 1005",
        "",
        "| Na | Detector | ARL ± SE | Delay ± SE | Coherent / detector delay | Confirmed |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for r in rows:
        if r["change_at"] == 1005:
            lines.append(
                f"| {r['Na']:g} | {r['label']} | {r['arl']:,.0f} ± {r['arl_se']:,.0f} | {r['delay']:.3f} ± {r['delay_se']:.3f} | {r['gain']:.3f} | {r['confirmed']} |"
            )
    confirmed_rows = [r for r in rows if r["change_at"] == 1005 and r["confirmed"]]
    if confirmed_rows:
        best = max(confirmed_rows, key=lambda r: r["gain"])
        lines += [
            "",
            f"At change pulse 1005, the largest confirmed speedup on this sampled grid is {best['gain']:.3f} times ({best['label']}, Na={best['Na']:g}; delay {best['delay']:.3f} pulses). This is a grid result, not a global optimum.",
        ]
    flagged = [m for m in report["models"] if not m["confirmed"]]
    if flagged:
        lines += ["", "Unconfirmed settings and their validation screens:"]
        for m in flagged:
            lines.append(
                f"- {m['label']}, Na={m['Na']:g}: ARL match={m['arl_match']['passed']}; repeat/weight screen={m['repeat_passed']}; moderate direct crosscheck={m['moderate_check']['passed']}."
            )
    lines += [
        "",
        "## Phase sensitivity for n=8",
        "",
        "Ranges below are descriptive minima/maxima over the eight sampled change times, not worst-case Lorden delay or confidence intervals for extrema. Phase estimates share a burn-in and are correlated.",
        "",
        "| Na | Minimum delay | Maximum delay | Range (pulses) | Confirmed |",
        "|---:|---:|---:|---:|---|",
    ]
    for a, values in zip(nas, matrix):
        confirmed = next(m["confirmed"] for m in report["models"] if m["Na"] == a and m["n"] == 8)
        lines.append(
            f"| {a:g} | {values.min():.3f} | {values.max():.3f} | {np.ptp(values):.3f} | {confirmed} |"
        )
    lines += [
        "",
        "## Calibration, validation and control rules",
        "",
        "Each positive-energy setting gets a fresh calibration, two independent large-ARL holdouts, and a fresh direct no-change check at h=4.5. Four log-ARL calibration corrections are followed by one independent refinement using twice the per-evaluation sample budget. All thresholds freeze before holdouts; no holdout can retune them. The direct check must agree with a separate importance-sampling estimate within four combined SE, and that weighted estimate must pass the diagnostic screen.",
        "",
        "The default budgets are 500,000 cycles per law per calibration evaluation, 1,000,000 per law per holdout repetition, and 500,000 per law for the moderate check. All these budgets are multiplied by three for n>=4. Direct moderate checks use 10,000 first-alarm runs with horizon 100,000; changed-run phase branches use 20,000 trials with horizon 10,000. All settings and seeds are in report.json.",
        "",
        "Confirmation requires full pointwise approximate 95% ARL intervals within ±5% of target and ARL-ratio intervals within ±5% of the same-energy coherent baseline. Two holdouts must agree within four combined SE. Weighted ESS>=1000, largest weight fraction<=1%, and relative SE<=2% must hold in both repetitions, pooled holdout and moderate estimate. A detector also requires its coherent comparator to pass repeat and moderate checks. Failing settings remain in tables, with crosses in the energy plot and asterisks in the phase heatmap. Their delay gains are exploratory, not confirmed fair comparisons.",
        "",
        "At Na=0 the covariance is I/4 and every model has the identical coherent pulse law. One scalar calibration/ARL validation is therefore shared explicitly across all five labels. No independent-ratio SE is applied to this shared ARL. Each label still gets an independent changed-run experiment. Every zero-control pointwise approximate 95% delay-ratio interval must fit inside ±3% of one.",
        "",
        "The batched estimator pools the numerator and denominator first/second moments, then forms their ratio; it never averages batch ARL ratios or normalizes likelihood weights. Each batch retains its seed, full moment/weight diagnostics and counts. One 100,000-cycle audit batch per holdout retains raw lengths, raw LLR sums, Page scores, success/completion flags for both sampling laws. The remaining cycle arrays can be regenerated exactly from recorded seeds, sizes and pinned dependencies. Changed stopping times and event flags are all retained. This reduces storage without discarding failed or extreme weighted batches.",
        "",
        "Phase branches share a common pre-change prefix ending at a block boundary, retaining Page evidence and first-alarm flags. Each continuation has an independent random stream and the exact mixed-loss within-block covariance. No early alarm is restarted. Delay is E[T−nu | T>=nu]; early alarms are reported separately. Relevant censoring refuses unrestricted estimates. Pulse 1005 is phase 1 for n=2,4 and phase 5 for n=8 under the explicit one-based convention.",
        "",
        "Intervals describe Monte Carlo sampling uncertainty conditional on fixed thresholds, not simultaneous coverage across all settings or optical-model uncertainty. Importance-sampling diagnostics do not prove finite variance or exclude unseen tail contributions. Zero-control labels share ARL evidence, and phase branches share burn-in; they are not independent ARL confirmations or independent phase measurements.",
        "",
        f"The supplemental Zenodo dataset linked by the paper could not be retrieved in this session. No numerical agreement with its raw Figure 7 data is claimed. The selected Na grid spans [{min(nas):g},{max(nas):g}] with {len(nas)} samples; it does not resolve the entire continuous curve. The model remains ideal homodyne with known losses. Next: resolve any unconfirmed settings and compare directly with the authors’ supplemental data when accessible.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("results_conditional_sweep"))
    p.add_argument("--energies", type=float, nargs="+", default=[0.0, 0.01, 0.1, 1.0, 5.0])
    p.add_argument("--seed", type=int, default=2026102201)
    p.add_argument("--target-arl", type=float, default=2e6)
    p.add_argument("--calibration-cycles", type=int, default=500000)
    p.add_argument("--validation-cycles", type=int, default=1000000)
    p.add_argument("--moderate-cycles", type=int, default=500000)
    p.add_argument("--batch-size", type=int, default=100000)
    p.add_argument("--direct-trials", type=int, default=10000)
    p.add_argument("--delay-trials", type=int, default=20000)
    args = p.parse_args(argv)
    if (
        any(not np.isfinite(a) or a < 0 for a in args.energies)
        or sorted(set(args.energies)) != args.energies
    ):
        p.error("energies must be finite, nonnegative, unique and increasing")
    if (
        not np.isfinite(args.target_arl)
        or args.target_arl <= 1
        or min(args.direct_trials, args.delay_trials, args.calibration_cycles) < 2
    ):
        p.error("invalid target ARL or ensemble sizes")
    if args.batch_size < 2 or any(
        c < 2 or c % args.batch_size for c in [args.validation_cycles, args.moderate_cycles]
    ):
        p.error("validation and moderate cycle counts must be positive multiples of batch size")
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "settings": vars(args)
        | {
            "output": str(args.output),
            "burnin": 1000,
            "offsets": list(range(1, 9)),
            "horizon": 10000,
            "direct_horizon": 100000,
            "N": 100.0,
            "eta1": 0.9,
            "eta2": 0.85,
            "tolerance": 0.05,
        },
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {key: importlib.metadata.version(key) for key in ["numpy", "matplotlib"]},
        "models": [],
    }
    streams = np.random.SeedSequence(args.seed).spawn(len(args.energies) * 5 * 7)
    for ei, a in enumerate(args.energies):
        for i, (kind, n, label) in enumerate(CHOICES):
            folder = args.output / f"energy_{ei:02d}_model_{i:02d}"
            folder.mkdir()
            seeds = [
                int(child.generate_state(1, dtype=np.uint64)[0])
                for child in streams[(ei * 5 + i) * 7 : (ei * 5 + i + 1) * 7]
            ]
            factor = 3 if n >= 4 else 1
            print(f"Na={a:g} {label}: calibration / independent validation", flush=True)
            shared = a == 0 and i > 0
            baseline = next(
                (m for m in report["models"] if m["Na"] == a and m["kind"] == "coherent"), None
            )
            if shared:
                calibration = baseline["calibration"]
                validation = baseline["validation"]
                moderate = baseline["moderate_check"]
                repeat_z = baseline["repeat_z"]
                repeat_pass = baseline["repeat_passed"]
                repetitions = []
                provenance = "same scalar law at Na=0: shared coherent calibration and validation"
            else:
                calibration = calibrate_rare(
                    kind, n, args.target_arl, args.calibration_cycles * factor, seeds[0], Na=a
                )
                refinement, _, _ = estimate_conditional_arl(
                    kind,
                    n,
                    calibration["threshold"],
                    cycles=2 * args.calibration_cycles * factor,
                    seed=seeds[1],
                    Na=a,
                )
                calibration["refinement"] = refinement
                calibration["threshold"] = max(
                    0.0,
                    calibration["threshold"]
                    + float(np.clip(np.log(args.target_arl / refinement["arl"]), -2.0, 2.0)),
                )
                repetitions = []
                for j in range(2):
                    batch = batched_estimate(
                        kind,
                        n,
                        calibration["threshold"],
                        a,
                        cycles=args.validation_cycles * factor,
                        batch_size=args.batch_size,
                        seed=seeds[2 + j],
                        audit_path=folder / f"audit_cycles_{j}.npz",
                    )
                    repetitions.append(batch)
                    (folder / f"validation_{j}.json").write_text(json.dumps(batch, indent=2) + "\n")
                s0, s1 = [b["summary"] for b in repetitions]
                validation = pool_cycle_summaries([s0, s1])
                repeat_z = float(abs(s0["arl"] - s1["arl"]) / np.hypot(s0["arl_se"], s1["arl_se"]))
                repeat_pass = bool(
                    repeat_z <= 4 and s0["diagnostics_passed"] and s1["diagnostics_passed"]
                )
                direct = summarize_arl(
                    simulate_conditional(
                        kind,
                        n,
                        Na=a,
                        threshold=4.5,
                        trials=args.direct_trials,
                        horizon=100000,
                        seed=seeds[4],
                    )
                )
                if direct["arl"] is None:
                    raise RuntimeError("censored direct moderate check")
                weighted = batched_estimate(
                    kind,
                    n,
                    4.5,
                    a,
                    cycles=args.moderate_cycles * factor,
                    batch_size=args.batch_size,
                    seed=seeds[5],
                )
                (folder / "moderate_batches.json").write_text(json.dumps(weighted, indent=2) + "\n")
                ws = weighted["summary"]
                z = float(abs(direct["arl"] - ws["arl"]) / np.hypot(direct["arl_se"], ws["arl_se"]))
                moderate = {
                    "direct": direct,
                    "weighted": ws,
                    "z": z,
                    "passed": bool(z <= 4 and ws["diagnostics_passed"]),
                }
                provenance = "fresh calibration, holdouts and direct moderate check"
            if i == 0 or shared:
                match = {
                    "passed": target_match(validation, args.target_arl, 0.05),
                    "ratio": 1.0,
                    "approx_95_interval": [1.0, 1.0],
                }
            else:
                match = assess_match(validation, baseline["validation"], args.target_arl, 0.05)
            own_pass = bool(match["passed"] and moderate["passed"] and repeat_pass)
            confirmed = own_pass and (baseline is None or baseline["confirmed"])
            branches, phase_seeds = simulate_phases(
                kind,
                n,
                Na=a,
                threshold=calibration["threshold"],
                trials=args.delay_trials,
                seed=seeds[6],
            )
            delay = {str(nu): summarize_delay(run) for nu, run in branches.items()}
            if any(d["conditional_delay"] is None for d in delay.values()):
                raise RuntimeError("censored/no-survivor changed-run delay")
            np.savez_compressed(
                folder / "changed_runs.npz",
                **{
                    f"nu{nu}_{key}": getattr(run, key)
                    for nu, run in branches.items()
                    for key in ["stop", "event"]
                },
            )
            model = {
                "Na": a,
                "kind": kind,
                "n": n,
                "label": label,
                "seeds": seeds,
                "calibration": calibration,
                "validation": validation,
                "repeat_z": repeat_z,
                "repeat_passed": repeat_pass,
                "moderate_check": moderate,
                "arl_match": match,
                "confirmed": confirmed,
                "validation_provenance": provenance,
                "shared_arl": shared,
                "phase_seeds": phase_seeds,
                "delay": delay,
            }
            report["models"].append(model)
            (folder / "report.json").write_text(json.dumps(model, indent=2) + "\n")
            (args.output / "checkpoint.json").write_text(json.dumps(report, indent=2) + "\n")
            print(
                f"  ARL={validation['arl']:,.0f} +/- {validation['arl_se']:,.0f}; confirmed={confirmed}; moderate z={moderate['z']:.2f}",
                flush=True,
            )
    write_outputs(args.output, report)
    print(
        f"Completed: {report['confirmed_settings']}/{report['total_settings']} settings confirmed; zero control={report['zero_control']['passed']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
