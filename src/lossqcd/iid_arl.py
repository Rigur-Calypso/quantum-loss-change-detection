"""Figure-3-style delay study using grid-checked numerical ARL calibration.

This is a sampled computational reproduction with a different calibration
method from Appendix F; no independent rare-event Monte Carlo ARL validation
at two million pulses is claimed.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .models import Gaussian, coherent_pair, squeezed_pair
from .arl_solver import mean_stopping_time, calibrate_numerical
from .calibration import independent_ratio
from .simulation import simulate_stopping_times, summarize_delay


def crosscheck_existing_monte_carlo(path):
    """Compare saved direct moderate-ARL runs with grid-based stopping-time solutions."""
    data = json.loads(path.read_text())
    rows = []
    for index, condition in enumerate(data["conditions"]):
        for name, model in condition["models"].items():
            pre, post = Gaussian(**model["pre"]), Gaussian(**model["post"])
            h = model["calibration"]["threshold"]
            base = mean_stopping_time(pre, post, h, cells=512)
            fine = mean_stopping_time(pre, post, h, cells=1024)
            mc = model["validation"]
            if mc["arl"] is None or mc["censored"]:
                raise RuntimeError("historical validation contains censored runs")
            value = fine["mean_stopping_time"]
            relative_error = abs(mc["arl"] / value - 1)
            tolerance = max(0.02, 3 * mc["arl_se"] / value)
            rows.append(
                {
                    "condition_index": index,
                    "Na": condition["Na"],
                    "target_arl": condition["target_arl"],
                    "model": name,
                    "threshold": h,
                    "monte_carlo_arl": mc["arl"],
                    "monte_carlo_se": mc["arl_se"],
                    "numerical_arl": value,
                    "relative_difference": relative_error,
                    "grid_relative_change": abs(value / base["mean_stopping_time"] - 1),
                    "acceptance_tolerance": tolerance,
                    "passed": relative_error <= tolerance,
                }
            )
    return {
        "source": str(path),
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "criterion": "absolute relative deviation <= max(2%, 3 Monte Carlo SE / numerical ARL)",
        "all_passed": all(x["passed"] for x in rows),
        "rows": rows,
    }


def write_results(root, report):
    """Save numerical calibration and its grid checks without claiming independent rare-event validation."""
    rows = report["summary"]
    with (root / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    main = [x for x in rows if x["change_at"] == 1000]
    limit = max(x["kl_ratio"] for x in main) * 1.05
    axes[0].plot([1, limit], [1, limit], "k--", label="Asymptotic 1:1 reference")
    for row in main:
        axes[0].errorbar(
            row["kl_ratio"],
            row["delay_ratio"],
            yerr=1.96 * row["delay_ratio_se"],
            fmt="o",
            capsize=3,
            color="#c5413c" if row["Na"] else "#52677c",
        )
        axes[0].annotate(
            f"Na={row['Na']:g}",
            (row["kl_ratio"], row["delay_ratio"]),
            xytext=(4, -13),
            textcoords="offset points",
            fontsize=8,
        )
    axes[0].set(
        xlabel="Squeezed / coherent KL divergence",
        ylabel="Coherent / squeezed conditional delay",
        title="Sampled Figure-3-style comparison",
    )
    axes[0].legend(fontsize=9)
    for i, name in enumerate(["coherent", "squeezed"]):
        values = [
            100 * (x[f"{name}_fine_grid_arl"] / report["settings"]["target_arl"] - 1) for x in main
        ]
        axes[1].plot(range(len(main)), values, "o-", label=name)
    axes[1].axhline(0, color="black", ls="--", lw=1)
    axes[1].set(
        xticks=range(len(main)),
        xticklabels=[f"{x['Na']:g}" for x in main],
        xlabel="Tested Na (categorical spacing)",
        ylabel="Fine-grid ARL departure from target (%)",
        title="Numerical calibration: fine-grid check",
    )
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_axisbelow(True)
    fig.suptitle(
        f"Target ARL={report['settings']['target_arl']:,.0f}; N=100, eta1=0.9, eta2=0.85; MC bars: approx. 95%"
    )
    fig.tight_layout()
    fig.savefig(root / "large_arl_comparison.png", dpi=180)
    fig.savefig(root / "large_arl_comparison.svg")
    plt.close(fig)
    check = report["moderate_arl_crosscheck"]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for name, color in [("coherent", "#52677c"), ("squeezed", "#c5413c")]:
        data = [x for x in check["rows"] if x["model"] == name]
        ax.errorbar(
            [x["numerical_arl"] for x in data],
            [x["monte_carlo_arl"] for x in data],
            yerr=[1.96 * x["monte_carlo_se"] for x in data],
            fmt="o",
            capsize=3,
            label=name,
            color=color,
        )
    ax.plot([400, 6000], [400, 6000], "k--", lw=1)
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Numerical ARL (pulses)",
        ylabel="Existing Monte Carlo ARL (pulses)",
        title="Independent method check at moderate ARLs",
    )
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(root / "solver_validation.png", dpi=180)
    plt.close(fig)
    lines = [
        "# Large-ARL numerical calibration and detection experiment",
        "",
        "Source: [arXiv:2503.12276v2, Figure 3 and Appendix F](https://arxiv.org/html/2503.12276v2#S0.F3).",
        "",
        f"Target ARL={report['settings']['target_arl']:g} pulses, N=100, eta1=0.9, eta2=0.85. Na grid={report['settings']['Na']}. Equal total transmitted energy N+Na within each pair. Change at pulse 1000, horizon 5000; change at pulse 1 is an additional initialization check.",
        "",
        "## What was verified",
        "",
        f"- All {len(check['rows'])} comparisons with the prior moderate-ARL Monte Carlo validation passed the declared criterion: relative deviation <= max(2%, 3 Monte Carlo SE / numerical ARL). These were previously saved independent no-change ensembles at target ARLs 500, 1000 and 5000, not new two-million-pulse ensembles.",
        f"- High-ARL thresholds were calibrated on {report['settings']['cells']} positive cells, then evaluated at the SAME thresholds on half and twice as many cells. Every final base-to-fine change and fine-grid deviation from target was <= {100 * report['settings']['mesh_tolerance']:g}%. Grid refinement is an empirical numerical check, not a rigorous error bound.",
        "- Raw changed-run stopping times/event flags were saved. Censored or unavailable unrestricted delay estimates cause the experiment to stop, rather than silently discarding runs.",
        "",
        "| Na | Fine-grid ARL ordinary / squeezed | Ordinary / squeezed delay (nu=1000) | Delay ratio (approx. 95% interval) | KL ratio | Early alarms ordinary / squeezed |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in main:
        lines.append(
            f"| {row['Na']:g} | {row['coherent_fine_grid_arl']:.0f} / {row['squeezed_fine_grid_arl']:.0f} | {row['coherent_delay']:.2f} / {row['squeezed_delay']:.2f} | {row['delay_ratio']:.3f} [{row['ratio_95_low']:.3f}, {row['ratio_95_high']:.3f}] | {row['kl_ratio']:.3f} | {100 * row['coherent_false_alarm_probability']:.3f}% / {100 * row['squeezed_false_alarm_probability']:.3f}% |"
        )
    max_change = max(
        x["base_to_fine_relative_change"]
        for c in report["conditions"]
        for x in [c["models"]["coherent"]["calibration"], c["models"]["squeezed"]["calibration"]]
    )
    lines += [
        "",
        f"Largest base-to-fine ARL change: {100 * max_change:.5f}%.",
        "",
        "## Method and interpretation",
        "",
        "For Page statistic w, let g(w) be expected remaining pulses to the first strict threshold crossing and F the pre-change LLR CDF. First-step conditioning gives g(w)=1+F(-w)g(0)+integral_(0,h] g(y)dF(y-w). We retain the point mass at zero and integrate exact Gaussian/quadratic-Gaussian LLR probabilities across each positive cell, approximating g by its midpoint value. A linear solve gives g(0). Threshold search sets this numerical estimate to the desired ARL. The CDF handles unequal variances, so squeezed LLRs retain their quadratic term.",
        "",
        "This differs from Appendix F’s shared-trajectory Monte Carlo calibration. It is a sampled, Figure-3-style numerical reproduction with an explicitly different ARL method. The threshold results have been checked against moderate-ARL simulations and mesh refinement; the actual rare-event ARLs have NOT been independently validated by no-change Monte Carlo at two million pulses.",
        "",
        "Reported detection delays are E[T-nu | T>=nu], not worst-case delays. False alarms before the change are reported separately; no restart is performed. nu=1 consumes T pulses but has paper-convention delay T-1. Zero squeezing uses a common threshold for identical models and independent observation streams. Uncertainty intervals are pointwise approximate Monte Carlo intervals conditional on the selected thresholds and exclude discretization/model uncertainty. No multiple-comparison correction is applied.",
        "",
        "Finite-threshold delay gains need not equal KL ratios. A closer approach to the diagonal at large ARL supports the paper’s asymptotic relationship but does not prove the limit or establish identical numerical results to its published scatter plot. The experiment does not include entangled probes.",
        "",
        "## Evidence and next work",
        "",
        "report.json contains settings, versions, seeds, crosschecks, threshold-search histories and all grid diagnostics. Per-condition report.json and raw_runs.npz hold the corresponding delay ensembles. summary.csv includes both change times and uncertainties. The original moderate-ARL evidence remains unchanged.",
        "",
        "Next: an independent rare-event ARL check (for example a regenerative estimator) and the correlated n-mode entanglement model. All numerical ARL claims should retain their discretization qualification until independently verified.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("results_large_arl"))
    p.add_argument("--validation-source", type=Path, default=Path("results_sweep/report.json"))
    p.add_argument("--seed", type=int, default=20261024)
    p.add_argument("--Na", type=float, nargs="+", default=[0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0])
    p.add_argument("--target-arl", type=float, default=2e6)
    p.add_argument("--cells", type=int, default=1024)
    p.add_argument("--check-cells", type=int, default=2048)
    p.add_argument("--mesh-tolerance", type=float, default=0.001)
    p.add_argument("--trials", type=int, default=50000)
    args = p.parse_args(argv)
    if any(not np.isfinite(x) or x < 0 for x in args.Na) or len(set(args.Na)) != len(args.Na):
        p.error("Na values must be finite, nonnegative and unique")
    if args.trials < 2 or args.cells < 16 or args.check_cells <= args.cells:
        p.error("need trials>=2 and check-cells>cells>=16")
    if not np.isfinite(args.target_arl) or args.target_arl <= 1:
        p.error("target ARL must be finite and greater than 1")
    if not args.validation_source.is_file():
        p.error("validation-source must name the saved sweep report.json")
    if not np.isfinite(args.mesh_tolerance) or not 0 < args.mesh_tolerance < 1:
        p.error("mesh tolerance must lie in (0,1)")
    args.output.mkdir(parents=True, exist_ok=False)
    print("Checking numerical method against prior Monte Carlo ensembles ...", flush=True)
    crosscheck = crosscheck_existing_monte_carlo(args.validation_source)
    (args.output / "moderate_arl_crosscheck.json").write_text(
        json.dumps(crosscheck, indent=2) + "\n"
    )
    if not crosscheck["all_passed"]:
        raise RuntimeError("moderate-ARL method check failed; investigate before extrapolating")
    print(f"All {len(crosscheck['rows'])} historical Monte Carlo checks passed.", flush=True)
    settings = vars(args) | {
        "output": str(args.output),
        "validation_source": str(args.validation_source),
    }
    report = {
        "settings": settings,
        "moderate_arl_crosscheck": crosscheck,
        "conditions": [],
        "summary": [],
        "paper": "https://arxiv.org/html/2503.12276v2",
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {k: importlib.metadata.version(k) for k in ["numpy", "scipy", "matplotlib"]},
        "rare_event_monte_carlo_arl_validated": False,
    }
    streams = np.random.SeedSequence(args.seed).spawn(len(args.Na) * 4)
    for index, Na in enumerate(args.Na):
        print(f"Na={Na:g}: numerical calibration and delay simulations ...", flush=True)
        models = {"coherent": coherent_pair(Na=Na), "squeezed": squeezed_pair(Na=Na)}
        condition = {"Na": Na, "models": {}, "stream_seeds": []}
        raw = {}
        for i, (name, (pre, post)) in enumerate(models.items()):
            if Na == 0 and name == "squeezed":
                pre, post = models["coherent"]
                calibration = condition["models"]["coherent"]["calibration"] | {
                    "reused_from": "coherent"
                }
            else:
                calibration = calibrate_numerical(
                    pre, post, args.target_arl, args.cells, args.check_cells, args.mesh_tolerance
                )
            delay = {}
            for j, nu in enumerate([1, 1000]):
                seed = int(streams[4 * index + 2 * j + i].generate_state(1, dtype=np.uint64)[0])
                condition["stream_seeds"].append({"model": name, "nu": nu, "seed": seed})
                runs = simulate_stopping_times(
                    pre, post, calibration["threshold"], args.trials, 5000, nu, seed
                )
                summary = summarize_delay(runs)
                if summary["conditional_delay"] is None:
                    raise RuntimeError("censored/no-survivor delay ensemble")
                delay[str(nu)] = summary
                raw[f"{name}_nu{nu}_stop"], raw[f"{name}_nu{nu}_event"] = runs.stop, runs.event
            condition["models"][name] = {
                "calibration": calibration,
                "delay": delay,
                "kl": post.kl(pre),
                "pre": {"mean": float(pre.mean), "variance": float(pre.variance)},
                "post": {"mean": float(post.mean), "variance": float(post.variance)},
            }
        a, b = condition["models"]["coherent"], condition["models"]["squeezed"]
        for nu in ["1", "1000"]:
            x, y = a["delay"][nu], b["delay"][nu]
            gain = independent_ratio(
                x["conditional_delay"],
                x["conditional_delay_se"],
                y["conditional_delay"],
                y["conditional_delay_se"],
            )
            report["summary"].append(
                {
                    "Na": Na,
                    "change_at": int(nu),
                    "kl_ratio": b["kl"] / a["kl"],
                    "delay_ratio": gain["ratio"],
                    "delay_ratio_se": gain["se"],
                    "ratio_95_low": gain["approx_95_interval"][0],
                    "ratio_95_high": gain["approx_95_interval"][1],
                    "delay_reduction_percent": 100 * (1 - 1 / gain["ratio"]),
                    "coherent_delay": x["conditional_delay"],
                    "coherent_delay_se": x["conditional_delay_se"],
                    "squeezed_delay": y["conditional_delay"],
                    "squeezed_delay_se": y["conditional_delay_se"],
                    "coherent_false_alarm_probability": x["false_alarm_probability"],
                    "squeezed_false_alarm_probability": y["false_alarm_probability"],
                    "coherent_fine_grid_arl": a["calibration"]["grid_checks"][-1][
                        "mean_stopping_time"
                    ],
                    "squeezed_fine_grid_arl": b["calibration"]["grid_checks"][-1][
                        "mean_stopping_time"
                    ],
                }
            )
        report["conditions"].append(condition)
        folder = args.output / f"condition_{index:02d}"
        folder.mkdir()
        (folder / "report.json").write_text(json.dumps(condition, indent=2) + "\n")
        np.savez_compressed(folder / "raw_runs.npz", **raw)
        print(
            f"  nu=1000 delay ratio: {report['summary'][-1]['delay_ratio']:.3f}; grids passed",
            flush=True,
        )
    write_results(args.output, report)
    print(
        "Large-ARL experiment complete; rare-event ARLs are numerically calibrated, not MC-validated.",
        flush=True,
    )


if __name__ == "__main__":
    main()
