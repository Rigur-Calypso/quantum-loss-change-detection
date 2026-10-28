"""Independent rare-event Monte Carlo checks of existing CUSUM thresholds.

This runner never recalibrates thresholds and never imports the grid solver.
It compares a regenerative importance-sampling estimate with frozen historical
Monte Carlo and numerical estimates. Raw cycle evidence is retained.
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
from .models import Gaussian
from .renewal import estimate_regenerative_arl, cycle_arrays
from .calibration import independent_ratio

REFERENCE = "https://web.stanford.edu/~glynn/papers/2001/CalvinGNakayama01.pdf"


def compare_estimate(stats, reference, reference_se=0.0, tolerance=0.02):
    """Compare an independent renewal estimate with a frozen reference and its uncertainty."""
    combined = float(np.hypot(stats["arl_se"], reference_se))
    difference = abs(stats["arl"] / reference - 1)
    allowed = max(tolerance, 3 * combined / reference)
    return {
        "reference_arl": reference,
        "reference_se": reference_se,
        "combined_se": combined,
        "relative_difference": difference,
        "agreement_tolerance": allowed,
        "agreement_passed": difference <= allowed,
        "standardized_difference": float((stats["arl"] - reference) / combined)
        if combined > 0
        else None,
    }


def write_outputs(root, report):
    """Write direct and rare-event checks without retuning any detector thresholds."""
    rows = report["checks"]
    with (root / "summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    high = [r for r in rows if r["phase"] == "large"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    nas = list(dict.fromkeys(r["Na"] for r in high))
    target = report["large_target_arl"]
    for name, color, offset in [("coherent", "#52677c", -0.10), ("squeezed", "#c5413c", 0.10)]:
        data = [r for r in high if r["model"] == name]
        axes[0].errorbar(
            np.arange(len(data)) + offset,
            [r["arl"] / 1e6 for r in data],
            yerr=[1.96 * r["arl_se"] / 1e6 for r in data],
            fmt="o",
            capsize=3,
            color=color,
            label=name,
        )
    axes[0].axhline(target / 1e6, color="black", ls="--", lw=1, label="Target")
    axes[0].set(
        xticks=range(len(nas)),
        xticklabels=[f"{x:g}" for x in nas],
        xlabel="Tested Na (categorical spacing)",
        ylabel="Independent ARL estimate (million pulses)",
        title="Rare-event Monte Carlo validation",
    )
    handles, labels = axes[0].get_legend_handles_labels()
    low = [r for r in rows if r["phase"] == "moderate"]
    for name, color in [("coherent", "#52677c"), ("squeezed", "#c5413c")]:
        data = [r for r in low if r["model"] == name]
        axes[1].errorbar(
            [r["reference_arl"] for r in data],
            [r["arl"] for r in data],
            yerr=[1.96 * r["arl_se"] for r in data],
            xerr=[1.96 * r["reference_se"] for r in data],
            fmt="o",
            capsize=2,
            color=color,
            label=name,
        )
    axes[1].plot([400, 6000], [400, 6000], "k--", lw=1)
    axes[1].set(
        xscale="log",
        yscale="log",
        xlabel="Historical ordinary Monte Carlo ARL",
        ylabel="Regenerative importance-sampling ARL",
        title="Check against direct simulation",
    )
    axes[1].set_xticks([500, 1000, 5000], labels=["500", "1,000", "5,000"])
    axes[1].set_yticks([500, 1000, 5000], labels=["500", "1,000", "5,000"])
    axes[1].minorticks_off()
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_axisbelow(True)
    fig.suptitle(
        "Independent algorithm and streams · fixed prior thresholds · approximate 95% Monte Carlo bars"
    )
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(root / "rare_event_validation.png", dpi=180)
    fig.savefig(root / "rare_event_validation.svg")
    plt.close(fig)
    lines = [
        "# Independent rare-event ARL validation",
        "",
        "The existing numerical-calibration thresholds were frozen. This experiment estimates their no-change ARLs using an independent regenerative importance-sampling algorithm, with no grid solver and no retuning.",
        "",
        f"**Moderate ARL checks:** {sum(r['agreement_passed'] for r in low)}/{len(low)} passed against historical direct Monte Carlo. **Large ARL checks:** {sum(r['agreement_passed'] for r in high)}/{len(high)} passed against the prior fine-grid numerical estimates.",
        "",
        "Agreement means absolute relative difference <= max(2%, 3 combined Monte Carlo SE / reference). The numerical reference has no assigned Monte Carlo SE; its discretization diagnostics are in the preceding experiment. All settings and thresholds were retained, including any failures.",
        "",
        "| Na | Method | Fixed threshold | Independent ARL ± SE | Approx. 95% interval | Numerical reference | Agrees |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for r in high:
        lines.append(
            f"| {r['Na']:g} | {r['model']} | {r['threshold']:.6f} | {r['arl']:.0f} ± {r['arl_se']:.0f} | [{r['ci_low']:.0f}, {r['ci_high']:.0f}] | {r['reference_arl']:.0f} | {r['agreement_passed']} |"
        )
    lines += [
        "",
        "## False-alarm matching within each pair",
        "",
        "The entire pointwise approximate 95% ordinary/squeezed ARL-ratio interval must fit inside a prespecified ±3% band, and each independently estimated ARL must lie within 3% of target. This is approximate equivalence at finite simulation precision, not proof of exact equality.",
        "",
        "| Na | Independent ARL ratio | Approx. 95% interval | Equivalence passed |",
        "|---|---:|---:|---|",
    ]
    for e in report["pair_equivalence"]:
        lo, hi = e["approx_95_interval"]
        lines.append(f"| {e['Na']:g} | {e['ratio']:.5f} | [{lo:.5f}, {hi:.5f}] | {e['passed']} |")
    lines += [
        "",
        "## Exact identities and estimator",
        "",
        "Let ell=log(p2(X)/p1(X)). Starting at zero, run an unreflected sum S until it returns to/below zero or strictly crosses h. This is one cycle, with length tau and upper-crossing event A. A failed cycle resets Page CUSUM to zero, giving m=E_pre[tau]+P_pre(A-complement)*m; hence ARL=m=E_pre[tau]/P_pre(A).",
        "",
        "Estimate E_pre[tau] with ordinary pre-change cycles. Independently sample cycles under the post-change law. For a stopped path the original/proposal density ratio is product(p1/p2)=exp(-S_tau), so P_pre(A)=E_post[exp(-S_tau) 1_A]. Every failed proposal cycle contributes ZERO, and the estimate divides by ALL proposal cycles. It is ordinary importance sampling, not self-normalized weighting. The terminal overshoot S_tau-h is retained. Successful weights are bounded by exp(-h). The two sample means are unbiased individually; their ratio is consistent but not exactly unbiased at finite sample size.",
        "",
        f"This is a specialization of measure-specific regenerative importance sampling; general background: [Calvin, Glynn & Nakayama (2001), section 3.2]({REFERENCE}). The source discusses estimating a regenerative numerator under the original measure and the rare-event denominator under an independent changed measure.",
        "",
        "## Evidence, uncertainty and scope",
        "",
        f"Per model/setting: {report['settings']['moderate_cycles']:,} pre and post cycles for each moderate-ARL check, and {report['settings']['large_cycles']:,} each for the large-ARL checks. Seeds are independently spawned for every setting and measure. No cycles were unfinished; the estimator refuses to discard unfinished cycles or count them as failed alarms.",
        "",
        "Per-check JSON files contain streams, cycle moments, proposal-success counts, effective weight sample size, actual simulated pulse counts and uncertainty. Compressed raw length/terminal-sum/upper-event/completion arrays are saved for both measures, so every estimate can be recomputed. Source report hashes identify the exact earlier thresholds. Code uses no new dependencies.",
        "",
        "Uncertainties use an independent-means delta-method SE and a positive log-scale approximate 95% interval. They do not include physical-model uncertainty or provide simultaneous coverage across all settings. Agreement checks are diagnostics, not exact hypothesis-free guarantees. This is weighted Monte Carlo validation; it does not mean that full no-change sequences were simulated until millions of pulses per alarm.",
        "",
        "The preceding detection-delay estimates are unchanged. The new check addresses the previously unverified rare-event false-alarm calibration under the same ideal Gaussian models. Entangled probes remain unimplemented; the next model extension is correlated n-mode blocks.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("results_rare_validation"))
    p.add_argument("--moderate-source", type=Path, default=Path("results_sweep/report.json"))
    p.add_argument("--large-source", type=Path, default=Path("results_large_arl/report.json"))
    p.add_argument("--seed", type=int, default=20261025)
    p.add_argument("--moderate-cycles", type=int, default=100000)
    p.add_argument("--large-cycles", type=int, default=500000)
    p.add_argument("--max-steps", type=int, default=100000)
    args = p.parse_args(argv)
    if min(args.moderate_cycles, args.large_cycles, args.max_steps) < 2:
        p.error("cycle counts and max-steps must be >=2")
    if not args.moderate_source.is_file() or not args.large_source.is_file():
        p.error("both historical source reports must exist")
    moderate = json.loads(args.moderate_source.read_text())
    large = json.loads(args.large_source.read_text())
    jobs = []
    for i, c in enumerate(moderate["conditions"]):
        for name, m in c["models"].items():
            if m["validation"]["arl"] is None or m["validation"]["censored"]:
                p.error("historical moderate source has censored ARL estimates")
            jobs.append(
                ("moderate", i, c["Na"], name, m, m["validation"]["arl"], m["validation"]["arl_se"])
            )
    for i, c in enumerate(large["conditions"]):
        for name, m in c["models"].items():
            jobs.append(
                (
                    "large",
                    i,
                    c["Na"],
                    name,
                    m,
                    m["calibration"]["grid_checks"][-1]["mean_stopping_time"],
                    0.0,
                )
            )
    args.output.mkdir(parents=True, exist_ok=False)
    settings = vars(args) | {
        k: str(getattr(args, k)) for k in ["output", "moderate_source", "large_source"]
    }
    report = {
        "settings": settings,
        "checks": [],
        "pair_equivalence": [],
        "large_target_arl": large["settings"]["target_arl"],
        "source_hashes": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [args.moderate_source, args.large_source]
        },
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {k: importlib.metadata.version(k) for k in ["numpy", "matplotlib"]},
        "method_reference": REFERENCE,
    }
    streams = np.random.SeedSequence(args.seed).spawn(len(jobs))
    for k, (phase, index, Na, name, m, reference, reference_se) in enumerate(jobs):
        count = args.moderate_cycles if phase == "moderate" else args.large_cycles
        seed = int(streams[k].generate_state(1, dtype=np.uint64)[0])
        pre, post = Gaussian(**m["pre"]), Gaussian(**m["post"])
        h = m["calibration"]["threshold"]
        stats, a, b = estimate_regenerative_arl(pre, post, h, count, seed, args.max_steps)
        agreement = compare_estimate(stats, reference, reference_se)
        identifier = f"{phase}_{index:02d}_{name}"
        detail = {
            "identifier": identifier,
            "phase": phase,
            "source_condition": index,
            "Na": Na,
            "model": name,
            "pre": m["pre"],
            "post": m["post"],
            "estimate": stats,
            "comparison": agreement,
        }
        folder = args.output / identifier
        folder.mkdir()
        (folder / "report.json").write_text(json.dumps(detail, indent=2) + "\n")
        np.savez_compressed(folder / "raw_cycles.npz", **cycle_arrays(a, b))
        row = {
            "phase": phase,
            "source_condition": index,
            "Na": Na,
            "model": name,
            "threshold": h,
            "arl": stats["arl"],
            "arl_se": stats["arl_se"],
            "ci_low": stats["approx_95_interval"][0],
            "ci_high": stats["approx_95_interval"][1],
            **agreement,
            "relative_se": stats["relative_se"],
            "cycle_alarm_probability": stats["cycle_alarm_probability"],
            "success_weight_ess": stats["success_weight_ess"],
            "pre_cycles": count,
            "post_cycles": count,
        }
        report["checks"].append(row)
        print(
            f"{k + 1}/{len(jobs)} {identifier} Na={Na:g}: ARL={stats['arl']:.0f} +/- {stats['arl_se']:.0f}; agreement={agreement['agreement_passed']}",
            flush=True,
        )
    high = [r for r in report["checks"] if r["phase"] == "large"]
    for Na in dict.fromkeys(r["Na"] for r in high):
        a = next(r for r in high if r["Na"] == Na and r["model"] == "coherent")
        b = next(r for r in high if r["Na"] == Na and r["model"] == "squeezed")
        ratio = independent_ratio(a["arl"], a["arl_se"], b["arl"], b["arl_se"])
        lo, hi = ratio["approx_95_interval"]
        passed = (
            0.97 <= lo
            and hi <= 1.03
            and all(abs(r["arl"] / report["large_target_arl"] - 1) <= 0.03 for r in [a, b])
        )
        report["pair_equivalence"].append({"Na": Na, **ratio, "tolerance": 0.03, "passed": passed})
    report["all_agreement_checks_passed"] = all(r["agreement_passed"] for r in report["checks"])
    report["all_high_arl_pairs_equivalent"] = all(r["passed"] for r in report["pair_equivalence"])
    write_outputs(args.output, report)
    print(
        f"All agreement checks passed: {report['all_agreement_checks_passed']}; all high-ARL pairs equivalent: {report['all_high_arl_pairs_equivalent']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
