"""Independent fixed-threshold follow-up of the two flagged sweep settings.

An evidence-dependent Gaussian mixture changes only ARL sampling, not the
detector or physical model. Original thresholds and changed runs are retained.
Original failed diagnostics remain in their original result directory.
"""

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .batch_renewal import batched_estimate, pool_cycle_summaries
from .corr_arl import assess_match
from .corr_cusum import simulate_conditional
from .simulation import summarize_arl
from .renewal import _mean_se
from .corr_grid import CHOICES, make_rows


def audit_mass(path):
    """Independent path-measure identities, including failed cycles."""
    with np.load(path) as raw:
        weight = np.exp(raw["post_log_weight"])
        mass, mass_se = _mean_se(weight)
        length, length_se = _mean_se(weight * raw["post_length"])
        pre, pre_se = _mean_se(raw["pre_length"])
        assert np.all(raw["pre_complete"]) and np.all(raw["post_complete"])
    mass_z = abs(mass - 1) / mass_se if mass_se > 0 else (0.0 if mass == 1 else float("inf"))
    combined_se = float(np.hypot(length_se, pre_se))
    length_z = (
        abs(length - pre) / combined_se
        if combined_se > 0
        else (0.0 if length == pre else float("inf"))
    )
    return {
        "mass": mass,
        "mass_se": mass_se,
        "mass_z": float(mass_z),
        "weighted_cycle_length": length,
        "weighted_cycle_length_se": length_se,
        "pre_cycle_length": pre,
        "pre_cycle_length_se": pre_se,
        "length_z": float(length_z),
        "passed": bool(mass_z <= 4 and length_z <= 4),
    }


def combined_rows(source, followup):
    """Overlay fresh ARL evidence on a copy; never mutate historical records."""
    models = copy.deepcopy(source["models"])
    for updated in followup["models"]:
        m = next(m for m in models if m["Na"] == updated["Na"] and m["n"] == 8)
        assert m["calibration"]["threshold"] == updated["threshold"]
        m["validation"] = updated["validation"]
        m["arl_match"] = updated["arl_match"]
        m["confirmed"] = updated["confirmed"]
    rows = make_rows(models)
    for row in rows:
        row["validation_source"] = (
            "adaptive_followup"
            if any(row["Na"] == m["Na"] and row["n"] == 8 for m in followup["models"])
            else "original_sweep"
        )
    return rows


def write_combined(root, source, report):
    """Plot the earlier sweep with only the two prespecified ARL follow-ups overlaid."""
    rows = combined_rows(source, report)
    with (root / "combined_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(9, 5.3), layout="constrained")
    for i, (_, _, label) in enumerate(CHOICES):
        data = [r for r in rows if r["label"] == label and r["change_at"] == 1005]
        ax.errorbar(
            [r["Na"] for r in data],
            [r["delay"] if r["confirmed"] else np.nan for r in data],
            yerr=[1.96 * r["delay_se"] for r in data],
            fmt="o--" if i == 0 else "o-",
            capsize=3,
            color=plt.cm.tab10(i),
            label=label,
        )
    marked = [
        r for r in rows if r["validation_source"] == "adaptive_followup" and r["change_at"] == 1005
    ]
    ax.scatter(
        [r["Na"] for r in marked],
        [r["delay"] for r in marked],
        marker="s",
        s=100,
        facecolors="none",
        edgecolors="black",
        linewidths=1.0,
        label="ARL checked by new adaptive follow-up",
    )
    nas = source["settings"]["energies"]
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xticks(nas, [f"{a:g}" for a in nas])
    ax.set(
        xlabel="Augmentation energy Na (photons/pulse)",
        ylabel="Conditional detection delay (pulses)",
        title=f"Combined evidence: {report['combined_supported_settings']}/{report['original_total_settings']} settings supported\nChange on pulse 1005; target ARL=2 million pulses; bars: approximate 95% MC intervals",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8.5)
    fig.savefig(root / "combined_energy_sweep.png", dpi=180)
    fig.savefig(root / "combined_energy_sweep.svg")
    plt.close(fig)


def write_outputs(root, report):
    """Save fixed-threshold follow-up evidence and its combined confirmation table."""
    report["followup_confirmed"] = sum(m["confirmed"] for m in report["models"])
    report["combined_supported_settings"] = (
        report["original_confirmed_settings"] + report["followup_confirmed"]
    )
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    index = np.arange(len(report["models"]))
    axes[0].errorbar(
        index,
        [m["validation"]["arl"] / 1e6 for m in report["models"]],
        yerr=[1.96 * m["validation"]["arl_se"] / 1e6 for m in report["models"]],
        fmt="o",
        capsize=4,
    )
    target = report["settings"]["target_arl"]
    axes[0].axhline(target / 1e6, color="black", ls="--", label="Original target")
    axes[0].axhspan(
        0.95 * target / 1e6, 1.05 * target / 1e6, alpha=0.1, color="green", label="±5% band"
    )
    axes[0].set(
        xticks=index,
        xticklabels=[f"Na={m['Na']:g}" for m in report["models"]],
        ylabel="No-change ARL (million pulses)",
        title="Fresh adaptive-mixture holdouts",
    )
    axes[0].legend(fontsize=9)
    for offset, key, label in [
        (-0.16, "original_validation", "Original post proposal"),
        (0.16, "validation", "New adaptive proposal"),
    ]:
        axes[1].scatter(
            index + offset,
            [100 * m[key]["largest_weight_fraction"] for m in report["models"]],
            s=65,
            label=label,
        )
    axes[1].axhline(1.0, color="black", ls="--", label="Screen: at most 1%")
    axes[1].set(
        xticks=index,
        xticklabels=[f"Na={m['Na']:g}" for m in report["models"]],
        yscale="log",
        ylabel="Largest contribution to pooled success-weight sum (%)",
        title="Tail concentration (different sample budgets)",
    )
    axes[1].legend(fontsize=8)
    fig.suptitle("n=8 fixed-threshold follow-up; approximate 95% Monte Carlo intervals")
    fig.savefig(root / "adaptive_check.png", dpi=180)
    fig.savefig(root / "adaptive_check.svg")
    plt.close(fig)
    lines = [
        "# Adaptive-mixture validation at the original thresholds",
        "",
        f"**{report['followup_confirmed']}/2 fixed-threshold follow-ups confirmed.** Combined with the original {report['original_confirmed_settings']} confirmations, {report['combined_supported_settings']}/{report['original_total_settings']} grid settings now have supporting evidence. The original 23/25 report and its failed repetitions remain unchanged.",
        "",
        "## Method and derivation",
        "",
        "At each pulse let p1(x|prefix) and p2(x|prefix) be the two conditional Gaussian densities, W the Page evidence BEFORE the update, and ell=log(p2/p1). With a fixed b=2, sample from q=(p1+exp(W−b)*p2)/(1+exp(W−b)). Equivalently choose the post component with probability sigmoid(W−b), then sample that conditional Gaussian. The component choice is auxiliary: the likelihood correction always uses the ENTIRE mixture density.",
        "",
        "The exact one-step correction is log(p1/q)=log(1+exp(W−b))−log(1+exp(W+ell−b)). Stable logaddexp computes it. Accumulate these corrections across every pulse of the cycle, including Page reflections and partial final blocks. Using the selected component density, exp(−sum ell), or the reflected W at the end would produce incorrect weights for this proposal. The raw sum is retained for auditing but is not used as the adaptive path weight.",
        "",
        "The mixture is predictable from observed history and has full support. Multiplying these conditional density ratios gives the stopped path density ratio. Cycles still renew ONLY when W=0 at a block boundary, or end at the first strict pulse alarm W>h. Thus ARL remains E_pre[cycle length] / E_q[exp(path log weight)*cycle-alarm indicator]. The numerator and denominator use independent ensembles; the denominator averages over all cycles, including failures with zero success weight. Neither self-normalization nor clipping is used.",
        "",
        "The general idea of retaining the target density in a proposal mixture is discussed in [Hesterberg (1995)](https://doi.org/10.1080/00401706.1995.10484303). This state-dependent sequential proposal and its Page-based schedule are our implementation, not a method claimed from the optical paper or a direct application of a static mixture variance theorem. There is no claimed global bound on stopped-path weights and no proof of finite variance.",
        "",
        "## Frozen protocol and independent evidence",
        "",
        f"The exploratory pilots are in {report['settings']['pilot_source']}. Offset b=2 was fixed before these follow-up holdouts to reduce tail concentration while keeping successes common. No follow-up threshold is tuned: both are copied exactly from the earlier sweep. Each setting uses two fresh repetitions of {report['settings']['cycles']:,} cycles per sampling law in batches of {report['settings']['batch_size']:,}. Two fresh moderate checks use h=2 and h=4.5, each with {report['settings']['direct_trials']:,} direct runs and {report['settings']['moderate_cycles']:,} cycles per law.",
        "",
        "Confirmation keeps the earlier criteria: each target and coherent ARL-ratio approximate 95% interval within ±5%; both repetitions and their pool pass ESS>=1000, largest contribution<=1%, relative SE<=2%; repeats agree within four combined SE; direct moderate comparisons agree within four combined SE and pass weighted diagnostics. Each first audit batch must also satisfy total path-weight mass one and weighted proposal cycle length equal to the independent pre-change cycle mean within four SE. These are pointwise approximate Monte Carlo checks, not rigorous bounds or simultaneous coverage.",
        "",
        "| Na | Original threshold | New ARL ± SE | Repeat z | Matched to coherent | Confirmed |",
        "|---:|---:|---:|---:|---|---|",
    ]
    for m in report["models"]:
        lines.append(
            f"| {m['Na']:g} | {m['threshold']:.9f} | {m['validation']['arl']:,.0f} ± {m['validation']['arl_se']:,.0f} | {m['repeat_z']:.2f} | {m['arl_match']['passed']} | {m['confirmed']} |"
        )
    lines += [
        "",
        "## Sampling diagnostics",
        "",
        "| Na | Method | Cycles per law (pooled) | Relative SE | Success-weight ESS | Largest weight contribution |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for m in report["models"]:
        for key, label in [
            ("original_validation", "Original post proposal"),
            ("validation", "New adaptive mixture"),
        ]:
            s = m[key]
            lines.append(
                f"| {m['Na']:g} | {label} | {s['post_cycles']:,} | {100 * s['relative_se']:.3f}% | {s['success_weight_ess']:,.0f} | {100 * s['largest_weight_fraction']:.5f}% |"
            )
    lines += [
        "",
        "Original failures occurred in an individual repetition, not necessarily in the pooled value in this table. The old and new sample sizes differ; this is a descriptive diagnostic comparison, not a controlled runtime benchmark.",
        "",
        "## Direct moderate-ARL and likelihood-mass checks",
        "",
        "| Na | h | Direct ARL ± SE | Adaptive ARL ± SE | Difference / combined SE | Passed |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for m in report["models"]:
        for c in m["moderate_checks"]:
            lines.append(
                f"| {m['Na']:g} | {c['threshold']:g} | {c['direct']['arl']:.3f} ± {c['direct']['arl_se']:.3f} | {c['weighted']['arl']:.3f} ± {c['weighted']['arl_se']:.3f} | {c['z']:.2f} | {c['passed']} |"
            )
    lines += [
        "",
        "For each holdout audit batch, mass and weighted-cycle-length z scores are in report.json. Raw audit arrays retain the adaptive path log weight in addition to the raw conditional LLR sum, Page score, cycle length and success/completion flags. Batch summaries and random seeds retain every weighted batch, including extreme observations.",
        "",
        "## Detection delays carried forward",
        "",
        "Only the ARL estimator changed. The detector, its frozen thresholds and its Gaussian observation model did not change, so all eight phase-delay ensembles are carried forward from the sweep and copied into this result folder with source hashes. No new detection-speed improvement is inferred from changing the sampler.",
        "",
        "| Na | Change pulse | Original n=8 delay ± SE | Coherent / n=8 delay | New ARL confirmation |",
        "|---:|---:|---:|---:|---|",
    ]
    for m in report["models"]:
        d = m["delay"]["1005"]
        base = m["baseline_delay"]["1005"]
        lines.append(
            f"| {m['Na']:g} | 1005 | {d['conditional_delay']:.3f} ± {d['conditional_delay_se']:.3f} | {base['conditional_delay'] / d['conditional_delay']:.3f} | {m['confirmed']} |"
        )
    lines += [
        "",
        "[Combined energy plot](combined_energy_sweep.png) and [combined CSV](combined_summary.csv) overlay the new ARL evidence on copies of the original changed-run summaries. The CSV validation_source column distinguishes original_sweep from adaptive_followup. Squares in the plot mark the two follow-up points. Original results are not edited.",
        "",
        "Delay is E[T−nu | T>=nu]; first alarms before the change remain excluded and reported separately. The phase branches share burn-in, so their ranges are descriptive. This follow-up supports the two previously flagged comparisons; it does not rewrite the historical sweep as having passed its original protocol. The combined grid is still five Na samples, not a direct reproduction against the authors’ supplemental data. Next: compare the numerical curves with that supplemental dataset, then test sensitivity to imperfect or uncertain channel parameters.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path("results_conditional_sweep"))
    p.add_argument("--output", type=Path, default=Path("results_adaptive_validation"))
    p.add_argument("--pilot-source", type=Path, default=Path("results_adaptive_pilots/pilots.json"))
    p.add_argument("--seed", type=int, default=2026102301)
    p.add_argument("--cycles", type=int, default=1000000)
    p.add_argument("--moderate-cycles", type=int, default=500000)
    p.add_argument("--batch-size", type=int, default=100000)
    p.add_argument("--direct-trials", type=int, default=50000)
    args = p.parse_args(argv)
    if (
        args.batch_size < 2
        or any(c < 2 or c % args.batch_size for c in [args.cycles, args.moderate_cycles])
        or args.direct_trials < 2
    ):
        p.error("cycle counts must be positive multiples of batch size; direct trials >=2")
    source_path = args.source / "report.json"
    source_bytes = source_path.read_bytes()
    source = json.loads(source_bytes)
    models = [next(m for m in source["models"] if m["Na"] == a and m["n"] == 8) for a in [0.1, 5.0]]
    if any(m["confirmed"] for m in models) or any(m["kind"] != "entangled" for m in models):
        p.error("source must retain the two original unconfirmed entangled settings")
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "settings": vars(args)
        | {
            "source": str(args.source),
            "output": str(args.output),
            "pilot_source": str(args.pilot_source),
            "proposal_offset": 2.0,
            "target_arl": source["settings"]["target_arl"],
        },
        "source_report_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "original_confirmed_settings": source["confirmed_settings"],
        "original_total_settings": source["total_settings"],
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {k: importlib.metadata.version(k) for k in ["numpy", "matplotlib"]},
        "models": [],
    }
    children = np.random.SeedSequence(args.seed).spawn(12)
    for i, m in enumerate(models):
        a = m["Na"]
        h = m["calibration"]["threshold"]
        folder = args.output / f"model_{i:02d}"
        folder.mkdir()
        seeds = [
            int(s.generate_state(1, dtype=np.uint64)[0]) for s in children[i * 6 : (i + 1) * 6]
        ]
        repetitions = []
        mass_checks = []
        print(f"Na={a:g}, n=8: frozen threshold {h:.9f}", flush=True)
        for j in range(2):
            path = folder / f"audit_cycles_{j}.npz"
            batch = batched_estimate(
                "entangled",
                8,
                h,
                a,
                cycles=args.cycles,
                batch_size=args.batch_size,
                seed=seeds[j],
                audit_path=path,
                proposal="adaptive",
                proposal_offset=2.0,
            )
            (folder / f"validation_{j}.json").write_text(json.dumps(batch, indent=2) + "\n")
            repetitions.append(batch["summary"])
            mass_checks.append(audit_mass(path))
            print(
                f"  repeat {j + 1}: ARL={batch['summary']['arl']:,.0f} +/- {batch['summary']['arl_se']:,.0f}; mass z={mass_checks[-1]['mass_z']:.2f}",
                flush=True,
            )
        combined = pool_cycle_summaries(repetitions)
        repeat_z = float(
            abs(repetitions[0]["arl"] - repetitions[1]["arl"])
            / np.hypot(repetitions[0]["arl_se"], repetitions[1]["arl_se"])
        )
        checks = []
        for j, moderate_h in enumerate([2.0, 4.5]):
            runs = simulate_conditional(
                n=8,
                Na=a,
                threshold=moderate_h,
                trials=args.direct_trials,
                horizon=100000,
                seed=seeds[2 + 2 * j],
            )
            direct = summarize_arl(runs)
            if direct["arl"] is None:
                raise RuntimeError("censored direct moderate check")
            np.savez_compressed(folder / f"direct_{j}.npz", stop=runs.stop, event=runs.event)
            batch = batched_estimate(
                "entangled",
                8,
                moderate_h,
                a,
                cycles=args.moderate_cycles,
                batch_size=args.batch_size,
                seed=seeds[3 + 2 * j],
                proposal="adaptive",
                proposal_offset=2.0,
            )
            (folder / f"moderate_{j}.json").write_text(json.dumps(batch, indent=2) + "\n")
            weighted = batch["summary"]
            z = float(
                abs(direct["arl"] - weighted["arl"])
                / np.hypot(direct["arl_se"], weighted["arl_se"])
            )
            checks.append(
                {
                    "threshold": moderate_h,
                    "direct": direct,
                    "weighted": weighted,
                    "z": z,
                    "passed": bool(z <= 4 and weighted["diagnostics_passed"]),
                }
            )
        base = next(b for b in source["models"] if b["Na"] == a and b["kind"] == "coherent")
        match = assess_match(combined, base["validation"], source["settings"]["target_arl"], 0.05)
        passed = bool(
            match["passed"]
            and repeat_z <= 4
            and all(s["diagnostics_passed"] for s in repetitions)
            and all(c["passed"] for c in checks)
            and all(c["passed"] for c in mass_checks)
            and base["confirmed"]
        )
        original_index = source["models"].index(m)
        changed_path = (
            args.source
            / f"energy_{original_index // 5:02d}_model_{original_index % 5:02d}"
            / "changed_runs.npz"
        )
        shutil.copyfile(changed_path, folder / "changed_runs.npz")
        model = {
            "Na": a,
            "kind": "entangled",
            "n": 8,
            "threshold": h,
            "seeds": seeds,
            "original_validation": m["validation"],
            "validation": combined,
            "validation_repetitions": repetitions,
            "repeat_z": repeat_z,
            "mass_checks": mass_checks,
            "moderate_checks": checks,
            "arl_match": match,
            "confirmed": passed,
            "baseline_validation": base["validation"],
            "delay": m["delay"],
            "baseline_delay": base["delay"],
            "source_changed_runs": str(changed_path),
            "changed_runs_sha256": hashlib.sha256(changed_path.read_bytes()).hexdigest(),
        }
        report["models"].append(model)
        (folder / "report.json").write_text(json.dumps(model, indent=2) + "\n")
        print(f"  confirmed={passed}; moderate z={[round(c['z'], 2) for c in checks]}", flush=True)
    write_outputs(args.output, report)
    write_combined(args.output, source, report)
    if hashlib.sha256(source_path.read_bytes()).hexdigest() != report["source_report_sha256"]:
        raise RuntimeError("source report changed during follow-up")
    print(
        f"Follow-up confirmed {report['followup_confirmed']}/2; combined support {report['combined_supported_settings']}/{report['original_total_settings']}."
    )


if __name__ == "__main__":
    main()
