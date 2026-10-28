"""Compare independent simulations with checksummed author data, without fitting.

Imports Octave data as text; the supplied plotting script is never executed.
All our thresholds calibrate to ARL, never to the authors' detection delays.
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
from .octave_data import load_octave
from .models import coherent_kl, squeezed_kl
from .corr_renewal import estimate_conditional_arl
from .batch_renewal import batched_estimate, pool_cycle_summaries
from .corr_cusum import simulate_conditional
from .corr_arl import assess_match
from .corr_grid import CHOICES, target_match
from .simulation import summarize_delay
from .calibration import independent_ratio


def author_path(path):
    """Resolve an absent historical author-data path; existing files take priority.

    Released study records retain their original directory name and hashes.
    Only that directory is mapped, and bytes are still verified by the caller.
    """
    path = Path(path)
    folder = path if path.name == "zenodo_17460928" else path.parent
    if not path.exists() and folder.name == "zenodo_17460928":
        return folder.with_name("authors") / path.relative_to(folder)
    return path


def author_data(root):
    """Verify published file hashes and map the author arrays to the Figure 3/7 series.

    The supplemental squeezing range fixes Figure 3 energy coordinates. The
    saved n=5000 is a horizon, not an inferred Monte Carlo sample count.
    """
    root = author_path(root)
    manifest = json.loads((root / "manifest.json").read_text())
    for entry in manifest["files"]:
        content = (root / entry["name"]).read_bytes()
        if (
            hashlib.md5(content).hexdigest() != entry["md5"]
            or hashlib.sha256(content).hexdigest() != entry["sha256"]
        ):
            raise ValueError(f"author file checksum failed: {entry['name']}")
    f3 = load_octave(root / "prl_fig3_data")
    f7 = load_octave(root / "prl_fig7_data")
    supp = load_octave(root / "prl_supplem_fig4b_data")
    if (
        f3["N"] != 100
        or f7["N"] != 100
        or f3["cp"] != 1000
        or f7["cp"] != 1005
        or f3["ARL0"] != 2e6
        or f7["ARL0"] != 2e6
    ):
        raise ValueError("unexpected author experiment settings")
    if list(f7["blocksizes"].ravel()) != [1, 2, 4, 8] or len(f7["mean_sqentlat"]) != 4:
        raise ValueError("unexpected Figure 7 block mapping")
    energies = f7["Na_range"]
    n = len(energies)
    means = [f7["mean_clat"].ravel(), *[v.ravel() for v in f7["mean_sqentlat"]]]
    if any(len(v) != n or np.any(v <= 0) for v in means):
        raise ValueError("invalid Figure 7 means")
    rs = supp["r_range"]
    nas = np.sinh(rs) ** 2
    x = f3["ratio_SS_CS"].ravel()
    if len(rs) != len(x):
        raise ValueError("Figure 3 / supplemental squeezing grids differ")
    computed = np.array(
        [squeezed_kl(100, a, 0.9, 0.85) / coherent_kl(100, a, 0.9, 0.85) for a in nas]
    )
    error = float(np.max(np.abs(computed - x)))
    if error > 1e-10:
        raise ValueError("inferred Figure 3 energy mapping does not reproduce author KL ratios")
    return {
        "manifest": manifest,
        "figure7": {
            "energies": energies.tolist(),
            "means": [v.tolist() for v in means],
            "N": f7["N"],
            "cp": int(f7["cp"]),
            "saved_n": int(f7["n"]),
            "ARL0": f7["ARL0"],
        },
        "figure3": {
            "r": rs.tolist(),
            "energies": nas.tolist(),
            "kl_ratios": x.tolist(),
            "kl_ratios_computed": computed.tolist(),
            "max_absolute_kl_error": error,
            "unmodulated_ratios": f3["ratio_det_latency"].ravel().tolist(),
            "bpsk_ratios": f3["ratio_det_latency_bpsk"].ravel().tolist(),
            "N": f3["N"],
            "cp": int(f3["cp"]),
            "saved_n": int(f3["n"]),
            "ARL0": f3["ARL0"],
        },
        "supplement_raw_latency_negative_count": int(np.sum(supp["lat"] < 0)),
        "supplement_raw_latency_count": int(supp["lat"].size),
    }


def simulate_case(folder, kind, n, Na, nu, seeds, args):
    """Calibrate a fixed-budget detector, validate twice, then save independent changed runs."""
    h = 12.0
    history = []
    for child in np.random.SeedSequence(seeds[0]).spawn(5):
        seed = int(child.generate_state(1, dtype=np.uint64)[0])
        stats, _, _ = estimate_conditional_arl(
            kind,
            n,
            h,
            Na=Na,
            cycles=args.calibration_cycles,
            seed=seed,
            proposal="adaptive",
            proposal_offset=2.0,
        )
        history.append(stats)
        h = max(0.0, h + float(np.clip(np.log(2e6 / stats["arl"]), -2.0, 2.0)))
    repetitions = []
    for j in range(2):
        batch = batched_estimate(
            kind,
            n,
            h,
            Na,
            cycles=args.validation_cycles,
            batch_size=args.batch_size,
            seed=seeds[1 + j],
            audit_path=folder / f"audit_{j}.npz",
            proposal="adaptive",
            proposal_offset=2.0,
        )
        (folder / f"validation_{j}.json").write_text(json.dumps(batch, indent=2) + "\n")
        repetitions.append(batch["summary"])
    combined = pool_cycle_summaries(repetitions)
    z = float(
        abs(repetitions[0]["arl"] - repetitions[1]["arl"])
        / np.hypot(repetitions[0]["arl_se"], repetitions[1]["arl_se"])
    )
    repeat_pass = bool(z <= 4 and all(s["diagnostics_passed"] for s in repetitions))
    # Saved n=5000 is the run horizon in the Figure 3 caption, not an assumed MC trial count.
    runs = simulate_conditional(
        kind,
        n,
        Na=Na,
        threshold=h,
        trials=args.delay_trials,
        horizon=5000,
        change_at=nu,
        seed=seeds[3],
    )
    delay = summarize_delay(runs)
    if delay["conditional_delay"] is None:
        raise RuntimeError("censored/no-survivor delay; unrestricted mean unavailable")
    np.savez_compressed(folder / "changed_runs.npz", stop=runs.stop, event=runs.event)
    model = {
        "kind": kind,
        "n": n,
        "Na": float(Na),
        "change_at": nu,
        "seeds": seeds,
        "calibration": {
            "threshold": h,
            "history": history,
            "seed": seeds[0],
            "method": "five fixed log-ARL updates using adaptive proposal; no latency fitting",
        },
        "validation": combined,
        "repeat_z": z,
        "repeat_passed": repeat_pass,
        "delay": delay,
    }
    return model


def write_outputs(root, report):
    """Write author overlays and differences while retaining unresolved latency gaps."""
    rows7 = report["figure7_rows"]
    rows3 = report["figure3_rows"]
    for name, rows in [("figure7_comparison.csv", rows7), ("figure3_comparison.csv", rows3)]:
        with (root / name).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    report["figure7_metrics"] = {}
    for _, _, label in CHOICES:
        rows = [r for r in rows7 if r["label"] == label]
        differences = np.array([r["relative_difference"] for r in rows])
        report["figure7_metrics"][label] = {
            "mean_relative_difference": float(differences.mean()),
            "max_absolute_relative_difference": float(np.max(abs(differences))),
            "points_within_5_percent": int(np.sum(abs(differences) <= 0.05)),
            "points": len(rows),
            "confirmed_arl_settings": sum(r["confirmed"] for r in rows),
        }
    report["figure3_metrics"] = {
        "max_absolute_kl_error": report["author_data"]["figure3"]["max_absolute_kl_error"],
        "mean_absolute_relative_ratio_difference": float(
            np.mean([abs(r["relative_difference"]) for r in rows3])
        ),
        "points_within_5_percent": sum(abs(r["relative_difference"]) <= 0.05 for r in rows3),
        "points": len(rows3),
        "confirmed_pairs": sum(r["confirmed"] for r in rows3),
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), layout="constrained")
    for i, (_, _, label) in enumerate(CHOICES):
        data = [r for r in rows7 if r["label"] == label]
        color = plt.cm.tab10(i)
        axes[0].plot(
            [r["Na"] for r in data],
            [r["author_delay"] for r in data],
            "--",
            color=color,
            alpha=0.75,
        )
        axes[0].errorbar(
            [r["Na"] for r in data],
            [r["our_delay"] for r in data],
            yerr=[1.96 * r["our_delay_se"] for r in data],
            fmt="o-",
            capsize=2,
            color=color,
            label=label,
        )
        axes[1].plot(
            [r["Na"] for r in data],
            [100 * r["relative_difference"] for r in data],
            "o-",
            color=color,
            label=label,
        )
    axes[0].plot([], [], "k--", label="Authors: dashed; ours: solid")
    axes[0].set(
        xlabel="Na (photons/pulse)",
        ylabel="Mean detection delay (pulses)",
        title="Figure 7: exact author energy grid",
    )
    axes[0].legend(fontsize=8)
    axes[1].axhline(0, color="black", ls="--")
    axes[1].axhspan(-5, 5, color="gray", alpha=0.15)
    axes[1].set(
        xlabel="Na (photons/pulse)",
        ylabel="100 × (ours / author − 1)",
        title="Numerical differences; author uncertainty unavailable",
    )
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle("Independent ARL-calibrated simulation versus authors’ published means")
    fig.savefig(root / "figure7_overlay.png", dpi=180)
    fig.savefig(root / "figure7_overlay.svg")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    x = [r["kl_ratio"] for r in rows3]
    ax.plot(x, [r["author_ratio"] for r in rows3], "o--", label="Authors: unmodulated")
    ax.plot(
        x,
        report["author_data"]["figure3"]["bpsk_ratios"],
        "x:",
        color="gray",
        label="Authors: BPSK (not simulated here)",
    )
    ax.errorbar(
        x,
        [r["our_ratio"] for r in rows3],
        yerr=[1.96 * r["our_ratio_se"] for r in rows3],
        fmt="o-",
        capsize=3,
        label="Our unmodulated simulation",
    )
    ax.plot([1, 4], [1, 4], "k:", label="Asymptotic 1:1 reference")
    ax.set(
        xlabel="Squeezed / coherent KL ratio",
        ylabel="Coherent / squeezed delay ratio",
        title="Figure 3 comparison: all eight author KL coordinates reproduced",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.savefig(root / "figure3_overlay.png", dpi=180)
    fig.savefig(root / "figure3_overlay.svg")
    plt.close(fig)
    (root / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Direct comparison with the authors’ supplemental dataset",
        "",
        "Source: Guha, Saikat and Basu, Prithwish, [Zenodo record 17460928](https://doi.org/10.5281/zenodo.17460928), published 2025-10-28, version v1, CC BY 4.0. The four downloaded ASCII/Octave plotting files match the published MD5 checksums; SHA-256 hashes and original metadata are preserved in data/authors. The MATLAB/Octave script is read as source material, never executed.",
        "",
        "## What matches and what does not",
        "",
        f"All eight published Figure 3 analytical KL ratios match our Gaussian formulas to maximum absolute error **{report['figure3_metrics']['max_absolute_kl_error']:.3g}**. This independently verifies the optical/statistical formulas at the published coordinates. It is not a validation of detection delays.",
        "",
        "The published latency means and our independently calibrated simulations differ. Values below are descriptive numerical differences, not combined statistical tests: author per-point standard errors and raw Figure 3/7 stopping times are not supplied. We do not tune thresholds, optical parameters or time offsets to match the reported means.",
        "",
        "| Figure 7 detector | Mean relative difference (ours/author−1) | Maximum absolute difference | Points within ±5% | Our confirmed ARL settings |",
        "|---|---:|---:|---:|---:|",
    ]
    for label, s in report["figure7_metrics"].items():
        lines.append(
            f"| {label} | {100 * s['mean_relative_difference']:.2f}% | {100 * s['max_absolute_relative_difference']:.2f}% | {s['points_within_5_percent']}/{s['points']} | {s['confirmed_arl_settings']}/{s['points']} |"
        )
    lines += [
        "",
        f"Figure 3 unmodulated delay ratios: {report['figure3_metrics']['points_within_5_percent']}/8 within ±5% numerical difference; mean absolute relative difference {100 * report['figure3_metrics']['mean_absolute_relative_ratio_difference']:.2f}%; {report['figure3_metrics']['confirmed_pairs']}/8 of our pairs pass ARL confirmation.",
        "",
        "## Figure 7 comparison at Na=0.1",
        "",
        "| Detector | Authors’ delay | Our conditional delay ± SE | Relative difference |",
        "|---|---:|---:|---:|",
    ]
    for r in rows7:
        if np.isclose(r["Na"], 0.1):
            lines.append(
                f"| {r['label']} | {r['author_delay']:.4f} | {r['our_delay']:.4f} ± {r['our_delay_se']:.4f} | {100 * r['relative_difference']:.2f}% |"
            )
    lines += [
        "",
        "## Matched inputs and explicit assumptions",
        "",
        "- Figure 7 Na=0.1:0.05:0.5 is read directly from the Octave double_range; the column-to-block mapping comes from blocksizes=[1,2,4,8] and the supplied plot script. N=100, change pulse 1005 and ARL0=2,000,000 are in the data file.",
        "- Figure 3 contains eight KL/latency ratios, N=100, r_min=0.1, N_a=1, cp=1000 and ARL0=2,000,000. Its r coordinates are inferred from the shared supplemental r_range (0.1:0.1:asinh(sqrt(1)), ending at 0.8). This mapping is checked by reproducing every supplied KL coordinate to machine precision; Figure 3 does not explicitly save its energy array.",
        "- eta1=0.9, eta2=0.85 are explicitly stated for Figure 3 in the paper and used throughout our comparison. They are not saved in the Figure 7 file. Figure 7 uses these same paper-context values as an explicit assumption, not independently verified per-file metadata.",
        "- The saved n=5000 is treated as an observation horizon, consistent with the Figure 3 paper caption and plot titles; it is not assumed to be the number of Monte Carlo trials. Our new simulations use horizon 5000 and report their separate trial count.",
        "- Our detector makes strict one-based Page crossings, retains prefix history through reflections, and generates the exact mixed-loss covariance for a straddling block. Delay is E[T−nu | T>=nu]; early first alarms are excluded and counted separately. All unrestricted estimates refuse relevant censoring.",
        "- BPSK data are retained and shown for context, but this implementation compares unmodulated homodyne models. No BPSK reproduction is claimed. The Supplemental Figure 4b experiment is archived but is not simulated in this milestone.",
        "",
        "## Independent simulation protocol",
        "",
        f"Each of the 45 Figure 7 settings and 16 Figure 3 scalar settings receives five fixed log-ARL calibration updates using {report['settings']['calibration_cycles']:,} cycles per law. Two fresh holdout repetitions each use {report['settings']['validation_cycles']:,} cycles per law, batched at {report['settings']['batch_size']:,}. Adaptive proposal offset b=2 is fixed. Each changed-run ensemble has {report['settings']['delay_trials']:,} trials, a separate stream, and horizon 5000. Thresholds are frozen before holdouts and are never fitted to author latency values.",
        "",
        "Confirmation requires both holdouts and their pool to pass ESS>=1000, largest weighted contribution<=1%, relative SE<=2%, repeat agreement within four combined SE, plus approximate target and same-energy coherent ARL-ratio 95% intervals fully inside ±5%. Any unconfirmed cases remain visible in CSV and JSON. Intervals are pointwise Monte Carlo uncertainty conditional on thresholds; they do not include model uncertainty or uncertainty in published author means. All raw changed-run arrays, audit cycle arrays, calibration histories, weighted batch moments, seeds and versions are retained.",
        "",
        "## Unresolved differences and information needed",
        "",
        "The dataset supplies plotting data and a plotting script, not the simulator or ARL lookup tables. It does not provide the detector thresholds, per-point ARL validation, raw first-alarm/event flags, seeds, author per-point SEs, or explicit Figure 7 transmissivities. These missing details prevent attribution of the latency gap to one cause. Agreement of KL formulas and our ARL diagnostics does not prove that either implementation is wrong.",
        "A one-sample indexing convention changes a delay by one pulse; the observed coherent discrepancy is larger, and cannot be explained solely by block phase because that detector has block size one. No index correction is applied to reduce the gap.",
        f"The separate Supplemental Figure 4b raw latency vector contains {report['author_data']['supplement_raw_latency_negative_count']} negative entries among {report['author_data']['supplement_raw_latency_count']} values, and its plotting script applies mean(lat). That is a different experiment; it does not establish how Figure 3 or 7 early alarms were processed. We record it only as a reason to request explicit first-alarm handling rather than infer it from summary means.",
        "Next research step: obtain or reconstruct the original per-setting thresholds and simulator conventions, then test candidate explanations with independent ARL checks. Contacting the authors is not performed by this script.",
        "",
    ]
    (root / "RESULTS.md").write_text("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=Path("data/authors"))
    p.add_argument("--output", type=Path, default=Path("results_author_comparison"))
    p.add_argument("--seed", type=int, default=2026102401)
    p.add_argument("--calibration-cycles", type=int, default=200000)
    p.add_argument("--validation-cycles", type=int, default=500000)
    p.add_argument("--batch-size", type=int, default=100000)
    p.add_argument("--delay-trials", type=int, default=30000)
    args = p.parse_args(argv)
    if (
        min(args.calibration_cycles, args.batch_size, args.delay_trials) < 2
        or args.validation_cycles < 2
        or args.validation_cycles % args.batch_size
    ):
        p.error("invalid budgets; validation cycles must be a multiple of batch size")
    authors = author_data(args.data)
    args.output.mkdir(parents=True, exist_ok=False)
    report = {
        "settings": vars(args)
        | {
            "data": str(args.data),
            "output": str(args.output),
            "N": 100.0,
            "eta1": 0.9,
            "eta2": 0.85,
            "horizon": 5000,
            "target_arl": 2e6,
            "proposal": "adaptive",
            "proposal_offset": 2.0,
        },
        "author_data": authors,
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {k: importlib.metadata.version(k) for k in ["numpy", "matplotlib"]},
        "models": [],
        "figure7_rows": [],
        "figure3_rows": [],
    }
    streams = np.random.SeedSequence(args.seed).spawn(61 * 4)
    counter = 0
    for figure, energies in [
        (7, authors["figure7"]["energies"]),
        (3, authors["figure3"]["energies"]),
    ]:
        for ei, energy in enumerate(energies):
            base = None
            pair = []
            for mi, (kind, n, label) in enumerate(CHOICES if figure == 7 else CHOICES[:2]):
                folder = args.output / f"f{figure}_e{ei:02d}_m{mi}"
                folder.mkdir()
                seeds = [
                    int(c.generate_state(1, dtype=np.uint64)[0])
                    for c in streams[counter * 4 : (counter + 1) * 4]
                ]
                counter += 1
                print(f"Figure {figure}, Na={energy:.6g}, {label}", flush=True)
                model = simulate_case(
                    folder, kind, n, energy, 1005 if figure == 7 else 1000, seeds, args
                )
                s = model["validation"]
                match = (
                    {"passed": target_match(s, 2e6, 0.05), "ratio": 1.0}
                    if base is None
                    else assess_match(s, base["validation"], 2e6, 0.05)
                )
                model.update(
                    {
                        "figure": figure,
                        "label": label,
                        "folder": str(folder),
                        "arl_match": match,
                        "confirmed": bool(
                            match["passed"]
                            and model["repeat_passed"]
                            and (base is None or base["confirmed"])
                        ),
                    }
                )
                if base is None:
                    base = model
                pair.append(model)
                report["models"].append(model)
                (folder / "report.json").write_text(json.dumps(model, indent=2) + "\n")
                d = model["delay"]
                print(
                    f"  ARL {s['arl']:,.0f}; delay {d['conditional_delay']:.3f}; confirmed={model['confirmed']}",
                    flush=True,
                )
                if figure == 7:
                    author_delay = authors["figure7"]["means"][mi][ei]
                    report["figure7_rows"].append(
                        {
                            "Na": float(energy),
                            "label": label,
                            "n": n,
                            "author_delay": author_delay,
                            "our_delay": d["conditional_delay"],
                            "our_delay_se": d["conditional_delay_se"],
                            "relative_difference": d["conditional_delay"] / author_delay - 1,
                            "arl": s["arl"],
                            "arl_se": s["arl_se"],
                            "confirmed": model["confirmed"],
                            "false_alarms": d["false_alarms"],
                            "survivors": d["survivors"],
                        }
                    )
            if figure == 3:
                c, q = [m["delay"] for m in pair]
                ratio = independent_ratio(
                    c["conditional_delay"],
                    c["conditional_delay_se"],
                    q["conditional_delay"],
                    q["conditional_delay_se"],
                )
                author = authors["figure3"]["unmodulated_ratios"][ei]
                report["figure3_rows"].append(
                    {
                        "r": authors["figure3"]["r"][ei],
                        "Na": float(energy),
                        "kl_ratio": authors["figure3"]["kl_ratios"][ei],
                        "author_ratio": author,
                        "our_ratio": ratio["ratio"],
                        "our_ratio_se": ratio["se"],
                        "relative_difference": ratio["ratio"] / author - 1,
                        "confirmed": all(m["confirmed"] for m in pair),
                    }
                )
            (args.output / "checkpoint.json").write_text(json.dumps(report, indent=2) + "\n")
    write_outputs(args.output, report)
    print("Comparison completed without fitting author latency data.", flush=True)


if __name__ == "__main__":
    main()
