"""Paper-focused additions: BPSK, single-photon augmentation and appendices.

Existing main-paper results are read-only. Formula-derived plots are explicitly
separated from comparisons with authors' numeric data and from LO sweeps whose
parameter grid is chosen here. Nothing fits detection latency means.
"""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from .models import coherent_pair, squeezed_pair, coherent_kl, squeezed_kl
from .entangled import entangled_kl_per_pulse, entangled_pair
from .paper_models import (
    BPSK,
    bpsk_pair,
    bpsk_capacity,
    DisplacedPhoton,
    Bernoulli,
    kennedy_residual,
    green_machine_covariance,
    augmentation_crossover,
    squeezed_energy_derivative,
)
from .receivers import PhaseHomodyne, PhaseKennedy
from .iid_detection import estimate_iid_arl, simulate_iid, simulate_threshold_grid
from .mismatch_renewal import cycle_arrays
from .simulation import summarize_delay, summarize_arl
from .channel_diagnostics import sha256
from .author_comparison import author_data
from .octave_data import load_octave


def build_pair(spec):
    """Construct the BPSK or photon pre/post observation laws described by a case specification."""
    if spec["family"] == "photon":
        return DisplacedPhoton(10, 0.9), DisplacedPhoton(10, 0.85)
    if spec["family"] == "bpsk_coherent":
        return bpsk_pair(*coherent_pair(Na=spec["Na"]))
    if spec["family"] == "bpsk_squeezed":
        return bpsk_pair(*squeezed_pair(Na=spec["Na"]))
    raise ValueError("unknown observation family")


def analytical_results():
    """Evaluate appendix formulas, receiver sweeps, mutual information and numerical controls."""
    fixed = []
    for db in np.linspace(0, 20, 81):
        s = db * np.log(10) / 20
        for n in (2, 4, 8, 16, 32, 64, 128, 256):
            na = float(np.sinh(s) ** 2 / n)
            fixed.append(
                {
                    "seed_db": float(db),
                    "n": n,
                    "Na": na,
                    "kl_per_pulse": entangled_kl_per_pulse(100, na, 0.9, 0.85, n),
                }
            )
    dv = DisplacedPhoton(10, 0.85).kl(DisplacedPhoton(10, 0.9))
    root = brentq(lambda r: squeezed_kl(100, float(np.sinh(r) ** 2), 0.9, 0.85) - dv, 0.001, 1.0)
    photon_curves = []
    weak = []
    receiver = []
    for loss in np.linspace(0.05, 10, 41):
        eta1 = float(10 ** (-loss / 10))
        eta2 = 0.944 * eta1
        photon_curves.append(
            {
                "loss_db": float(loss),
                "eta1": eta1,
                "eta2": eta2,
                "coherent": coherent_kl(100, 1, eta1, eta2),
                "squeezed": squeezed_kl(100, 1, eta1, eta2),
                "entangled_n8": entangled_kl_per_pulse(100, 1, eta1, eta2, 8),
                "photon": DisplacedPhoton(10, eta2).kl(DisplacedPhoton(10, eta1)),
            }
        )
        eta2 = 0.9 * eta1
        weak.append(
            {
                "loss_db": float(loss),
                "eta1": eta1,
                "eta2": eta2,
                "coherent_homodyne": coherent_kl(0, 1, eta1, eta2),
                "squeezed_vacuum": squeezed_kl(0, 1, eta1, eta2),
                "entangled_n8": entangled_kl_per_pulse(0, 1, eta1, eta2, 8),
                "photon_homodyne": DisplacedPhoton(0, eta2).kl(DisplacedPhoton(0, eta1)),
                "photon_counting": Bernoulli(eta2).kl(Bernoulli(eta1)),
                "kennedy_residual_1e4": kennedy_residual(0, 1, eta1, eta2, 1e-4),
                "kennedy_residual_1e5": kennedy_residual(0, 1, eta1, eta2, 1e-5),
            }
        )
        # Eq.10-11 sweep: source Figure 5 parameters not fully saved numerically.
        # Chosen explicit setup N=100, Na=1, tap=.9, fractional LO error=.01.
        row = {"loss_db": float(loss), "eta1": eta1, "eta2": eta2}
        for sigma in (0.0, 0.01, 0.03):
            a, b = PhaseHomodyne(100, 1, eta1, sigma), PhaseHomodyne(100, 1, eta2, sigma)
            c, d = (
                PhaseKennedy(101, eta1, eta1, 0.01, sigma),
                PhaseKennedy(101, eta2, eta1, 0.01, sigma),
            )
            row[f"homodyne_sigma_{sigma:g}"] = b.kl(a)
            row[f"kennedy_pnr_sigma_{sigma:g}"] = d.kl(c)
        receiver.append(row)
    capacity = []
    for n in np.linspace(0, 2, 81):
        a, _ = coherent_pair(N=n, Na=0.1)
        b, _ = squeezed_pair(N=n, Na=0.1)
        capacity.append(
            {
                "N": float(n),
                "coherent_bits": bpsk_capacity(a.mean, a.variance),
                "squeezed_bits": bpsk_capacity(b.mean, b.variance),
            }
        )
    refinement = []
    for loss in (0.05, 0.5, 2.0, 10.0):
        eta1 = float(10 ** (-loss / 10))
        eta2 = 0.9 * eta1
        coarse = PhaseHomodyne(100, 1, eta2, 0.03).kl(PhaseHomodyne(100, 1, eta1, 0.03))
        fine = PhaseHomodyne(100, 1, eta2, 0.03, order=160).kl(
            PhaseHomodyne(100, 1, eta1, 0.03, order=160)
        )
        ck = PhaseKennedy(101, eta2, eta1, 0.01, 0.03).kl(PhaseKennedy(101, eta1, eta1, 0.01, 0.03))
        fk = PhaseKennedy(101, eta2, eta1, 0.01, 0.03, 160).kl(
            PhaseKennedy(101, eta1, eta1, 0.01, 0.03, 160)
        )
        refinement.append(
            {
                "loss_db": loss,
                "homodyne_relative_change": abs(fine / coarse - 1),
                "kennedy_relative_change": abs(fk / ck - 1),
            }
        )
    covariance_error = max(
        float(
            np.max(
                abs(
                    green_machine_covariance(n, 0.1, 0.9)
                    - entangled_pair(Na=0.1, n=n)[0].covariance
                )
            )
        )
        for n in (2, 4, 8, 16, 32, 64, 128, 256)
    )
    return {
        "fixed_seed": fixed,
        "figure11a": photon_curves,
        "figure12": weak,
        "receiver_sweep": receiver,
        "phase_refinement": refinement,
        "capacity": capacity,
        "benchmarks": {
            "dv_kl": dv,
            "coherent_kl": coherent_kl(100, 1, 0.9, 0.85),
            "dv_squeezing_crossover_r": root,
            "dv_squeezing_crossover_db": float(20 * root / np.log(10)),
            "dv_squeezing_crossover_Na": float(np.sinh(root) ** 2),
            "energy_crossover": augmentation_crossover(),
            "green_machine_covariance_max_error": covariance_error,
            "squeezed_derivatives": [
                {"Na": a, "derivative": squeezed_energy_derivative(100, a, 0.9, 0.85)}
                for a in (1e-6, 1e-4, 0.01, 0.1)
            ],
        },
        "settings": {
            "fixed_seed_db": [0, 20],
            "figure11a_tap": 0.944,
            "figure12_tap": 0.9,
            "receiver_sweep_N": 100,
            "receiver_sweep_Na": 1,
            "receiver_sweep_tap": 0.9,
            "receiver_sweep_fractional_amplitude_error": 0.01,
            "phase_sds": [0, 0.01, 0.03],
            "phase_order": 96,
            "receiver_sweep_is_exact_figure5_reproduction": False,
            "figure12_entanglement_block": 8,
            "kennedy_residual_convention": "PNR; absolute positive residual amplitude, Eq.7",
        },
    }


def run_case(folder, spec, h, nu, seeds, cycles, trials):
    """Validate a frozen threshold on two independent ARL streams and one changed-run stream."""
    folder.mkdir()
    pre, post = build_pair(spec)
    repeats = []
    for j in range(2):
        stats, a, b = estimate_iid_arl(pre, post, h, cycles, seeds[j])
        np.savez_compressed(folder / f"cycles_{j}.npz", **cycle_arrays(a, b))
        repeats.append(stats)
    # Separate independent repetitions are retained; validate each against target.
    z = float(
        abs(repeats[0]["arl"] - repeats[1]["arl"])
        / np.hypot(repeats[0]["arl_se"], repeats[1]["arl_se"])
    )
    passed = z <= 4 and all(
        r["diagnostics_passed"]
        and 1.9e6 <= r["approx_95_interval"][0]
        and r["approx_95_interval"][1] <= 2.1e6
        for r in repeats
    )
    runs = simulate_iid(pre, post, h, trials=trials, horizon=5000, change_at=nu, seed=seeds[2])
    np.savez_compressed(folder / "changed.npz", stop=runs.stop, event=runs.event)
    delay = summarize_delay(runs)
    return {
        "spec": spec,
        "threshold": h,
        "change_at": nu,
        "seeds": seeds,
        "folder": folder.name,
        "repeats": repeats,
        "repeat_z": z,
        "arl_confirmed": bool(passed),
        "delay": delay,
        "kl": post.kl(pre),
        "passed": bool(passed and delay["conditional_delay"] is not None),
    }


def shared_grid_demo(output, seed):
    """Check several thresholds on shared, unreset null paths against independent renewal estimates."""
    a, b = BPSK(0.8, 0.25), BPSK(0.5, 0.3)
    thresholds = [0.5, 1.0, 2.0]
    seeds = [
        int(c.generate_state(1, dtype=np.uint64)[0]) for c in np.random.SeedSequence(seed).spawn(4)
    ]
    runs = simulate_threshold_grid(a, b, thresholds, trials=4000, horizon=100000, seed=seeds[0])
    arrays = {
        f"{field}_{i}": getattr(run, field)
        for i, run in enumerate(runs)
        for field in ("stop", "event")
    }
    np.savez_compressed(output / "shared_grid.npz", **arrays)
    rows = []
    for i, (h, run) in enumerate(zip(thresholds, runs)):
        weighted, num, den = estimate_iid_arl(a, b, h, 100000, seeds[i + 1])
        np.savez_compressed(output / f"shared_grid_cycles_{i}.npz", **cycle_arrays(num, den))
        direct = summarize_arl(run)
        z = (
            None
            if direct["arl"] is None
            else float(
                abs(direct["arl"] - weighted["arl"])
                / np.hypot(direct["arl_se"], weighted["arl_se"])
            )
        )
        rows.append(
            {
                "threshold": h,
                "direct": direct,
                "weighted": weighted,
                "z": z,
                "passed": bool(z is not None and z <= 4 and weighted["diagnostics_passed"]),
            }
        )
    return {
        "seed": seed,
        "seeds": seeds,
        "rows": rows,
        "passed": all(r["passed"] for r in rows),
        "note": "Shared paths make ARL estimates across thresholds correlated. Do not use independent difference SEs.",
    }


def run(source_root, output, cycles=250000, trials=30000, seed=2026110101):
    """Reuse frozen author-grid thresholds for BPSK and independently calibrate the photon case."""
    source_root, output = Path(source_root).resolve(), Path(output).resolve()
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 2 for x in (cycles, trials)):
        raise ValueError("integer counts >=2 required")
    source = json.loads((source_root / "report.json").read_text())
    author = author_data(source["settings"]["data"])
    models = [m for m in source["models"] if m["figure"] == 3]
    if len(models) != 16 or not all(m["confirmed"] for m in models):
        raise ValueError("expected 16 confirmed scalar Figure 3 settings")
    hashes = {str(source_root / "report.json"): sha256(source_root / "report.json")}
    for m in models:
        path = source_root / Path(m["folder"]).name / "changed_runs.npz"
        hashes[str(path)] = sha256(path)
    for name in ("prl_fig3_data", "prl_supplem_fig4b_data", "manifest.json"):
        path = Path(source["settings"]["data"]).resolve() / name
        hashes[str(path)] = sha256(path)
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "source_root": str(source_root),
        "source_hashes": hashes,
        "settings": {"cycles": cycles, "trials": trials, "seed": seed, "horizon": 5000},
        "source_url": "https://arxiv.org/html/2503.12276v2",
        "author_figure3": author["figure3"],
        "analytic": analytical_results(),
        "cases": [],
        "bpsk_rows": [],
        "controls": [],
    }
    children = np.random.SeedSequence(seed).spawn(61)
    seed_at = lambda i: int(children[i].generate_state(1, dtype=np.uint64)[0])
    for i, m in enumerate(models):
        spec = {
            "family": "bpsk_coherent" if m["kind"] == "coherent" else "bpsk_squeezed",
            "Na": m["Na"],
        }
        case = run_case(
            output / f"bpsk_{i:02d}",
            spec,
            m["calibration"]["threshold"],
            m["change_at"],
            [seed_at(3 * i + j) for j in range(3)],
            cycles,
            trials,
        )
        report["cases"].append(case)
        print(f"BPSK {i + 1}/16: {spec}, ARL confirmed={case['arl_confirmed']}", flush=True)
    for i, na in enumerate(author["figure3"]["energies"]):
        a, b = report["cases"][2 * i : 2 * i + 2]
        ratio = a["delay"]["conditional_delay"] / b["delay"]["conditional_delay"]
        ratio_se = ratio * np.hypot(
            a["delay"]["conditional_delay_se"] / a["delay"]["conditional_delay"],
            b["delay"]["conditional_delay_se"] / b["delay"]["conditional_delay"],
        )
        published = author["figure3"]["bpsk_ratios"][i]
        report["bpsk_rows"].append(
            {
                "Na": na,
                "kl_ratio": b["kl"] / a["kl"],
                "our_ratio": ratio,
                "our_ratio_se": float(ratio_se),
                "author_ratio": published,
                "relative_difference": ratio / published - 1,
                "unmodulated_our_ratio": source["figure3_rows"][i]["our_ratio"],
            }
        )
    # Independent ARL calibration of the non-Gaussian single-photon probe.
    pre, post = build_pair({"family": "photon"})
    h = 12.0
    history = []
    for i in range(4):
        stats, _, _ = estimate_iid_arl(pre, post, h, cycles, seed_at(48 + i))
        history.append({"threshold": h, "summary": stats})
        h = max(0.0, h + np.log(2e6 / stats["arl"]))
    dv = run_case(
        output / "photon",
        {"family": "photon"},
        float(h),
        1000,
        [seed_at(i) for i in (52, 53, 54)],
        cycles,
        trials,
    )
    dv["calibration_history"] = history
    report["cases"].append(dv)
    for i, pair in enumerate(
        ((BPSK(0.8, 0.25), BPSK(0.5, 0.3)), (DisplacedPhoton(2, 0.9), DisplacedPhoton(2, 0.7)))
    ):
        a, b = pair
        stats, num, den = estimate_iid_arl(a, b, 2, 100000, seed_at(55 + 2 * i))
        direct = simulate_iid(a, b, 2, trials=10000, horizon=100000, seed=seed_at(56 + 2 * i))
        summary = summarize_arl(direct)
        z = (
            None
            if summary["arl"] is None
            else float(
                abs(stats["arl"] - summary["arl"]) / np.hypot(stats["arl_se"], summary["arl_se"])
            )
        )
        filename = f"control_{i}.npz"
        np.savez_compressed(
            output / filename, stop=direct.stop, event=direct.event, **cycle_arrays(num, den)
        )
        report["controls"].append(
            {
                "family": "low_snr_bpsk" if i == 0 else "photon",
                "estimate": stats,
                "direct": summary,
                "z": z,
                "file": filename,
                "estimate_seed": seed_at(55 + 2 * i),
                "direct_seed": seed_at(56 + 2 * i),
                "passed": bool(z is not None and z <= 4 and stats["diagnostics_passed"]),
            }
        )
    report["shared_grid"] = shared_grid_demo(output, seed_at(59))
    report["source_unchanged"] = all(sha256(p) == v for p, v in hashes.items())
    report["all_detection_checks_passed"] = (
        report["source_unchanged"]
        and all(c["passed"] for c in report["cases"] + report["controls"])
        and report["shared_grid"]["passed"]
    )
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    write_outputs(report, output)
    return report


def write_outputs(report, output):
    """Write appendix calculations, author comparisons and explicitly separate chosen receiver sweeps."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for name, rows in report["analytic"].items():
        if isinstance(rows, list):
            with (output / f"{name}.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    with (output / "bpsk_comparison.csv").open("w", newline="") as handle:
        rows = report["bpsk_rows"]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    rows = report["bpsk_rows"]
    x = [r["kl_ratio"] for r in rows]
    axes[0].plot(x, [r["author_ratio"] for r in rows], "o", label="Published BPSK ratios")
    axes[0].errorbar(
        x,
        [r["our_ratio"] for r in rows],
        yerr=[1.96 * r["our_ratio_se"] for r in rows],
        fmt="s",
        capsize=3,
        label="Our BPSK ±1.96 SE",
    )
    axes[0].plot(
        x, [r["unmodulated_our_ratio"] for r in rows], "+", label="Our earlier unmodulated ratios"
    )
    axes[0].plot([1, 4], [1, 4], "k:", label="Asymptotic 1:1 reference")
    axes[0].set(
        xlabel="Squeezed/coherent KL ratio",
        ylabel="Coherent/squeezed delay ratio",
        title="Figure 3: BPSK extension at independently checked ARL",
    )
    fixed = report["analytic"]["fixed_seed"]
    for n in (2, 4, 8, 16, 32, 64, 128, 256):
        part = [r for r in fixed if r["n"] == n]
        axes[1].plot(
            [r["seed_db"] for r in part], [r["kl_per_pulse"] for r in part], label=f"n={n}"
        )
    axes[1].set(
        xlabel="Seed squeezing (dB)",
        ylabel="KL per pulse (nats)",
        title="Figure 4: fixed seed energy, Na=sinh²(s)/n",
    )
    for ax in axes:
        ax.legend(fontsize=7)
        ax.grid(alpha=0.2)
    fig.savefig(output / "bpsk_and_fixed_seed.png", dpi=180)
    fig.savefig(output / "bpsk_and_fixed_seed.svg")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    for ax, name, title in (
        (axes[0], "figure11a", "Figure 11(a): bright coherent pulses, tap=0.944"),
        (axes[1], "figure12", "Figure 12 models: one photon/pulse, tap=0.9"),
    ):
        rows = report["analytic"][name]
        for key in rows[0]:
            if key not in ("loss_db", "eta1", "eta2"):
                ax.semilogy(
                    [r["loss_db"] for r in rows],
                    [r[key] for r in rows],
                    label=key.replace("_", " "),
                )
        ax.set(xlabel="Pre-change loss (dB)", ylabel="KL per pulse (nats)", title=title)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.2)
    fig.savefig(output / "photon_probe_comparisons.png", dpi=180)
    fig.savefig(output / "photon_probe_comparisons.svg")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    rows = report["analytic"]["receiver_sweep"]
    for key in rows[0]:
        if key not in ("loss_db", "eta1", "eta2"):
            axes[0].semilogy(
                [r["loss_db"] for r in rows], [r[key] for r in rows], label=key.replace("_", " ")
            )
    axes[0].set(
        xlabel="Pre-change loss (dB)",
        ylabel="KL (nats)",
        title="Eqs.10–11 LO sensitivity: explicitly chosen parameter sweep",
    )
    rows = report["analytic"]["capacity"]
    for key in ("coherent_bits", "squeezed_bits"):
        axes[1].plot([r["N"] for r in rows], [r[key] for r in rows], label=key)
    axes[1].set(
        xlabel="Coherent energy N",
        ylabel="BPSK mutual information (bits/pulse)",
        title="Joint communications model: Na=0.1, eta=0.9",
    )
    for ax in axes:
        ax.legend(fontsize=7)
        ax.grid(alpha=0.2)
    fig.savefig(output / "receiver_and_communications.png", dpi=180)
    fig.savefig(output / "receiver_and_communications.svg")
    plt.close(fig)
    source = json.loads((Path(report["source_root"]) / "report.json").read_text())
    old = [m for m in source["models"] if m["figure"] == 3]
    supplement = load_octave(Path(source["settings"]["data"]) / "prl_supplem_fig4b_data")
    db = 20 * np.asarray(report["author_figure3"]["r"]) / np.log(10)
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for offset, label, color in ((0, "Coherent", "C0"), (1, "Squeezed", "C1")):
        unmod = old[offset::2]
        bpsk = report["cases"][offset:16:2]
        ax.errorbar(
            db,
            [m["delay"]["conditional_delay"] for m in unmod],
            yerr=[1.96 * m["delay"]["conditional_delay_se"] for m in unmod],
            fmt="o-",
            color=color,
            capsize=2,
            label=label + " unmodulated",
        )
        ax.plot(
            db,
            [m["delay"]["conditional_delay"] for m in bpsk],
            "x--",
            color=color,
            label=label + " BPSK",
        )
    dv = report["cases"][-1]
    ax.axhline(
        dv["delay"]["conditional_delay"],
        color="black",
        label="Our DV conditional first-alarm delay",
    )
    ax.axhline(
        float(supplement["lat"].mean()),
        color="gray",
        ls=":",
        label="Author DV signed mean (different early-alarm convention)",
    )
    ax.set(
        xlabel="Per-pulse squeezing (dB)",
        ylabel="Mean delay / signed latency statistic (pulses)",
        title="Figure 11(b) model comparison; original source results preserved",
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    fig.savefig(output / "figure11b_comparison.png", dpi=180)
    fig.savefig(output / "figure11b_comparison.svg")
    plt.close(fig)
    dv = report["cases"][-1]
    bench = report["analytic"]["benchmarks"]
    error = np.mean([abs(r["relative_difference"]) for r in report["bpsk_rows"]])
    (output / "RESULTS.md").write_text(f"""# Paper-specific implementation additions

Source: [arXiv:2503.12276v2](https://arxiv.org/html/2503.12276v2).
The previous robustness extensions and original source arrays remain unchanged.

- BPSK Gaussian-mixture densities, stable exact LLRs, sampling, KL and uniform-input
  mutual information are implemented. Unknown signs are marginalized, not supplied
  by an oracle. Physical raw homodyne units are used consistently before/after loss.
- Figure 3 now includes 16 new BPSK runs on the authors' eight energy points.
  Source scalar thresholds are fixed and newly checked using two independent
  {report["settings"]["cycles"]:,}-cycle ensembles per law. No threshold is fitted to
  latency. All non-Gaussian detection checks pass: **{report["all_detection_checks_passed"]}**.
  Mean absolute relative difference from the eight published BPSK delay ratios:
  **{100 * error:.2f}%**. Author errors and thresholds are unknown, so disagreement
  remains descriptive; this does not establish exact published-latency reproduction.
- Figure 4's fixed-squeezer analytical family is generated for n=2 through 256.
  The chosen seed range is 0–20 dB; these are formula-derived curves, not digitized data.
- Eq.22 one-photon homodyne densities, Eq.23 CDFs and exact mixture sampling are
  implemented. The sampler mixes vacuum and a signed chi-square(3) radius; it does
  not need numerical inverse-CDF interpolation. Independent tests check its CDF.
- The DV benchmark KL is **{bench["dv_kl"]:.6f} nats**, matching the paper's rounded
  0.2063 value. Its squeezing crossover is **{bench["dv_squeezing_crossover_db"]:.3f} dB**
  at eta1=0.9/eta2=0.85. The independently calibrated DV mean conditional delay is
  **{dv["delay"]["conditional_delay"]:.3f} ± {dv["delay"]["conditional_delay_se"]:.3f} pulses**
  (±1 Monte Carlo SE), with source-style change at 1000 and horizon 5000.
- Figure 11(a) and Figure 12 model families are evaluated. Figure 11(a) uses the
  stated rounded tap=0.944; the separate benchmark uses exactly 0.9/0.85. The weak
  entangled example uses n=8; Kennedy residuals use PNR Eq.7 and an explicit residual
  amplitude convention. These choices are saved rather than silently equated.
- Appendix A Eqs.10–11 phase-mixture homodyne and PNR Kennedy receivers are
  implemented, including phase quadrature and thermal occupation in Eq.11. The
  saved LO sweep uses chosen N=100, Na=1, tap=0.9, fractional amplitude error=0.01,
  and phase SDs 0,0.01,0.03. It is not an exact Figure 5 parameter reproduction.
  Zero-noise limits and quadrature refinement are tested independently.
- Figure 11(b) now has a model-comparison plot using the earlier unmodulated
  simulations, the new BPSK simulations and independently calibrated DV result.
  The author DV reference is the plotting script's signed mean; its two negative
  latencies are retained. It is explicitly distinguished from our conditional
  first-alarm delay. Exact numerical agreement is not claimed.
- Appendix C's Hadamard-network covariance is propagated independently and
  agrees with the analytic entangled covariance to machine precision. Appendix
  B's energy derivative is evaluated by differentiating Eq.3 directly, and its
  high-energy squeezing/classical crossover is computed without dropping the
  small variance-entropy term. The approximate Appendix B.3 threshold is retained
  separately. Appendix F's shared-path multiple-threshold first-alarm estimator
  is implemented and tested, including same-pulse overshoots and censoring.
- Bernoulli photon counting and Poisson likelihoods support the discrete receiver
  models. Ideal Kennedy nulling has infinite KL and infinite no-change ARL; its
  first post-change click still has a geometric waiting time, not zero expected
  latency at finite pulse energy.

The weak-probe and receiver-noise curves are formula/model evaluations with
explicit conventions and chosen grids. Their qualitative ranking is not imposed
in code; it can depend on loss range, residual-error convention and PNR versus
on/off observations. They are not presented as exact digitizations of Figures
5 or 12. The open coding/communications tradeoff problems are not solved here.

Every changed simulation stops at its first strict Page crossing; early alarms
and censoring are explicit. Independent moderate-ARL direct controls check the
generic non-Gaussian renewal engine. Successful iid weights are bounded by exp(-h),
with overshoot retained. No iid engine is applied to entangled correlated increments.

See the root PAPER_IMPLEMENTATION.md for complete figure-by-figure coverage and
remaining reproduction work. Full optical hardware and capacity-achieving codes
are schematics/discussion in the source, not claimed software implementations.
""")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", default="results_author_comparison")
    p.add_argument("--output", default="results_paper_completion")
    p.add_argument("--cycles", type=int, default=250000)
    p.add_argument("--trials", type=int, default=30000)
    p.add_argument("--seed", type=int, default=2026110101)
    a = p.parse_args()
    run(a.source, a.output, a.cycles, a.trials, a.seed)


if __name__ == "__main__":
    main()
