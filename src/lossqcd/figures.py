"""PDF-verified settings; equation evaluation is distinct from plotted agreement.

No author latency or threshold is fitted. Figure 5's black curve is manually
digitized at four unobscured points, with explicit pixel uncertainty.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .models import coherent_kl, squeezed_kl
from .entangled import entangled_kl_per_pulse
from .probes import DisplacedPhoton, Bernoulli, kennedy_residual
from .receivers import PhaseHomodyne, PhaseKennedy


def figure5(loss_db, order=96):
    """Evaluate the bright-probe curves at the verified Figure 5 phase/amplitude settings."""
    eta1 = float(10 ** (-loss_db / 10))
    eta2 = 0.9 * eta1
    # Energy-matched classical probe has 101 photons; squeezed displacement 100.
    a, b = (
        PhaseKennedy(101, eta1, eta1, 0.01, 0.1, order),
        PhaseKennedy(101, eta2, eta1, 0.01, 0.1, order),
    )
    return {
        "loss_db": float(loss_db),
        "coherent_ideal": coherent_kl(100, 1, eta1, eta2),
        "squeezed_ideal": squeezed_kl(100, 1, eta1, eta2),
        "entangled_n8_ideal": entangled_kl_per_pulse(100, 1, eta1, eta2, 8),
        "coherent_phase": PhaseHomodyne(101, 0, eta2, 0.1, order=order).kl(
            PhaseHomodyne(101, 0, eta1, 0.1, order=order)
        ),
        "squeezed_phase": PhaseHomodyne(100, 1, eta2, 0.1, order=order).kl(
            PhaseHomodyne(100, 1, eta1, 0.1, order=order)
        ),
        "kennedy_pnr": b.kl(a),
        "kennedy_onoff": b.onoff_kl(a),
    }


def figure12(loss_db):
    """Evaluate weak-probe receiver KLs at the verified Figure 12 settings."""
    eta1 = float(10 ** (-loss_db / 10))
    eta2 = 0.9 * eta1
    return {
        "loss_db": float(loss_db),
        "coherent": coherent_kl(0, 1, eta1, eta2),
        "squeezed_vacuum": squeezed_kl(0, 1, eta1, eta2),
        "entangled_n8": entangled_kl_per_pulse(0, 1, eta1, eta2, 8),
        "photon_homodyne": DisplacedPhoton(0, eta2).kl(DisplacedPhoton(0, eta1)),
        "photon_counting": Bernoulli(eta2).kl(Bernoulli(eta1)),
        "kennedy_pnr_1e4": kennedy_residual(0, 1, eta1, eta2, 1e-4),
        "kennedy_pnr_1e5": kennedy_residual(0, 1, eta1, eta2, 1e-5),
        "kennedy_onoff_1e4": kennedy_residual(0, 1, eta1, eta2, 1e-4, mode="onoff"),
        "kennedy_onoff_1e5": kennedy_residual(0, 1, eta1, eta2, 1e-5, mode="onoff"),
    }


def fixed_seed(db, n):
    """Convert seed squeezing in dB to per-pulse energy and independently evaluate A/n+B."""
    s = float(db * np.log(10) / 20)
    v1 = (0.1 + 0.9 * np.exp(-2 * s)) / 4
    v2 = (0.15 + 0.85 * np.exp(-2 * s)) / 4
    q = v2 / v1
    a = 0.5 * (q - 1 - np.log(q))
    b = 100 * (np.sqrt(0.85) - np.sqrt(0.9)) ** 2 / (2 * v1)
    return {
        "seed_db": float(db),
        "n": n,
        "Na": float(np.sinh(s) ** 2 / n),
        "kl_per_pulse": entangled_kl_per_pulse(100, float(np.sinh(s) ** 2 / n), 0.9, 0.85, n),
        "independent_A_over_n_plus_B": float(a / n + b),
    }


def run(output, pdf=None):
    """Save equation curves, refinement checks and manually digitized comparison anchors."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    five = [figure5(x) for x in np.linspace(0, 6, 61)]
    # At zero initial loss photon counting KL is infinite: save separately.
    twelve = [figure12(x) for x in np.linspace(0.05, 5, 100)]
    fixed = [
        fixed_seed(x, n) for x in np.linspace(0, 18, 73) for n in (1, 2, 4, 8, 16, 32, 64, 128, 256)
    ]
    checks = []
    for loss in (0, 2, 4, 6):
        coarse = figure5(loss)
        fine = figure5(loss, 160)
        for key in coarse:
            if key != "loss_db":
                checks.append(
                    {
                        "loss_db": loss,
                        "curve": key,
                        "relative_change": abs(fine[key] / coarse[key] - 1),
                    }
                )
    # Source raster 971x734: y=118 at 10 nats, y=625 at .001 nats.
    # Manual centre reads of BLACK Kennedy '+' markers; no model-guided fitting.
    anchors = []
    for loss, pixel in ((0, 538), (2, 566), (4, 594), (6, 619)):
        published = 10 ** (-3 + (625 - pixel) * 4 / 507)
        row = figure5(loss)
        anchors.append(
            {
                "loss_db": loss,
                "source_y_pixel": pixel,
                "approx_published_kennedy": published,
                "pixel_uncertainty": 3,
                "multiplicative_uncertainty": 10 ** (12 / 507),
                "pnr": row["kennedy_pnr"],
                "onoff": row["kennedy_onoff"],
                "pnr_to_digitized_ratio": row["kennedy_pnr"] / published,
            }
        )
    # Figure 12 embedded raster: 1226x976, y=54 at .1, y=835 at .001.
    # Read at x=350, corresponding to 1 dB (axis x=152..1141, 0..5 dB).
    weak_anchors = []
    predicted = figure12(1.0)
    # Red curve center at y=478 confirmed by RGB scan of source column x=350.
    for key, pixel in (
        ("coherent", 594),
        ("squeezed_vacuum", 478),
        ("entangled_n8", 657),
        ("photon_homodyne", 549),
        ("photon_counting", 345),
        ("kennedy_pnr_1e4", 491),
        ("kennedy_pnr_1e5", 435),
    ):
        published = 10 ** (-1 - (pixel - 54) * 2 / 781)
        weak_anchors.append(
            {
                "curve": key,
                "loss_db": 1.0,
                "source_y_pixel": pixel,
                "approx_published": published,
                "pixel_uncertainty": 5,
                "multiplicative_uncertainty": 10 ** (10 / 781),
                "equation_value": predicted[key],
                "relative_difference": predicted[key] / published - 1,
            }
        )
    report = {
        "source_url": "https://arxiv.org/pdf/2503.12276v2",
        "source_sha256": hashlib.sha256(Path(pdf).read_bytes()).hexdigest() if pdf else None,
        "settings": {
            "figure5": {
                "N": 100,
                "Na": 1,
                "tap": 0.9,
                "phase_variance": 0.01,
                "phase_sd": 0.1,
                "fractional_amplitude_error": 0.01,
                "loss_range_db": [0, 6],
                "thermal": 0,
                "n": 8,
            },
            "figure12": {
                "N": 0,
                "Na": 1,
                "tap": 0.9,
                "n": 8,
                "loss_plot_range_db": [0, 5],
                "finite_grid_db": [0.05, 5],
                "residual_photons": [0.0001, 0.00001],
            },
            "grid_note": "Axis ranges and legend parameters verified from PDF; original sample grids unavailable.",
        },
        "figure12_zero_loss_photon_counting": "infinite (pre-change survival probability 1)",
        "refinement": checks,
        "figure5_kennedy_digitization": anchors,
        "figure12_digitization": weak_anchors,
        "fixed_seed_max_identity_error": max(
            abs(r["kl_per_pulse"] - r["independent_A_over_n_plus_B"]) for r in fixed
        ),
        "figure4_issue": "At fixed seed, A/n+B with A>=0 implies n=1 >= n=2 >= n=4. The published inset places n=2 above n=1. Equation and plotted ordering are unresolved; no curve-matching adjustment applied.",
        "digitization_note": "Approximate manual raster readings: Figure 5 (971x734) black-marker centers +/-3 pixels; Figure 12 (1226x976), seven curves at 1 dB +/-5 pixels. These are not author numeric data or statistical error bars.",
    }
    for name, rows in [
        ("figure5", five),
        ("figure12", twelve),
        ("figure4_fixed_seed", fixed),
        ("phase_refinement", checks),
        ("figure5_kennedy_anchors", anchors),
        ("figure12_anchors", weak_anchors),
    ]:
        with (output / f"{name}.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    plot(output, five, twelve, fixed, anchors)
    return report


def plot(output, five, twelve, fixed, anchors):
    """Plot calculated receiver curves separately from approximate published pixel readings."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for ax, rows, title in (
        (axes[0], five, "Figure 5 settings: phase SD=0.1, amplitude error=0.01"),
        (axes[1], twelve, "Figure 12 settings: N=0, Na=1, tap=0.9"),
    ):
        for key in rows[0]:
            if key != "loss_db":
                ax.semilogy(
                    [r["loss_db"] for r in rows],
                    [r[key] for r in rows],
                    label=key.replace("_", " "),
                )
        ax.set(xlabel="Initial loss (dB)", ylabel="KL per pulse (nats)", title=title)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)
    axes[0].errorbar(
        [r["loss_db"] for r in anchors],
        [r["approx_published_kennedy"] for r in anchors],
        yerr=np.array(
            [
                [
                    r["approx_published_kennedy"] * (1 - 1 / r["multiplicative_uncertainty"])
                    for r in anchors
                ],
                [
                    r["approx_published_kennedy"] * (r["multiplicative_uncertainty"] - 1)
                    for r in anchors
                ],
            ]
        ),
        fmt="k+",
        capsize=3,
        label="Published Kennedy (approx. pixels)",
    )
    axes[0].legend(fontsize=7)
    axes[1].set(ylim=(1e-3, 0.1), xlim=(0, 5))
    for ext in ("png", "svg"):
        fig.savefig(output / f"verified_receiver_curves.{ext}", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for n in (1, 2, 4, 8, 16, 32, 64, 128, 256):
        rows = [r for r in fixed if r["n"] == n]
        ax.plot([r["seed_db"] for r in rows], [r["kl_per_pulse"] for r in rows], label=f"n={n}")
    x = np.linspace(0, 18, 73)
    ax.plot(
        x,
        [coherent_kl(100, float(np.sinh(t * np.log(10) / 20) ** 2), 0.9, 0.85) for t in x],
        "k--",
        label="Coherent: energy matched to n=1",
    )
    ax.set(
        xlabel="Seed squeezing before splitting (dB)",
        ylabel="KL per pulse (nats)",
        title="Figure 4 equations: published n=2 > n=1 ordering remains unresolved",
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    for ext in ("png", "svg"):
        fig.savefig(output / f"fixed_seed_equation_check.{ext}", dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", default="results_verified_figures")
    p.add_argument("--pdf")
    a = p.parse_args()
    r = run(a.output, a.pdf)
    print(
        json.dumps(
            {
                "max_refinement_change": max(x["relative_change"] for x in r["refinement"]),
                "fixed_seed_identity_error": r["fixed_seed_max_identity_error"],
                "kennedy_ratios": [
                    x["pnr_to_digitized_ratio"] for x in r["figure5_kennedy_digitization"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
