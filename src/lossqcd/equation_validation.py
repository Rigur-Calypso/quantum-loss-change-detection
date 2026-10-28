"""Independent integration and sampling of unresolved Figure 4/5 models.

Ten-photon Kennedy curves are a counterfactual diagnostic, not a correction of
the published 100-photon setting. Original studies are read-only inputs.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm, poisson
from .receivers import PhaseKennedy
from .verified_figures import fixed_seed


def integrated_counts(photons, eta, reference_eta, epsilon=0.01, sigma=0.1, cutoff=128):
    """Eq.10 with adaptive real integration, independent of phase_quadrature.

    Integrate in standard-normal phase coordinates. |Z|>12 has mass <4e-33;
    this is numerical truncation, separate from the physical [-pi,pi] cut.
    """
    if photons <= 0 or sigma <= 0 or cutoff < 1:
        raise ValueError("positive energy, phase SD and count cutoff required")
    ref = np.sqrt(reference_eta) * (1 + epsilon)
    signal = np.sqrt(eta)
    bound = min(12.0, np.pi / sigma)
    normalization = norm.cdf(np.pi / sigma) - norm.cdf(-np.pi / sigma)

    def rate(z):
        return photons * abs(signal - ref * np.exp(1j * sigma * z)) ** 2

    def integrate(f):
        return quad(
            lambda z: norm.pdf(z) * f(rate(z)) / normalization,
            -bound,
            bound,
            epsabs=1e-20,
            epsrel=1e-10,
            limit=300,
        )[0]

    probs = np.array([integrate(lambda lam: poisson.pmf(k, lam)) for k in range(cutoff + 1)])
    tail = integrate(lambda lam: poisson.sf(cutoff, lam))
    return probs, float(tail)


def count_moments(photons, eta, reference_eta, epsilon=0.01, sigma=0.1):
    """Untruncated Gaussian phase moments (truncation negligible at SD=.1)."""
    a = photons * (eta + reference_eta * (1 + epsilon) ** 2)
    b = 2 * photons * np.sqrt(eta * reference_eta) * (1 + epsilon)
    mean = a - b * np.exp(-sigma * sigma / 2)
    variance = mean + b * b * np.expm1(-sigma * sigma) ** 2 / 2
    return float(mean), float(variance)


def sample_counts(photons, eta, ref, rng, trials):
    # Direct complex-amplitude rate; no receiver sampler or quadrature calls.
    """Sample the complex-amplitude Kennedy rate independently of receiver quadrature code."""
    theta = rng.normal(0, 0.1, trials)
    if np.any(abs(theta) > np.pi):
        raise RuntimeError("phase sample outside physical truncation")
    rate = photons * abs(np.sqrt(eta) - 1.01 * np.sqrt(ref) * np.exp(1j * theta)) ** 2
    return rng.poisson(rate)


def summary(x):
    """Return the empirical mean, Monte Carlo SE and sample count."""
    return {
        "mean": float(np.mean(x)),
        "se": float(np.std(x, ddof=1) / np.sqrt(len(x))),
        "n": len(x),
    }


def dense_figure4(db, n, trials, seed):
    """Build Eqs.4-6 directly; sample with dense Cholesky and score via solves."""
    s = db * np.log(10) / 20
    one = np.ones((n, n))
    eye = np.eye(n)
    c1 = (eye + 0.9 * (np.exp(-2 * s) - 1) * one / n) / 4
    c2 = (eye + 0.85 * (np.exp(-2 * s) - 1) * one / n) / 4
    m1 = np.full(n, np.sqrt(90))
    m2 = np.full(n, np.sqrt(85))
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(trials, n)) @ np.linalg.cholesky(c2).T + m2
    d1 = x - m1
    d2 = x - m2
    score = (
        np.linalg.slogdet(c1)[1]
        - np.linalg.slogdet(c2)[1]
        + np.sum(d1 * np.linalg.solve(c1, d1.T).T, axis=1)
        - np.sum(d2 * np.linalg.solve(c2, d2.T).T, axis=1)
    ) / (2 * n)
    return score


def run(output, source="results_verified_figures_v2", trials=500000, seed=2026102401):
    """Crosscheck Figures 4/5 with independent integration and dense Gaussian sampling.

    The energy fit uses one digitized point and reserves the other three as
    held out; it does not replace the paper's stated probe energy.
    """
    source = Path(source).resolve()
    output = Path(output)
    original = source / "report.json"
    source_hash = hashlib.sha256(original.read_bytes()).hexdigest()
    prior = json.loads(original.read_text())
    anchors = prior["figure5_kennedy_digitization"]
    output.mkdir(parents=True, exist_ok=False)
    children = np.random.SeedSequence(seed).spawn(10)
    seeds = [int(c.generate_state(1, dtype=np.uint64)[0]) for c in children]
    rows = []
    arrays = {}
    diagnostic = []
    # Fit one approximate point ONLY, reserving remaining three as held out.
    fitted = brentq(
        lambda e: (
            PhaseKennedy(e, 0.9, 1, 0.01, 0.1).kl(PhaseKennedy(e, 1, 1, 0.01, 0.1))
            - anchors[0]["approx_published_kennedy"]
        ),
        1,
        101,
    )
    for i, anchor in enumerate(anchors):
        loss = anchor["loss_db"]
        eta1 = 10 ** (-loss / 10)
        eta2 = 0.9 * eta1
        p, pt = integrated_counts(101, eta1, eta1)
        q, qt = integrated_counts(101, eta2, eta1)
        llr = np.log(q) - np.log(p)
        kl = float(q @ llr)
        kl64 = float(q[:65] @ llr[:65])
        receiver = PhaseKennedy(101, eta2, eta1, 0.01, 0.1)
        pre = PhaseKennedy(101, eta1, eta1, 0.01, 0.1)
        counts = sample_counts(101, eta2, eta1, np.random.default_rng(seeds[i]), trials)
        if counts.max() >= len(llr):
            raise RuntimeError("extend integrated count table before scoring samples")
        scores = llr[counts]
        mc = summary(scores)
        mean, var = count_moments(101, eta2, eta1)
        row = {
            "loss_db": loss,
            "integrated_kl": kl,
            "production_kl": receiver.kl(pre),
            "relative_integration_difference": abs(receiver.kl(pre) / kl - 1),
            "pre_tail_mass_above_128": pt,
            "post_tail_mass_above_128": qt,
            "kl_cutoff_64_to_128_change": abs(kl - kl64),
            "pre_mass": float(p.sum()),
            "post_mass": float(q.sum()),
            "mc": mc,
            "mc_z": float(abs(mc["mean"] - kl) / mc["se"]),
            "sample_count_mean": float(counts.mean()),
            "theoretical_count_mean": mean,
            "theoretical_count_variance": var,
            "count_mean_z": float(abs(counts.mean() - mean) / np.sqrt(var / trials)),
            "seed": seeds[i],
            "arrays_prefix": f"kennedy_{i}",
        }
        row["passed"] = bool(
            row["relative_integration_difference"] < 1e-7
            and mc["se"] > 0
            and row["mc_z"] < 4
            and row["count_mean_z"] < 4
            and max(pt, qt) < 1e-10
            and abs(p.sum() - 1) < 1e-10
            and abs(q.sum() - 1) < 1e-10
            and abs(kl - kl64) < 1e-10
        )
        rows.append(row)
        arrays[f"kennedy_{i}_counts"] = counts
        arrays[f"kennedy_{i}_pre"] = p
        arrays[f"kennedy_{i}_post"] = q
        for photons, label in (
            (101, "stated_energy"),
            (10, "fixed_ten_photon_counterfactual"),
            (fitted, "one_point_energy_fit"),
        ):
            prediction = PhaseKennedy(photons, eta2, eta1, 0.01, 0.1).kl(
                PhaseKennedy(photons, eta1, eta1, 0.01, 0.1)
            )
            diagnostic.append(
                {
                    "loss_db": loss,
                    "case": label,
                    "photons": photons,
                    "predicted_kl": prediction,
                    "approx_published_kl": anchor["approx_published_kennedy"],
                    "relative_difference": prediction / anchor["approx_published_kennedy"] - 1,
                    "used_for_fit": label == "one_point_energy_fit" and i == 0,
                }
            )
        print(f"Kennedy loss={loss}: independent integration/MC passed={row['passed']}", flush=True)
    dense = []
    for j, (db, n) in enumerate(((db, n) for db in (6.0, 18.0) for n in (1, 2, 8))):
        scores = dense_figure4(db, n, trials, seeds[4 + j])
        expected = fixed_seed(db, n)["kl_per_pulse"]
        mc = summary(scores)
        z = float(abs(mc["mean"] - expected) / mc["se"])
        prefix = f"dense_{j}"
        arrays[prefix + "_scores"] = scores
        dense.append(
            {
                "db": db,
                "n": n,
                "equation_kl": expected,
                "mc": mc,
                "z": z,
                "seed": seeds[4 + j],
                "arrays_prefix": prefix,
                "passed": bool(z < 4),
            }
        )
    np.savez_compressed(output / "independent_samples.npz", **arrays)
    report = {
        "source": str(original),
        "source_sha256": source_hash,
        "seed": seed,
        "trials": trials,
        "kennedy": rows,
        "figure4_dense": dense,
        "energy_diagnostic": diagnostic,
        "fitted_energy": fitted,
        "interpretation": "Counterfactual energy agreement does not identify author energy or justify changing the stated model.",
        "source_unchanged": hashlib.sha256(original.read_bytes()).hexdigest() == source_hash,
        "passed": all(r["passed"] for r in rows + dense),
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (output / "energy_diagnostic.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(diagnostic[0]))
        w.writeheader()
        w.writerows(diagnostic)
    plot(output, report)
    return report


def plot(output, report):
    """Contrast independent equation checks with the explicitly counterfactual energy diagnostic."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    for label in ("stated_energy", "fixed_ten_photon_counterfactual", "one_point_energy_fit"):
        rows = [r for r in report["energy_diagnostic"] if r["case"] == label]
        axes[0].semilogy(
            [r["loss_db"] for r in rows],
            [r["predicted_kl"] for r in rows],
            "-o",
            label=label.replace("_", " "),
        )
    axes[0].semilogy(
        [r["loss_db"] for r in rows],
        [r["approx_published_kl"] for r in rows],
        "k+",
        ms=10,
        label="Approx. published pixels",
    )
    axes[0].set(
        xlabel="Initial loss (dB)",
        ylabel="Kennedy KL (nats)",
        title="Energy diagnostic; stated setting remains 101 photons",
    )
    for n in (1, 2, 8):
        rows = [r for r in report["figure4_dense"] if r["n"] == n]
        axes[1].errorbar(
            [r["db"] for r in rows],
            [r["mc"]["mean"] for r in rows],
            yerr=[1.96 * r["mc"]["se"] for r in rows],
            fmt="o",
            capsize=3,
            label=f"Dense simulation n={n}",
        )
        axes[1].plot(
            [r["db"] for r in rows], [r["equation_kl"] for r in rows], ":", label=f"Equation n={n}"
        )
    axes[1].set(
        xlabel="Seed squeezing (dB)",
        ylabel="KL per pulse (nats)",
        title="Independent dense Gaussian simulation ±1.96 SE",
    )
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)
    for ext in ("png", "svg"):
        fig.savefig(output / f"independent_checks.{ext}", dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", default="results_equation_validation")
    p.add_argument("--source", default="results_verified_figures_v2")
    p.add_argument("--trials", type=int, default=500000)
    a = p.parse_args()
    if a.trials < 2:
        p.error("at least two trials required")
    r = run(a.output, a.source, a.trials)
    print("All independent checks passed:", r["passed"])


if __name__ == "__main__":
    main()
