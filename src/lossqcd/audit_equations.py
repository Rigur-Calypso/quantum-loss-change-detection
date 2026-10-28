"""Audit independently integrated receiver laws, raw samples and holdout fits."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .equations import (
    integrated_counts,
    count_moments,
    sample_counts,
    dense_figure4,
    summary,
)
from .figures import fixed_seed
from .receivers import PhaseKennedy


def audit(root, regenerate=False):
    """Check independent quadrature, dense sampling and held-out energy-fit diagnostics.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    r = json.loads((root / "report.json").read_text())
    checks = []

    def check(name, value):
        checks.append({"name": name, "passed": bool(value)})

    check(
        "source unchanged",
        hashlib.sha256(Path(r["source"]).read_bytes()).hexdigest() == r["source_sha256"],
    )
    prior = json.loads(Path(r["source"]).read_text())
    anchors = prior["figure5_kennedy_digitization"]
    seeds = [
        int(c.generate_state(1, dtype=np.uint64)[0])
        for c in np.random.SeedSequence(r["seed"]).spawn(10)
    ]
    raw = np.load(root / "independent_samples.npz", allow_pickle=False)
    check(
        "four independent receiver cases and six dense cases",
        len(r["kennedy"]) == 4 and len(r["figure4_dense"]) == 6,
    )
    for i, row in enumerate(r["kennedy"]):
        prefix = row["arrays_prefix"]
        p = raw[prefix + "_pre"]
        q = raw[prefix + "_post"]
        counts = raw[prefix + "_counts"]
        e = 10 ** (-row["loss_db"] / 10)
        llr = np.log(q) - np.log(p)
        mc = summary(llr[counts])
        kl = float(q @ llr)
        check(f"case {i} independent seed", row["seed"] == seeds[i])
        check(f"case {i} sample size", len(counts) == r["trials"])
        check(f"case {i} regenerated summary", mc == row["mc"])
        check(f"case {i} integrated KL", kl == row["integrated_kl"])
        check(f"case {i} normalized laws", abs(p.sum() - 1) < 1e-10 and abs(q.sum() - 1) < 1e-10)
        check(
            f"case {i} production agreement",
            np.isclose(
                PhaseKennedy(101, 0.9 * e, e, 0.01, 0.1).kl(PhaseKennedy(101, e, e, 0.01, 0.1)),
                kl,
                rtol=1e-7,
            ),
        )
        check(
            f"case {i} MC agrees with integral",
            row["passed"] and abs(mc["mean"] - kl) < 4 * mc["se"],
        )
        mean, var = count_moments(101, 0.9 * e, e)
        check(
            f"case {i} count moments",
            abs(counts.mean() - mean) < 4 * np.sqrt(var / len(counts))
            and np.isclose(var, row["theoretical_count_variance"]),
        )
        if regenerate:
            check(
                f"case {i} exact raw count regeneration",
                np.array_equal(
                    counts,
                    sample_counts(101, 0.9 * e, e, np.random.default_rng(seeds[i]), r["trials"]),
                ),
            )
            pp, pt = integrated_counts(101, e, e)
            qq, qt = integrated_counts(101, 0.9 * e, e)
            check(
                f"case {i} independent law regeneration",
                np.array_equal(pp, p)
                and np.array_equal(qq, q)
                and pt == row["pre_tail_mass_above_128"]
                and qt == row["post_tail_mass_above_128"],
            )
    for j, row in enumerate(r["figure4_dense"]):
        x = raw[row["arrays_prefix"] + "_scores"]
        mc = summary(x)
        expected = fixed_seed(row["db"], row["n"])["kl_per_pulse"]
        check(f"dense {j} seed", row["seed"] == seeds[4 + j])
        check(f"dense {j} sample size", len(x) == r["trials"])
        check(f"dense {j} raw summary", mc == row["mc"])
        check(
            f"dense {j} Gaussian entropy agreement",
            row["passed"] and abs(mc["mean"] - expected) < 4 * mc["se"],
        )
        if regenerate:
            check(
                f"dense {j} exact sample regeneration",
                np.array_equal(x, dense_figure4(row["db"], row["n"], r["trials"], row["seed"])),
            )
    for i, row in enumerate(r["energy_diagnostic"]):
        e = 10 ** (-row["loss_db"] / 10)
        prediction = PhaseKennedy(row["photons"], 0.9 * e, e, 0.01, 0.1).kl(
            PhaseKennedy(row["photons"], e, e, 0.01, 0.1)
        )
        a = next(x for x in anchors if x["loss_db"] == row["loss_db"])
        check(f"energy {i} recompute", np.isclose(prediction, row["predicted_kl"], rtol=1e-12))
        check(
            f"energy {i} source anchor", a["approx_published_kennedy"] == row["approx_published_kl"]
        )
        check(
            f"energy {i} fit versus holdout",
            row["used_for_fit"] == (row["case"] == "one_point_energy_fit" and row["loss_db"] == 0),
        )
    fit = [x for x in r["energy_diagnostic"] if x["used_for_fit"]]
    check("one point only fitted", len(fit) == 1 and abs(fit[0]["relative_difference"]) < 1e-10)
    check("distinct seeds", len(set(seeds)) == 10)
    with (root / "energy_diagnostic.csv").open() as f:
        rows = list(csv.DictReader(f))
    check(
        "CSV matches report",
        rows == [{k: str(v) for k, v in row.items()} for row in r["energy_diagnostic"]],
    )
    check("reported pass", r["passed"] and r["source_unchanged"])
    for ext in ("png", "svg"):
        check("plot " + ext, (root / f"independent_checks.{ext}").stat().st_size > 1000)
    result = {
        "passed": all(x["passed"] for x in checks),
        "count": len(checks),
        "regenerated": regenerate,
        "checks": checks,
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root")
    p.add_argument("--regenerate", action="store_true")
    a = p.parse_args()
    r = audit(a.root, a.regenerate)
    print(f"{sum(x['passed'] for x in r['checks'])}/{r['count']} checks passed")
    raise SystemExit(0 if r["passed"] else 1)


if __name__ == "__main__":
    main()
