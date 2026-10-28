"""Recalculate saved figure tables and check numerical invariants."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .verified_figures import figure5, figure12, fixed_seed


def audit(root, pdf=None):
    """Recalculate verified-setting curves, refinement checks and sparse digitization comparisons.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    checks = []

    def check(name, value):
        checks.append({"name": name, "passed": bool(value)})

    if pdf:
        check(
            "original PDF checksum",
            hashlib.sha256(Path(pdf).read_bytes()).hexdigest() == report["source_sha256"],
        )
    tables = {}
    for name, expected in (("figure5", 61), ("figure12", 100), ("figure4_fixed_seed", 657)):
        with (root / f"{name}.csv").open() as f:
            rows = list(csv.DictReader(f))
        tables[name] = rows
        check(name + " row count", len(rows) == expected)
        for i, row in enumerate(rows):
            if name == "figure5":
                fresh = figure5(float(row["loss_db"]))
            elif name == "figure12":
                fresh = figure12(float(row["loss_db"]))
            else:
                fresh = fixed_seed(float(row["seed_db"]), int(row["n"]))
            check(
                f"{name} recomputed row {i}",
                row.keys() == fresh.keys()
                and all(
                    np.isclose(float(row[k]), v, rtol=1e-10, atol=1e-13) for k, v in fresh.items()
                ),
            )
    for i, row in enumerate(report["refinement"]):
        fresh = figure5(row["loss_db"], 160)
        base = figure5(row["loss_db"])
        change = abs(fresh[row["curve"]] / base[row["curve"]] - 1)
        check(
            f"phase refinement {i}",
            change < 1e-7 and np.isclose(change, row["relative_change"], atol=1e-13),
        )
    check("fixed seed identity", report["fixed_seed_max_identity_error"] < 1e-12)
    for name, key in (
        ("figure5_kennedy_anchors", "figure5_kennedy_digitization"),
        ("figure12_anchors", "figure12_digitization"),
        ("phase_refinement", "refinement"),
    ):
        with (root / f"{name}.csv").open() as f:
            rows = list(csv.DictReader(f))
        check(
            name + " CSV matches report",
            rows == [{k: str(v) for k, v in row.items()} for row in report[key]],
        )
    for i, row in enumerate(tables["figure12"]):
        check(
            f"weak photon data processing {i}",
            float(row["photon_counting"]) >= float(row["photon_homodyne"]),
        )
    for name in ("verified_receiver_curves", "fixed_seed_equation_check"):
        for ext in ("png", "svg"):
            check(name + "." + ext, (root / f"{name}.{ext}").stat().st_size > 1000)
    result = {"passed": all(r["passed"] for r in checks), "checks": checks, "count": len(checks)}
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root")
    p.add_argument("--pdf")
    a = p.parse_args()
    r = audit(a.root, a.pdf)
    print(f"{sum(c['passed'] for c in r['checks'])}/{r['count']} checks passed")
    raise SystemExit(0 if r["passed"] else 1)


if __name__ == "__main__":
    main()
