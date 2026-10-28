"""Audit actual-null correlated-cycle weights, summaries and seeded samples."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel import sha256
from .block_null import models_for_study, parameters, read_cycle_pair, exported_rows
from .corr_renewal import (
    summarize_conditional_cycles,
    estimate_conditional_arl,
    conditional_cycle_arrays,
)
from .corr_cusum import simulate_conditional
from .batch_renewal import pool_cycle_summaries
from .simulation import Runs, summarize_delay, summarize_arl


def audit(root, regenerate=False):
    """Check actual-null weights, independent repetitions and direct correlated controls.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    s = report["settings"]
    source = json.loads((Path(report["source_root"]) / "report.json").read_text())
    models = models_for_study(source)
    checks = []

    def check(name, passed):
        checks.append({"name": name, "passed": bool(passed)})

    for path, digest in report["source_hashes"].items():
        check("source " + path, sha256(path) == digest)
    check(
        "coverage", len(report["rows"]) == 9 and len(report["controls"]) == 6 and len(models) == 3
    )
    check(
        "fixed scope",
        s["true_eta1_grid"] == [0.899, 0.9, 0.901]
        and (s["N"], s["Na"], s["assumed_eta1"], s["eta2"]) == (100, 0.1, 0.9, 0.85)
        and s["proposal"] == "adaptive"
        and s["proposal_offset"] == 2,
    )
    children = np.random.SeedSequence(s["seed"]).spawn(39)
    seed = lambda i: int(children[i].generate_state(1, dtype=np.uint64)[0])
    for i, row in enumerate(report["rows"]):
        m = models[i // 3]
        eta = s["true_eta1_grid"][i % 3]
        prefix = f"row {i}: "
        check(
            prefix + "model mapping",
            row["n"] == m["n"]
            and row["label"] == m["label"]
            and row["true_eta1"] == eta
            and row["change_at"] == m["change_at"]
            and row["threshold"] == m["calibration"]["threshold"],
        )
        check(
            prefix + "independent seeds",
            row["seeds"] == [seed(3 * i), seed(3 * i + 1)]
            and row["changed_seed"] == seed(3 * i + 2),
        )
        for repeat in range(2):
            saved = row["repeats"][repeat]
            with np.load(root / row["cycle_files"][repeat]) as data:
                pre, post = read_cycle_pair(data, row["n"], eta, row["threshold"])
                stats = summarize_conditional_cycles(pre, post)
                check(
                    prefix + f"raw cycle counts {repeat}",
                    len(pre.length) == len(post.length) == s["cycles_per_law_per_repeat"],
                )
                check(
                    prefix + f"raw summaries {repeat}", all(saved[k] == v for k, v in stats.items())
                )
                check(
                    prefix + f"correct actual-null metadata {repeat}",
                    saved["parameters"]
                    == {
                        k: v
                        for k, v in parameters(row["n"], eta).items()
                        if k not in ("kind", "n") and (k != "sampling_eta1" or eta != 0.9)
                    },
                )
                expected = [
                    int(c.generate_state(1, dtype=np.uint64)[0])
                    for c in np.random.SeedSequence(row["seeds"][repeat]).spawn(2)
                ]
                check(
                    prefix + f"independent law streams {repeat}",
                    saved["seed"] == row["seeds"][repeat]
                    and saved["stream_seeds"] == expected
                    and expected[0] != expected[1],
                )
                check(prefix + f"weight diagnostics {repeat}", stats["diagnostics_passed"])
                if regenerate:
                    _, a, b = estimate_conditional_arl(
                        threshold=row["threshold"],
                        cycles=s["cycles_per_law_per_repeat"],
                        seed=row["seeds"][repeat],
                        proposal="adaptive",
                        proposal_offset=2.0,
                        **parameters(row["n"], eta),
                    )
                    check(
                        prefix + f"exact cycle regeneration {repeat}",
                        all(
                            np.array_equal(data[k], v)
                            for k, v in conditional_cycle_arrays(a, b).items()
                        ),
                    )
        pooled = pool_cycle_summaries(row["repeats"])
        check(prefix + "pooled moments", pooled == row["pooled"])
        z = abs(row["repeats"][0]["arl"] - row["repeats"][1]["arl"]) / np.hypot(
            row["repeats"][0]["arl_se"], row["repeats"][1]["arl_se"]
        )
        check(prefix + "repeatability", z == row["repeat_z"] and z <= 4)
        with np.load(root / row["changed_file"]) as data:
            runs = Runs(data["stop"], data["event"], s["horizon"], row["change_at"])
        changed = summarize_delay(runs)
        check(
            prefix + "raw changed runs",
            len(runs.stop) == s["changed_trials"] and changed == row["changed"],
        )
        check(prefix + "no censored survivor", changed["conditional_delay"] is not None)
        check(prefix + "classification", row["checks_passed"] and pooled["diagnostics_passed"])
        if regenerate:
            fresh = simulate_conditional(
                threshold=row["threshold"],
                trials=s["changed_trials"],
                horizon=s["horizon"],
                change_at=row["change_at"],
                seed=row["changed_seed"],
                **parameters(row["n"], eta),
            )
            check(
                prefix + "exact changed regeneration",
                np.array_equal(fresh.stop, runs.stop) and np.array_equal(fresh.event, runs.event),
            )
    for i, row in enumerate(report["controls"]):
        prefix = f"control {i}: "
        eta = (0.899, 0.901)[i % 2]
        n = models[i // 2]["n"]
        check(
            prefix + "mapping and seeds",
            row["n"] == n
            and row["true_eta1"] == eta
            and row["threshold"] == 3
            and row["direct_seed"] == seed(27 + 2 * i)
            and row["estimate_seed"] == seed(28 + 2 * i),
        )
        with np.load(root / row["file"]) as data:
            pre, post = read_cycle_pair(data, n, eta, 3.0)
            stats = summarize_conditional_cycles(pre, post)
            runs = Runs(data["stop"], data["event"], s["control_horizon"], None)
            check(
                prefix + "cycle counts", len(pre.length) == len(post.length) == s["control_cycles"]
            )
            check(
                prefix + "raw weighted summaries",
                all(row["estimate"][k] == v for k, v in stats.items()),
            )
            if regenerate:
                _, a, b = estimate_conditional_arl(
                    threshold=3,
                    cycles=s["control_cycles"],
                    seed=row["estimate_seed"],
                    proposal="adaptive",
                    proposal_offset=2.0,
                    **parameters(n, eta),
                )
                check(
                    prefix + "exact cycle regeneration",
                    all(
                        np.array_equal(data[k], v)
                        for k, v in conditional_cycle_arrays(a, b).items()
                    ),
                )
        direct = summarize_arl(runs)
        check(
            prefix + "raw direct runs",
            len(runs.stop) == s["control_trials"] and direct == row["direct"],
        )
        z = (
            None
            if direct["arl"] is None
            else float(
                abs(direct["arl"] - stats["arl"]) / np.hypot(direct["arl_se"], stats["arl_se"])
            )
        )
        check(
            prefix + "independent direct agreement",
            z is not None
            and z <= 4
            and z == row["comparison_z"]
            and stats["diagnostics_passed"]
            and row["checks_passed"],
        )
        if regenerate:
            fresh = simulate_conditional(
                threshold=3,
                trials=s["control_trials"],
                horizon=s["control_horizon"],
                seed=row["direct_seed"],
                **parameters(n, eta),
            )
            check(
                prefix + "exact direct regeneration",
                np.array_equal(fresh.stop, runs.stop) and np.array_equal(fresh.event, runs.event),
            )
    with (root / "sensitivity.csv").open() as handle:
        exported = list(csv.DictReader(handle))
    check("CSV", exported == [{k: str(v) for k, v in row.items()} for row in exported_rows(report)])
    check("status", report["source_unchanged"] and report["all_checks_passed"])
    result = {
        "passed": all(c["passed"] for c in checks),
        "regenerated": regenerate,
        "checks": checks,
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    if not result["passed"]:
        raise AssertionError([c["name"] for c in checks if not c["passed"]])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", nargs="?", default="results_entangled_null_diagnostics")
    p.add_argument("--regenerate", action="store_true")
    a = p.parse_args()
    audit(a.root, a.regenerate)


if __name__ == "__main__":
    main()
