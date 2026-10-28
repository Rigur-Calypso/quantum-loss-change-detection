"""Audit quantum-state samples, first-alarm clocks and dual-rail ARL cycles."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .dual_rail import DualRailRuns, simulate_heralds, summarize_dual_rail
from .rail_study import state_checks
from .probes import Bernoulli
from .iid_cusum import estimate_iid_arl
from .mismatch import read_cycles, summarize_mismatch_cycles, cycle_arrays
from .simulation import Runs


def audit(root, regenerate=False):
    """Check qubit states, pair/bin timing, calibration holdouts and saved heralding trajectories.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    r = json.loads((root / "report.json").read_text())
    s = r["settings"]
    checks = []

    def check(name, value):
        checks.append({"name": name, "passed": bool(value)})

    seeds = [
        int(c.generate_state(1, dtype=np.uint64)[0])
        for c in np.random.SeedSequence(s["seed"]).spawn(19)
    ]
    pre, post = Bernoulli(s["eta1"]), Bernoulli(s["eta2"])
    h = 8.0
    check(
        "physical target uses two time bins per qubit",
        s["target_arl_physical_bins"] == 2 * s["target_arl_pairs"],
    )
    check("six calibration rounds", len(r["calibration"]) == 6)
    for i, entry in enumerate(r["calibration"]):
        check(
            f"calibration {i} threshold and seed",
            entry["threshold"] == h and entry["summary"]["seed"] == seeds[i],
        )
        if regenerate:
            fresh, _, _ = estimate_iid_arl(pre, post, h, s["cycles"], seeds[i])
            check(f"calibration {i} exact regeneration", fresh == entry["summary"])
        h = max(0.0, h + np.log(s["target_arl_pairs"] / entry["summary"]["arl"]))
    check("frozen threshold", h == r["threshold"])
    for i, saved in enumerate(r["validation"]):
        with np.load(root / f"validation_{i}.npz", allow_pickle=False) as raw:
            a, b = [read_cycles(raw, p) for p in ("null", "proposal")]
            fresh = summarize_mismatch_cycles(a, b)
            check(
                f"holdout {i} complete independent cycles",
                len(a.length) == s["cycles"]
                and len(b.length) == s["cycles"]
                and np.all(a.complete)
                and np.all(b.complete),
            )
            check(f"holdout {i} raw summary", all(fresh[k] == saved[k] for k in fresh))
            check(f"holdout {i} seed", saved["seed"] == seeds[6 + i])
            check(
                f"holdout {i} ARL match",
                saved["diagnostics_passed"]
                and 0.95e6 <= saved["approx_95_interval"][0]
                and saved["approx_95_interval"][1] <= 1.05e6,
            )
            if regenerate:
                again, c, d = estimate_iid_arl(pre, post, h, s["cycles"], seeds[6 + i])
                check(
                    f"holdout {i} exact samples",
                    again == saved
                    and all(np.array_equal(raw[k], v) for k, v in cycle_arrays(c, d).items()),
                )
    check("holdout repetition and status", r["arl_confirmed"] and r["validation_repeat_z"] <= 4)
    control = r["control"]
    with np.load(root / "control.npz", allow_pickle=False) as raw:
        a, b = [read_cycles(raw, p) for p in ("null", "proposal")]
        fresh = summarize_mismatch_cycles(a, b)
        run = DualRailRuns(Runs(raw["stop"], raw["event"], 200000, None), None)
        check(
            "control raw weighted summary", all(fresh[k] == control["weighted"][k] for k in fresh)
        )
        check("control raw physical summary", summarize_dual_rail(run) == control["direct"])
        check("control direct comparison", control["passed"] and control["z"] <= 4)
        if regenerate:
            again, c, d = estimate_iid_arl(pre, post, 2, s["cycles"], seeds[8])
            direct = simulate_heralds(
                0.9, 0.81, 2, trials=8000, horizon_pairs=200000, seed=seeds[9]
            )
            check(
                "control exact weighted regeneration",
                again == control["weighted"]
                and all(np.array_equal(raw[k], v) for k, v in cycle_arrays(c, d).items()),
            )
            check(
                "control exact direct regeneration",
                np.array_equal(direct.pairs.stop, run.pairs.stop)
                and np.array_equal(direct.pairs.event, run.pairs.event),
            )
    check(
        "complete encoding/phase grid",
        [(c["change_pulse"], c["encoding"]) for c in r["changed"]]
        == [(cp, e) for cp in (2001, 2002) for e in ("early", "late", "plus", "haar")],
    )
    for i, c in enumerate(r["changed"]):
        cp = c["change_pulse"]
        with np.load(root / c["file"], allow_pickle=False) as raw:
            run = DualRailRuns(
                Runs(raw["stop"], raw["event"], s["horizon_pairs"], (cp + 1) // 2), cp
            )
            summary = summarize_dual_rail(run)
            check(f"case {i} seed", c["seed"] == seeds[10 + i])
            check(f"case {i} raw summary", summary == c["summary"])
            false = run.pairs.event & (run.readout_pulses < cp)
            survivors = ~false
            check(
                f"case {i} physical false-alarm classification",
                int(false.sum()) == summary["pair_statistics"]["false_alarms"],
            )
            check(
                f"case {i} physical delay",
                float(np.mean(run.readout_pulses[survivors] - cp))
                == summary["physical_time_bin_statistics"]["conditional_delay"],
            )
            check(
                f"case {i} complete and passed",
                c["passed"] and summary["pair_statistics"]["post_change_censored"] == 0,
            )
            if regenerate:
                fresh = simulate_heralds(
                    0.9, 0.81, h, s["trials"], s["horizon_pairs"], cp, c["encoding"], c["seed"]
                )
                check(
                    f"case {i} exact samples",
                    np.array_equal(fresh.pairs.stop, run.pairs.stop)
                    and np.array_equal(fresh.pairs.event, run.pairs.event),
                )
    with np.load(root / "states.npz", allow_pickle=False) as raw:
        fresh = state_checks(seeds[18], s["state_count"])
        check(
            "all state samples regenerate", all(np.array_equal(raw[k], v) for k, v in fresh.items())
        )
        check(
            "state survival and fidelity",
            np.max(abs(raw["uniform_probability"] - 0.9)) < 1e-12
            and np.max(abs(raw["uniform_fidelity"] - 1)) < 1e-12,
        )
    with (root / "delays.csv").open() as f:
        rows = list(csv.DictReader(f))
    check("delay CSV size", len(rows) == 8)
    for i, (row, c) in enumerate(zip(rows, r["changed"])):
        check(
            f"CSV physical delay {i}",
            float(row["conditional_delay_physical_bins"])
            == c["summary"]["physical_time_bin_statistics"]["conditional_delay"],
        )
    check("overall status", r["passed"])
    result = {
        "passed": all(c["passed"] for c in checks),
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
    print(f"{sum(c['passed'] for c in r['checks'])}/{r['count']} checks passed")
    raise SystemExit(0 if r["passed"] else 1)


if __name__ == "__main__":
    main()
